"""
output.py
Result streaming to JSON, CSV, and Parquet formats.
"""
import csv
import json
import logging
import queue
import sys
import threading
from pathlib import Path
from typing import List, Dict, Any

from . import schema

logger = logging.getLogger("xpid.output")

# Columns included in non-verbose output by default. P-slab columns are added
# only when explicitly requested by the CLI/API caller.
BASE_SIMPLE_COLS = schema.SIMPLE_NAMES

P_GEOMETRY_SIMPLE_COLS = [c for c in schema.P_SLAB_NAMES if c != 'is_p_slab']

P_SLAB_SIMPLE_COLS = schema.P_SLAB_NAMES

CANDIDATE_SIMPLE_COLS = schema.CANDIDATE_NAMES

COORDINATE_SIMPLE_COLS = [c for c in schema.COORDS_NAMES if c in schema._FIELD_MAP]

SASA_SIMPLE_COLS = schema.SASA_NAMES
COOP_SIMPLE_COLS = schema.COOP_NAMES
HBOND_SIMPLE_COLS = schema.HBOND_NAMES


SIMPLE_COLS = BASE_SIMPLE_COLS
P_SLAB_LABEL_KEYS = {'is_p_slab'}
P_GEOMETRY_OUTPUT_KEYS = set(P_GEOMETRY_SIMPLE_COLS)
P_SLAB_OUTPUT_KEYS = set(P_SLAB_SIMPLE_COLS)
CANDIDATE_OUTPUT_KEYS = set(CANDIDATE_SIMPLE_COLS)
COORDINATE_OUTPUT_KEYS = set(COORDINATE_SIMPLE_COLS)

FLOAT_COLS = {
    'resolution', 'dist_X_centroid', 'dist_X_Pi', 'proj_dist', 'theta',
    'angle_XPCN', 'angle_XH_Pi', 'P_radius', 'P_slab_half_thickness',
    'h_proj_dist', 'H_ray_t', 'H_ray_entry_dist', 'h_plane_proj_dist',
    'H_plane_t', 'H_plane_entry_dist', 'delta_h_proj_dist',
    'xh_centroid_cos', 'xh_lateral_inward_score',
    'pi_avg_b', 'pi_center_x', 'pi_center_y', 'pi_center_z',
    'pi_normal_x', 'pi_normal_y', 'pi_normal_z',
    'X_b', 'X_xyz_x', 'X_xyz_y', 'X_xyz_z',
    'H_xyz_x', 'H_xyz_y', 'H_xyz_z',
    'pi_sasa_avg', 'X_sasa', 'H_sasa',
    'hbond_HA_dist', 'hbond_DHA_angle', 'hbond_vs_xhpi_score',
    'pi_sasa_avg', 'X_sasa', 'H_sasa',
}

INT_COLS = set(schema._INT_NAMES)


class ResultStreamer:
    """Context-managed streaming writer for detection results.

    Supports JSON, CSV, and Parquet output formats. Writes results
    incrementally to avoid holding full datasets in memory.
    """

    def __init__(self, output_path: Path, file_type: str, verbose: bool,
                 include_p_slab: bool = False,
                 include_xh_candidates: bool = False,
                 include_coordinates: bool = False,
                 include_sasa: bool = False,
                 include_cooperativity: bool = False,
                 include_hbond_comp: bool = False):
        self.output_path = output_path
        self.file_type = file_type.lower()
        self.verbose = verbose
        self.include_p_slab = include_p_slab
        self.include_xh_candidates = include_xh_candidates
        self.include_coordinates = include_coordinates
        self.include_sasa = include_sasa
        self.include_cooperativity = include_cooperativity
        self.include_sasa = include_sasa
        self.include_cooperativity = include_cooperativity
        self.include_hbond_comp = include_hbond_comp
        self.file_handle = None
        self.csv_writer = None
        self.parquet_writer = None
        self.parquet_schema = None
        self.parquet_columns = None
        self.is_first_chunk = True
        self._write_queue: "queue.Queue[List[Dict[str, Any]]]" = queue.Queue(maxsize=200)
        self._writer_thread: threading.Thread | None = None
        self._writer_error: BaseException | None = None
        self._sentinel = object()

        # Validate parquet dependencies early
        if self.file_type == 'parquet':
            try:
                import pandas as pd
                import pyarrow as pa
                import pyarrow.parquet as pq

                self.pd = pd
                self.pa = pa
                self.pq = pq
            except ImportError:
                logger.error("To use --file-type parquet, install 'pandas' and 'pyarrow'.")
                logger.error("Try: pip install 'xpid[parquet]'")
                sys.exit(1)

    def __enter__(self):
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

        if self.file_type in ('json', 'csv'):
            self.file_handle = open(self.output_path, 'w', newline='', encoding='utf-8')
            if self.file_type == 'json':
                self.file_handle.write('[\n')

        # The writer thread is *not* started here.  It is lazily created
        # on the first write_chunk call, which happens after the
        # multiprocessing Pool has been forked.  This keeps the process
        # free of threads at fork time.
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._write_queue.put(self._sentinel)
        if self._writer_thread is not None:
            self._writer_thread.join(timeout=120)

        if self.file_type == 'json' and self.file_handle:
            self.file_handle.write('\n]')

        if self.file_handle:
            self.file_handle.close()

        if self.parquet_writer:
            self.parquet_writer.close()

    def _ensure_writer(self) -> None:
        """Start the background writer thread if it has not been started yet.

        Must be called after the multiprocessing Pool has been forked
        (i.e. from the main loop, not from ``__enter__``).
        """
        if self._writer_thread is not None:
            return
        self._writer_thread = threading.Thread(
            target=self._run_writer, daemon=True, name="xpid-writer")
        self._writer_thread.start()

    def write_chunk(self, results: List[Dict[str, Any]]):
        """Enqueue a batch of result dicts for background writing.

        Raises :exc:`RuntimeError` if the background writer thread has
        died, turning what would otherwise be a silent deadlock into a
        visible error.
        """
        if not results:
            return
        self._ensure_writer()
        if self._writer_error is not None:
            raise RuntimeError(
                "Background writer thread has failed; "
                "output is incomplete."
            ) from self._writer_error
        # Use a timeout so a full queue never blocks the main thread
        # indefinitely — the writer thread may have died without setting
        # _writer_error (e.g. the process is being killed).
        while True:
            try:
                self._write_queue.put(results, timeout=5.0)
                break
            except queue.Full:
                if self._writer_error is not None:
                    raise RuntimeError(
                        "Background writer thread has failed; "
                        "output is incomplete."
                    ) from self._writer_error
                if (self._writer_thread is not None
                        and not self._writer_thread.is_alive()):
                    raise RuntimeError(
                        "Background writer thread died unexpectedly; "
                        "output is incomplete."
                    )

    def _run_writer(self) -> None:
        """Background thread: drain the write queue and persist batches."""
        try:
            while True:
                item = self._write_queue.get()
                if item is self._sentinel:
                    break
                self._write_batch(item)
        except BaseException as exc:
            self._writer_error = exc
            logger.error("Background writer thread failed: %s", exc)
            # Drain remaining items so the main thread never blocks on put().
            while True:
                try:
                    item = self._write_queue.get(timeout=0.1)
                except queue.Empty:
                    break
                if item is self._sentinel:
                    break

    def _write_batch(self, results: List[Dict[str, Any]]) -> None:
        """Write a single batch of result dicts to the output file.

        This is called exclusively by the background writer thread.
        """
        if not results:
            return

        if self.file_type == 'csv':
            if self.is_first_chunk:
                headers = self._headers(results[0])
                self.csv_writer = csv.DictWriter(self.file_handle, fieldnames=headers)
                self.csv_writer.writeheader()
                self.is_first_chunk = False
            rows = [self._row_for_output(r) for r in results]
            self.csv_writer.writerows(rows)

        elif self.file_type == 'json':
            serialized = ",\n".join(
                json.dumps(self._row_for_output(r),
                           indent=2 if self.verbose else None)
                for r in results
            )
            if serialized:
                prefix = "" if self.is_first_chunk else ",\n"
                self.file_handle.write(prefix + serialized)
            self.is_first_chunk = False

        elif self.file_type == 'parquet':

            df = self._dataframe_for_parquet(results)
            if self.is_first_chunk:
                self.parquet_columns = list(df.columns)
                self.parquet_schema = self._arrow_schema(
                    self.parquet_columns)
            else:
                unexpected = set(df.columns) - set(self.parquet_columns)
                if unexpected:
                    raise ValueError(
                        "Parquet output columns changed between chunks: "
                        + ", ".join(sorted(unexpected)))
                df = df.reindex(columns=self.parquet_columns)
            table = self.pa.Table.from_pandas(
                df, schema=self.parquet_schema, preserve_index=False)
            if self.is_first_chunk:
                self.parquet_writer = self.pq.ParquetWriter(self.output_path, table.schema)
                self.is_first_chunk = False
            self.parquet_writer.write_table(table)

    def _dataframe_for_parquet(self, results: List[Dict[str, Any]]):
        df = self.pd.DataFrame([self._row_for_output(r) for r in results])

        for col in FLOAT_COLS.intersection(df.columns):
            df[col] = self.pd.to_numeric(df[col], errors='coerce')
        for col in INT_COLS.intersection(df.columns):
            df[col] = self.pd.to_numeric(df[col], errors='coerce').astype('Int64')

        return df

    def _arrow_schema(self, columns: List[str]):
        """Build a stable Arrow schema from the central output registry."""
        arrow_fields = []
        arrow_types = {
            "str": self.pa.large_string(),
            "float": self.pa.float64(),
            "int": self.pa.int64(),
        }
        for column in columns:
            descriptor = schema.field(column)
            if descriptor is None:
                raise ValueError(
                    f"No declared output type for Parquet column {column!r}")
            arrow_fields.append(
                self.pa.field(column, arrow_types[descriptor.dtype]))
        return self.pa.schema(arrow_fields)

    def _simple_cols(self) -> List[str]:
        def with_coordinates(cols: List[str]) -> List[str]:
            if not self.include_coordinates:
                return cols
            return cols + [col for col in COORDINATE_SIMPLE_COLS if col not in cols]

        after_x_chain = BASE_SIMPLE_COLS.index("X_chain") + 1
        if self.include_xh_candidates:
            cols = (
                BASE_SIMPLE_COLS[:after_x_chain] +
                CANDIDATE_SIMPLE_COLS[:1] +
                (['is_p_slab'] if self.include_p_slab else []) +
                BASE_SIMPLE_COLS[after_x_chain:] +
                CANDIDATE_SIMPLE_COLS[1:] +
                P_GEOMETRY_SIMPLE_COLS
            )
            return with_coordinates(cols)
        if self.include_p_slab:
            cols = (
                BASE_SIMPLE_COLS[:after_x_chain] +
                P_SLAB_SIMPLE_COLS[:1] +
                BASE_SIMPLE_COLS[after_x_chain:] +
                P_SLAB_SIMPLE_COLS[1:]
            )
            return with_coordinates(cols)
        cols = list(BASE_SIMPLE_COLS)
        if self.include_sasa:
            cols.extend(SASA_SIMPLE_COLS)
        if self.include_cooperativity:
            cols.extend(COOP_SIMPLE_COLS)
        if self.verbose:
            cols.extend(HBOND_SIMPLE_COLS)
        return with_coordinates(cols)

    def _headers(self, first_result: Dict[str, Any]):
        if not self.verbose:
            return self._simple_cols()
        return self._row_for_output(first_result).keys()

    def _row_for_output(self, row: Dict[str, Any]) -> Dict[str, Any]:
        if self.verbose:
            if self.include_sasa:
                return row
            if self.include_cooperativity:
                return row
            if self.include_hbond_comp:
                return row
            if self.include_p_slab:
                return row
            if self.include_xh_candidates:
                return {k: v for k, v in row.items() if k not in P_SLAB_LABEL_KEYS}
            return {k: v for k, v in row.items() if k not in P_SLAB_OUTPUT_KEYS}
        return {k: row.get(k) for k in self._simple_cols()}
