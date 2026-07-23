"""
output.py
Result streaming to JSON, CSV, and Parquet formats.
"""
import csv
import json
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any

logger = logging.getLogger("xpid.output")

# Columns included in non-verbose output by default. P-slab columns are added
# only when explicitly requested by the CLI/API caller.
BASE_SIMPLE_COLS = [
    'pdb', 'resolution',
    'pi_chain', 'pi_res', 'pi_id',
    'X_chain', 'X_res', 'X_id', 'X_atom', 'H_atom',
    'H_source',
    'is_hudson', 'is_plevin',
    'dist_X_centroid', 'dist_X_Pi',
    'proj_dist', 'theta', 'angle_XPCN', 'angle_XH_Pi',
    'is_trp_5ring_acceptor', 'is_pi_pi_tshaped', 'sym_op'
]

P_GEOMETRY_SIMPLE_COLS = [
    'P_radius', 'P_slab_half_thickness', 'h_proj_dist', 'H_ray_t',
    'H_ray_entry_dist', 'h_plane_proj_dist', 'H_plane_t',
    'H_plane_entry_dist', 'delta_h_proj_dist',
]

P_SLAB_SIMPLE_COLS = ['is_p_slab'] + P_GEOMETRY_SIMPLE_COLS

CANDIDATE_SIMPLE_COLS = [
    'is_xh_candidate', 'is_hudson_spatial', 'is_plevin_spatial',
    'hudson_dist_ok', 'hudson_proj_ok', 'hudson_direction_ok',
    'plevin_dist_ok', 'plevin_xpcn_ok', 'plevin_direction_ok',
    'xh_centroid_cos', 'xh_lateral_inward_score',
]

COORDINATE_SIMPLE_COLS = [
    'pi_center_x', 'pi_center_y', 'pi_center_z',
    'pi_normal_x', 'pi_normal_y', 'pi_normal_z',
    'X_xyz_x', 'X_xyz_y', 'X_xyz_z',
    'H_xyz_x', 'H_xyz_y', 'H_xyz_z',
    'pi_sasa_avg', 'X_sasa', 'H_sasa',
    'X_side_of_pi',
]

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
}

INT_COLS = {
    'is_hudson', 'is_plevin', 'is_p_slab', 'is_trp_5ring_acceptor',
    'is_pi_pi_tshaped', 'seq_sep', 'sym_op', 'is_xh_candidate',
    'is_hudson_spatial', 'is_plevin_spatial', 'hudson_dist_ok',
    'hudson_proj_ok', 'hudson_direction_ok', 'plevin_dist_ok',
    'plevin_xpcn_ok', 'plevin_direction_ok', 'X_side_of_pi',
}


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
                 include_cooperativity: bool = False):
        self.output_path = output_path
        self.file_type = file_type.lower()
        self.verbose = verbose
        self.include_p_slab = include_p_slab
        self.include_xh_candidates = include_xh_candidates
        self.include_coordinates = include_coordinates
        self.include_sasa = include_sasa
        self.include_cooperativity = include_cooperativity
        self.file_handle = None
        self.csv_writer = None
        self.parquet_writer = None
        self.is_first_chunk = True

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

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_type == 'json' and self.file_handle:
            self.file_handle.write('\n]')

        if self.file_handle:
            self.file_handle.close()

        if self.parquet_writer:
            self.parquet_writer.close()

    def write_chunk(self, results: List[Dict[str, Any]]):
        """Write a batch of result dicts to the output file."""
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
            comma = '' if self.is_first_chunk else ',\n'
            for r in results:
                self.file_handle.write(comma + json.dumps(
                    self._row_for_output(r), indent=2 if self.verbose else None))
                comma = ',\n'
            self.is_first_chunk = False

        elif self.file_type == 'parquet':

            df = self._dataframe_for_parquet(results)
            table = self.pa.Table.from_pandas(df, preserve_index=False)
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

    def _simple_cols(self) -> List[str]:
        def with_coordinates(cols: List[str]) -> List[str]:
            if not self.include_coordinates:
                return cols
            return cols + [col for col in COORDINATE_SIMPLE_COLS if col not in cols]

        if self.include_xh_candidates:
            cols = (
                BASE_SIMPLE_COLS[:11] +
                CANDIDATE_SIMPLE_COLS[:1] +
                (['is_p_slab'] if self.include_p_slab else []) +
                BASE_SIMPLE_COLS[11:] +
                CANDIDATE_SIMPLE_COLS[1:] +
                P_GEOMETRY_SIMPLE_COLS
            )
            return with_coordinates(cols)
        if self.include_p_slab:
            cols = BASE_SIMPLE_COLS[:11] + P_SLAB_SIMPLE_COLS[:1] + BASE_SIMPLE_COLS[11:] + P_SLAB_SIMPLE_COLS[1:]
            return with_coordinates(cols)
        if self.include_sasa:
            cols = cols + SASA_SIMPLE_COLS
        if self.include_cooperativity:
            cols = cols + COOP_SIMPLE_COLS
        return with_coordinates(BASE_SIMPLE_COLS)

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
            if self.include_p_slab:
                return row
            if self.include_xh_candidates:
                return {k: v for k, v in row.items() if k not in P_SLAB_LABEL_KEYS}
            return {k: v for k, v in row.items() if k not in P_SLAB_OUTPUT_KEYS}
        return {k: row.get(k) for k in self._simple_cols()}
