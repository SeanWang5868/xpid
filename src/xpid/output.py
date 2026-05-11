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

# Columns included in non-verbose output
SIMPLE_COLS = [
    'pdb', 'resolution',
    'pi_chain', 'pi_res', 'pi_id',
    'X_chain', 'X_res', 'X_id', 'X_atom', 'H_atom',
    'H_source', 'dist_X_Pi', 'is_plevin', 'is_hudson',
    'is_trp_5ring_acceptor', 'is_pi_pi_tshaped', 'sym_op'
]


class ResultStreamer:
    """Context-managed streaming writer for detection results.

    Supports JSON, CSV, and Parquet output formats. Writes results
    incrementally to avoid holding full datasets in memory.
    """

    def __init__(self, output_path: Path, file_type: str, verbose: bool):
        self.output_path = output_path
        self.file_type = file_type.lower()
        self.verbose = verbose
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
                headers = results[0].keys() if self.verbose else SIMPLE_COLS
                self.csv_writer = csv.DictWriter(self.file_handle, fieldnames=headers)
                self.csv_writer.writeheader()
                self.is_first_chunk = False
            rows = results if self.verbose else [{k: r[k] for k in SIMPLE_COLS} for r in results]
            self.csv_writer.writerows(rows)

        elif self.file_type == 'json':
            comma = '' if self.is_first_chunk else ',\n'
            for r in results:
                self.file_handle.write(comma + json.dumps(r, indent=2 if self.verbose else None))
                comma = ',\n'
            self.is_first_chunk = False

        elif self.file_type == 'parquet':

            df = self.pd.DataFrame(results)
            table = self.pa.Table.from_pandas(df if self.verbose else df[SIMPLE_COLS])
            if self.is_first_chunk:
                self.parquet_writer = self.pq.ParquetWriter(self.output_path, table.schema)
                self.is_first_chunk = False
            self.parquet_writer.write_table(table)
