"""
cli.py
Command Line Interface for xpid.
Handles argument parsing, validation, execution flow, and streaming output.
"""
import argparse
import sys
import re
import multiprocessing
import logging
import gemmi
import json
import csv
import os
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional, Set

# Ensure package is accessible
try:
    from xpid import prep, core, config
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from xpid import prep, core, config

# --- Constants ---
H_MODE_MAP = {
    0: "NoChange", 1: "Shift", 2: "Remove", 
    3: "ReAdd", 4: "ReAddButWater", 5: "ReAddKnown"
}

SIMPLE_COLS = [
    'pdb', 'resolution', 
    'pi_chain', 'pi_res', 'pi_id', 
    'X_chain', 'X_res', 'X_id', 'X_atom', 'H_atom', 
    'dist_X_Pi', 'is_plevin', 'is_hudson', 'remark'
]

# --- Logging ---
logger = logging.getLogger('xpid')

def setup_logging(log_file: Path):
    """Configures logging to file and stdout."""
    if log_file.parent:
        log_file.parent.mkdir(parents=True, exist_ok=True)

    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file, mode='w', encoding='utf-8')
    ]

    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%H:%M:%S',
        handlers=handlers,
        force=True 
    )

# --- Helper Classes (Streaming) ---

class ResultStreamer:
    """
    Handles streaming output to CSV or JSON files to prevent Memory OOM.
    Writes data incrementally as it becomes available.
    """
    def __init__(self, output_path: Path, file_type: str, verbose: bool):
        self.output_path = output_path
        self.file_type = file_type.lower()
        self.verbose = verbose
        self.file_handle = None
        self.csv_writer = None
        self.is_first_chunk = True
        self.has_written_data = False

    def __enter__(self):
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.file_handle = open(self.output_path, 'w', newline='', encoding='utf-8')
        
        if self.file_type == 'json':
            self.file_handle.write('[\n')
        
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_type == 'json':
            self.file_handle.write('\n]')
        
        if self.file_handle:
            self.file_handle.close()

    def _get_fieldnames(self, first_record: Dict[str, Any]) -> List[str]:
        if self.verbose:
            return list(first_record.keys())
        else:
            keys = list(first_record.keys())
            return [k for k in SIMPLE_COLS if k in keys]

    def write_chunk(self, results: List[Dict[str, Any]]):
        if not results:
            return

        self.has_written_data = True

        if self.file_type == 'csv':
            if self.is_first_chunk:
                fieldnames = self._get_fieldnames(results[0])
                self.csv_writer = csv.DictWriter(self.file_handle, fieldnames=fieldnames, extrasaction='ignore')
                self.csv_writer.writeheader()
                self.is_first_chunk = False
            
            self.csv_writer.writerows(results)

        elif self.file_type == 'json':
            keys_to_keep = None
            if not self.verbose and self.is_first_chunk:
                sample_keys = results[0].keys()
                keys_to_keep = set([k for k in SIMPLE_COLS if k in sample_keys])

            for item in results:
                if not self.verbose:
                    if keys_to_keep is None:
                        sample_keys = item.keys()
                        keys_to_keep = set([k for k in SIMPLE_COLS if k in sample_keys])
                    clean_item = {k: v for k, v in item.items() if k in keys_to_keep}
                else:
                    clean_item = item

                if not self.is_first_chunk:
                    self.file_handle.write(',\n')
                else:
                    self.is_first_chunk = False
                
                json.dump(clean_item, self.file_handle, indent=2)

# --- Helpers ---

def find_files(inputs: List[str]) -> List[Path]:
    file_list = set()
    pattern = re.compile(r'^[a-zA-Z0-9]{4}\.(cif|pdb)(\.gz)?$', re.IGNORECASE)
    for inp in inputs:
        path = Path(inp)
        if path.is_file():
            file_list.add(path.resolve())
        elif path.is_dir():
            for p in path.rglob("*"):
                if p.is_file() and pattern.match(p.name):
                    file_list.add(p.resolve())
    return sorted(list(file_list))

def save_single_file_results(results: List[Dict[str, Any]], output_path: Path, file_type: str, verbose: bool) -> None:
    """Helper for saving individual files in separate mode (non-streaming context)."""
    with ResultStreamer(output_path, file_type, verbose) as streamer:
        streamer.write_chunk(results)

def process_one_file(args_packet: Tuple) -> Tuple[Optional[str], int, Optional[List[Dict[str, Any]]], Optional[str]]:
    filepath, mon_lib, ftype, hmode, output_dir, separate_mode, filters, verbose, model_mode = args_packet
    
    if re.match(r'^[a-zA-Z0-9]{4}\.', filepath.name):
        pdb_name = filepath.name[:4]
    else:
        pdb_name = filepath.stem.replace('.cif', '').replace('.pdb', '')
    
    try:
        try:
            structure = gemmi.read_structure(str(filepath))
        except Exception as e:
            return (f"Read Error ({pdb_name}): {e}", 0, None, None)

        structure = prep.add_hydrogens_memory(structure, mon_lib, h_change_val=hmode)
        if not structure:
            return (f"AddH Failed ({pdb_name})", 0, None, None)

        results = core.detect_interactions_in_structure(
            structure, 
            pdb_name, 
            filter_pi=filters['pi'], 
            filter_donor=filters['donor'],
            filter_donor_atom=filters['donor_atom'],
            model_mode=model_mode
        )

        count = len(results)

        if count > 0:
            if separate_mode:
                filename = f"{pdb_name}_xpid.{ftype}"
                out_path = Path(output_dir) / filename
                save_single_file_results(results, out_path, ftype, verbose)
                return (None, count, None, str(out_path.parent))
            else:
                return (None, count, results, None)
        else:
            return (None, 0, None, None)
            
    except Exception as e:
        return (f"Critical Error ({pdb_name}): {e}", 0, None, None)

def get_parser():
    h_mode_help = (
        "Hydrogen handling mode (Default: 4):\n"
        "  0: NoChange       - Keep existing H (Recommended for Neutron)\n"
        "  1: Shift          - Move H to standard distances\n"
        "  2: Remove         - Remove all H\n"
        "  3: ReAdd          - Remove and Re-add all H\n"
        "  4: ReAddButWater  - Re-add all H, skip waters (Standard)\n"
        "  5: ReAddKnown     - Only add H to known residues"
    )

    parser = argparse.ArgumentParser(
        prog="xpid",
        description="xpid: XH-/pi interaction detector.",
        formatter_class=argparse.RawTextHelpFormatter,
        usage="xpid [options] <input>"  # [修改点] 简洁的基础用法提示
    )

    parser.add_argument('inputs', nargs='*', metavar='PATH', help="Input file(s) or directory.")

    io_group = parser.add_argument_group("Input/Output Options")
    io_group.add_argument('--separate', action='store_true', 
                          help="Save separate result files for each PDB (default: merge into one file).")
    io_group.add_argument('--out-dir', type=str, metavar='DIR', 
                          help="Override default output directory.")
    io_group.add_argument('--output-name', type=str, default='xpid_results', metavar='NAME',
                          help="Base name for the merged result file (default: 'xpid_results').")
    io_group.add_argument('--file-type', default='json', choices=['json', 'csv'], 
                          help="Output format. Default: json.")
    io_group.add_argument('-v', '--verbose', action='store_true', help="Output detailed metrics.")
    io_group.add_argument('--log', action='store_true', help="Save log file.")

    set_group = parser.add_argument_group("Configuration")
    set_group.add_argument('--h-mode', type=int, default=4, choices=[0,1,2,3,4,5], metavar='N', help=h_mode_help)
    set_group.add_argument('--jobs', type=int, default=1, metavar='N', help="CPU cores.")
    set_group.add_argument('--model', type=str, default="0", metavar='ID', help="Model index (0,1...) or 'all'. Default: 0.")
    set_group.add_argument('--mon-lib', type=str, metavar='DIR', help="Custom Monomer Library path.")
    set_group.add_argument('--set-mon-lib', type=str, metavar='DIR', help="Set default Monomer Library.")

    filt_group = parser.add_argument_group("Filters")
    filt_group.add_argument('--pi-res', type=str, help="Limit Pi residues (TRP,TYR).")
    filt_group.add_argument('--donor-res', type=str, help="Limit Donor residues (HIS,ARG).")
    filt_group.add_argument('--donor-atom', type=str, help="Limit Donor elements (C,N,O,S).")

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()

    # 1. Config Mode
    if args.set_mon_lib:
        if os.path.isdir(args.set_mon_lib):
            path = config.save_mon_lib_path(args.set_mon_lib)
            print(f"Configuration saved: {path}")
            sys.exit(0)
        else:
            print(f"Error: Invalid path {args.set_mon_lib}")
            sys.exit(1)

    if not args.inputs:
        parser.print_help()
        sys.exit(0)

    # 2. File Finding
    files = find_files(args.inputs)
    if not files:
        print("ERROR: No valid structure files found.")
        sys.exit(1)

    # 3. Output Dir & Log Setup
    if args.out_dir:
        output_dir = Path(args.out_dir)
    else:
        output_dir = files[0].parent / "xpid_output"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = output_dir / "xpid_run.log"
    if args.log:
        setup_logging(log_file)
    else:
        logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%H:%M:%S', force=True)

    logger.info("Scanning files...")

    # 4. Validation
    mon_lib_path = args.mon_lib if args.mon_lib else config.DEFAULT_MON_LIB_PATH
    if mon_lib_path and not os.path.isdir(mon_lib_path):
        logger.error(f"Monomer Library path not found: {mon_lib_path}")
        sys.exit(1)

    # 5. Filters
    filters = {
        'pi': [x.strip().upper() for x in args.pi_res.split(',')] if args.pi_res else None,
        'donor': [x.strip().upper() for x in args.donor_res.split(',')] if args.donor_res else None,
        'donor_atom': None
    }
    if args.donor_atom:
        valid_elements = {'C', 'N', 'O', 'S'}
        inputs = [x.strip().upper() for x in args.donor_atom.split(',')]
        for inp in inputs:
            if inp not in valid_elements:
                logger.error(f"Invalid donor element: '{inp}'. Allowed: C, N, O, S")
                sys.exit(1)
        filters['donor_atom'] = inputs

    # 6. Status Log
    separate_mode = args.separate
    ftype = args.file_type.lower()
    
    h_mode_desc = H_MODE_MAP.get(args.h_mode, "Unknown")
    model_desc = "First Model" if args.model == "0" else ("All Models" if args.model == "all" else "Specific Index")
    output_mode_desc = "Detailed" if args.verbose else "Simple"

    logger.info("--- Xpid Initialization ---")
    logger.info(f"Targets    : {len(files)} files")
    logger.info(f"Output Dir : {output_dir.resolve()}")
    logger.info(f"Format     : {ftype.upper()} ({'Separate Files' if separate_mode else 'Merged File'})")
    logger.info(f"H-Mode     : {args.h_mode} ({h_mode_desc})")
    logger.info(f"Model      : {args.model} ({model_desc})")
    logger.info(f"Filters    : {filters if any(filters.values()) else 'None'}")
    logger.info(f"Columns    : {output_mode_desc}")

    # 7. Execution
    tasks = [(f, mon_lib_path, ftype, args.h_mode, str(output_dir), separate_mode, filters, args.verbose, args.model) for f in files]
    
    error_logs = []
    total_found = 0
    
    try:
        merge_file_path = None
        streamer = None
        
        if not separate_mode:
            merge_filename = f"{args.output_name}.{ftype}"
            merge_file_path = output_dir / merge_filename
            streamer = ResultStreamer(merge_file_path, ftype, args.verbose)
            streamer.__enter__()

        with multiprocessing.Pool(args.jobs, maxtasksperchild=100) as pool:
            for i, (err, count, data, out_path) in enumerate(pool.imap_unordered(process_one_file, tasks), 1):
                if err: 
                    error_logs.append(err)
                    logger.warning(err)
                
                total_found += count
                
                if not separate_mode and data:
                    streamer.write_chunk(data)
                
                msg = f"[INFO] Progress   : {i}/{len(files)} files processed..."
                sys.stdout.write(f"\r{msg}")
                sys.stdout.flush()

        if streamer:
            streamer.__exit__(None, None, None)

        print("")

        # 8. Summary
        print("-" * 60)
        print(f"[SUMMARY] Total XH-pi interactions detected: {total_found}")
        
        if error_logs:
            print(f"[WARNING] {len(error_logs)} files failed processing. See log for details.")

        if total_found > 0:
            if not separate_mode:
                print(f"[OUTPUT] Merged result saved to:\n  -> {merge_file_path.resolve()}")
            else:
                print(f"[OUTPUT] Separate result files saved in:\n  -> {output_dir.resolve()}")
            
            if args.log:
                print(f"[OUTPUT] Log saved to:\n  -> {log_file.resolve()}")
        else:
            print("[OUTPUT] No interactions found.")

    except KeyboardInterrupt:
        print("\n[Aborted by user]")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()