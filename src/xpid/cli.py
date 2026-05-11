"""
cli.py
Command-line interface for xpid — the XH-π interaction detector.
"""
import argparse
import logging
import multiprocessing
import os
import sys
from pathlib import Path
from typing import List, Dict, Any, NamedTuple

import gemmi

try:
    from xpid import prep, core, config
    from xpid.output import ResultStreamer
    from xpid.resolver import gather_inputs
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from xpid import prep, core, config
    from xpid.output import ResultStreamer
    from xpid.resolver import gather_inputs


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

H_MODE_MAP = {
    0: "NoChange", 1: "Shift", 2: "Remove",
    3: "ReAdd", 4: "ReAddButWater", 5: "ReAddKnown",
}

logger = logging.getLogger("xpid")


class TaskPacket(NamedTuple):
    """All parameters needed to process one structure file."""
    filepath: Path
    ftype_arg: str
    h_mode: int
    output_dir_str: str
    separate: bool
    filters: dict
    verbose: bool
    model_mode: str
    use_cone: bool
    min_occ: float
    sym_contacts: bool
    include_water: bool
    max_b: float


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

def setup_logging(log_file: Path):
    """Configure logging to both stdout and a file."""
    if log_file.parent:
        log_file.parent.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file, mode='w', encoding='utf-8'),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%H:%M:%S',
        handlers=handlers,
        force=True,
    )


# ---------------------------------------------------------------------------
# Per-file worker (runs in multiprocessing pool)
# ---------------------------------------------------------------------------

def process_one_file(task: TaskPacket):
    """Process a single structure file. Returns (error, count, results, path)."""
    output_dir = Path(task.output_dir_str)
    pdb_code = task.filepath.stem.split('.')[0].lower()

    try:
        structure = gemmi.read_structure(str(task.filepath))

        if not structure or len(structure) == 0:
            return f"Empty or invalid structure: {task.filepath}", 0, [], None

        if task.h_mode > 0:
            structure = prep.add_hydrogens_memory(
                structure, h_change_val=task.h_mode)
            if structure is None:
                return f"Hydrogen addition failed: {task.filepath}", 0, [], None

        results = core.detect_interactions_in_structure(
            structure,
            pdb_name=pdb_code,
            filter_pi=task.filters.get('pi'),
            filter_donor=task.filters.get('donor'),
            filter_donor_atom=task.filters.get('donor_atom'),
            model_mode=task.model_mode,
            use_cone=task.use_cone,
            min_occ=task.min_occ,
            sym_contacts=task.sym_contacts,
            include_water=task.include_water,
            max_b=task.max_b,
        )

        count = len(results)

        if task.separate:
            out_path = output_dir / f"{pdb_code}.{task.ftype_arg}"
            with ResultStreamer(out_path, task.ftype_arg, task.verbose) as streamer:
                streamer.write_chunk(results)
            return None, count, [], str(out_path)
        else:
            return None, count, results, None

    except Exception as e:
        import traceback
        return f"{task.filepath}: {e}\n{traceback.format_exc()}", 0, [], None


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="XH-π interaction detector")
    parser.add_argument('inputs', nargs='*', help="PDB/CIF files or directories")
    parser.add_argument('--pdb-list', type=str, help="Text file with PDB codes")
    parser.add_argument('--pdb-mirror', type=str, help="Local PDB mirror root")
    parser.add_argument('--redo-mirror', type=str,
                        help="Local PDB-REDO mirror root (prioritized over standard PDB)")

    out = parser.add_argument_group("Output Options")
    out.add_argument('--separate', action='store_true',
                     help="Write separate output files for each PDB.")
    out.add_argument('--out-dir', type=str, help="Directory for output files.")
    out.add_argument('--output-name', type=str, default='xpid_results',
                     help="Filename for merged output.")
    out.add_argument('--file-type', default='json',
                     choices=['json', 'csv', 'parquet'], help="Output format.")
    out.add_argument('-v', '--verbose', action='store_true',
                     help="Include detailed geometric columns.")
    out.add_argument('--log', action='store_true', help="Save run log to file.")

    proc = parser.add_argument_group("Processing Options")
    proc.add_argument('--h-mode', type=int, default=4,
                      help="Hydrogen handling mode (0-5). Default: 4.")
    proc.add_argument('--jobs', type=int, default=1,
                      help="Number of CPU cores to use.")
    proc.add_argument('--model', type=str, default="0",
                      help="Model index to analyze (or 'all').")
    proc.add_argument('--cone', dest='use_cone', action='store_true', default=True,
                      help="Enable implicit Cone logic for rotatable groups (default).")
    proc.add_argument('--no-cone', dest='use_cone', action='store_false',
                      help="Disable implicit Cone logic for rotatable groups.")
    proc.add_argument('--sym-contacts', action='store_true',
                      help="Detect XH-π interactions across crystallographic symmetry mates.")
    proc.add_argument('--include-water', action='store_true',
                      help="Include water molecules as potential XH-π donors.")
    proc.add_argument('--max-b', type=float, default=0.0,
                      help="Maximum B-factor for any π-ring atom or donor X atom (0 = no filter).")

    filt = parser.add_argument_group("Filters & Config")
    filt.add_argument('--pi-res', type=str,
                      help="Filter: π residues (e.g. TRP,TYR).")
    filt.add_argument('--donor-res', type=str,
                      help="Filter: Donor residues (e.g. LYS,ARG).")
    filt.add_argument('--donor-atom', type=str,
                      help="Filter: Donor element symbols or exact atom names (e.g. N,O,C or OG,NZ).")
    filt.add_argument('--min-occ', type=float, default=0.0,
                      help="Minimum combined occupancy to report (default: 0.0).")

    return parser


def _default_output_dir(inputs: List[str], files: List[Path]) -> Path:
    """Default to the scanned folder, or to a single structure's parent folder."""
    if len(inputs) == 1:
        inp = Path(inputs[0])
        if inp.is_dir():
            return inp.resolve()
        if inp.is_file():
            return inp.resolve().parent

    if len(files) == 1:
        return files[0].resolve().parent

    return Path.cwd()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = _build_parser()
    args = parser.parse_args()

    # Step 1: Gather input files
    files = gather_inputs(args.inputs, args.pdb_list, args.pdb_mirror, args.redo_mirror)

    if not files:
        print("No valid input files found. Please check inputs or list/mirror paths.")
        sys.exit(1)

    # Output directory
    output_dir = Path(args.out_dir).resolve() if args.out_dir else _default_output_dir(args.inputs, files)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Logging
    log_file = output_dir / "xpid_run.log"
    if args.log:
        setup_logging(log_file)
    else:
        logging.basicConfig(
            level=logging.INFO, format='%(message)s',
            handlers=[logging.StreamHandler(sys.stdout)])

    logger.info("--- Xpid Initialization ---")

    # Step 2: Build configuration
    filters = {
        'pi': [x.strip().upper() for x in args.pi_res.split(',')] if args.pi_res else None,
        'donor': [x.strip().upper() for x in args.donor_res.split(',')] if args.donor_res else None,
        'donor_atom': [x.strip().upper() for x in args.donor_atom.split(',')] if args.donor_atom else None,
    }

    h_mode_desc = H_MODE_MAP.get(args.h_mode, "Unknown")
    cone_status = "Enabled (default)" if args.use_cone else "Disabled"
    sym_status = "Enabled" if args.sym_contacts else "Disabled"
    water_status = "Included" if args.include_water else "Excluded (default)"
    max_b_status = f"{args.max_b:.1f} \u00c5\u00b2" if args.max_b > 0 else "No filter"

    logger.info(f"Targets     : {len(files)} unique structures")
    logger.info(f"Output Dir  : {output_dir.resolve()}")
    logger.info(f"Format      : {args.file_type.upper()} "
                f"({'Separate Files' if args.separate else 'Merged File'})")
    logger.info(f"H-Mode      : {args.h_mode} ({h_mode_desc})")
    logger.info(f"Monomer Lib : {config.DEFAULT_MON_LIB_PATH}")
    logger.info(f"Cone Logic  : {cone_status}")
    logger.info(f"Sym Contacts: {sym_status}")
    logger.info(f"Water       : {water_status}")
    logger.info(f"Max B-factor: {max_b_status}")
    print("")

    # Step 3: Execute
    ftype_arg = args.file_type.lower()
    tasks = [
        TaskPacket(f, ftype_arg, args.h_mode, str(output_dir),
                   args.separate, filters, args.verbose, args.model, args.use_cone,
                   args.min_occ, args.sym_contacts, args.include_water, args.max_b)
        for f in files
    ]

    error_logs: List[str] = []
    total_found = 0

    try:
        merge_file_path = None
        streamer = None

        if not args.separate:
            merge_file_path = output_dir / f"{args.output_name}.{ftype_arg}"
            streamer = ResultStreamer(merge_file_path, ftype_arg, args.verbose)
            streamer.__enter__()

        with multiprocessing.Pool(args.jobs, maxtasksperchild=100) as pool:
            for i, (err, count, data, out_path) in enumerate(
                    pool.imap_unordered(process_one_file, tasks), 1):
                if err:
                    error_logs.append(err)
                    logging.warning(f"\n[WARN] {err}")

                total_found += count

                if not args.separate and data:
                    streamer.write_chunk(data)

                percent = (i / len(files)) * 100
                sys.stdout.write(
                    f"\r[INFO] Progress   : {i}/{len(files)} ({percent:.1f}%)")
                sys.stdout.flush()

        if streamer:
            streamer.__exit__(None, None, None)

        print("\n")

        # Step 4: Summary
        print("-" * 60)
        print(f"[SUMMARY] Total XH-π interactions detected: {total_found}")

        if error_logs:
            print(f"[WARNING] {len(error_logs)} files failed processing. Check log for details.")

        if total_found > 0:
            if not args.separate:
                print(f"[OUTPUT] Merged result saved to:\n  -> {merge_file_path.resolve()}")
            else:
                print(f"[OUTPUT] Separate result files saved in:\n  -> {output_dir.resolve()}")
        else:
            print("[OUTPUT] No interactions found.")

    except Exception as e:
        logger.error(f"\nExecution failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
