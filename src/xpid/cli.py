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
from typing import List, Dict, Any, NamedTuple, Optional

try:
    from xpid import prep, core, config, structure_io
    from xpid.output import ResultStreamer
    from xpid.resolver import gather_inputs
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from xpid import prep, core, config, structure_io
    from xpid.output import ResultStreamer
    from xpid.resolver import gather_inputs
from xpid import provenance


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

H_MODE_MAP = {
    0: "NoChange", 1: "Shift", 2: "Remove",
    3: "ReAdd", 4: "ReAddButWater", 5: "ReAddKnown",
}

logger = logging.getLogger("xpid")


BASE_SYSTEM_SUMMARY_KEYS = (
    "hudson",
    "plevin",
    "hudson_plevin_union",
    "hudson_plevin",
)

P_SLAB_SYSTEM_SUMMARY_KEYS = (
    "p_slab",
    "all_three",
)


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
    cone_mode: str  # "auto" | "none"
    no_hbond_gate: bool
    include_p_slab: bool
    report_xh_candidates: bool
    include_coordinates: bool
    include_sasa: bool
    include_cooperativity: bool
    residue_pair: Optional[tuple[str, str]]
    allow_remote_recovery: bool
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


def empty_system_summary(include_p_slab: bool = False) -> Dict[str, int]:
    keys = BASE_SYSTEM_SUMMARY_KEYS
    if include_p_slab:
        keys = keys + P_SLAB_SYSTEM_SUMMARY_KEYS
    return {key: 0 for key in keys}


def summarize_systems(results: List[Dict[str, Any]], include_p_slab: bool = False) -> Dict[str, int]:
    summary = empty_system_summary(include_p_slab)
    for row in results:
        is_hudson = int(row.get("is_hudson", 0))
        is_plevin = int(row.get("is_plevin", 0))

        summary["hudson"] += is_hudson
        summary["plevin"] += is_plevin
        summary["hudson_plevin_union"] += int(bool(is_hudson or is_plevin))
        summary["hudson_plevin"] += int(bool(is_hudson and is_plevin))
        if include_p_slab:
            is_p_slab = int(row.get("is_p_slab", 0))
            summary["p_slab"] += is_p_slab
            summary["all_three"] += int(bool(is_hudson and is_plevin and is_p_slab))
    return summary


def add_system_summary(total: Dict[str, int], chunk: Dict[str, int]):
    for key in total:
        total[key] += chunk.get(key, 0)


# ---------------------------------------------------------------------------
# Per-file worker (runs in multiprocessing pool)
# ---------------------------------------------------------------------------

def process_one_file(task: TaskPacket):
    """Process a single structure file. Returns (error, count, results, path, system_summary)."""
    output_dir = Path(task.output_dir_str)
    pdb_code = task.filepath.stem.split('.')[0].lower()

    try:
        structure = structure_io.read_structure(
            task.filepath, allow_remote_recovery=task.allow_remote_recovery)

        if not structure or len(structure) == 0:
            return f"Empty or invalid structure: {task.filepath}", 0, [], None, empty_system_summary(task.include_p_slab)

        if task.h_mode > 0:
            structure = prep.add_hydrogens_memory(
                structure, h_change_val=task.h_mode)
            if structure is None:
                return f"Hydrogen addition failed: {task.filepath}", 0, [], None, empty_system_summary(task.include_p_slab)

        results = core.detect_interactions_in_structure(
            structure,
            pdb_name=pdb_code,
            filter_pi=task.filters.get('pi'),
            filter_donor=task.filters.get('donor'),
            filter_donor_atom=task.filters.get('donor_atom'),
            model_mode=task.model_mode,
            cone_mode=task.cone_mode,
            no_hbond_gate=task.no_hbond_gate,
            include_p_slab=task.include_p_slab,
            report_xh_candidates=task.report_xh_candidates,
            include_coordinates=task.include_coordinates,
            compute_sasa=task.include_sasa,
            annotate_cooperativity=task.include_cooperativity,
            residue_pair=task.residue_pair,
            min_occ=task.min_occ,
            sym_contacts=task.sym_contacts,
            include_water=task.include_water,
            max_b=task.max_b,
        )

        count = len(results)
        system_summary = summarize_systems(results, include_p_slab=task.include_p_slab)

        if task.separate:
            out_path = output_dir / f"{pdb_code}.{task.ftype_arg}"
            with ResultStreamer(
                out_path, task.ftype_arg, task.verbose,
                include_p_slab=task.include_p_slab,
                include_xh_candidates=task.report_xh_candidates,
                include_coordinates=task.include_coordinates,
                include_sasa=task.include_sasa,
                include_cooperativity=task.include_cooperativity) as streamer:
                streamer.write_chunk(results)
            return None, count, [], str(out_path), system_summary
        else:
            return None, count, results, None, system_summary

    except Exception as e:
        import traceback
        return f"{task.filepath}: {e}\n{traceback.format_exc()}", 0, [], None, empty_system_summary(task.include_p_slab)


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
    out.add_argument('--include-coordinates', action='store_true',
                     help="Include absolute pi-center, X, H coordinates, pi normal, and X-side columns.")
    out.add_argument('--sasa', action='store_true',
                     help='Include solvent accessible surface area columns (pi_sasa_avg, X_sasa, H_sasa).')
    out.add_argument('--cooperativity', action='store_true', default=True,
                     help='Annotate hits with cooperativity metrics (default: on).')
    out.add_argument('--no-cooperativity', action='store_false', dest='cooperativity',
                     help='Disable cooperativity annotation.')
    out.add_argument('--log', action='store_true', help="Save run log to file.")
    out.add_argument('--provenance', action='store_true',
                     help='Write a _metadata.json companion file recording all run parameters.')

    proc = parser.add_argument_group("Processing Options")
    proc.add_argument('--h-mode', type=int, default=4,
                      help="Hydrogen handling mode (0-5). Default: 4.")
    proc.add_argument('--jobs', type=int, default=1,
                      help="Number of CPU cores to use.")
    proc.add_argument('--model', type=str, default="0",
                      help="Model index to analyze (or 'all').")
    proc.add_argument('--no-cone', action="store_true",
                      help="Disable cone rescue. Use explicit hydrogen positions only.")
    proc.add_argument('--no-hbond-gate', action="store_true",
                      help="Disable H-bond competition gate in cone rescue. All sterically-valid conformers are evaluated.")
    proc.add_argument('--include-p-slab', '--p-slab', dest='include_p_slab',
                      action='store_true', default=False,
                      help="Include the optional P-slab system and P-slab output columns.")
    proc.add_argument('--xh-candidates', dest='report_xh_candidates',
                      action='store_true', default=False,
                      help="Export all X-H bonds passing Hudson/Plevin X-position filters, including direction-failed candidates.")
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
    filt.add_argument('--residue-pair', nargs=2, metavar=('SEL1', 'SEL2'),
                      help="Restrict detection to XH-pi interactions between two residue selections, e.g. //A/12 //A/18.")
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


def _allow_remote_recovery(args: argparse.Namespace, files: List[Path]) -> bool:
    """Allow RCSB recovery only for one explicitly provided structure file."""
    if args.pdb_list:
        return False
    if len(files) != 1 or len(args.inputs) != 1:
        return False

    inp = Path(args.inputs[0])
    return inp.is_file()


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
    cone_mode = "none" if args.no_cone else "auto"
    cone_status = "Auto (cone for rotatable)" if cone_mode == "auto" else "Disabled (explicit H only)"
    p_slab_status = "Included" if args.include_p_slab else "Excluded (default)"
    candidate_status = "Enabled" if args.report_xh_candidates else "Disabled (default)"
    coordinate_status = "Included" if args.include_coordinates else "Excluded (default)"
    sym_status = "Enabled" if args.sym_contacts else "Disabled"
    water_status = "Included" if args.include_water else "Excluded (default)"
    max_b_status = f"{args.max_b:.1f} \u00c5\u00b2" if args.max_b > 0 else "No filter"

    logger.info(f"Targets     : {len(files)} unique structures")
    logger.info(f"Output Dir  : {output_dir.resolve()}")
    logger.info(f"Format      : {args.file_type.upper()} "
                f"({'Separate Files' if args.separate else 'Merged File'})")
    logger.info(f"H-Mode      : {args.h_mode} ({h_mode_desc})")
    logger.info(f"Monomer Lib : {config.DEFAULT_MON_LIB_PATH}")
    logger.info(f"Cone Mode   : {cone_mode} ({cone_status})")
    logger.info(f"P-slab      : {p_slab_status}")
    logger.info(f"X-H Export  : {candidate_status}")
    logger.info(f"Coordinates : {coordinate_status}")
    logger.info(f"Sym Contacts: {sym_status}")
    logger.info(f"Water       : {water_status}")
    logger.info(f"Max B-factor: {max_b_status}")
    if args.residue_pair:
        logger.info(f"Residue Pair: {args.residue_pair[0]} <-> {args.residue_pair[1]}")
    allow_remote_recovery = _allow_remote_recovery(args, files)
    recovery_status = "Enabled for single direct file" if allow_remote_recovery else "Disabled for batch input"
    logger.info(f"RCSB Rescue : {recovery_status}")
    if args.provenance:
        meta = provenance.build_metadata(args, output_dir, args.output_name, len(files))
        provenance.write_metadata(meta, output_dir, args.output_name)
        _meta_path = output_dir / f"{args.output_name}_metadata.json"
        logger.info(f"[PROVENANCE] Written to {_meta_path}")
    print("")

    # Step 3: Execute
    ftype_arg = args.file_type.lower()
    tasks = [
        TaskPacket(f, ftype_arg, args.h_mode, str(output_dir),
                   args.separate, filters, args.verbose, args.model,
                   cone_mode, args.include_p_slab, args.report_xh_candidates,
                   args.include_coordinates,
                   args.sasa,
                   args.cooperativity,
                   tuple(args.residue_pair) if args.residue_pair else None,
                   allow_remote_recovery, args.min_occ,
                   args.sym_contacts, args.include_water, args.max_b)
        for f in files
    ]

    error_logs: List[str] = []
    total_found = 0
    system_totals = empty_system_summary(args.include_p_slab)

    try:
        merge_file_path = None
        streamer = None

        if not args.separate:
            merge_file_path = output_dir / f"{args.output_name}.{ftype_arg}"
            streamer = ResultStreamer(
                merge_file_path, ftype_arg, args.verbose,
                include_p_slab=args.include_p_slab,
                include_xh_candidates=args.report_xh_candidates,
                include_coordinates=args.include_coordinates,
                include_sasa=args.sasa,
                include_cooperativity=args.cooperativity)
            streamer.__enter__()

        with multiprocessing.Pool(args.jobs, maxtasksperchild=100) as pool:
            for i, (err, count, data, out_path, system_summary) in enumerate(
                    pool.imap_unordered(process_one_file, tasks), 1):
                if err:
                    error_logs.append(err)
                    logging.warning(f"\n[WARN] {err}")

                total_found += count
                add_system_summary(system_totals, system_summary)

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
        if args.report_xh_candidates:
            print(f"[SUMMARY] Total X-H candidates exported       : {total_found}")
            print(f"[SUMMARY] Hudson-positive among candidates   : {system_totals['hudson']}")
            print(f"[SUMMARY] Plevin-positive among candidates   : {system_totals['plevin']}")
            print(f"[SUMMARY] Hudson/Plevin union among candidates: {system_totals['hudson_plevin_union']}")
            print(f"[SUMMARY] Hudson+Plevin overlap among candidates: {system_totals['hudson_plevin']}")
            print(f"[SUMMARY] Direction-filter negative candidates: {total_found - system_totals['hudson_plevin_union']}")
        else:
            print(f"[SUMMARY] Total XH-π interactions detected: {total_found}")
            print(f"[SUMMARY] Hudson-positive       : {system_totals['hudson']}")
            print(f"[SUMMARY] Plevin-positive       : {system_totals['plevin']}")
            print(f"[SUMMARY] Hudson/Plevin union   : {system_totals['hudson_plevin_union']}")
            print(f"[SUMMARY] Hudson+Plevin overlap : {system_totals['hudson_plevin']}")
        if args.include_p_slab:
            print(f"[SUMMARY] P-slab-positive       : {system_totals['p_slab']}")
            print(f"[SUMMARY] All-three overlap     : {system_totals['all_three']}")

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
