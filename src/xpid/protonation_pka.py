"""pH provenance and optional, isolated PROPKA calculations.

No file is modified and PROPKA's protonated model is never used for detection.
Only its site pKa estimates are imported. Coordinates remain under Gemmi/CCP4
control. Reference values are PROPKA 3.5 model-compound pKas, not predictions.
"""
from __future__ import annotations

import gzip
import io
import logging
import math
import re
import string
import threading
from importlib.metadata import version
from pathlib import Path

import gemmi

REFERENCE_PKA = {"ASP": 3.8, "GLU": 4.5, "HIS": 6.5, "CYS": 9.0,
                 "TYR": 10.0, "LYS": 10.5, "ARG": 12.5}
_BACKEND_LOCK = threading.RLock()


class _WarningCapture(logging.Handler):
    def __init__(self, messages):
        super().__init__(logging.WARNING)
        self.messages = messages

    def emit(self, record):
        self.messages.append(record.getMessage())


def _label_warnings(messages, start, mapping):
    """Keep backend text and append original identities to recognisable labels."""
    for i in range(start, len(messages)):
        labels = [(chain, int(number), name) for name, number, chain in
                  re.findall(r"\b([A-Z0-9]{3})\s+(\d+)\s+([A-Za-z0-9])\b", messages[i])]
        labels += [(chain, int(number), name) for number, name, chain in
                   re.findall(r"\b(\d+)-([A-Z0-9]{3})\s+\(([A-Za-z0-9])\)", messages[i])]
        identities = []
        for chain, number, name in labels:
            key = mapping.get((chain, number))
            if key is not None and key[3] == name:
                identities.append(f"model {key[0]}, chain {key[1]!r}, residue {key[2]} {key[3]}, altloc {key[4]!r}")
        if identities:
            messages[i] += " [XPID input: " + "; ".join(sorted(set(identities))) + "]"


def validate_options(ph, method, h_mode):
    if method not in {"propka", "reference"}:
        raise ValueError("pKa method must be propka or reference")
    if ph is None:
        return
    if h_mode not in {3, 4, 5}:
        raise ValueError("--ph requires a rebuilding H mode (3, 4 or 5)")
    if str(ph).lower() != "auto":
        numeric_ph(ph)
    if method == "propka":
        try:
            import propka.run  # noqa: F401
        except ImportError as exc:
            raise ValueError('pH prediction requires PROPKA: pip install "xpid[ph]"; '
                             'use --pka-method reference only for an explicit approximation') from exc


def numeric_ph(value):
    # Never interpret a range, standard uncertainty, or free text as a point pH.
    text = str(value).strip().strip("'\"")
    if not re.fullmatch(r"[+-]?\d+(?:\.\d+)?", text):
        raise ValueError(f"Expected a single numeric pH, got {value!r}")
    result = float(text)
    if not math.isfinite(result) or not 0 <= result <= 14:
        raise ValueError("Supported pH range is 0–14")
    return result


def resolve_ph(path, ph):
    if str(ph).lower() != "auto":
        return numeric_ph(ph), "user"
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        text = handle.read()
    values = []
    suffix = Path(path.stem).suffix.lower() if path.suffix.lower() == ".gz" else path.suffix.lower()
    if suffix in {".cif", ".mmcif"} or text.lstrip().startswith("data_"):
        block = gemmi.cif.read_string(text).sole_block()
        methods = " ".join(block.find_values("_exptl.method")).upper()
        tags = []
        if "DIFFRACTION" in methods:
            tags.append("_exptl_crystal_grow.pH")
        if "MICROSCOPY" in methods:
            tags.append("_em_buffer.pH")
        if "NMR" in methods:
            tags.append("_pdbx_nmr_exptl_sample_conditions.pH")
        for tag in tags:
            for value in block.find_values(tag):
                if value not in {"?", ".", ""}:
                    values.append((numeric_ph(value), tag))
    else:
        for line in text.splitlines():
            # Only the experimental-condition record, never arbitrary REMARKs.
            if line.startswith("REMARK 200") and "PH" in line.upper():
                match = re.search(r"\bPH\s*:\s*([^;]+)", line, re.I)
                if match and match[1].strip() not in {"NULL", "N/A"}:
                    values.append((numeric_ph(match[1]), "PDB REMARK 200"))
    if not values or len({value for value, _ in values}) != 1:
        raise ValueError("No unique experimental pH; specify --ph explicitly (no default is assumed)")
    return values[0][0], "; ".join(sorted({source for _, source in values}))


def alt(atom):
    return "" if atom.altloc in {"", "\0"} else atom.altloc


def residue_key(model_index, chain, residue, label=""):
    return (model_index, chain, str(residue.seqid).strip(), residue.name, label)


def predict_model(model, model_index, messages=None):
    """Return mapped site predictions for one occupied conformer per residue.

    Long chain names, insertion codes and author numbering are mapped reversibly.
    Alternate conformers not used in the pKa calculation receive no prediction.
    The site layer handles them as unresolved, never by copying another pKa.
    """
    from propka.run import single

    alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
    if len(model) > len(alphabet):
        raise ValueError("PROPKA PDB adapter supports at most 62 chains per model")
    st = gemmi.Structure()
    st.add_model(gemmi.Model("1"))
    mapping = {}
    for ci, chain in enumerate(model):
        if len(chain) > 9999:
            raise ValueError("PROPKA PDB adapter supports at most 9999 residues per chain")
        target = gemmi.Chain(alphabet[ci])
        for ri, residue in enumerate(chain, 1):
            labels = {alt(a) for a in residue if not a.is_hydrogen() and alt(a)}
            def occupancy(label):
                values = [a.occ for a in residue if not a.is_hydrogen() and alt(a) == label]
                return sum(values) / len(values)
            chosen = min(labels, key=lambda label: (-occupancy(label), label)) if labels else ""
            copy = residue.clone()
            for index in range(len(copy) - 1, -1, -1):
                atom = copy[index]
                if (atom.is_hydrogen() or atom.occ < 0.1
                        or alt(atom) not in {"", chosen}):
                    del copy[index]
                else:
                    atom.altloc = "\0"
            copy.seqid = gemmi.SeqId(ri, " ")
            if len({a.name for a in copy}) != len(copy):
                raise ValueError("Ambiguous duplicate atom names in PROPKA input conformer")
            mapping[(alphabet[ci], ri)] = residue_key(model_index, chain.name, residue, chosen)
            target.add_residue(copy)
        st[0].add_chain(target)
    # Keep large PROPKA tables off the progress display, but retain scientific
    # warnings in diagnostics. Restore the caller's logging configuration.
    with _BACKEND_LOCK:
        backend_logger = logging.getLogger("propka")
        old_level, old_propagate = backend_logger.level, backend_logger.propagate
        messages = messages if messages is not None else []
        start = len(messages)
        capture = _WarningCapture(messages)
        backend_logger.addHandler(capture)
        backend_logger.setLevel(logging.WARNING)
        backend_logger.propagate = False
        try:
            result = single("xpid.pdb", stream=io.StringIO(st.make_pdb_string()), write_pka=False)
        finally:
            backend_logger.removeHandler(capture)
            backend_logger.setLevel(old_level)
            backend_logger.propagate = old_propagate
            _label_warnings(messages, start, mapping)
    predictions = {}
    # Inspect the calculated conformer, not averaged groups that can lose coupling flags.
    conformer = result.conformations[result.conformation_names[0]]
    for group in conformer.groups:
        if not group.titratable or group.atom.res_name.strip() not in REFERENCE_PKA:
            continue
        if group.type in {"N+", "C-"}:
            continue
        key = mapping.get((group.atom.chain_id, group.atom.res_num))
        if key is None:
            continue
        coupled = bool(group.covalently_coupled_groups or group.non_covalently_coupled_groups
                       or group.coupled_titrating_group)
        predictions[key] = {"pka": float(group.pka_value), "coupled": coupled}
    return predictions, version("propka")
