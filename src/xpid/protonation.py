"""Opt-in, auditable protein-sidechain protonation on fixed heavy atoms.

Population estimates are NOT atomic occupancies or confidence intervals.
Uncertain sites are not unioned into unconditional XH-pi positives. Charged
His rings are excluded from the neutral-pi analysis, like cationic donors.
"""
from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path

import gemmi

from . import config, covalent_state, hbond_acceptors, hydrogen_prep, monomer_bonds
from . import protonation_geometry as geometry
from . import protonation_pka as pka

logger = logging.getLogger("xpid")
# Dominant-state policy, not an uncertainty estimate for PROPKA.
MIN_POPULATION = 0.9
SITE_H = {
    "HIS": {"ND1": ("HD1",), "NE2": ("HE2",)},
    "ASP": {"OD1": ("HD1",), "OD2": ("HD2",)},
    "GLU": {"OE1": ("HE1",), "OE2": ("HE2",)},
    "CYS": {"SG": ("HG",)}, "TYR": {"OH": ("HH",)},
    "LYS": {"NZ": ("HZ1", "HZ2", "HZ3")},
    "ARG": {"NE": ("HE",), "NH1": ("HH11", "HH12"), "NH2": ("HH21", "HH22")},
}
METALS = {"Zn", "Fe", "Cu", "Mn", "Mg", "Co", "Ni", "Ca", "Cd", "Hg"}
STATE_NAMES = {
    "HIS": {"delta", "epsilon", "doubly_protonated"},
    "CYS": {"protonated", "deprotonated"}, "TYR": {"protonated", "deprotonated"},
    "LYS": {"protonated", "deprotonated"}, "ARG": {"protonated", "deprotonated"},
    "ASP": {f"protonated_OD{i}_{direction}" for i in (1, 2) for direction in ("cis", "trans")} | {"deprotonated"},
    "GLU": {f"protonated_OE{i}_{direction}" for i in (1, 2) for direction in ("cis", "trans")} | {"deprotonated"},
}


@dataclass
class Site:
    key: tuple
    state: str = "unresolved"
    status: str = "unresolved"
    pka: float | None = None
    protonated_fraction: float | None = None
    charge: int | None = None
    source: str = ""
    reason: str = ""
    donor_h: dict = field(default_factory=dict)
    acceptors: dict = field(default_factory=dict)
    candidates: list = field(default_factory=list)
    site_type: str = "sidechain"


@dataclass
class ProtonationContext:
    ph: float
    ph_source: str
    method: str
    sites: list = field(default_factory=list)
    backend_version: str | None = None
    warnings: list = field(default_factory=list)
    _atoms: dict = field(default_factory=dict, repr=False)
    _alternates: dict = field(default_factory=dict, repr=False)

    def add(self, site):
        self.sites.append(site)
        for atom_name in site.donor_h:
            self._atoms[(site.key, atom_name)] = site
            self._alternates.setdefault((site.key[:4], atom_name), []).append(site)

    def lookup(self, model_index, chain, residue, atom, label=None):
        label = pka.alt(atom) if label is None else label
        key = pka.residue_key(model_index, chain, residue, label)
        match = self._atoms.get((key, atom.name)) or self._atoms.get((key[:-1] + ("",), atom.name))
        if match is None and not label:
            alternatives = self._alternates.get((key[:4], atom.name), ())
            if alternatives:
                charges = {s.charge for s in alternatives}
                names = set.intersection(*(set(s.donor_h[atom.name]) for s in alternatives))
                return Site(key, state="shared_altloc_consensus",
                            reason="shared_atom_with_alternate_states",
                            charge=next(iter(charges)) if len(charges) == 1 else None,
                            donor_h={atom.name: sorted(names)},
                            acceptors={atom.name: all(s.acceptors.get(atom.name, False) for s in alternatives)})
        return match

    def hydrogen_names(self, model_index, chain, residue, atom, label=None):
        site = self.lookup(model_index, chain, residue, atom, label)
        if site is not None and residue.name == "HIS" and site.charge == 1:
            return set()  # no imidazolium N-H/cation-pi donors
        return None if site is None else set(site.donor_h.get(atom.name, ()))

    def shared_alternates(self, model_index, chain, residue, atom):
        key = pka.residue_key(model_index, chain, residue)
        return not pka.alt(atom) and any(
            s.key[-1] for s in self._alternates.get((key[:4], atom.name), ()))

    def is_acceptor(self, model_index, chain, residue, atom, label=None):
        site = self.lookup(model_index, chain, residue, atom, label)
        if site is not None:
            return bool(site.acceptors.get(atom.name, False))
        return hbond_acceptors.is_hbond_acceptor(residue, atom, label)

    def pi_allowed(self, model_index, chain, residue, label):
        if residue.name != "HIS":
            return True
        atom = geometry.atom_for(residue, "ND1", label)
        if atom is None:
            return False
        site = self.lookup(model_index, chain, residue, atom, label)
        return site is not None and site.charge == 0

    def hit_fields(self, model_index, chain, residue, atom, pi_chain, pi_residue, pi_alt, h_alt):
        site = self.lookup(model_index, chain, residue, atom, h_alt or pka.alt(atom))
        pi_site = self._atoms.get((pka.residue_key(model_index, pi_chain, pi_residue, pi_alt), "ND1"))
        return {"pH": self.ph, "protonation_method": self.method,
                "X_protonation_state": site.state if site else "not_titrated",
                "X_protonation_source": site.source if site else "monomer_library",
                "X_site_pKa": site.pka if site else None,
                "X_protonated_fraction": site.protonated_fraction if site else None,
                "pi_protonation_state": pi_site.state if pi_site else "not_titrated"}

    def to_dict(self):
        return {"pH": self.ph, "pH_source": self.ph_source, "method": self.method,
                "backend_version": self.backend_version, "minimum_dominant_population": MIN_POPULATION,
                "population_is_occupancy": False, "warnings": self.warnings,
                "sites": [asdict(site) for site in self.sites],
                "resolved_sites": sum(s.status == "resolved" for s in self.sites),
                "unresolved_sites": sum(s.status != "resolved" for s in self.sites)}


def read_overrides(path):
    """Strict JSON records, with explicit model, residue, altloc and state."""
    if path is None:
        return {}
    with open(Path(path), encoding="utf-8") as handle:
        records = json.load(handle)
    if not isinstance(records, list):
        raise ValueError("Protonation overrides must be a JSON list")
    result = {}
    fields = {"model_index", "chain", "residue_id", "residue_name", "altloc", "state"}
    for record in records:
        if not isinstance(record, dict) or set(record) != fields:
            raise ValueError(f"Each protonation override must contain exactly {sorted(fields)}")
        key = (record["model_index"], record["chain"], str(record["residue_id"]),
               record["residue_name"], record["altloc"])
        if (not isinstance(key[0], int) or isinstance(key[0], bool) or key[0] < 0
                or any(not isinstance(value, str) for value in key[1:])
                or len(key[4]) > 1 or key in result):
            raise ValueError(f"Invalid or duplicate override identity: {key}")
        if record["state"] not in STATE_NAMES.get(key[3], set()):
            raise ValueError(f"Unsupported state {record['state']!r} for {key[3]}")
        result[key] = record["state"]
    return result


def _special_state(chain_name, residue, label, model, ns, connections):
    for name in SITE_H[residue.name]:
        atom = geometry.atom_for(residue, name, label)
        if atom is None or atom.occ < config.MIN_ATOM_OCCUPANCY:
            return "missing_or_unoccupied_site_atom"
        bonded = monomer_bonds.get_bonded_heavy_atom_names(residue.name, name)
        if not bonded or any(geometry.atom_for(residue, parent, label) is None for parent in bonded):
            return "missing_dictionary_or_covalent_parent"
        key = (chain_name, str(residue.seqid).strip(), residue.name, name)
        if any(hbond_acceptors.altlocs_compatible(label, source.altloc)
               and hbond_acceptors.altlocs_compatible(label, target.altloc)
               for source, target in connections.get(key, ())):
            return "covalent_link_or_coordination"
        for mark in ns.find_atoms(atom.pos, alt=label or "\0", radius=2.8):
            cra = mark.to_cra(model)
            other = cra.atom
            if other.is_hydrogen() or other.occ < config.MIN_ATOM_OCCUPANCY:
                continue
            if not hbond_acceptors.altlocs_compatible(label, other.altloc):
                continue
            pos = ns.grid_cell.find_nearest_pbc_position(atom.pos, other.pos, mark.image_idx)
            distance = atom.pos.dist(pos)
            if distance < 0.01:
                continue
            if other.element.name in METALS and distance < 2.8:
                return "metal_coordination_requires_specialist_assignment"
            same = (cra.chain.name == chain_name and cra.residue.seqid == residue.seqid
                    and cra.residue.name == residue.name and mark.image_idx == 0)
            if not same and distance < atom.element.covalent_r + other.element.covalent_r + 0.25:
                return "possible_covalent_link"
    return ""


def _assign(site, residue, label, ph, override, monlib, chain_name, model, ns):
    name = residue.name
    candidates = []
    fraction = site.protonated_fraction
    if override:
        protonated = override not in {"delta", "epsilon", "deprotonated"}
    elif fraction is None:
        site.reason = "no_site_pKa"
        return {}
    elif 1 - MIN_POPULATION < fraction < MIN_POPULATION:
        site.reason = "mixed_charge_states_near_pKa"
        return {}
    else:
        protonated = fraction >= MIN_POPULATION
    site.charge = (int(protonated) if name in {"HIS", "LYS", "ARG"} else int(protonated) - 1)
    if not protonated and name in {"ASP", "GLU", "CYS", "TYR"}:
        site.state, site.status = "deprotonated", "resolved"
        site.acceptors = {atom: True for atom in SITE_H[name]}
        return {}
    if name == "HIS":
        hd = geometry.hydrogen_copy(residue, "ND1", "HD1", label, monlib)
        he = geometry.hydrogen_copy(residue, "NE2", "HE2", label, monlib)
        if hd is None or he is None:
            site.reason = "missing_dictionary_hydrogen_geometry"
            return {}
        candidates = [("doubly_protonated", {"ND1": [hd], "NE2": [he]})] if protonated else [
            ("delta", {"ND1": [hd]}), ("epsilon", {"NE2": [he]})]
    elif name in {"ASP", "GLU"}:
        candidates = geometry.acid_candidates(residue, label, monlib)
    elif name in {"CYS", "TYR", "LYS", "ARG"}:
        if not protonated:
            # Neutral lysine/arginine donor geometries are not part of the
            # production detector. Do not arbitrarily remove one of N's Hs.
            site.state = "deprotonated"
            site.reason = "neutral_amine_or_guanidine_geometry_not_supported"
            if name == "LYS":
                site.acceptors = {"NZ": True}
            return {}
        hs = {parent: [geometry.hydrogen_copy(residue, parent, h, label, monlib) for h in names]
              for parent, names in SITE_H[name].items()}
        if any(h is None for values in hs.values() for h in values):
            site.reason = "missing_dictionary_hydrogen_geometry"
            return {}
        # Rotatable groups retain existing riding geometry; their alternative
        # orientations are evaluated by the existing complete Cone detector.
        site.state, site.status = "protonated", "resolved"
        site.donor_h.update({parent: [h.name for h in values] for parent, values in hs.items()})
        site.acceptors = {atom: name == "TYR" for atom in SITE_H[name]}
        return hs
    selected, site.candidates = geometry.choose_candidate(
        candidates, residue, chain_name, label, model, ns, override)
    if selected is None:
        site.reason = "ambiguous_or_clashing_microstate" if candidates else "missing_dictionary_geometry"
        return {}
    site.state, hs = selected
    site.status = "resolved"
    site.donor_h.update({parent: [h.name for h in values] for parent, values in hs.items()})
    site.acceptors = {atom: atom not in hs for atom in SITE_H[name]}
    return hs


def prepare(structure, path, ph, method="propka", overrides=None, excluded_residues=None,
            model_mode="all"):
    """Apply after Gemmi H rebuilding, before any neighbour search/detection.

    Original atom names, heavy coordinates, occupancies and formal charges are
    preserved. Chemical state lives in this explicit context, not global caches.
    """
    value, source = pka.resolve_ph(path, ph)
    context = ProtonationContext(value, source, method)
    overrides = read_overrides(overrides) if not isinstance(overrides, dict) else overrides
    pending = set(overrides)
    excluded_residues = excluded_residues or set()
    monlib = hydrogen_prep._get_shared_monlib(set(SITE_H))
    for mi, model in enumerate(structure):
        if model_mode != "all" and mi != int(model_mode):
            continue
        predictions = {}
        if method == "propka":
            messages = []
            try:
                predictions, context.backend_version = pka.predict_model(model, mi, messages)
            except Exception as exc:
                message = f"Model {mi}: PROPKA failed: {exc}; no reference-pKa fallback"
                context.warnings.append(message)
                logger.warning("%s", message)
            if messages:
                context.warnings.extend(f"PROPKA model {mi}: {message}" for message in messages)
                logger.warning("PROPKA model %d emitted %d warning(s); see protonation diagnostics", mi, len(messages))
        ns = gemmi.NeighborSearch(model, structure.cell, 4.0).populate(include_h=False)
        connections = covalent_state.connection_index(structure, model)
        changes = []
        for chain in model:
            for residue in chain:
                if residue.name not in SITE_H:
                    continue
                labels = sorted({pka.alt(a) for a in residue if not a.is_hydrogen() and pka.alt(a)}) or [""]
                replacement = []
                for label in labels:
                    key = pka.residue_key(mi, chain.name, residue, label)
                    site = Site(key, source=method, donor_h={name: [] for name in SITE_H[residue.name]})
                    override = overrides.get(key)
                    pending.discard(key)
                    coupled = False
                    if method == "reference":
                        site.pka = pka.REFERENCE_PKA[residue.name]
                    else:
                        prediction = predictions.get(key)
                        if prediction and math.isfinite(prediction["pka"]):
                            site.pka = prediction["pka"]
                            coupled = prediction["coupled"]
                            if coupled and not override:
                                site.reason = "coupled_sites_require_microstate_assignment"
                    if site.pka is not None and not coupled:
                        site.protonated_fraction = 1 / (1 + 10 ** max(-100, min(100, value - site.pka)))
                    if override:
                        site.source = "user_override"
                    if key[:4] in excluded_residues:
                        site.reason = "hydrogen_preparation_quarantined"
                    else:
                        site.reason = _special_state(chain.name, residue, label, model, ns, connections) or site.reason
                    if not site.reason:
                        hs = _assign(site, residue, label, value, override, monlib, chain.name, model, ns)
                        replacement.extend(h for values in hs.values() for h in values)
                    if override and site.status != "resolved":
                        context.warnings.append(f"Override {key} not applied: {site.reason}")
                    context.add(site)
                changes.append((residue, replacement))
        # Do not mutate atom arrays while a neighbour-search index is in use.
        for residue, replacement in changes:
            names = {h for values in SITE_H[residue.name].values() for h in values}
            for i in range(len(residue) - 1, -1, -1):
                atom = residue[i]
                if atom.is_hydrogen() and monomer_bonds.canonical_hydrogen_name(atom.name) in names:
                    del residue[i]
            for hydrogen in replacement:
                residue.add_atom(hydrogen)
        # Terminal Hs can be generated at a chain break as well as a genuine
        # terminus. Do not let an unmodelled terminal amine masquerade as a
        # neutral backbone NH. Terminal pKas/capping are outside this release.
        terminal_ns = gemmi.NeighborSearch(model, structure.cell, 2.0).populate(include_h=False)
        for chain in model:
            for residue in chain:
                if not gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                    continue
                labels = sorted({pka.alt(a) for a in residue if not a.is_hydrogen() and pka.alt(a)}) or [""]
                for label in labels:
                    names = {monomer_bonds.canonical_hydrogen_name(a.name) for a in residue
                             if a.is_hydrogen() and hbond_acceptors.altlocs_compatible(label, a.altloc)}
                    nitrogen = geometry.atom_for(residue, "N", label)
                    acyl_link = False
                    if nitrogen is not None:
                        for mark in terminal_ns.find_atoms(nitrogen.pos, alt=label or "\0", radius=1.8):
                            cra = mark.to_cra(model)
                            if cra.atom.element.name != "C" or cra.atom.occ < config.MIN_ATOM_OCCUPANCY:
                                continue
                            if cra.chain.name == chain.name and cra.residue.seqid == residue.seqid and mark.image_idx == 0:
                                continue
                            if not hbond_acceptors.altlocs_compatible(label, cra.atom.altloc):
                                continue
                            pos = terminal_ns.grid_cell.find_nearest_pbc_position(nitrogen.pos, cra.atom.pos, mark.image_idx)
                            if 1.1 < nitrogen.pos.dist(pos) < 1.8:
                                acyl_link = True
                                break
                    if nitrogen is not None and (names & {"H2", "H3"} or not acyl_link):
                        context.add(Site(pka.residue_key(mi, chain.name, residue, label),
                                         source="unmodelled", reason="terminal_amine_or_chain_break",
                                         donor_h={"N": []}, site_type="N_terminus"))
                    if geometry.atom_for(residue, "OXT", label) is not None:
                        context.add(Site(pka.residue_key(mi, chain.name, residue, label),
                                         source="unmodelled", reason="terminal_carboxyl_protonation_not_modelled",
                                         donor_h={"O": [], "OXT": []}, site_type="C_terminus"))
    if pending:
        raise ValueError(f"Protonation override(s) match no site: {sorted(pending)}")
    context.warnings.append("Terminal amines/carboxyls are unresolved and excluded at their affected atoms; "
                            "ligands/nonstandard residues retain dictionary preparation without pH assignment. "
                            "LYS/ARG N donors remain excluded.")
    unresolved = sum(s.status != "resolved" for s in context.sites)
    logger.info("pH %.2f (%s): %d resolved / %d unresolved sites; see diagnostics",
                value, method, len(context.sites) - unresolved, unresolved)
    return context
