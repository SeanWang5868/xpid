"""Fixed-heavy-atom candidate geometry and conservative local state selection.

The score below is a geometric discriminator, NOT an energy or probability.
It never considers an aromatic centroid or an XH–pi detection result.
"""
from __future__ import annotations

import math
import gemmi
import numpy as np

from . import config, hbond_acceptors, monomer_bonds
from .protonation_pka import alt


def atom_for(residue, name, label):
    exact = [a for a in residue if a.name == name and alt(a) == label]
    shared = [a for a in residue if a.name == name and not alt(a)]
    candidates = exact or shared
    return candidates[0] if len(candidates) == 1 else None


def vector(atom):
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z])


def conformer_occupancy(residue, parent, label):
    # A shared X (occupancy 1) must not give each A/B hydrogen occupancy 1.
    labelled = [a.occ for a in residue if not a.is_hydrogen() and label and alt(a) == label]
    return min([parent.occ] + labelled)


def hydrogen_copy(residue, parent_name, h_name, label, monlib):
    parent = atom_for(residue, parent_name, label)
    hydrogen = atom_for(residue, h_name, label)
    if parent is None or hydrogen is None:
        return None
    if label and not alt(hydrogen):
        heavy_names = monomer_bonds.get_bonded_heavy_atom_names(residue.name, parent_name) or set()
        neighbours = [atom_for(residue, name, label) for name in heavy_names]
        if alt(parent) or any(a is not None and alt(a) for a in neighbours):
            # A shared H direction cannot be assumed valid for two different
            # parent frameworks. Gemmi must supply a conformer-specific H.
            return None
    distance = parent.pos.dist(hydrogen.pos)
    maximum = config.COVALENT_XH_MAX.get(parent.element.name.upper(), 1.2)
    if not math.isfinite(distance) or not config.MIN_COVALENT_XH < distance <= maximum:
        return None
    if residue.name not in monlib.monomers:
        return None
    restraint = next((bond for bond in monlib.monomers[residue.name].rt.bonds
                      if {bond.id1.atom, bond.id2.atom} == {parent_name, h_name}), None)
    if (restraint is None or not math.isfinite(restraint.value_nucleus)
            or not config.MIN_COVALENT_XH < restraint.value_nucleus <= maximum):
        return None
    h = hydrogen.clone()
    # Generated pH-state H represents a nucleus, not an X-ray density maximum.
    h.pos = gemmi.Position(*(vector(parent) + (vector(hydrogen) - vector(parent))
                            * restraint.value_nucleus / distance))
    h.altloc = label or "\0"
    h.occ = conformer_occupancy(residue, parent, label)
    h.flag = "P"
    h.element = gemmi.Element("H")
    return h


def acid_candidates(residue, label, monlib):
    """Use CCP4 ASPprot/GLUprot restraints, not tetrahedral O-H guesses.

    Both carboxyl oxygens and both planar OH directions are considered. The
    dictionary's O-H nuclear distance and C-O-H angle are used explicitly.
    """
    mod_id = residue.name + "prot"
    if mod_id not in monlib.modifications:
        return []
    mod = monlib.modifications[mod_id]
    carbon_name, oxygens, names = ("CG", ("OD1", "OD2"), ("HD1", "HD2")) if residue.name == "ASP" else (
        "CD", ("OE1", "OE2"), ("HE1", "HE2"))
    restraint = next((b for b in mod.rt.bonds if {b.id1.atom, b.id2.atom} == {oxygens[1], names[1]}), None)
    angle = next((a for a in mod.rt.angles if (a.id1.atom, a.id2.atom, a.id3.atom) == (
        carbon_name, oxygens[1], names[1])), None)
    if (restraint is None or angle is None
            or not math.isfinite(restraint.value_nucleus)
            or not config.MIN_COVALENT_XH < restraint.value_nucleus <= config.COVALENT_XH_MAX["O"]
            or not math.isfinite(angle.value) or not 90 <= angle.value <= 130):
        return []
    carbon = atom_for(residue, carbon_name, label)
    oxygen_atoms = [atom_for(residue, name, label) for name in oxygens]
    if carbon is None or any(a is None for a in oxygen_atoms):
        return []
    candidates = []
    for index, oxygen in enumerate(oxygen_atoms):
        axis = vector(carbon) - vector(oxygen)
        norm = np.linalg.norm(axis)
        if not 1.05 < norm < 1.65:
            return []
        axis /= norm
        transverse = vector(oxygen_atoms[1 - index]) - vector(carbon)
        transverse -= np.dot(transverse, axis) * axis
        norm = np.linalg.norm(transverse)
        if norm < 0.1:
            return []
        transverse /= norm
        for sign in (1, -1):
            direction = math.cos(math.radians(angle.value)) * axis + sign * math.sin(
                math.radians(angle.value)) * transverse
            h = gemmi.Atom()
            h.name, h.element, h.flag = names[index], gemmi.Element("H"), "P"
            h.altloc = label or "\0"
            h.occ, h.b_iso = conformer_occupancy(residue, oxygen, label), oxygen.b_iso
            h.pos = gemmi.Position(*(vector(oxygen) + restraint.value_nucleus * direction))
            candidates.append((f"protonated_{oxygens[index]}_{'cis' if sign == 1 else 'trans'}",
                               {oxygens[index]: [h]}))
    return candidates


def choose_candidate(candidates, residue, chain_name, label, model, ns, override=None):
    """Require a nonclashing, distinguishable local H-bond preference.

    Overrides choose a chemical microstate but do not bypass a severe clash.
    Unresolved neighbouring titratable side chains do not vote for a tautomer.
    """
    scores = []
    titratable = {"HIS": {"ND1", "NE2"}, "ASP": {"OD1", "OD2"},
                  "GLU": {"OE1", "OE2"}, "LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"},
                  "CYS": {"SG"}, "TYR": {"OH"}}
    for name, hydrogens in candidates:
        clash, score = False, 0.0
        for parent_name, hs in hydrogens.items():
            parent = atom_for(residue, parent_name, label)
            covalent = monomer_bonds.get_bonded_heavy_atom_names(residue.name, parent_name) or set()
            for h in hs:
                for mark in ns.find_atoms(h.pos, alt=label or "\0", radius=3.0):
                    cra = mark.to_cra(model)
                    a = cra.atom
                    if a.is_hydrogen() or a.occ < config.MIN_ATOM_OCCUPANCY:
                        continue
                    if not hbond_acceptors.altlocs_compatible(label, a.altloc):
                        continue
                    pos = ns.grid_cell.find_nearest_pbc_position(h.pos, a.pos, mark.image_idx)
                    same = (cra.chain.name == chain_name and cra.residue.seqid == residue.seqid
                            and cra.residue.name == residue.name and mark.image_idx == 0)
                    if same and (a.name == parent_name or a.name in covalent):
                        continue
                    distance = h.pos.dist(pos)
                    if distance < 1.45:
                        clash = True
                    if a.name in titratable.get(cra.residue.name, set()):
                        continue
                    if same or not hbond_acceptors.is_hbond_acceptor(cra.residue, a, label):
                        continue
                    if 1.5 <= distance <= 2.6:
                        hd = vector(parent) - vector(h)
                        ha = np.array([pos.x, pos.y, pos.z]) - vector(h)
                        cosine = np.dot(hd, ha) / (np.linalg.norm(hd) * distance)
                        angle = math.degrees(math.acos(np.clip(cosine, -1, 1)))
                        if angle >= 140:
                            score += (2.6 - distance) * (angle - 120) / 60
        scores.append({"state": name, "clash": clash, "score": round(score, 6)})
    valid = sorted((s for s in scores if not s["clash"]), key=lambda s: (-s["score"], s["state"]))
    selected = None
    if override:
        selected = next((s["state"] for s in valid if s["state"] == override), None)
    elif len(valid) == 1 or (len(valid) > 1 and valid[0]["score"] - valid[1]["score"] >= 0.25):
        selected = valid[0]["state"]
    return next((item for item in candidates if item[0] == selected), None), scores
