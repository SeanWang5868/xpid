"""Read-only Cone/H-bond compatibility audit; never an occupancy or veto model.

Production steric masks and positive directions are observed, not recomputed.
Stricter acceptor typing is used ONLY here, never to change detection.
"""
from __future__ import annotations

import numpy as np
import gemmi

from . import config, cone, hbond_acceptors as acceptors, monomer_bonds

MAX_PAIRS = 10000
METALS = {"Zn", "Fe", "Cu", "Mn", "Mg", "Co", "Ni", "Ca", "Cd", "Hg", "Na", "K"}


def new_audit():
    return {"schema_version": 1, "mode": "annotation_only", "affects_detection": False,
            "scope": "Cone-eligible donor/ring conformers after spatial and chemical prefilters",
            "altloc_policy": "production-compatible labels; cross-residue correlations are not established",
            "geometry": {"HA_min": cone.ABSOLUTE_MIN_HA,
                         "HA_max": cone.HBOND_CONTACT_HA_MAX,
                         "XHA_min": cone.HBOND_CONTACT_DHA_MIN},
            "mask_encoding": "numpy.packbits bitorder=big; conformer-major, H-minor; padded with zero bits",
            "population_or_energy": False, "max_pairs": MAX_PAIRS,
            "pairs_seen": 0, "pairs_omitted": 0, "pairs": []}


def classify_acceptor(residue, atom, label="", site=None, special=False):
    """Return (status, reason, source). 'yes' means typed, not experimentally proven."""
    element, name = atom.element.name, atom.name
    if element not in {"N", "O", "S", "Se"}:
        return "no", "not_a_supported_acceptor_element", "element"
    standard = gemmi.find_tabulated_residue(residue.name).is_amino_acid()
    # Ordinary peptide links are inter-residue covalent connections, not a
    # reason to promote an amide N into the uncertain-acceptor bucket. These
    # nitrogens lack an available lone pair in the supported standard chemistry.
    # Leave titratable amines (e.g. Lys NZ) to the assigned microstate below.
    if element == "N":
        if (standard and name == "N") or (residue.name, name) in {
                ("ASN", "ND2"), ("GLN", "NE2")}:
            return "no", "non_acceptor_nitrogen", "standard_residue"
        if (residue.name, name) == ("TRP", "NE1"):
            return "no", "pyrrole_nitrogen", "standard_residue"
    if special:
        return "unknown", "covalent_modification_or_metal_environment", "structure"
    if site is not None:
        if site.status != "resolved":
            return "unknown", site.reason or "unresolved_protonation", site.source
        return ("yes" if site.acceptors.get(name, False) else "no"), "assigned_microstate", site.source
    if atom.charge > 0:
        return "unknown", "positive_formal_charge", "structure"
    if standard and name == "O":
        return "yes", "backbone_carbonyl", "standard_residue"
    if standard and name == "OXT":
        return "unknown", "terminal_protonation_not_established", "structure"
    if (residue.name, name) in acceptors.NON_ACCEPTOR_NITROGEN or (standard and name == "N"):
        return "no", "non_acceptor_nitrogen", "standard_residue"
    if (residue.name, name) == ("TRP", "NE1"):
        return "no", "pyrrole_nitrogen", "standard_residue"
    if residue.is_water() and element == "O":
        return "yes", "water_oxygen_orientation_not_known", "standard_component"
    if (residue.name, name) in acceptors.HYDROXYL_ACCEPTORS:
        return "yes", "hydroxyl_oxygen", "standard_residue"
    if residue.name in {"HIS", "ASP", "GLU", "CYS"}:
        # A missing or zero-occupancy riding H is NOT proof of deprotonation.
        if residue.name == "HIS" and name in {"ND1", "NE2"}:
            if acceptors._has_bonded_hydrogen(residue, atom, label):
                return "no", "occupied_ring_nitrogen", "explicit_hydrogen"
        return "unknown", "protonation_not_established", "structure"
    if name in acceptors.SIDECHAIN_ACCEPTORS.get(residue.name, set()):
        if acceptors._has_bonded_hydrogen(residue, atom, label):
            return "unknown", "unexpected_bonded_hydrogen", "structure"
        return "yes", "standard_acceptor", "standard_residue"
    return "unknown", "unsupported_component_or_atom_chemistry", "structure"


def pack(mask):
    return np.packbits(np.asarray(mask, dtype=bool).ravel()).tobytes().hex()


def relations(pi, hb):
    """Independent relations over all allowed conformers, not one winning H."""
    pi, hb = np.asarray(pi, dtype=bool), np.asarray(hb, dtype=bool)
    pc, hc = np.any(pi, axis=1), np.any(hb, axis=1)
    other = np.any(pi & ((hb.sum(axis=1)[:, None] - hb.astype(int)) > 0), axis=1)
    return {"hbond_conformers": int(hc.sum()),
            "same_hydrogen_conformers": int(np.any(pi & hb, axis=1).sum()),
            "same_conformer_other_hydrogen": int(other.sum()),
            "hbond_without_pi_conformers": int((hc & ~pc).sum()),
            "pi_without_hbond_conformers": int((pc & ~hc).sum()),
            "disjoint_only": bool(pc.any() and hc.any() and not (pc & hc).any())}


def contact_mask(x_pos, positions, acceptor_position, valid):
    ha = np.asarray(acceptor_position) - positions
    hx = np.asarray(x_pos) - positions
    distances = np.linalg.norm(ha, axis=-1)
    denom = distances * np.linalg.norm(hx, axis=-1)
    cosine = np.divide(np.sum(ha * hx, axis=-1), denom,
                       out=np.ones_like(distances), where=denom > 0)
    angles = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    mask = (valid[:, None] & (distances >= cone.ABSOLUTE_MIN_HA)
            & (distances <= cone.HBOND_CONTACT_HA_MAX)
            & (angles >= cone.HBOND_CONTACT_DHA_MIN))
    return mask, distances, angles


def summarize(raw, x_pos, contacts, identity, donor_element):
    conformers, valid, pi = raw["conformers"], raw["valid_mask"], raw["pi_mask"]
    known, unknown = np.zeros_like(pi), np.zeros_like(pi)
    record = {**identity, "donor_class": "weak_CH_contact_only" if donor_element == "C" else "polar_XH",
              "rotation_step_degrees": 1, "mask_shape": list(pi.shape),
              "total_conformers": len(conformers), "sterically_allowed": int(valid.sum()),
              "pi_positive_conformers": int(np.any(pi, axis=1).sum()),
              "allowed_mask": pack(valid), "pi_mask": pack(pi), "acceptors": []}
    if conformers:
        positions = np.asarray([c.hydrogen_positions for c in conformers])
        for candidate in contacts:
            if candidate["status"] == "no":
                continue
            mask, distances, angles = contact_mask(x_pos, positions, candidate["position"], valid)
            if not mask.any():
                continue
            if candidate["status"] == "yes":
                known |= mask
            else:
                unknown |= mask
            ci, hi = np.argwhere(mask)[0]
            record["acceptors"].append({**candidate, **relations(pi, mask), "contact_mask": pack(mask),
                                       "example": {"phi": conformers[ci].phi, "hydrogen_index": int(hi),
                                                   "HA": float(distances[ci, hi]), "XHA": float(angles[ci, hi])}})
    record.update(typed_contacts=relations(pi, known), uncertain_contacts=relations(pi, unknown))
    return record


def gather_contacts(rctx, donor, x_pos, parent_pos, definition, combined_occ):
    """Separate audit search: does not alter the production 4 A environment."""
    reference = gemmi.Position(*x_pos)
    radius = definition.xh_bond_length + cone.HBOND_CONTACT_HA_MAX + 1e-6
    label = donor.active_altloc or rctx.pi_alt
    result, seen = [], set()
    one_three = monomer_bonds.get_bonded_heavy_atom_names(
        donor.residue_name, donor.parent_atom.name) or set()
    for mark in rctx.ns.find_atoms(reference, radius=radius):
        cra = mark.to_cra(rctx.model)
        atom, residue = cra.atom, cra.residue
        if atom.element.name not in {"O", "N", "S", "Se"}:
            continue
        if (atom.occ < config.MIN_ATOM_OCCUPANCY
                or min(combined_occ, donor.occupancy, atom.occ) < rctx.min_occ
                or (rctx.max_b > 0 and atom.b_iso > rctx.max_b)):
            continue
        if any(not acceptors.altlocs_compatible(a, atom.altloc)
               for a in (donor.active_altloc, rctx.pi_alt)):
            continue
        position = rctx.ns.grid_cell.find_nearest_pbc_position(reference, atom.pos, mark.image_idx)
        image = rctx.ns.grid_cell.find_nearest_pbc_image(reference, atom.pos, mark.image_idx)
        same_res = (cra.chain.name == donor.chain_name and str(residue.seqid).strip() == donor.residue_seqid
                    and residue.name == donor.residue_name)
        if same_res:
            if (atom.name == donor.x_atom.name and atom.altloc == donor.x_atom.altloc
                    and position.dist(reference) < 0.01):
                continue
            parent_distance = position.dist(gemmi.Position(*parent_pos))
            if (atom.name == donor.parent_atom.name and atom.altloc == donor.parent_atom.altloc
                    and parent_distance < 0.01):
                continue
            if atom.name in one_three and parent_distance < (
                    atom.element.covalent_r + donor.parent_atom.element.covalent_r + .25):
                continue
        key = (cra.chain.name, str(residue.seqid).strip(), residue.name, atom.name)
        dedup = key + (atom.altloc, image.symmetry_code())
        if dedup in seen:
            continue
        seen.add(dedup)
        special = any(acceptors.altlocs_compatible(label, source.altloc)
                      and acceptors.altlocs_compatible(label, target.altloc)
                      for source, target in rctx.connections.get(key, ()))
        if not special:
            bonded = monomer_bonds.get_bonded_heavy_atom_names(residue.name, atom.name)
            for neighbour in rctx.ns.find_atoms(atom.pos, radius=2.8):
                other_cra = neighbour.to_cra(rctx.model)
                other = other_cra.atom
                if (other.is_hydrogen() or other.occ < config.MIN_ATOM_OCCUPANCY
                        or not acceptors.altlocs_compatible(label, other.altloc)
                        or not acceptors.altlocs_compatible(atom.altloc, other.altloc)):
                    continue
                near = rctx.ns.grid_cell.find_nearest_pbc_position(atom.pos, other.pos, neighbour.image_idx)
                distance = atom.pos.dist(near)
                if distance < .01:
                    continue
                expected_bond = (bonded is not None and other.name in bonded
                                 and other_cra.chain.name == cra.chain.name
                                 and other_cra.residue.seqid == residue.seqid
                                 and other_cra.residue.name == residue.name
                                 and near.dist(other.pos) < .01)
                if (other.element.name in METALS or
                        (not expected_bond and distance < atom.element.covalent_r + other.element.covalent_r + .25)):
                    special = True
                    break
        site = (rctx.protonation.lookup(rctx.model_index, cra.chain.name, residue, atom, label)
                if rctx.protonation is not None else None)
        status, reason, source = classify_acceptor(residue, atom, label, site, special)
        if not same_res and (atom.altloc not in {"", "\0"}) and donor.active_altloc:
            if status == "yes":
                status, reason = "unknown", "cross_residue_altloc_correlation_unestablished"
        result.append({"chain": key[0], "residue_id": key[1], "residue_name": key[2],
                       "atom": atom.name, "altloc": atom.altloc.strip("\0"),
                       "symmetry_code": image.symmetry_code(), "position": position.tolist(),
                       "occupancy": atom.occ, "status": status, "reason": reason, "source": source})
    return sorted(result, key=lambda a: (a["chain"], a["residue_id"], a["atom"], a["altloc"], a["symmetry_code"]))
