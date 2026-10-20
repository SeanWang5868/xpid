"""Actual covalent-state checks for terminal Cone donor groups."""
from __future__ import annotations

import gemmi

from . import config, hbond_acceptors, monomer_bonds, rotatable_groups


def _address_key(address: gemmi.AtomAddress) -> tuple:
    return (address.chain_name, str(address.res_id.seqid).strip(),
            address.res_id.name, address.atom_name)


def connection_index(structure: gemmi.Structure, model: gemmi.Model) -> dict:
    """Keep declared cross-links even when H preparation used a recovery clone."""
    index = {}
    bond_types = {gemmi.ConnectionType.Covale, gemmi.ConnectionType.Disulf,
                  gemmi.ConnectionType.MetalC}
    for connection in structure.connections:
        if connection.type not in bond_types:
            continue
        for source, target in ((connection.partner1, connection.partner2),
                               (connection.partner2, connection.partner1)):
            partner = model.find_cra(target)
            if partner.atom is not None and partner.atom.element.name.upper() in {"H", "D"}:
                continue
            index.setdefault(_address_key(source), []).append((source, target))
    return index


def allows_cone(residue: gemmi.Residue, chain_name: str,
                conformer: rotatable_groups.DonorConformer,
                definition: rotatable_groups.RotatableGroupDefinition,
                connections: dict, ns: gemmi.NeighborSearch,
                model: gemmi.Model, formal_charge=None) -> bool:
    """A terminal group needs one heavy parent and an available X–H valence.

    Declared covalent/coordination links take precedence over distance guesses.
    A short, compatible heavy-atom contact also prevents virtual protonation of
    an unrecorded ether, disulfide or substituted methyl. This is deliberately
    conservative and does not assign protonation from proximity alone.
    """
    x, parent = conformer.x_atom, conformer.parent_atom
    altloc = conformer.active_altloc
    charge = x.charge if formal_charge is None else formal_charge
    if x.element.name.upper() != definition.element or charge != 0:
        return False
    heavy_names = monomer_bonds.get_bonded_heavy_atom_names(residue.name, x.name)
    if heavy_names is not None and heavy_names != {definition.parent_atom_name}:
        return False
    h_names = monomer_bonds.get_bonded_hydrogen_names(residue.name, x.name)
    if h_names is not None and len(h_names) != definition.hydrogen_count:
        return False
    key = (chain_name, str(residue.seqid).strip(), residue.name, x.name)
    for source, target in connections.get(key, ()):
        if not hbond_acceptors.altlocs_compatible(altloc, source.altloc):
            continue
        if not hbond_acceptors.altlocs_compatible(altloc, target.altloc):
            continue
        if (_address_key(target) != key[:-1] + (parent.name,)):
            return False
    # The parent must itself form a plausible covalent heavy-atom bond.
    parent_limit = x.element.covalent_r + parent.element.covalent_r + 0.25
    if not (config.MIN_COVALENT_XH < x.pos.dist(parent.pos) <= parent_limit):
        return False
    for mark in ns.find_atoms(x.pos, radius=2.5):
        cra = mark.to_cra(model)
        atom = cra.atom
        if atom.element.name.upper() in {"H", "D"}:
            continue
        if atom.occ < config.MIN_ATOM_OCCUPANCY:
            continue
        if not hbond_acceptors.altlocs_compatible(altloc, atom.altloc):
            continue
        pos = ns.grid_cell.find_nearest_pbc_position(x.pos, atom.pos, mark.image_idx)
        distance = x.pos.dist(pos)
        if distance < 0.01:
            continue
        if (cra.chain.name == chain_name and cra.residue.seqid == residue.seqid
                and cra.residue.name == residue.name and atom.name == parent.name
                and atom.altloc == parent.altloc and pos.dist(parent.pos) < 0.01):
            continue
        if distance <= x.element.covalent_r + atom.element.covalent_r + 0.25:
            return False
    return True
