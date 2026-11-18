"""Histidine states without ``--ph``: explicit per-site overrides and opt-in policies.

Policies (``his_policy`` in ``api.detect``):

- ``legacy`` (default): the 2.1.13 behaviour. With ``his_overrides``, only the named His
  sites change; every other His keeps the legacy preparation.
- ``neutral``: every His neutral; the tautomer from the local discriminator.
- ``default``: neutral His, with only high-confidence His+ (PROPKA, pKa - pH >= 1.5) cationic.

Overrides map a residue key ``(model_index, chain, seqid, "HIS", altloc)`` to
``"delta"`` (HID), ``"epsilon"`` (HIE) or ``"doubly_protonated"`` (HIP: the ring is not a
pi acceptor and has no N-H donor). An override chooses the microstate but does not bypass
a severe clash of the placed H; such a site stays neutral without an N-H donor.

Unlike the opt-in ``--ph`` model, this touches histidine only: other titratable
residues and termini keep the dictionary preparation, a missing pH is not an
error, and an undecided His stays a neutral pi acceptor without an N-H donor.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

import gemmi

from . import covalent_state, hbond_acceptors, hydrogen_prep, monomer_bonds, protonation
from . import protonation_geometry as geometry
from . import protonation_pka as pka

logger = logging.getLogger("xpid")
POLICIES = ("default", "neutral", "legacy")
# pKa - pH needed for a His+ call: the margin of the XH-pi pH sensitivity study.
HIP_MARGIN = 1.5
HIS_H = protonation.SITE_H["HIS"]


@dataclass
class HisContext(protonation.ProtonationContext):
    """Only sites of known charge +1 lose their ring; His without a site keeps it."""

    def pi_allowed(self, model_index, chain, residue, label):
        if residue.name != "HIS":
            return True
        atom = geometry.atom_for(residue, "ND1", label)
        site = None if atom is None else self.lookup(model_index, chain, residue, atom, label)
        return site is None or site.charge != 1

    def hit_fields(self, *args, **kwargs):
        return {}


def validate_options(policy, his_ph=None):
    if policy not in POLICIES:
        raise ValueError(f"His policy must be one of {', '.join(POLICIES)}")
    if his_ph is not None and str(his_ph).lower() != "auto":
        pka.numeric_ph(his_ph)
    if policy == "default":
        try:
            import propka.run  # noqa: F401
        except ImportError as exc:
            raise ValueError('The default His policy requires PROPKA; his_policy="neutral" '
                             'needs neither PROPKA nor a pH') from exc


def _resolve_ph(path, his_ph, context):
    try:
        return pka.resolve_ph(path, "auto" if his_ph is None else his_ph)
    except ValueError as exc:
        context.warnings.append(f"No pH for His charge: {exc}; all His neutral")
        return None, ""


def _site(key, label, residue, chain, model, ns, connections, excluded, monlib, prediction, ph, override):
    site = protonation.Site(key, source="his_default", charge=0, donor_h={name: [] for name in HIS_H})
    legacy_acceptors = {}
    for name in HIS_H:
        atom = geometry.atom_for(residue, name, label)
        legacy_acceptors[name] = atom is not None and hbond_acceptors.is_hbond_acceptor(residue, atom, label)
    site.acceptors = legacy_acceptors
    if prediction:
        site.pka = prediction["pka"]
    if override:
        site.source = "user_override"
    if key[:4] in excluded:
        site.reason = "hydrogen_preparation_quarantined"
    else:
        site.reason = protonation._special_state(chain.name, residue, label, model, ns, connections) or ""
    if site.reason:
        site.state = "neutral_unresolved"
        return site, []
    coupled = bool(prediction and prediction["coupled"])
    if override == "doubly_protonated" or (
            not override and ph is not None and site.pka is not None
            and not coupled and site.pka - ph >= HIP_MARGIN):
        site.state, site.status, site.charge = "doubly_protonated", "resolved", 1
        return site, []
    hd = geometry.hydrogen_copy(residue, "ND1", "HD1", label, monlib)
    he = geometry.hydrogen_copy(residue, "NE2", "HE2", label, monlib)
    if hd is None or he is None:
        site.state, site.reason = "neutral_unresolved", "missing_dictionary_hydrogen_geometry"
        return site, []
    selected, site.candidates = geometry.choose_candidate(
        [("delta", {"ND1": [hd]}), ("epsilon", {"NE2": [he]})],
        residue, chain.name, label, model, ns, override)
    if selected is None:
        site.state, site.reason = "neutral_unresolved", "ambiguous_or_clashing_microstate"
        return site, []
    site.state, hydrogens = selected
    site.status = "resolved"
    site.donor_h.update({parent: [h.name for h in values] for parent, values in hydrogens.items()})
    site.acceptors = {name: name not in hydrogens for name in HIS_H}
    return site, [h for values in hydrogens.values() for h in values]


def prepare(structure, path, policy="default", his_ph=None, overrides=None,
            excluded_residues=None, model_mode=0):
    """Apply after Gemmi H rebuilding, before any neighbour search/detection."""
    context = HisContext(None, "", f"his_{policy}")
    if policy == "default":
        context.ph, context.ph_source = _resolve_ph(path, his_ph, context)
    overrides = dict(overrides or {})
    excluded_residues = excluded_residues or set()
    monlib = hydrogen_prep._get_shared_monlib({"HIS"})
    for mi, model in enumerate(structure):
        if model_mode != "all" and mi != int(model_mode):
            continue
        predictions = {}
        if context.ph is not None:
            messages = []
            try:
                predictions, context.backend_version = pka.predict_model(model, mi, messages)
            except Exception as exc:
                context.warnings.append(f"Model {mi}: PROPKA failed: {exc}; all His neutral")
                logger.warning("Model %d: PROPKA failed; all His neutral", mi)
            context.warnings.extend(f"PROPKA model {mi}: {message}" for message in messages)
        ns = gemmi.NeighborSearch(model, structure.cell, 4.0).populate(include_h=False)
        connections = covalent_state.connection_index(structure, model)
        changes = []
        for chain in model:
            for residue in chain:
                if residue.name != "HIS":
                    continue
                labels = sorted({pka.alt(a) for a in residue if not a.is_hydrogen() and pka.alt(a)}) or [""]
                if policy == "legacy" and not any(
                        pka.residue_key(mi, chain.name, residue, label) in overrides for label in labels):
                    continue   # legacy: only the overridden sites change
                replacement = []
                for label in labels:
                    key = pka.residue_key(mi, chain.name, residue, label)
                    site, hydrogens = _site(key, label, residue, chain, model, ns, connections,
                                            excluded_residues, monlib, predictions.get(key),
                                            context.ph, overrides.pop(key, None))
                    replacement.extend(hydrogens)
                    context.add(site)
                changes.append((residue, replacement))
        for residue, replacement in changes:
            for i in range(len(residue) - 1, -1, -1):
                atom = residue[i]
                if atom.is_hydrogen() and monomer_bonds.canonical_hydrogen_name(atom.name) in ("HD1", "HE2"):
                    del residue[i]
            for hydrogen in replacement:
                residue.add_atom(hydrogen)
    if overrides:
        raise ValueError(f"His override(s) match no site: {sorted(overrides)}")
    return context
