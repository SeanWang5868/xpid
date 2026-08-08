#!/usr/bin/env python3
"""Generate the complete Section-02 chemical-landscape evidence package.

The input must be a schema-v8 ``generate_background.py`` output.  The master
contact universe contains primary neutral standard-protein donor-ring pairs for
which either endpoint is on a PISCES-selected chain.  Endpoint propensities use
the corresponding selected donor or selected acceptor risk set.

All uncertainty intervals use the PDB structure as the resampling unit.  No
contact-level p values are produced because contacts within a structure are not
independent observations.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from background_core import _canonical_ring_pair_id
from xpid import geometry
from analyze_residue_donor import (
    DONOR_CLASS_ORDER,
    PRIMARY_DONOR_CLASSES,
    STANDARD_AMINO_ACIDS,
    is_primary_donor,
)


ACCEPTOR_ORDER = ("PHE-6", "TYR-6", "TRP-6", "TRP-5", "HIS-5")
ANGLE_SCENARIOS = (
    ("strict_5deg", 35.0, 125.0),
    ("primary", 40.0, 120.0),
    ("threshold_envelope_4p82deg", 44.82, 115.18),
    ("relaxed_5deg", 45.0, 115.0),
)
DISTANCE_CONTRASTS = (
    ("Amide side-chain N-H", "Methylene C-H"),
    ("Aromatic side-chain N-H", "Methylene C-H"),
    ("Backbone N-H", "Methylene C-H"),
    ("Hydroxyl O-H", "Methylene C-H"),
    ("Aromatic C-H", "Methylene C-H"),
    ("Thiol S-H", "Hydroxyl O-H"),
    ("Methine C-H", "Methylene C-H"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260801)
    parser.add_argument("--max-structures", type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def text(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acceptor_type(residue: Any, ring_size: Any) -> str:
    residue = text(residue).upper()
    size = int(ring_size) if text(ring_size) else 0
    label = f"{residue}-{size}"
    return label if label in ACCEPTOR_ORDER else "Other"


def read_columns(path: Path, columns: Sequence[str]) -> pd.DataFrame:
    names = set(pq.read_schema(path).names)
    missing = [name for name in columns if name not in names]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    return pq.read_table(path, columns=list(columns)).to_pandas()


def read_columns_or_empty(path: Path, columns: Sequence[str]) -> pd.DataFrame:
    """Read a partition, allowing generator-omitted zero-row partitions."""
    if not path.exists():
        return pd.DataFrame({name: pd.Series(dtype="object") for name in columns})
    return read_columns(path, columns)


def collapse_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse alternate conformers using existence semantics."""
    if frame.empty:
        return frame.copy()
    frame = frame.copy()
    occupancy = pd.to_numeric(frame["combined_occupancy"], errors="coerce").fillna(0)
    representative = frame.loc[
        occupancy.groupby(frame["spatial_pair_id"]).idxmax()
    ].copy()
    representative = representative.set_index("spatial_pair_id", drop=False)
    binary = (
        "within_element_distance_cutoff", "within_common_4p5_distance_cutoff",
        "ring_face_spatial", "hudson_spatial", "plevin_spatial",
        "is_positive_pair", "is_hudson_positive_pair",
        "is_plevin_positive_pair", "is_pi_pi_tshaped_spatial",
    )
    for field in binary:
        if field in frame.columns:
            values = pd.to_numeric(frame[field], errors="coerce")
            grouped = values.groupby(frame["spatial_pair_id"]).max()
            representative[field] = grouped.reindex(representative.index).values
    # Prefer a conformer with reconstructed ring or N-H geometry when present.
    for status_field, value_fields in (
        ("pi_pi_geometry_status", (
            "donor_aromatic_ring_id", "donor_aromatic_ring_site_id",
            "donor_aromatic_ring_image_id",
            "aromatic_ring_pair_id", "pi_pi_centroid_distance",
            "pi_pi_interplanar_angle", "pi_pi_lateral_offset",
            "is_pi_pi_tshaped_spatial",
        )),
        ("backbone_nh_direction_status", (
            "backbone_nh_hydrogen_count", "backbone_nh_theta_deg",
            "backbone_nh_xh_pi_angle_deg",
        )),
    ):
        if status_field not in frame.columns:
            continue
        ok = frame[frame[status_field] == "ok"]
        if ok.empty:
            continue
        ok_occ = pd.to_numeric(ok["combined_occupancy"], errors="coerce").fillna(0)
        chosen = ok.loc[ok_occ.groupby(ok["spatial_pair_id"]).idxmax()]
        chosen = chosen.set_index("spatial_pair_id")
        common = representative.index.intersection(chosen.index)
        representative.loc[common, status_field] = "ok"
        for field in value_fields:
            if field in chosen.columns:
                representative.loc[common, field] = chosen.loc[common, field]
    return representative.reset_index(drop=True)


def primary_mask(frame: pd.DataFrame) -> pd.Series:
    return pd.Series([
        is_primary_donor(donor_class, residue, assignment)
        for donor_class, residue, assignment in zip(
            frame["donor_class"], frame["x_residue_name"],
            frame["cys_sg_thiol_assignment"],
        )
    ], index=frame.index)


def q025(values: Sequence[float]) -> float:
    return float(np.quantile(values, 0.025)) if values else float("nan")


def q975(values: Sequence[float]) -> float:
    return float(np.quantile(values, 0.975)) if values else float("nan")


def ratio_bootstrap(
    per_pdb: pd.DataFrame,
    numerator: str,
    denominator: str,
    pdb_ids: Sequence[str],
    rng: np.random.Generator,
    replicates: int,
) -> Tuple[float, float, float]:
    table = per_pdb.set_index("pdb_id").reindex(pdb_ids, fill_value=0)
    num = table[numerator].to_numpy(dtype=float)
    den = table[denominator].to_numpy(dtype=float)
    estimate = float(num.sum() / den.sum()) if den.sum() else float("nan")
    samples: List[float] = []
    n = len(pdb_ids)
    for _ in range(replicates):
        weights = rng.multinomial(n, np.full(n, 1.0 / n))
        sampled_den = float(np.dot(weights, den))
        if sampled_den:
            samples.append(float(np.dot(weights, num) / sampled_den))
    return estimate, q025(samples), q975(samples)


def ratio_of_ratios_bootstrap(
    per_pdb: pd.DataFrame,
    a_num: str, a_den: str, b_num: str, b_den: str,
    pdb_ids: Sequence[str], rng: np.random.Generator, replicates: int,
) -> Tuple[float, float, float]:
    table = per_pdb.set_index("pdb_id").reindex(pdb_ids, fill_value=0)
    arrays = [table[name].to_numpy(dtype=float) for name in (a_num, a_den, b_num, b_den)]
    def statistic(weights: np.ndarray) -> float:
        an, ad, bn, bd = (float(np.dot(weights, value)) for value in arrays)
        return (an / ad) / (bn / bd) if ad and bn and bd else float("nan")
    point = statistic(np.ones(len(pdb_ids)))
    samples: List[float] = []
    n = len(pdb_ids)
    for _ in range(replicates):
        value = statistic(rng.multinomial(n, np.full(n, 1.0 / n)))
        if math.isfinite(value):
            samples.append(value)
    return point, q025(samples), q975(samples)


def conditional_enrichment_bootstrap(
    positives: pd.DataFrame,
    ring_face: pd.DataFrame,
    donor_order: Sequence[str],
    acceptor_order: Sequence[str],
    pdb_ids: Sequence[str],
    rng: np.random.Generator,
    replicates: int,
) -> Dict[Tuple[str, str], Tuple[float, float, float]]:
    """Bootstrap acceptor enrichment conditional on donor class.

    For each donor class, the acceptor distribution among positive pairs is
    divided by the acceptor distribution in the H-independent ring-face risk
    set.  PDB structures, rather than contacts, are resampled.
    """
    columns = pd.MultiIndex.from_product(
        [list(donor_order), list(acceptor_order)],
        names=["donor_class", "acceptor_type"],
    )

    def count_matrix(frame: pd.DataFrame) -> np.ndarray:
        if frame.empty:
            return np.zeros((len(pdb_ids), len(columns)), dtype=float)
        table = frame.groupby(
            ["pdb_id", "donor_class", "acceptor_type"]
        ).size().unstack(["donor_class", "acceptor_type"], fill_value=0)
        table = table.reindex(index=pdb_ids, columns=columns, fill_value=0)
        return table.to_numpy(dtype=float)

    positive_matrix = count_matrix(positives)
    spatial_matrix = count_matrix(ring_face)
    n_donors, n_acceptors = len(donor_order), len(acceptor_order)

    def enrichments(
        positive_counts: np.ndarray,
        spatial_counts: np.ndarray,
    ) -> np.ndarray:
        shape = positive_counts.shape[:-1] + (n_donors, n_acceptors)
        p = positive_counts.reshape(shape)
        s = spatial_counts.reshape(shape)
        with np.errstate(divide="ignore", invalid="ignore"):
            p_share = p / p.sum(axis=-1, keepdims=True)
            s_share = s / s.sum(axis=-1, keepdims=True)
            values = p_share / s_share
        return values.reshape(positive_counts.shape)

    point = enrichments(
        positive_matrix.sum(axis=0, keepdims=True),
        spatial_matrix.sum(axis=0, keepdims=True),
    )[0]
    samples: List[np.ndarray] = []
    n = len(pdb_ids)
    probabilities = np.full(n, 1.0 / n)
    for start in range(0, replicates, 128):
        size = min(128, replicates - start)
        weights = rng.multinomial(n, probabilities, size=size)
        samples.append(enrichments(
            weights @ positive_matrix,
            weights @ spatial_matrix,
        ))
    sample_matrix = np.concatenate(samples, axis=0)
    low = np.nanquantile(sample_matrix, 0.025, axis=0)
    high = np.nanquantile(sample_matrix, 0.975, axis=0)
    return {
        tuple(key): (float(point[index]), float(low[index]), float(high[index]))
        for index, key in enumerate(columns)
    }


def median_cluster_bootstrap(
    arrays_by_pdb: Mapping[str, Mapping[str, np.ndarray]],
    group_a: str, group_b: str, pdb_ids: Sequence[str],
    rng: np.random.Generator, replicates: int,
) -> Tuple[float, float, float, float, float]:
    def pooled(group: str, sampled: Iterable[str]) -> np.ndarray:
        values = [arrays_by_pdb.get(pdb, {}).get(group, np.empty(0)) for pdb in sampled]
        values = [value for value in values if len(value)]
        return np.concatenate(values) if values else np.empty(0)
    a = pooled(group_a, pdb_ids)
    b = pooled(group_b, pdb_ids)
    median_a = float(np.median(a)) if len(a) else float("nan")
    median_b = float(np.median(b)) if len(b) else float("nan")
    point = median_a - median_b
    samples: List[float] = []
    n = len(pdb_ids)
    for _ in range(replicates):
        sampled = rng.choice(pdb_ids, size=n, replace=True)
        aa, bb = pooled(group_a, sampled), pooled(group_b, sampled)
        if len(aa) and len(bb):
            samples.append(float(np.median(aa) - np.median(bb)))
    return median_a, median_b, point, q025(samples), q975(samples)


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    pd.DataFrame(list(rows)).to_csv(path, index=False)


def build_asu_aromatic_ring_background(rings: pd.DataFrame) -> pd.DataFrame:
    """Build the complete non-symmetry ring-pair universe from stored rings."""
    rows: Dict[str, Dict[str, Any]] = {}
    scores: Dict[str, Tuple[float, float]] = {}
    identity = np.eye(3)
    zero = np.zeros(3)
    for pdb_id, group in rings.groupby("pdb_id", sort=False):
        records = group.to_dict("records")
        for i, first in enumerate(records):
            first_center = np.array([
                first["ring_center_x"], first["ring_center_y"],
                first["ring_center_z"],
            ], dtype=float)
            first_normal = np.array([
                first["ring_normal_x"], first["ring_normal_y"],
                first["ring_normal_z"],
            ], dtype=float)
            for second in records[i + 1:]:
                if first["ring_site_id"] == second["ring_site_id"]:
                    continue
                if not (
                    int(first["is_standard_protein_acceptor"]) == 1
                    or int(second["is_standard_protein_acceptor"]) == 1
                ):
                    continue
                if not (
                    int(first["is_selected_chain"]) == 1
                    or int(second["is_selected_chain"]) == 1
                ):
                    continue
                if (
                    text(first["altloc"]) and text(second["altloc"])
                    and text(first["altloc"]) != text(second["altloc"])
                ):
                    continue
                second_center = np.array([
                    second["ring_center_x"], second["ring_center_y"],
                    second["ring_center_z"],
                ], dtype=float)
                second_normal = np.array([
                    second["ring_normal_x"], second["ring_normal_y"],
                    second["ring_normal_z"],
                ], dtype=float)
                distance, angle, offset = geometry.calculate_pi_pi_geometry(
                    first_center, first_normal, second_center, second_normal)
                if not 3.0 <= distance <= 5.5:
                    continue
                pair_id = _canonical_ring_pair_id(
                    first["ring_site_id"], second["ring_site_id"],
                    identity, zero)
                combined_occ = min(
                    float(first["mean_occupancy"]),
                    float(second["mean_occupancy"]),
                )
                is_tshaped = int(angle >= 60.0)
                row = {
                    "pdb_id": pdb_id,
                    "aromatic_ring_pair_id": pair_id,
                    "ring_a_site_id": min(
                        first["ring_site_id"], second["ring_site_id"]),
                    "ring_b_site_id": max(
                        first["ring_site_id"], second["ring_site_id"]),
                    "either_endpoint_on_selected_chain": 1,
                    "is_symmetry_contact": 0,
                    "combined_occupancy": combined_occ,
                    "pi_pi_centroid_distance": float(distance),
                    "pi_pi_interplanar_angle": float(angle),
                    "pi_pi_lateral_offset": float(offset),
                    "is_pi_pi_tshaped_spatial": is_tshaped,
                }
                score = (combined_occ, -float(distance))
                if pair_id not in rows or score > scores[pair_id]:
                    previous = int(rows.get(pair_id, {}).get(
                        "is_pi_pi_tshaped_spatial", 0))
                    row["is_pi_pi_tshaped_spatial"] = max(
                        previous, is_tshaped)
                    rows[pair_id] = row
                    scores[pair_id] = score
                elif is_tshaped:
                    rows[pair_id]["is_pi_pi_tshaped_spatial"] = 1
    return pd.DataFrame(list(rows.values()))


def main() -> int:
    args = parse_args()
    root = args.background.resolve()
    metadata_path = root / "background_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    if int(metadata.get("output_schema_version", 0)) < 8:
        raise ValueError("Schema v8 is required; rerun generate_background.py with the updated code.")
    output = (args.output or root / "analysis" / "chemical_landscape").resolve()
    expected_outputs = (
        "global_summary.csv", "global_structure_summary.csv",
        "main_contact_composition.csv", "acceptor_composition.csv",
        "acceptor_propensity.csv", "donor_acceptor_matrix.csv",
        "tshape_summary.csv", "gly_backbone_nh_funnel.csv",
        "backbone_nh_angle_sensitivity.csv", "distance_by_donor_class.csv",
        "backbone_nh_primary_mismatches.csv", "bootstrap_estimates.csv",
        "qc_report.json", "analysis_metadata.json",
    )
    if output.exists() and not args.force and any((output / name).exists() for name in expected_outputs):
        raise FileExistsError(f"{output} already contains analysis files; use --force or a new --output.")
    output.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(root / "structure_manifest.csv")
    pdb_ids = manifest.loc[manifest["status"] == "ok", "pdb_id"].astype(str).tolist()
    if args.max_structures:
        pdb_ids = pdb_ids[:args.max_structures]

    # Full-model headline statistics.  The pair summary contains unique
    # physical pairs, whereas positive_rows_matched contains raw XPID rows;
    # both are retained explicitly so the two units cannot be conflated.
    manifest_ok = manifest[manifest["pdb_id"].astype(str).isin(pdb_ids)].copy()
    pair_summary = pd.read_csv(root / "summaries" / "pair_summary.csv")
    pair_summary["pdb_id"] = pair_summary["pdb_id"].astype(str)
    pair_summary = pair_summary[pair_summary["pdb_id"].isin(pdb_ids)]
    pair_by_pdb = pair_summary.groupby("pdb_id")["positive_pair_count"].sum()
    global_structure = manifest_ok[[
        "pdb_id", "standard_amino_acid_position_count", "positive_rows_matched"
    ]].copy()
    global_structure["unique_positive_pair_count"] = (
        global_structure["pdb_id"].map(pair_by_pdb).fillna(0).astype(int)
    )
    global_structure["positive_pairs_per_1000_standard_residues"] = (
        1000.0 * global_structure["unique_positive_pair_count"]
        / global_structure["standard_amino_acid_position_count"]
    )
    density = global_structure[
        "positive_pairs_per_1000_standard_residues"
    ].to_numpy(dtype=float)
    unique_positive_total = int(global_structure["unique_positive_pair_count"].sum())
    raw_positive_total = int(global_structure["positive_rows_matched"].sum())
    residue_total = int(global_structure["standard_amino_acid_position_count"].sum())
    structures_with_positive = int(
        (global_structure["unique_positive_pair_count"] > 0).sum())
    global_summary_rows = [{
        "structure_count": len(global_structure),
        "structures_with_at_least_one_positive_pair": structures_with_positive,
        "structure_coverage_fraction": structures_with_positive / len(global_structure),
        "structure_coverage_percent": 100.0 * structures_with_positive / len(global_structure),
        "unique_physical_positive_pair_count": unique_positive_total,
        "raw_xpid_positive_row_count": raw_positive_total,
        "standard_amino_acid_residue_count": residue_total,
        "pooled_positive_pairs_per_1000_residues": 1000.0 * unique_positive_total / residue_total,
        "median_positive_pairs_per_1000_residues": float(np.median(density)),
        "q1_positive_pairs_per_1000_residues": float(np.quantile(density, 0.25)),
        "q3_positive_pairs_per_1000_residues": float(np.quantile(density, 0.75)),
        "residues_per_contact_from_median_structure_density": float(
            1000.0 / np.median(density)
        ),
        "pair_unit": "unique physical donor-ring spatial_pair_id",
        "residue_unit": "unique standard-amino-acid residue position in the full coordinate model",
    }]

    main_positive_parts: List[pd.DataFrame] = []
    spatial_matrix_parts: List[pd.DataFrame] = []
    distance_matrix_parts: List[pd.DataFrame] = []
    acceptor_ring_parts: List[pd.DataFrame] = []
    aromatic_ring_source_parts: List[pd.DataFrame] = []
    aromatic_ring_background_parts: List[pd.DataFrame] = []
    gly_donor_parts: List[pd.DataFrame] = []
    gly_distance_parts: List[pd.DataFrame] = []
    gly_spatial_conformer_parts: List[pd.DataFrame] = []
    failures: List[Dict[str, str]] = []
    absent_partition_counts = defaultdict(int)

    pair_columns = (
        "pdb_id", "spatial_pair_id", "ring_site_id", "donor_site_id",
        "pi_residue_name", "pi_ring_size", "x_residue_name", "x_atom_name",
        "pi_altloc", "x_altloc",
        "x_element", "donor_class", "donor_treatment", "cys_sg_thiol_assignment",
        "donor_on_selected_chain", "acceptor_on_selected_chain",
        "either_endpoint_on_selected_chain", "is_symmetry_contact",
        "combined_occupancy", "dist_X_plane", "dist_X_centroid",
        "projection_distance", "angle_XPCN", "within_element_distance_cutoff",
        "within_common_4p5_distance_cutoff", "ring_face_spatial",
        "hudson_spatial", "plevin_spatial", "is_positive_pair",
        "is_hudson_positive_pair", "is_plevin_positive_pair",
        "donor_aromatic_ring_site_id", "donor_aromatic_ring_image_id",
        "aromatic_ring_pair_id",
        "pi_pi_centroid_distance", "pi_pi_interplanar_angle",
        "pi_pi_lateral_offset", "is_pi_pi_tshaped_spatial",
        "pi_pi_geometry_status", "backbone_nh_hydrogen_count",
        "backbone_nh_theta_deg", "backbone_nh_xh_pi_angle_deg",
        "backbone_nh_direction_status",
    )
    ring_columns = (
        "pdb_id", "ring_site_id", "residue_name", "ring_id", "ring_size",
        "is_standard_protein_acceptor", "is_selected_chain", "eligible",
        "mean_occupancy", "altloc", "ring_center_x", "ring_center_y",
        "ring_center_z", "ring_normal_x", "ring_normal_y", "ring_normal_z",
    )
    donor_columns = (
        "pdb_id", "donor_site_id", "residue_name", "donor_class",
        "cys_sg_thiol_assignment", "is_selected_chain", "eligible", "occupancy",
    )
    aromatic_ring_pair_columns = (
        "pdb_id", "aromatic_ring_pair_id", "ring_a_site_id",
        "ring_b_site_id", "ring_a_residue_name", "ring_a_size",
        "ring_b_residue_name", "ring_b_size",
        "either_endpoint_on_selected_chain", "is_symmetry_contact",
        "combined_occupancy", "pi_pi_centroid_distance",
        "pi_pi_interplanar_angle", "pi_pi_lateral_offset",
        "is_pi_pi_tshaped_spatial",
    )

    for index, pdb_id in enumerate(pdb_ids, 1):
        try:
            spatial_path = root / "datasets" / "spatial_pair_background" / f"{pdb_id}.parquet"
            distance_path = root / "datasets" / "distance_pair_background" / f"{pdb_id}.parquet"
            ring_path = root / "datasets" / "acceptor_ring_background" / f"{pdb_id}.parquet"
            donor_path = root / "datasets" / "donor_site_background" / f"{pdb_id}.parquet"
            aromatic_ring_path = root / "datasets" / "aromatic_ring_pair_background" / f"{pdb_id}.parquet"
            for label, path in (
                ("spatial_pair_background", spatial_path),
                ("distance_pair_background", distance_path),
                ("acceptor_ring_background", ring_path),
                ("donor_site_background", donor_path),
                ("aromatic_ring_pair_background", aromatic_ring_path),
            ):
                if not path.exists():
                    absent_partition_counts[label] += 1
            spatial_raw = read_columns_or_empty(spatial_path, pair_columns)
            distance_raw = read_columns_or_empty(distance_path, pair_columns)
            spatial = collapse_pairs(spatial_raw)
            distance = collapse_pairs(distance_raw)
            for frame in (spatial, distance):
                frame["acceptor_type"] = [
                    acceptor_type(residue, size)
                    for residue, size in zip(frame["pi_residue_name"], frame["pi_ring_size"])
                ]
                frame["primary_donor"] = primary_mask(frame)
            main_positive = spatial[
                spatial["primary_donor"]
                & (spatial["either_endpoint_on_selected_chain"] == 1)
                & (spatial["is_positive_pair"] == 1)
            ].copy()
            main_positive_parts.append(main_positive)
            spatial_matrix_parts.append(spatial[
                spatial["primary_donor"]
                & (spatial["either_endpoint_on_selected_chain"] == 1)
                & spatial["acceptor_type"].isin(ACCEPTOR_ORDER)
            ].copy())
            distance_matrix_parts.append(distance[
                distance["primary_donor"]
                & (distance["either_endpoint_on_selected_chain"] == 1)
                & (distance["within_element_distance_cutoff"] == 1)
                & distance["acceptor_type"].isin(ACCEPTOR_ORDER)
            ].copy())

            rings = read_columns_or_empty(ring_path, ring_columns)
            aromatic_ring_source_parts.append(
                rings[rings["eligible"] == 1].copy())
            rings = rings.sort_values("mean_occupancy", ascending=False).drop_duplicates("ring_site_id")
            acceptor_ring_parts.append(rings[
                (rings["eligible"] == 1)
                & (rings["is_selected_chain"] == 1)
                & (rings["is_standard_protein_acceptor"] == 1)
            ].copy())

            donors = read_columns_or_empty(donor_path, donor_columns)
            donors = donors.sort_values("occupancy", ascending=False).drop_duplicates("donor_site_id")
            gly_donor_parts.append(donors[
                (donors["eligible"] == 1)
                & (donors["is_selected_chain"] == 1)
                & (donors["donor_class"] == "Backbone N-H")
                & donors["residue_name"].isin(STANDARD_AMINO_ACIDS)
                & (donors["residue_name"] != "PRO")
            ].copy())
            gly_distance_parts.append(distance[
                distance["primary_donor"]
                & (distance["donor_on_selected_chain"] == 1)
                & (distance["donor_class"] == "Backbone N-H")
                & (distance["within_element_distance_cutoff"] == 1)
            ].copy())
            gly_spatial_conformer_parts.append(spatial_raw[
                primary_mask(spatial_raw)
                & (spatial_raw["donor_on_selected_chain"] == 1)
                & (spatial_raw["donor_class"] == "Backbone N-H")
            ].copy())
            aromatic_ring_background = read_columns_or_empty(
                aromatic_ring_path, aromatic_ring_pair_columns)
            aromatic_ring_background_parts.append(
                aromatic_ring_background[
                    aromatic_ring_background[
                        "either_endpoint_on_selected_chain"] == 1
                ].copy())
        except Exception as exc:
            failures.append({"pdb_id": pdb_id, "error": repr(exc)})
        if index % 250 == 0:
            print(f"Loaded {index}/{len(pdb_ids)} structures", flush=True)
    if failures:
        write_csv(output / "processing_failures.csv", failures)
        raise RuntimeError(f"{len(failures)} structures failed; see processing_failures.csv")

    positives = pd.concat(main_positive_parts, ignore_index=True)
    spatial = pd.concat(spatial_matrix_parts, ignore_index=True)
    distance = pd.concat(distance_matrix_parts, ignore_index=True)
    rings = pd.concat(acceptor_ring_parts, ignore_index=True)
    donors = pd.concat(gly_donor_parts, ignore_index=True)
    gly_distance = pd.concat(gly_distance_parts, ignore_index=True)
    gly_spatial_raw = pd.concat(gly_spatial_conformer_parts, ignore_index=True)
    aromatic_ring_background = pd.concat(
        aromatic_ring_background_parts, ignore_index=True)
    aromatic_ring_source = pd.concat(
        aromatic_ring_source_parts, ignore_index=True)
    asu_aromatic_ring_background = build_asu_aromatic_ring_background(
        aromatic_ring_source)
    symmetry_aromatic_ring_background = aromatic_ring_background[
        aromatic_ring_background["is_symmetry_contact"] == 1]
    aromatic_ring_background = pd.concat(
        [asu_aromatic_ring_background, symmetry_aromatic_ring_background],
        ignore_index=True,
    ).sort_values(
        "combined_occupancy", ascending=False
    ).drop_duplicates("aromatic_ring_pair_id")

    # Main contact composition and donor-class geometry.
    total_positive = len(positives)
    composition_rows: List[Dict[str, Any]] = []
    distance_rows: List[Dict[str, Any]] = []
    distance_arrays: Dict[str, Dict[str, np.ndarray]] = defaultdict(dict)
    for donor_class in DONOR_CLASS_ORDER:
        current = positives[positives["donor_class"] == donor_class]
        if current.empty:
            continue
        composition_rows.append({
            "donor_class": donor_class,
            "element": text(current["x_element"].iloc[0]),
            "positive_pair_count": len(current),
            "composition_percent": 100.0 * len(current) / total_positive,
            "structure_count": current["pdb_id"].nunique(),
        })
        for metric in ("dist_X_plane", "dist_X_centroid", "projection_distance", "angle_XPCN"):
            values = pd.to_numeric(current[metric], errors="coerce").dropna().to_numpy()
            distance_rows.append({
                "donor_class": donor_class, "metric": metric,
                "n_pairs": len(values), "n_structures": current["pdb_id"].nunique(),
                "median": float(np.median(values)), "q1": float(np.quantile(values, 0.25)),
                "q3": float(np.quantile(values, 0.75)),
            })
        for pdb_id, group in current.groupby("pdb_id"):
            distance_arrays[pdb_id][donor_class] = pd.to_numeric(
                group["dist_X_plane"], errors="coerce").dropna().to_numpy()

    # Acceptor composition and selected-acceptor propensity.
    positives["acceptor_type"] = [
        acceptor_type(residue, size)
        for residue, size in zip(positives["pi_residue_name"], positives["pi_ring_size"])
    ]
    acceptor_composition_rows = []
    for kind, group in positives.groupby("acceptor_type", dropna=False):
        acceptor_composition_rows.append({
            "acceptor_type": kind, "positive_pair_count": len(group),
            "composition_percent": 100.0 * len(group) / len(positives),
            "structure_count": group["pdb_id"].nunique(),
        })
    selected_acceptor_positive = positives[
        (positives["acceptor_on_selected_chain"] == 1)
        & positives["acceptor_type"].isin(ACCEPTOR_ORDER)
    ]
    participating_rings = selected_acceptor_positive.groupby("ring_site_id").size().rename("positive_pairs")
    rings["acceptor_type"] = [
        acceptor_type(residue, size) for residue, size in zip(rings["residue_name"], rings["ring_size"])
    ]
    rings = rings.join(participating_rings, on="ring_site_id")
    rings["positive_pairs"] = rings["positive_pairs"].fillna(0).astype(int)
    rings["participating"] = (rings["positive_pairs"] > 0).astype(int)
    acceptor_propensity_rows = []
    for kind in ACCEPTOR_ORDER:
        group = rings[rings["acceptor_type"] == kind]
        participating = int(group["participating"].sum())
        positive_pairs = int(group["positive_pairs"].sum())
        acceptor_propensity_rows.append({
            "acceptor_type": kind, "eligible_ring_count": len(group),
            "participating_ring_count": participating,
            "participating_ring_percent": 100.0 * participating / len(group) if len(group) else np.nan,
            "positive_pair_count": positive_pairs,
            "positive_pairs_per_1000_eligible_rings": 1000.0 * positive_pairs / len(group) if len(group) else np.nan,
            "positive_pairs_per_participating_ring": positive_pairs / participating if participating else np.nan,
        })

    # Donor x acceptor matrices at distance, ring-face and positive levels.
    matrix_rows = []
    distance_counts = distance.groupby(["donor_class", "acceptor_type"]).size()
    spatial_counts = spatial.groupby(["donor_class", "acceptor_type"]).size()
    positive_standard = positives[positives["acceptor_type"].isin(ACCEPTOR_ORDER)]
    positive_counts = positive_standard.groupby(["donor_class", "acceptor_type"]).size()
    distance_by_donor = distance_counts.groupby(level=0).sum()
    spatial_by_donor = spatial_counts.groupby(level=0).sum()
    positive_by_donor = positive_counts.groupby(level=0).sum()
    for donor_class in DONOR_CLASS_ORDER:
        for kind in ACCEPTOR_ORDER:
            key = (donor_class, kind)
            dc, sc, pc = int(distance_counts.get(key, 0)), int(spatial_counts.get(key, 0)), int(positive_counts.get(key, 0))
            donor_dc = int(distance_by_donor.get(donor_class, 0))
            donor_sc = int(spatial_by_donor.get(donor_class, 0))
            donor_pc = int(positive_by_donor.get(donor_class, 0))
            spatial_share = sc / donor_sc if donor_sc else np.nan
            positive_share = pc / donor_pc if donor_pc else np.nan
            enrichment = positive_share / spatial_share if spatial_share else np.nan
            matrix_rows.append({
                "donor_class": donor_class, "acceptor_type": kind,
                "element_cutoff_distance_pair_count": dc,
                "ring_face_spatial_pair_count": sc,
                "positive_pair_count": pc,
                "positive_fraction_of_ring_face": pc / sc if sc else np.nan,
                "donor_class_distance_pair_count": donor_dc,
                "donor_class_ring_face_pair_count": donor_sc,
                "donor_class_positive_pair_count": donor_pc,
                "positive_share_within_donor_class": positive_share,
                "ring_face_share_within_donor_class": spatial_share,
                "acceptor_enrichment_within_donor_class": enrichment,
                "log2_acceptor_enrichment_within_donor_class": (
                    math.log2(enrichment)
                    if enrichment and enrichment > 0 else np.nan
                ),
            })

    # Exact contact fraction, then ring-pair enrichment against an entirely
    # H-independent nearby aromatic-ring universe.
    aromatic_positive = positives[
        (positives["donor_class"] == "Aromatic C-H")
        & (positives["pi_pi_geometry_status"] == "ok")
        & (positives["aromatic_ring_pair_id"].astype(str) != "")
    ].copy()
    positive_ring_pairs = aromatic_positive.groupby(
        "aromatic_ring_pair_id", as_index=False).agg(
        pdb_id=("pdb_id", "first"),
        is_tshaped=("is_pi_pi_tshaped_spatial", "max"),
        is_symmetry_contact=("is_symmetry_contact", "max"),
        centroid_distance=("pi_pi_centroid_distance", "median"),
        positive_contact_count=("spatial_pair_id", "nunique"),
    )
    nearby_positive_ring_pairs = positive_ring_pairs[
        positive_ring_pairs["centroid_distance"].between(
            3.0, 5.5, inclusive="both")
        & (positive_ring_pairs["is_symmetry_contact"] == 0)]
    aromatic_ring_background_primary = aromatic_ring_background[
        aromatic_ring_background["is_symmetry_contact"] == 0]
    background_t_fraction = float(
        aromatic_ring_background_primary["is_pi_pi_tshaped_spatial"].mean())
    positive_near_t_fraction = float(
        nearby_positive_ring_pairs["is_tshaped"].mean())
    tshape_rows = [{
        "analysis_unit": "aromatic_CH_positive_contact",
        "background_count": np.nan, "background_tshaped_count": np.nan,
        "background_tshaped_fraction": np.nan,
        "positive_count": len(aromatic_positive),
        "positive_tshaped_count": int(
            aromatic_positive["is_pi_pi_tshaped_spatial"].sum()),
        "positive_tshaped_fraction": float(
            aromatic_positive["is_pi_pi_tshaped_spatial"].mean()),
        "tshaped_fraction_enrichment": np.nan,
        "scope": "exact master-set Aromatic C-H positives",
    }, {
        "analysis_unit": "unique_nearby_aromatic_ring_image_pair",
        "background_count": len(aromatic_ring_background_primary),
        "background_tshaped_count": int(
            aromatic_ring_background_primary[
                "is_pi_pi_tshaped_spatial"].sum()),
        "background_tshaped_fraction": background_t_fraction,
        "positive_count": len(nearby_positive_ring_pairs),
        "positive_tshaped_count": int(
            nearby_positive_ring_pairs["is_tshaped"].sum()),
        "positive_tshaped_fraction": positive_near_t_fraction,
        "tshaped_fraction_enrichment": (
            positive_near_t_fraction / background_t_fraction
            if background_t_fraction else np.nan),
        "scope": "3.0-5.5 A centroid distance; H-independent non-symmetry background; either endpoint on selected chain",
    }]

    # Gly/non-Gly four-stage funnel, reported as donor sites and physical pairs.
    stage_frames = {
        "eligible_backbone_NH": None,
        "within_N_distance_cutoff": gly_distance,
        "ring_face_spatial": collapse_pairs(gly_spatial_raw),
        "XH_direction_positive": collapse_pairs(gly_spatial_raw)[lambda x: x["is_positive_pair"] == 1],
    }
    funnel_rows = []
    groups = ("GLY", "NON_GLY", "ALL")
    for label in groups:
        if label == "GLY":
            donor_subset = donors[donors["residue_name"] == "GLY"]
        elif label == "NON_GLY":
            donor_subset = donors[donors["residue_name"] != "GLY"]
        else:
            donor_subset = donors
        eligible_sites = set(donor_subset["donor_site_id"])
        previous_pairs = None
        for stage, frame in stage_frames.items():
            if frame is None:
                pair_count, site_count = 0, len(eligible_sites)
            else:
                current = frame[frame["donor_site_id"].isin(eligible_sites)]
                pair_count = current["spatial_pair_id"].nunique()
                site_count = current["donor_site_id"].nunique()
            funnel_rows.append({
                "residue_group": label, "stage": stage,
                "eligible_donor_site_count": len(eligible_sites),
                "donor_sites_reaching_stage": site_count,
                "physical_pair_count": pair_count,
                "donor_site_reach_percent": 100.0 * site_count / len(eligible_sites) if eligible_sites else np.nan,
                "pairs_per_1000_eligible_sites": 1000.0 * pair_count / len(eligible_sites) if eligible_sites else np.nan,
                "pair_retention_from_previous_stage": pair_count / previous_pairs if previous_pairs else np.nan,
            })
            if frame is not None:
                previous_pairs = pair_count

    # Fixed-H threshold sensitivity. Any compatible conformer may pass.
    angle_rows = []
    angle_scenario_collapsed: Dict[str, pd.DataFrame] = {}
    primary_angle_mismatch = 0
    primary_mismatch_frame = pd.DataFrame()
    angle_geometry_ok = (
        gly_spatial_raw["backbone_nh_direction_status"] == "ok")
    angle_altloc_unambiguous = (
        gly_spatial_raw["pi_altloc"].fillna("").astype(str).eq("")
        & gly_spatial_raw["x_altloc"].fillna("").astype(str).eq(""))
    angle_source = gly_spatial_raw[
        angle_geometry_ok & angle_altloc_unambiguous].copy()
    angle_altloc_excluded_pair_count = gly_spatial_raw[
        angle_geometry_ok & ~angle_altloc_unambiguous
    ]["spatial_pair_id"].nunique()
    for scenario, theta_max, xh_pi_min in ANGLE_SCENARIOS:
        current = angle_source.copy()
        theta = pd.to_numeric(current["backbone_nh_theta_deg"], errors="coerce")
        xh_pi = pd.to_numeric(current["backbone_nh_xh_pi_angle_deg"], errors="coerce")
        current["sensitivity_positive"] = (
            ((current["hudson_spatial"] == 1) & theta.notna() & (theta <= theta_max))
            | ((current["plevin_spatial"] == 1) & xh_pi.notna() & (xh_pi >= xh_pi_min))
        ).astype(int)
        collapsed = current.groupby("spatial_pair_id", as_index=False).agg(
            pdb_id=("pdb_id", "first"), donor_site_id=("donor_site_id", "first"),
            x_residue_name=("x_residue_name", "first"),
            donor_treatment=("donor_treatment", "first"),
            is_symmetry_contact=("is_symmetry_contact", "max"),
            backbone_nh_hydrogen_count=("backbone_nh_hydrogen_count", "max"),
            min_theta_deg=("backbone_nh_theta_deg", "min"),
            max_xh_pi_angle_deg=("backbone_nh_xh_pi_angle_deg", "max"),
            sensitivity_positive=("sensitivity_positive", "max"),
            recorded_positive=("is_positive_pair", "max"),
        )
        angle_scenario_collapsed[scenario] = collapsed
        for residue_group, mask in (
            ("GLY", collapsed["x_residue_name"] == "GLY"),
            ("NON_GLY", collapsed["x_residue_name"] != "GLY"),
            ("ALL", pd.Series(True, index=collapsed.index)),
        ):
            group = collapsed[mask & (collapsed["sensitivity_positive"] == 1)]
            if residue_group == "GLY":
                eligible_count = int((donors["residue_name"] == "GLY").sum())
            elif residue_group == "NON_GLY":
                eligible_count = int((donors["residue_name"] != "GLY").sum())
            else:
                eligible_count = len(donors)
            angle_rows.append({
                "scenario": scenario, "hudson_theta_max_deg": theta_max,
                "plevin_XH_pi_min_deg": xh_pi_min, "residue_group": residue_group,
                "positive_pair_count": len(group),
                "eligible_donor_site_count": eligible_count,
                "positive_pairs_per_1000_eligible_sites": (
                    1000.0 * len(group) / eligible_count
                    if eligible_count else np.nan),
                "participating_donor_site_count": group["donor_site_id"].nunique(),
                "structure_count": group["pdb_id"].nunique(),
                "scope": "fixed Backbone N-H; blank donor and ring altloc",
            })
        if scenario == "primary":
            primary_mismatch_frame = collapsed[
                collapsed["sensitivity_positive"]
                != collapsed["recorded_positive"]
            ].copy()
            mismatch = len(primary_mismatch_frame)
            primary_angle_mismatch = mismatch

    # PDB-cluster bootstrap tables.
    rng = np.random.default_rng(args.seed)
    bootstrap_rows: List[Dict[str, Any]] = []
    donor_counts = positives.groupby(["pdb_id", "donor_class"]).size().unstack(fill_value=0)
    donor_counts["total"] = donor_counts.sum(axis=1)
    donor_counts = donor_counts.reset_index()
    for donor_class in DONOR_CLASS_ORDER:
        if donor_class not in donor_counts:
            continue
        estimate, low, high = ratio_bootstrap(
            donor_counts, donor_class, "total", pdb_ids, rng, args.bootstrap)
        bootstrap_rows.append({
            "analysis": "donor_class_composition", "contrast": donor_class,
            "estimate": estimate, "ci95_low": low, "ci95_high": high,
            "unit": "fraction", "independent_unit": "PDB structure",
        })
    ring_pdb = rings.groupby(["pdb_id", "acceptor_type"]).agg(
        eligible=("ring_site_id", "nunique"), participating=("participating", "sum")
    ).reset_index()
    for kind in ACCEPTOR_ORDER:
        table = ring_pdb[ring_pdb["acceptor_type"] == kind][["pdb_id", "eligible", "participating"]]
        estimate, low, high = ratio_bootstrap(table, "participating", "eligible", pdb_ids, rng, args.bootstrap)
        bootstrap_rows.append({
            "analysis": "acceptor_ring_participation", "contrast": kind,
            "estimate": estimate, "ci95_low": low, "ci95_high": high,
            "unit": "fraction", "independent_unit": "PDB structure",
        })
    acceptor_counts = positives.groupby(
        ["pdb_id", "acceptor_type"]).size().unstack(fill_value=0)
    acceptor_counts["total"] = acceptor_counts.sum(axis=1)
    acceptor_counts = acceptor_counts.reset_index()
    for kind in list(ACCEPTOR_ORDER) + ["Other"]:
        if kind not in acceptor_counts:
            continue
        estimate, low, high = ratio_bootstrap(
            acceptor_counts, kind, "total", pdb_ids, rng, args.bootstrap)
        bootstrap_rows.append({
            "analysis": "acceptor_contact_composition", "contrast": kind,
            "estimate": estimate, "ci95_low": low, "ci95_high": high,
            "unit": "fraction", "independent_unit": "PDB structure",
        })
    contact_t = aromatic_positive.groupby("pdb_id").agg(
        total=("spatial_pair_id", "size"),
        tshaped=("is_pi_pi_tshaped_spatial", "sum"),
    ).reset_index()
    estimate, low, high = ratio_bootstrap(
        contact_t, "tshaped", "total", pdb_ids, rng, args.bootstrap)
    bootstrap_rows.append({
        "analysis": "tshaped_positive_contact_fraction",
        "contrast": "Aromatic C-H contacts", "estimate": estimate,
        "ci95_low": low, "ci95_high": high, "unit": "fraction",
        "independent_unit": "PDB structure",
    })
    ring_background_pdb = aromatic_ring_background_primary.groupby("pdb_id").agg(
        background_total=("aromatic_ring_pair_id", "size"),
        background_tshaped=("is_pi_pi_tshaped_spatial", "sum"),
    )
    positive_ring_pdb = nearby_positive_ring_pairs.groupby("pdb_id").agg(
        positive_total=("aromatic_ring_pair_id", "size"),
        positive_tshaped=("is_tshaped", "sum"),
    )
    ring_boot = ring_background_pdb.join(
        positive_ring_pdb, how="outer").fillna(0).reset_index()
    estimate, low, high = ratio_of_ratios_bootstrap(
        ring_boot, "positive_tshaped", "positive_total",
        "background_tshaped", "background_total",
        pdb_ids, rng, args.bootstrap,
    )
    bootstrap_rows.append({
        "analysis": "tshaped_fraction_enrichment",
        "contrast": "positive / H-independent nearby ring pairs",
        "estimate": estimate, "ci95_low": low, "ci95_high": high,
        "unit": "risk ratio", "independent_unit": "PDB structure",
    })
    for group_a, group_b in DISTANCE_CONTRASTS:
        ma, mb, estimate, low, high = median_cluster_bootstrap(
            distance_arrays, group_a, group_b, pdb_ids, rng, args.bootstrap)
        bootstrap_rows.append({
            "analysis": "median_dist_X_plane_difference",
            "contrast": f"{group_a} minus {group_b}", "group_a_median": ma,
            "group_b_median": mb, "estimate": estimate,
            "ci95_low": low, "ci95_high": high, "unit": "angstrom",
            "independent_unit": "PDB structure",
        })
    # Gly/non-Gly pair-rate ratio.
    donor_pdb = donors.assign(group=np.where(donors["residue_name"] == "GLY", "GLY", "NON_GLY")).groupby(["pdb_id", "group"]).size().unstack(fill_value=0)
    positive_nh = collapse_pairs(gly_spatial_raw)
    positive_nh = positive_nh[positive_nh["is_positive_pair"] == 1].assign(group=lambda x: np.where(x["x_residue_name"] == "GLY", "GLY", "NON_GLY"))
    positive_pdb = positive_nh.groupby(["pdb_id", "group"]).size().unstack(fill_value=0)
    gly_boot = pd.DataFrame(index=pdb_ids)
    for name, source in (("eligible", donor_pdb), ("positive", positive_pdb)):
        gly_boot[f"gly_{name}"] = source.get("GLY", pd.Series(dtype=float)).reindex(pdb_ids, fill_value=0)
        gly_boot[f"nongly_{name}"] = source.get("NON_GLY", pd.Series(dtype=float)).reindex(pdb_ids, fill_value=0)
    gly_boot = gly_boot.reset_index(names="pdb_id")
    estimate, low, high = ratio_of_ratios_bootstrap(
        gly_boot, "gly_positive", "gly_eligible", "nongly_positive", "nongly_eligible",
        pdb_ids, rng, args.bootstrap)
    bootstrap_rows.append({
        "analysis": "backbone_NH_pair_rate_ratio", "contrast": "GLY / NON_GLY",
        "estimate": estimate, "ci95_low": low, "ci95_high": high,
        "unit": "rate ratio", "independent_unit": "PDB structure",
    })
    eligible_by_pdb = donors.assign(
        group=np.where(
            donors["residue_name"] == "GLY", "GLY", "NON_GLY")
    ).groupby(["pdb_id", "group"]).size().unstack(fill_value=0)
    for scenario, collapsed in angle_scenario_collapsed.items():
        positive = collapsed[
            collapsed["sensitivity_positive"] == 1
        ].assign(group=lambda x: np.where(
            x["x_residue_name"] == "GLY", "GLY", "NON_GLY"))
        positive_by_pdb = positive.groupby(
            ["pdb_id", "group"]).size().unstack(fill_value=0)
        table = pd.DataFrame(index=pdb_ids)
        for group, prefix in (("GLY", "gly"), ("NON_GLY", "nongly")):
            table[f"{prefix}_eligible"] = eligible_by_pdb.get(
                group, pd.Series(dtype=float)).reindex(pdb_ids, fill_value=0)
            table[f"{prefix}_positive"] = positive_by_pdb.get(
                group, pd.Series(dtype=float)).reindex(pdb_ids, fill_value=0)
        table = table.reset_index(names="pdb_id")
        estimate, low, high = ratio_of_ratios_bootstrap(
            table, "gly_positive", "gly_eligible",
            "nongly_positive", "nongly_eligible",
            pdb_ids, rng, args.bootstrap)
        bootstrap_rows.append({
            "analysis": "backbone_NH_angle_sensitivity_rate_ratio",
            "contrast": f"{scenario}: GLY / NON_GLY",
            "estimate": estimate, "ci95_low": low, "ci95_high": high,
            "unit": "rate ratio", "independent_unit": "PDB structure",
        })

    # Donor-conditional acceptor enrichment.  This removes the large baseline
    # differences in directional realization among donor classes and asks only
    # whether an acceptor is over- or under-represented relative to the
    # ring-face opportunities available to the same donor class.
    matrix_intervals = conditional_enrichment_bootstrap(
        positive_standard, spatial, DONOR_CLASS_ORDER, ACCEPTOR_ORDER,
        pdb_ids, rng, args.bootstrap,
    )
    for row in matrix_rows:
        key = (row["donor_class"], row["acceptor_type"])
        estimate, low, high = matrix_intervals[key]
        row_estimate = row["acceptor_enrichment_within_donor_class"]
        both_missing = not math.isfinite(estimate) and not math.isfinite(
            row_estimate
        )
        if not both_missing and not math.isclose(
            estimate, row_estimate, rel_tol=1e-12, abs_tol=1e-12
        ):
            raise AssertionError(f"Conditional enrichment mismatch for {key}")
        row["conditional_enrichment_ci95_low"] = low
        row["conditional_enrichment_ci95_high"] = high
        row["independent_unit"] = "PDB structure"
        bootstrap_rows.append({
            "analysis": "donor_conditional_acceptor_enrichment",
            "contrast": f"{row['donor_class']} | {row['acceptor_type']}",
            "estimate": estimate,
            "ci95_low": low,
            "ci95_high": high,
            "unit": "fold versus donor-matched ring-face acceptor distribution",
            "independent_unit": "PDB structure",
        })

    global_structure.to_csv(
        output / "global_structure_summary.csv", index=False)
    write_csv(output / "global_summary.csv", global_summary_rows)
    write_csv(output / "main_contact_composition.csv", composition_rows)
    write_csv(output / "acceptor_composition.csv", acceptor_composition_rows)
    write_csv(output / "acceptor_propensity.csv", acceptor_propensity_rows)
    write_csv(output / "donor_acceptor_matrix.csv", matrix_rows)
    write_csv(output / "tshape_summary.csv", tshape_rows)
    write_csv(output / "gly_backbone_nh_funnel.csv", funnel_rows)
    write_csv(output / "backbone_nh_angle_sensitivity.csv", angle_rows)
    primary_mismatch_frame.to_csv(
        output / "backbone_nh_primary_mismatches.csv", index=False)
    write_csv(output / "distance_by_donor_class.csv", distance_rows)
    write_csv(output / "bootstrap_estimates.csv", bootstrap_rows)

    expected_main = 173732 if args.max_structures is None else None
    aromatic_master_positive = positives[
        positives["donor_class"] == "Aromatic C-H"]
    aromatic_positive_geometry_resolved = aromatic_master_positive[
        aromatic_master_positive["pi_pi_geometry_status"] == "ok"]
    aromatic_spatial_candidate = spatial[
        spatial["donor_class"] == "Aromatic C-H"]
    aromatic_spatial_geometry_resolved = aromatic_spatial_candidate[
        aromatic_spatial_candidate["pi_pi_geometry_status"] == "ok"]
    unresolved_positive_aromatic = (
        len(aromatic_master_positive) - len(aromatic_positive_geometry_resolved))
    unresolved_spatial_aromatic = (
        len(aromatic_spatial_candidate) - len(aromatic_spatial_geometry_resolved))
    missing_positive_ring_pairs = sorted(
        set(nearby_positive_ring_pairs["aromatic_ring_pair_id"])
        - set(aromatic_ring_background["aromatic_ring_pair_id"])
    )
    qc = {
        "status": "pass",
        "background": str(root), "schema_version": metadata["output_schema_version"],
        "structures_analyzed": len(pdb_ids), "main_positive_pair_count": total_positive,
        "full_model_unique_positive_pair_count": unique_positive_total,
        "full_model_raw_xpid_positive_row_count": raw_positive_total,
        "full_model_structures_with_positive_pair": structures_with_positive,
        "expected_full_main_positive_pair_count": expected_main,
        "main_positive_count_matches_expected": expected_main is None or total_positive == expected_main,
        "primary_angle_reclassification_mismatch_count": primary_angle_mismatch,
        "backbone_NH_angle_sensitivity_altloc_ambiguous_pairs_excluded": int(
            angle_altloc_excluded_pair_count),
        "aromatic_CH_positive_count_in_master": len(aromatic_master_positive),
        "aromatic_CH_positive_geometry_resolved": len(aromatic_positive_geometry_resolved),
        "aromatic_CH_positive_geometry_unresolved": unresolved_positive_aromatic,
        "aromatic_CH_spatial_candidate_count": len(aromatic_spatial_candidate),
        "aromatic_CH_spatial_geometry_resolved": len(aromatic_spatial_geometry_resolved),
        "aromatic_CH_spatial_geometry_unresolved": unresolved_spatial_aromatic,
        "nearby_positive_ring_pairs_missing_from_H_independent_background": len(
            missing_positive_ring_pairs),
        "missing_positive_ring_pair_examples": missing_positive_ring_pairs[:20],
        "absent_zero_row_partition_counts": dict(absent_partition_counts),
        "notes": [
            "Global density uses unique physical pairs per full-model standard-amino-acid residue position; raw XPID rows are reported separately.",
            "Contact composition uses the either-endpoint representative-chain master set.",
            "Donor and acceptor propensities use endpoint-matched selected-chain risk sets.",
            "Donor-acceptor enrichment is conditioned within donor class against its H-independent ring-face acceptor distribution; it is not an affinity estimate.",
            "The exact Aromatic C-H contact fraction includes symmetry contacts; T-shape enrichment uses the fully matched non-symmetry ring-pair universe.",
            "Backbone N-H angular-threshold sensitivity is restricted to blank-altloc donor/ring pairs; the all-contact Gly funnel is unchanged.",
            "Bootstrap independent unit is the PDB structure; no contact-level p values are reported.",
        ],
    }
    if expected_main is not None and total_positive != expected_main:
        qc["status"] = "fail"
    if args.max_structures is None and unique_positive_total != 280843:
        qc["status"] = "fail"
    if args.max_structures is None and structures_with_positive != 4644:
        qc["status"] = "fail"
    if primary_angle_mismatch != 0:
        qc["status"] = "fail"
    if unresolved_positive_aromatic != 0 or unresolved_spatial_aromatic != 0:
        qc["status"] = "fail"
    if missing_positive_ring_pairs:
        qc["status"] = "fail"
    (output / "qc_report.json").write_text(json.dumps(qc, indent=2) + "\n")
    analysis_metadata = {
        "tool": "xpid-section02-chemical-landscape",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0], "platform": platform.platform(),
        "bootstrap_replicates": args.bootstrap, "random_seed": args.seed,
        "independent_unit": "PDB structure",
        "master_contact_scope": "primary neutral standard donor; either endpoint on PISCES-selected chain",
        "donor_propensity_scope": "donor on PISCES-selected chain",
        "acceptor_propensity_scope": "acceptor on PISCES-selected chain",
        "input_metadata_sha256": sha256(metadata_path),
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_directory": str(output),
    }
    (output / "analysis_metadata.json").write_text(json.dumps(analysis_metadata, indent=2) + "\n")
    print(json.dumps(qc, indent=2))
    print(f"Output: {output}")
    return 0 if qc["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
