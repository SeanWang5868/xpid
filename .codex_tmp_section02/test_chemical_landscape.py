#!/usr/bin/env python3
"""Small deterministic tests for the Section-02 landscape analysis."""
from __future__ import annotations

import unittest
from pathlib import Path
import tempfile

import gemmi
import numpy as np
import pandas as pd

from analyze_chemical_landscape import (
    acceptor_type,
    build_asu_aromatic_ring_background,
    collapse_pairs,
    conditional_enrichment_bootstrap,
    ratio_bootstrap,
    ratio_of_ratios_bootstrap,
    read_columns_or_empty,
)
from background_core import (
    _canonical_ring_pair_id,
    _transform_image_positions,
)
from analyze_residue_donor import cluster_ratio_interval


class AcceptorTypeTests(unittest.TestCase):
    def test_trp_rings_remain_distinct(self):
        self.assertEqual(acceptor_type("TRP", 6), "TRP-6")
        self.assertEqual(acceptor_type("TRP", 5), "TRP-5")
        self.assertEqual(acceptor_type("FAD", 6), "Other")


class PairCollapseTests(unittest.TestCase):
    def test_existence_flags_survive_lower_occupancy_conformer(self):
        rows = pd.DataFrame([
            {
                "spatial_pair_id": "p1", "combined_occupancy": 1.0,
                "is_positive_pair": 0, "ring_face_spatial": 1,
                "pi_pi_geometry_status": "not_aromatic_ch",
                "backbone_nh_direction_status": "not_backbone_nh",
            },
            {
                "spatial_pair_id": "p1", "combined_occupancy": 0.5,
                "is_positive_pair": 1, "ring_face_spatial": 1,
                "pi_pi_geometry_status": "not_aromatic_ch",
                "backbone_nh_direction_status": "not_backbone_nh",
            },
        ])
        collapsed = collapse_pairs(rows)
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(int(collapsed.iloc[0]["is_positive_pair"]), 1)

    def test_reconstructed_geometry_is_preferred(self):
        rows = pd.DataFrame([
            {
                "spatial_pair_id": "p1", "combined_occupancy": 1.0,
                "is_positive_pair": 0, "ring_face_spatial": 1,
                "pi_pi_geometry_status": "donor_ring_not_resolved",
                "pi_pi_centroid_distance": np.nan,
                "backbone_nh_direction_status": "not_backbone_nh",
            },
            {
                "spatial_pair_id": "p1", "combined_occupancy": 0.6,
                "is_positive_pair": 0, "ring_face_spatial": 1,
                "pi_pi_geometry_status": "ok",
                "pi_pi_centroid_distance": 4.7,
                "backbone_nh_direction_status": "not_backbone_nh",
            },
        ])
        collapsed = collapse_pairs(rows)
        self.assertEqual(collapsed.iloc[0]["pi_pi_geometry_status"], "ok")
        self.assertAlmostEqual(float(collapsed.iloc[0]["pi_pi_centroid_distance"]), 4.7)


class BootstrapTests(unittest.TestCase):
    def test_residue_rate_interval_uses_pdb_weights(self):
        numerator = np.array([1.0, 2.0, 3.0])
        denominator = np.array([2.0, 4.0, 6.0])
        weights = np.random.default_rng(4).multinomial(
            3, np.full(3, 1.0 / 3.0), size=100
        ).astype(float)
        estimate, low, high = cluster_ratio_interval(
            numerator, denominator, weights, 1000.0
        )
        self.assertAlmostEqual(estimate, 500.0)
        self.assertAlmostEqual(low, 500.0)
        self.assertAlmostEqual(high, 500.0)

    def test_ratio_uses_pdb_clusters(self):
        table = pd.DataFrame({
            "pdb_id": ["a", "b", "c"],
            "num": [1, 2, 3], "den": [2, 4, 6],
        })
        estimate, low, high = ratio_bootstrap(
            table, "num", "den", ["a", "b", "c"],
            np.random.default_rng(1), 100,
        )
        self.assertAlmostEqual(estimate, 0.5)
        self.assertAlmostEqual(low, 0.5)
        self.assertAlmostEqual(high, 0.5)

    def test_ratio_of_ratios(self):
        table = pd.DataFrame({
            "pdb_id": ["a", "b"],
            "an": [2, 2], "ad": [4, 4],
            "bn": [1, 1], "bd": [4, 4],
        })
        estimate, _, _ = ratio_of_ratios_bootstrap(
            table, "an", "ad", "bn", "bd", ["a", "b"],
            np.random.default_rng(2), 20,
        )
        self.assertAlmostEqual(estimate, 2.0)

    def test_acceptor_enrichment_is_conditioned_within_donor(self):
        positive = pd.DataFrame([
            {"pdb_id": "a", "donor_class": "D", "acceptor_type": "A"},
            {"pdb_id": "a", "donor_class": "D", "acceptor_type": "A"},
            {"pdb_id": "b", "donor_class": "D", "acceptor_type": "B"},
        ])
        spatial = pd.DataFrame([
            {"pdb_id": "a", "donor_class": "D", "acceptor_type": "A"},
            {"pdb_id": "b", "donor_class": "D", "acceptor_type": "B"},
            {"pdb_id": "b", "donor_class": "D", "acceptor_type": "B"},
        ])
        result = conditional_enrichment_bootstrap(
            positive, spatial, ["D"], ["A", "B"], ["a", "b"],
            np.random.default_rng(3), 100,
        )
        self.assertAlmostEqual(result[("D", "A")][0], 2.0)
        self.assertAlmostEqual(result[("D", "B")][0], 0.5)


class SchemaBoundaryTests(unittest.TestCase):
    def test_asu_ring_background_is_hydrogen_independent(self):
        common = {
            "pdb_id": "test", "eligible": 1, "altloc": "",
            "mean_occupancy": 1.0, "is_selected_chain": 1,
            "ring_normal_x": 0.0, "ring_normal_y": 0.0,
            "ring_normal_z": 1.0,
        }
        rings = pd.DataFrame([
            dict(common, ring_site_id="A", residue_name="PHE",
                 is_standard_protein_acceptor=1,
                 ring_center_x=0.0, ring_center_y=0.0, ring_center_z=0.0),
            dict(common, ring_site_id="B", residue_name="TYR",
                 is_standard_protein_acceptor=1,
                 ring_center_x=0.0, ring_center_y=0.0, ring_center_z=4.5,
                 ring_normal_x=1.0, ring_normal_z=0.0),
        ])
        result = build_asu_aromatic_ring_background(rings)
        self.assertEqual(len(result), 1)
        self.assertEqual(int(result.iloc[0]["is_pi_pi_tshaped_spatial"]), 1)

    def test_missing_zero_row_partition_is_allowed(self):
        with tempfile.TemporaryDirectory() as directory:
            frame = read_columns_or_empty(
                Path(directory) / "absent.parquet", ("a", "b"))
        self.assertTrue(frame.empty)
        self.assertEqual(list(frame.columns), ["a", "b"])

    def test_complete_group_uses_marked_symmetry_image(self):
        structure = gemmi.Structure()
        structure.cell = gemmi.UnitCell(10, 10, 10, 90, 90, 90)
        structure.spacegroup_hm = "P 21 21 21"
        structure.setup_cell_images()
        original = gemmi.Position(1.0, 2.0, 3.0)
        reference = gemmi.Position(8.0, 8.0, 8.0)
        marked = structure.cell.find_nearest_pbc_position(
            reference, original, 1)
        second = gemmi.Position(1.5, 2.0, 3.0)
        transformed = _transform_image_positions(
            structure.cell, 1, original, marked, [original, second])
        self.assertIsNotNone(transformed)
        self.assertLess(np.linalg.norm(
            transformed[0] - np.array([marked.x, marked.y, marked.z])), 1e-6)
        self.assertAlmostEqual(
            np.linalg.norm(transformed[1] - transformed[0]), 0.5, places=6)

    def test_identity_operation_retains_lattice_translation(self):
        cell = gemmi.UnitCell(10, 10, 10, 90, 90, 90)
        original = gemmi.Position(1.0, 2.0, 3.0)
        marked = gemmi.Position(11.0, 2.0, 3.0)
        transformed = _transform_image_positions(
            cell, 0, original, marked,
            [original, gemmi.Position(1.5, 2.0, 3.0)])
        self.assertIsNotNone(transformed)
        self.assertTrue(np.allclose(transformed[0], [11.0, 2.0, 3.0]))
        self.assertTrue(np.allclose(transformed[1], [11.5, 2.0, 3.0]))

    def test_ring_pair_id_is_invariant_to_contact_direction(self):
        rotation = np.diag([-1.0, -1.0, 1.0])
        vector = np.array([0.5, 0.0, 0.5])
        forward = _canonical_ring_pair_id("A", "B", rotation, vector)
        inverse_rotation = np.linalg.inv(rotation)
        inverse_vector = -(inverse_rotation @ vector)
        backward = _canonical_ring_pair_id(
            "B", "A", inverse_rotation, inverse_vector)
        self.assertEqual(forward, backward)


if __name__ == "__main__":
    unittest.main(verbosity=2)
