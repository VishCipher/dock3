"""
Tests for validate_metal_geometry.py -- the independent coordination-
geometry check computed from raw atomic coordinates.

Regression coverage: a real 3BH4 run showed ascorbic_acid at a genuine
2.15 A contact distance getting incorrectly flagged as zero coordinating
atoms, because the calcium lower distance bound was set to 2.2 A. These
tests lock in the corrected 2.0 A bound and the model-1-only pose reading.
"""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

from validate_metal_geometry import (  # noqa: E402
    distance,
    evaluate_ligand,
    load_top_pose_atoms,
    IDEAL_DISTANCE_RANGES,
)

FIXTURE_PDBQT = os.path.join(os.path.dirname(__file__), "fixtures", "mini_ligand_pose.pdbqt")


def test_distance_is_correct_euclidean():
    assert distance((0, 0, 0), (3, 4, 0)) == 5.0


def test_distance_is_symmetric():
    p1, p2 = (1.5, 2.5, 3.5), (4.0, 1.0, 0.0)
    assert distance(p1, p2) == distance(p2, p1)


def test_only_model_1_atoms_are_read():
    """The fixture has a model 2 atom at (50, 50, 50) -- it must be ignored."""
    atoms = load_top_pose_atoms(FIXTURE_PDBQT)
    positions = [pos for _, pos in atoms]
    assert (50.0, 50.0, 50.0) not in positions


def test_non_coordinating_atom_types_are_excluded():
    """Carbon and hydrogen (HD) should never count as coordinating atoms."""
    atoms = load_top_pose_atoms(FIXTURE_PDBQT)
    atom_types_found = [t for t, _ in atoms]
    assert "C" not in atom_types_found
    assert "HD" not in atom_types_found


def test_oxygen_and_nitrogen_atoms_are_included():
    atoms = load_top_pose_atoms(FIXTURE_PDBQT)
    atom_types_found = [t for t, _ in atoms]
    assert "OA" in atom_types_found
    assert "N" in atom_types_found


def test_calcium_lower_bound_is_2_0_not_2_2():
    """
    Direct regression test for the real bug: a 2.15 A contact was being
    missed because the lower bound was 2.2. Must stay at 2.0 or lower.
    """
    lower_bound, _ = IDEAL_DISTANCE_RANGES["Ca"]
    assert lower_bound <= 2.0


def test_a_2_15_angstrom_contact_counts_as_coordinating():
    """
    The exact real-world case: an oxygen at 2.15 A from calcium must be
    counted, not silently dropped by an overly strict distance cutoff.
    """
    metal_pos = (0.0, 0.0, 0.0)
    result = evaluate_ligand("Ca", metal_pos, FIXTURE_PDBQT)
    assert result is not None
    assert result["coordination_number"] >= 1
    assert result["closest_atom_distance"] == 2.15


def test_far_away_atom_is_excluded_from_coordination_count():
    """
    The fixture's N1 atom sits at 7.2 A -- far outside any calcium
    coordination range. Only the two oxygens (2.15 A, 2.5 A) should count,
    giving coordination_number == 2, not 3.
    """
    metal_pos = (0.0, 0.0, 0.0)
    result = evaluate_ligand("Ca", metal_pos, FIXTURE_PDBQT)
    assert result["coordination_number"] == 2


def test_flag_is_low_coordination_when_below_typical_range():
    """This synthetic fixture only has 2 coordinating atoms -- well below Ca's typical 6-8."""
    metal_pos = (0.0, 0.0, 0.0)
    result = evaluate_ligand("Ca", metal_pos, FIXTURE_PDBQT)
    assert result["flag"] == "LOW_COORDINATION"


def test_unknown_metal_falls_back_to_default_range():
    """A metal symbol not in IDEAL_DISTANCE_RANGES should use the generic default, not crash."""
    metal_pos = (0.0, 0.0, 0.0)
    result = evaluate_ligand("Xx", metal_pos, FIXTURE_PDBQT)
    assert result is not None  # must not raise/crash on an unrecognized metal
