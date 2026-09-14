"""
Tests for compute_admet.py -- RDKit-based drug-likeness descriptors.

Uses a real aspirin structure (generated and verified via RDKit itself
from its canonical SMILES, not hand-drawn coordinates) as a fixture with
well-known, checkable real-world values.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

from compute_admet import compute_properties  # noqa: E402

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "aspirin.sdf")


def test_aspirin_molecular_weight_matches_known_value():
    props = compute_properties(FIXTURE_PATH)
    assert props is not None
    # Aspirin's real MW is 180.16 g/mol -- allow tiny float tolerance
    assert abs(props["mw"] - 180.16) < 0.5


def test_aspirin_has_one_hydrogen_bond_donor():
    """Aspirin has exactly one H-bond donor: the carboxylic acid OH."""
    props = compute_properties(FIXTURE_PATH)
    assert props["hbd"] == 1


def test_aspirin_passes_lipinski_rule_of_five():
    """Aspirin is a textbook drug-like small molecule -- should have zero violations."""
    props = compute_properties(FIXTURE_PATH)
    assert props["lipinski_violations"] == 0
    assert props["violation_flags"] == "none"


def test_returns_none_for_unparseable_file(tmp_path):
    bad_file = tmp_path / "not_a_molecule.sdf"
    bad_file.write_text("this is not valid SDF content at all")
    assert compute_properties(str(bad_file)) is None
