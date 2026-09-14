"""
Tests for detect_cofactors.py -- the checkpoint that scans a receptor PDB
for functional metal cofactors.

Regression coverage: this pipeline's real 3BH4 run was broken for a long
time because sodium (a crystallization artifact) was originally included
in the metal whitelist and got treated as a genuine cofactor. These tests
lock in that sodium/potassium must NEVER be detected as a metal, and that
genuine functional metals (calcium, in this fixture) must always be found.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

from detect_cofactors import scan_receptor, METAL_CODES  # noqa: E402

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "mini_receptor_with_metal.pdb")


def test_sodium_and_potassium_are_never_in_the_whitelist():
    """Regression test: alkali ions caused a real bug (AutoGrid4 crash on
    an unsupported 'Na' atom type) when they were treated as cofactors."""
    assert "NA" not in METAL_CODES
    assert "K" not in METAL_CODES


def test_calcium_is_detected_in_fixture():
    metals = scan_receptor(FIXTURE_PATH)
    symbols_found = [m["symbol"] for m in metals]
    assert "Ca" in symbols_found


def test_sodium_in_fixture_is_not_detected():
    """The fixture deliberately includes a sodium HETATM -- it must be ignored."""
    metals = scan_receptor(FIXTURE_PATH)
    symbols_found = [m["symbol"] for m in metals]
    assert "Na" not in symbols_found


def test_water_is_never_detected_as_a_metal():
    metals = scan_receptor(FIXTURE_PATH)
    resnames_found = [m["resname"] for m in metals]
    assert "HOH" not in resnames_found


def test_detected_calcium_has_correct_coordinates():
    metals = scan_receptor(FIXTURE_PATH)
    calcium = next(m for m in metals if m["symbol"] == "Ca")
    assert calcium["x"] == 15.000
    assert calcium["y"] == 10.000
    assert calcium["z"] == 10.000
    assert calcium["chain"] == "A"
    assert calcium["resnum"] == "101"


def test_only_one_metal_found_in_fixture():
    """Exactly one real metal (calcium) should be found -- sodium and water excluded."""
    metals = scan_receptor(FIXTURE_PATH)
    assert len(metals) == 1
