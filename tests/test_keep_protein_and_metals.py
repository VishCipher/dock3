"""
Tests for keep_protein_and_metals.py -- cleans a receptor PDB down to
protein atoms + functional metal cofactors only.

Regression coverage: a real 3BH4 run failed AutoGrid4 with "Unknown
receptor type: Na" because a sodium ion (crystallization artifact) reached
the grid-generation step uncleaned. These tests lock in that sodium and
waters never survive cleaning, while a real functional metal always does.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

from keep_protein_and_metals import clean_receptor  # noqa: E402

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "mini_receptor_with_metal.pdb")


def test_protein_atom_lines_are_preserved(tmp_path):
    out_path = tmp_path / "cleaned.pdb"
    clean_receptor(FIXTURE_PATH, str(out_path))
    content = out_path.read_text()
    assert "ALA A   1" in content


def test_calcium_hetatm_is_kept(tmp_path):
    out_path = tmp_path / "cleaned.pdb"
    clean_receptor(FIXTURE_PATH, str(out_path))
    content = out_path.read_text()
    assert any(line.startswith("HETATM") and " CA " in line for line in content.splitlines())


def test_sodium_hetatm_is_removed(tmp_path):
    """The exact real bug this project hit: sodium must never reach downstream tools."""
    out_path = tmp_path / "cleaned.pdb"
    clean_receptor(FIXTURE_PATH, str(out_path))
    content = out_path.read_text()
    assert " NA " not in content
    assert "NA A 102" not in content


def test_water_is_removed(tmp_path):
    out_path = tmp_path / "cleaned.pdb"
    clean_receptor(FIXTURE_PATH, str(out_path))
    content = out_path.read_text()
    assert "HOH" not in content


def test_returns_correct_kept_count(tmp_path):
    out_path = tmp_path / "cleaned.pdb"
    kept = clean_receptor(FIXTURE_PATH, str(out_path))
    assert kept == 1  # only the calcium -- sodium and water don't count
