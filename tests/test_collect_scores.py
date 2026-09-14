"""
Tests for collect_scores.py's parsing functions (Vina result tables).

collect_smina_scores.py shares this exact parsing logic (smina is a Vina
fork with an identical log table format), so these tests implicitly cover
both -- see that script's own docstring for the explicit note.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

from collect_scores import parse_best_affinity, ligand_id_from_log_path  # noqa: E402

VINA_LOG_CONTENT = """\
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -11.5      0.000      0.000
   2       -10.7      1.395      4.253
   3       -10.6      2.025     11.047
"""

EMPTY_LOG_CONTENT = "ERROR: something went wrong before any poses were generated\n"


def test_parses_best_affinity_from_mode_1(tmp_path):
    log_path = tmp_path / "EDTA_log.txt"
    log_path.write_text(VINA_LOG_CONTENT)
    assert parse_best_affinity(str(log_path)) == -11.5


def test_ignores_later_modes_even_if_more_negative(tmp_path):
    """Mode 1 is always Vina's best pose -- must never accidentally pick mode 2/3."""
    content = VINA_LOG_CONTENT.replace("2       -10.7", "2       -99.9")
    log_path = tmp_path / "test_log.txt"
    log_path.write_text(content)
    assert parse_best_affinity(str(log_path)) == -11.5


def test_returns_none_for_a_failed_run(tmp_path):
    log_path = tmp_path / "failed_log.txt"
    log_path.write_text(EMPTY_LOG_CONTENT)
    assert parse_best_affinity(str(log_path)) is None


def test_ligand_id_extracted_from_log_filename():
    assert ligand_id_from_log_path("results/3BH4/docking/vina/EDTA_log.txt") == "EDTA"


def test_ligand_id_handles_underscored_names():
    assert ligand_id_from_log_path("/some/path/citric_acid_log.txt") == "citric_acid"
