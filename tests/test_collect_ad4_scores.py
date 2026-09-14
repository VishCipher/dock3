"""
Tests for collect_ad4_scores.py -- parses AutoDock4's .dlg output, which
repeats a binding-energy line once per docking run in the file.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))

from collect_ad4_scores import parse_best_energy, ligand_id_from_dlg_path  # noqa: E402

DLG_CONTENT = """\
DOCKED: USER    Run = 1
DOCKED: USER    Estimated Free Energy of Binding    =   -6.42 kcal/mol  [=(1)+(2)+(3)-(4)]
DOCKED: USER    Run = 2
DOCKED: USER    Estimated Free Energy of Binding    =   -9.18 kcal/mol  [=(1)+(2)+(3)-(4)]
DOCKED: USER    Run = 3
DOCKED: USER    Estimated Free Energy of Binding    =   -5.03 kcal/mol  [=(1)+(2)+(3)-(4)]
"""

NO_ENERGY_CONTENT = "DOCKED: USER    Run failed before completion\n"


def test_takes_the_most_negative_energy_across_all_runs(tmp_path):
    dlg_path = tmp_path / "EDTA.dlg"
    dlg_path.write_text(DLG_CONTENT)
    assert parse_best_energy(str(dlg_path)) == -9.18


def test_returns_none_when_no_energy_line_present(tmp_path):
    dlg_path = tmp_path / "failed.dlg"
    dlg_path.write_text(NO_ENERGY_CONTENT)
    assert parse_best_energy(str(dlg_path)) is None


def test_ligand_id_extracted_from_dlg_filename():
    assert ligand_id_from_dlg_path("results/3BH4/docking/ad4/citric_acid.dlg") == "citric_acid"
