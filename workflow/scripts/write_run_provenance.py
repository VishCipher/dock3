"""
Capture basic run provenance: when this run happened, which targets/ligands/
engine branches were involved, and (if available) which git commit of the
pipeline produced it.

Given how many config/file-sync bugs came up during this pipeline's
development, being able to answer "which exact version of the Snakefile
produced this results folder" is cheap insurance worth having by default.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import datetime
import json
import subprocess


def get_git_info():
    """Best-effort: returns None if the project isn't under git at all."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
        has_uncommitted_changes = subprocess.call(
            ["git", "diff", "--quiet"], stderr=subprocess.DEVNULL
        ) != 0
        return {"commit": commit, "uncommitted_changes": has_uncommitted_changes}
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def main():
    output_path = str(snakemake.output)  # noqa: F821
    targets = list(snakemake.params.targets)  # noqa: F821
    ligand_ids = list(snakemake.params.ligand_ids)  # noqa: F821

    provenance = {
        "run_timestamp": datetime.datetime.now().isoformat(),
        "targets": targets,
        "n_ligands": len(ligand_ids),
        "ligand_ids": ligand_ids,
        "engines_active": {
            "vina": True,
            "smina": True,
            "ad4": True,
            "pandadock_metal": "auto-detected per target (see each target's cofactor_report.json)",
        },
        "git": get_git_info(),
    }

    with open(output_path, "w") as f:
        json.dump(provenance, f, indent=2)

    print(f"[write_run_provenance] Wrote run metadata to {output_path}")


if __name__ == "__main__":
    main()
