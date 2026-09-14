"""
Parse AutoDock4 .dlg output files and write a ranked TSV summary.

AD4's .dlg format repeats a block per docking run, each containing a line
like:

    DOCKED: USER    Estimated Free Energy of Binding    =   -7.42 kcal/mol  [...]

We take the lowest (most negative) value across all runs in the file as the
ligand's best score. This is a reasonable first-pass approximation; the more
rigorous approach (best energy from the largest/most-populated cluster in
the file's "CLUSTERING HISTOGRAM" section) matters more once you're
comparing top candidates closely, not for an initial screen.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import os
import re

ENERGY_PATTERN = re.compile(
    r"Estimated Free Energy of Binding\s*=\s*(-?\d+\.?\d*)"
)


def parse_best_energy(dlg_path: str):
    best = None
    with open(dlg_path) as f:
        for line in f:
            match = ENERGY_PATTERN.search(line)
            if match:
                value = float(match.group(1))
                if best is None or value < best:
                    best = value
    return best


def ligand_id_from_dlg_path(dlg_path: str) -> str:
    return os.path.basename(dlg_path).replace(".dlg", "")


def main():
    dlg_paths = list(snakemake.input.dlgs)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    rows = []
    for dlg_path in dlg_paths:
        ligand_id = ligand_id_from_dlg_path(dlg_path)
        energy = parse_best_energy(dlg_path)
        if energy is None:
            print(f"[WARN] No binding energy found in {dlg_path} — check the run completed.")
        rows.append((ligand_id, energy))

    rows.sort(key=lambda row: (row[1] is None, row[1]))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as out:
        out.write("rank\tligand_id\tbest_energy_kcal_mol\n")
        for rank, (ligand_id, energy) in enumerate(rows, start=1):
            energy_str = f"{energy:.2f}" if energy is not None else "NA"
            out.write(f"{rank}\t{ligand_id}\t{energy_str}\n")

    print(f"Wrote {len(rows)} AD4 ligand scores to {output_path}")


if __name__ == "__main__":
    main()
