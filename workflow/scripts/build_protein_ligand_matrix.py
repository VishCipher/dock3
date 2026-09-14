"""
Concatenate every protein target's merged_all_engines.tsv into one master
long-format table: target, ligand_id, engine, score_kcal_mol.

This is the actual "screen a library across multiple proteins" output --
open it in pandas and pivot on `target` to get a wide ligand x protein
matrix per engine, e.g.:

    import pandas as pd
    df = pd.read_csv("protein_ligand_matrix.tsv", sep="\t")
    df[df.engine == "vina"].pivot(index="ligand_id", columns="target", values="score_kcal_mol")

Run only via Snakemake (uses the injected `snakemake` object).
"""

import csv
import os


def main():
    merged_paths = list(snakemake.input.merged)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    rows = []
    for path in merged_paths:
        # path looks like results/<target>/summary/merged_all_engines.tsv
        parts = os.path.normpath(path).split(os.sep)
        target = parts[1]

        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rows.append((target, row["ligand_id"], row["engine"], row["score_kcal_mol"]))

    rows.sort(key=lambda r: (r[1], r[0], r[2]))  # by ligand, then target, then engine

    with open(output_path, "w") as out:
        out.write("target\tligand_id\tengine\tscore_kcal_mol\n")
        for target, ligand_id, engine, score in rows:
            out.write(f"{target}\t{ligand_id}\t{engine}\t{score}\n")

    n_targets = len(set(r[0] for r in rows))
    n_ligands = len(set(r[1] for r in rows))
    print(
        f"Wrote protein x ligand matrix: {n_targets} target(s) x {n_ligands} "
        f"ligand(s), {len(rows)} total rows, to {output_path}"
    )


if __name__ == "__main__":
    main()
