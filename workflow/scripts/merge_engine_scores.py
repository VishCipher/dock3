"""
Merge Vina, smina, and AD4 per-engine score TSVs into one long-format table:
columns = ligand_id, engine, score. Makes cross-engine comparison a single
groupby in pandas/Excel instead of juggling three separate files.

Note: the three engines' scores are NOT directly comparable in absolute
terms (different scoring functions, different units of "goodness"), so this
table is for eyeballing relative agreement/disagreement between engines per
ligand, not for arithmetic like averaging the three scores together.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import csv


def read_scores(tsv_path, score_column):
    rows = []
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append((row["ligand_id"], row[score_column]))
    return rows


def main():
    vina_path = str(snakemake.input.vina)  # noqa: F821
    smina_path = str(snakemake.input.smina)  # noqa: F821
    ad4_path = str(snakemake.input.ad4)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    engine_files = [
        ("vina", vina_path, "best_affinity_kcal_mol"),
        ("smina", smina_path, "best_affinity_kcal_mol"),
        ("ad4", ad4_path, "best_energy_kcal_mol"),
    ]

    merged_rows = []
    for engine_name, path, score_column in engine_files:
        for ligand_id, score in read_scores(path, score_column):
            merged_rows.append((ligand_id, engine_name, score))

    # Sort for readability: group by ligand, then by engine name
    merged_rows.sort(key=lambda r: (r[0], r[1]))

    with open(output_path, "w") as out:
        out.write("ligand_id\tengine\tscore_kcal_mol\n")
        for ligand_id, engine_name, score in merged_rows:
            out.write(f"{ligand_id}\t{engine_name}\t{score}\n")

    print(f"Wrote merged comparison table ({len(merged_rows)} rows) to {output_path}")


if __name__ == "__main__":
    main()
