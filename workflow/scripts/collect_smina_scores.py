"""
Parse smina log files and write a ranked TSV summary.

smina shares Vina's result-table log format exactly (it's a Vina fork), so
this parsing logic is identical to collect_scores.py -- kept as a separate
file so each engine's output stays independently inspectable/debuggable.

Run only via Snakemake (uses the injected `snakemake` object for
input/output paths) — not meant to be run standalone.
"""

import os
import re

# Vina's result table looks like:
#    mode |   affinity | dist from best mode
#         | (kcal/mol) | rmsd l.b.| rmsd u.b.
# -----+------------+----------+----------
#    1        -7.2      0.000      0.000
# We want the affinity from mode 1 (the best pose).
MODE_1_PATTERN = re.compile(r"^\s*1\s+(-?\d+\.?\d*)")


def parse_best_affinity(log_path: str) -> float | None:
    with open(log_path) as f:
        for line in f:
            match = MODE_1_PATTERN.match(line)
            if match:
                return float(match.group(1))
    return None


def ligand_id_from_log_path(log_path: str) -> str:
    filename = os.path.basename(log_path)
    return filename.replace("_log.txt", "")


def main():
    log_paths = list(snakemake.input.logs)  # noqa: F821 (injected by Snakemake)
    output_path = str(snakemake.output)  # noqa: F821

    rows = []
    for log_path in log_paths:
        ligand_id = ligand_id_from_log_path(log_path)
        affinity = parse_best_affinity(log_path)
        if affinity is None:
            print(f"[WARN] Could not parse affinity from {log_path}")
        rows.append((ligand_id, affinity))

    # Rank best (most negative) first; unparsed scores sink to the bottom.
    rows.sort(key=lambda row: (row[1] is None, row[1]))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as out:
        out.write("rank\tligand_id\tbest_affinity_kcal_mol\n")
        for rank, (ligand_id, affinity) in enumerate(rows, start=1):
            affinity_str = f"{affinity:.2f}" if affinity is not None else "NA"
            out.write(f"{rank}\t{ligand_id}\t{affinity_str}\n")

    print(f"Wrote {len(rows)} ligand scores to {output_path}")


if __name__ == "__main__":
    main()
