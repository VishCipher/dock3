"""
Gather PandaDock metal-site docking results across all ligands into one
ranked TSV.

IMPORTANT FINDING: metal_docking_summary.json's own avg_metal_binding_score
and best_coordination_score fields come back as a flat 0.000 for EVERY
ligand regardless of actual chemistry, even after supplying a real
AutoDock4 parameter file (AD4.1_bound.dat) via --metal-params -- confirmed
across multiple real runs. That part of PandaDock's own reported summary
is not usable as-is.

However, its LOG output contains a real, differentiating value that IS
usable: a line like

    INFO:pandadock.docking.pandadock:Search complete: 30 runs, 12780 local
    minima, 1468109 energy evaluations, best -9.283 kcal/mol

This is a genuine physics-based (Vina-function-driven) search energy,
confirmed in real testing to correctly rank a known calcium chelator
(EDTA) as the most negative/best among the library. This script parses
THAT value from the log file as the primary usable PandaDock score,
alongside the (currently uninformative) JSON fields for reference.

Ranks by best_search_energy_kcal_mol ascending (more negative = better),
consistent with every other engine in this pipeline.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import json
import os
import re

SEARCH_ENERGY_PATTERN = re.compile(r"Search complete:.*best\s+(-?\d+\.?\d*)\s*kcal/mol")


def ligand_id_from_outdir(outdir_path: str) -> str:
    return os.path.basename(outdir_path.rstrip("/"))


def parse_best_search_energy(log_path):
    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        for line in f:
            match = SEARCH_ENERGY_PATTERN.search(line)
            if match:
                return float(match.group(1))
    return None


def main():
    outdirs = list(snakemake.input.outdirs)  # noqa: F821
    log_paths = list(snakemake.input.logs)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    # Match each outdir to its corresponding log by ligand_id
    log_by_ligand = {}
    for log_path in log_paths:
        # log filename looks like TARGET_dock_metal_LIGAND.log
        stem = os.path.basename(log_path).replace(".log", "")
        ligand_id = stem.split("_dock_metal_")[-1]
        log_by_ligand[ligand_id] = log_path

    rows = []
    for outdir in outdirs:
        ligand_id = ligand_id_from_outdir(outdir)
        summary_path = os.path.join(outdir, "metal_docking_summary.json")

        avg_binding, best_coord, viol_rate = None, None, None
        if os.path.exists(summary_path):
            with open(summary_path) as f:
                summary = json.load(f)
            avg_binding = summary.get("avg_metal_binding_score")
            best_coord = summary.get("best_coordination_score")
            viol_rate = summary.get("violation_rate")
        else:
            print(f"[collect_metal_scores] WARN: no summary JSON found in {outdir}")

        log_path = log_by_ligand.get(ligand_id)
        search_energy = parse_best_search_energy(log_path) if log_path else None
        if search_energy is None:
            print(f"[collect_metal_scores] WARN: no search energy found for {ligand_id}")

        rows.append((ligand_id, search_energy, avg_binding, best_coord, viol_rate))

    rows.sort(key=lambda r: (r[1] is None, r[1]))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as out:
        out.write(
            "rank\tligand_id\tbest_search_energy_kcal_mol\tavg_metal_binding_score\t"
            "best_coordination_score\tviolation_rate\n"
        )
        for rank, (ligand_id, energy, binding, coord, viol) in enumerate(rows, start=1):
            energy_str = f"{energy:.3f}" if energy is not None else "NA"
            binding_str = f"{binding:.3f}" if binding is not None else "NA"
            coord_str = f"{coord:.3f}" if coord is not None else "NA"
            viol_str = f"{viol:.2f}" if viol is not None else "NA"
            out.write(f"{rank}\t{ligand_id}\t{energy_str}\t{binding_str}\t{coord_str}\t{viol_str}\n")

    print(f"Wrote {len(rows)} metal-site ligand scores to {output_path}")


if __name__ == "__main__":
    main()
