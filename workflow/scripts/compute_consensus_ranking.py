"""
Combine every independently-computed score for this target's ligand
library into one consensus ranking.

REAL FINDING from the first version of this script: a flat average of
ranks across all metrics put aspirin in 1st place overall, ahead of both
EDTA and EGTA -- despite aspirin never coming near the calcium at all. The
4 generic docking/interaction metrics (Vina, smina, AD4, PLIP interaction
count) all reward general pocket fit, which aspirin's aromatic ring does
well regardless of the metal, and outvoted the 2 metal-specific metrics in
a flat average.

FIX: this is a GATED consensus for cofactor-containing targets, not a flat
average. Ligands are first split into two tiers based on real evidence of
metal engagement (coordination_number >= 1, i.e. at least one atom found
genuinely within bonding distance of the metal), then ranked by the
generic metrics only WITHIN each tier. A compound with zero evidence of
touching the metal can never outrank one that does, no matter how well it
fits the pocket generally -- but among compounds that DO show real metal
engagement, the generic scores still meaningfully differentiate a good
binder from a mediocre one.

Note: the gate uses coordination_number >= 1, NOT the 'flag' column from
metal_geometry_validation.tsv. That flag requires the FULL typical
coordination number (6-8 simultaneous contacts for Ca2+), a bar no single
rigid docking pose realistically achieves -- an earlier version of this
script gated on flag=='OK' and rejected every ligand in the library,
including confirmed real chelators EDTA/EGTA. "Touches the metal at all"
is the meaningful, achievable signal for a single docked pose.

For targets with no cofactor at all, this reduces to a flat average across
whichever generic metrics are available (there's no metal-specific tier to
gate on).

This is a simple, transparent heuristic, not a statistically validated
ensemble method -- treat it as "which compounds have real, corroborated
evidence of the binding mode that matters here," not a rigorously
calibrated combined score.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import csv


def read_rank_column(path, ligand_col="ligand_id", rank_col="rank"):
    """Read a TSV that already has its own precomputed rank column."""
    ranks = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rank_val = row.get(rank_col)
            if rank_val and rank_val != "NA":
                ranks[row[ligand_col]] = int(rank_val)
    return ranks


def read_geometry_evidence(path, ligand_col="ligand_id", distance_col="closest_atom_distance_A", coord_num_col="coordination_number"):
    """
    Returns (rank_by_distance, has_real_evidence) for every ligand.

    has_real_evidence means coordination_number >= 1 -- at least one atom
    genuinely within bonding distance of the metal. NOT the same as the
    'flag' column in metal_geometry_validation.tsv, which requires the FULL
    typical coordination number (6-8 simultaneous contacts for Ca2+) --
    an unrealistically strict bar that no single rigid docking pose
    achieves, chelator or not (confirmed: gating on flag=='OK' rejected
    every ligand in this library, including EDTA/EGTA). "Touches the
    metal at all" is the meaningful, achievable signal here.
    """
    rows = []
    evidence = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            dist = row.get(distance_col)
            if dist and dist != "NA":
                rows.append((row[ligand_col], float(dist)))
            coord_num = row.get(coord_num_col)
            evidence[row[ligand_col]] = coord_num not in (None, "NA") and int(coord_num) >= 1
    rows.sort(key=lambda r: r[1])
    rank_by_distance = {ligand_id: rank for rank, (ligand_id, _) in enumerate(rows, start=1)}
    return rank_by_distance, evidence


def main():
    input_ = snakemake.input  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    generic_ranks = {
        "vina": read_rank_column(str(input_.vina)),
        "smina": read_rank_column(str(input_.smina)),
        "ad4": read_rank_column(str(input_.ad4)),
        "plip": read_rank_column(str(input_.plip)),
    }

    has_metal = hasattr(input_, "metal_energy")
    metal_energy_rank = {}
    geometry_rank = {}
    real_metal_evidence = {}

    if has_metal:
        metal_energy_rank = read_rank_column(str(input_.metal_energy))
        geometry_rank, real_metal_evidence = read_geometry_evidence(str(input_.geometry))

    all_ligands = set()
    for ranks in generic_ranks.values():
        all_ligands.update(ranks.keys())
    all_ligands.update(metal_energy_rank.keys())

    rows = []
    for ligand_id in all_ligands:
        generic_available = {
            name: ranks[ligand_id] for name, ranks in generic_ranks.items() if ligand_id in ranks
        }
        avg_generic_rank = (
            sum(generic_available.values()) / len(generic_available) if generic_available else None
        )

        if has_metal:
            has_evidence = real_metal_evidence.get(ligand_id, False)
            tier = 1 if has_evidence else 2  # tier 1 = real metal evidence, sorts first
            m_energy = metal_energy_rank.get(ligand_id)
            g_rank = geometry_rank.get(ligand_id)
        else:
            has_evidence = None
            tier = None
            m_energy = None
            g_rank = None

        rows.append((ligand_id, tier, avg_generic_rank, generic_available, m_energy, g_rank, has_evidence))

    if has_metal:
        # Sort by tier first (real evidence beats no evidence, always),
        # then by generic-metric average rank as the tiebreaker within tier.
        rows.sort(key=lambda r: (r[1], r[2] if r[2] is not None else 999))
    else:
        rows.sort(key=lambda r: r[2] if r[2] is not None else 999)

    metric_names = list(generic_ranks.keys())
    with open(output_path, "w") as out:
        header = ["overall_rank", "ligand_id"]
        if has_metal:
            header += ["metal_evidence_tier", "has_real_metal_contact"]
        header += ["avg_generic_rank"] + [f"{m}_rank" for m in metric_names]
        if has_metal:
            header += ["metal_energy_rank", "geometry_rank"]
        out.write("\t".join(header) + "\n")

        for overall_rank, (ligand_id, tier, avg_rank, available, m_energy, g_rank, has_evidence) in enumerate(rows, start=1):
            row_values = [str(overall_rank), ligand_id]
            if has_metal:
                row_values += [str(tier), str(has_evidence)]
            row_values += [f"{avg_rank:.2f}" if avg_rank is not None else "NA"]
            row_values += [str(available.get(m, "NA")) for m in metric_names]
            if has_metal:
                row_values += [str(m_energy) if m_energy is not None else "NA"]
                row_values += [str(g_rank) if g_rank is not None else "NA"]
            out.write("\t".join(row_values) + "\n")

    print(
        f"[compute_consensus_ranking] Wrote {'gated' if has_metal else 'flat-average'} "
        f"consensus ranking for {len(rows)} ligands to {output_path}"
    )


if __name__ == "__main__":
    main()
