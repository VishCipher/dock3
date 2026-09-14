"""
Independently verify metal coordination geometry for each ligand's
top-ranked Vina pose, rather than trusting any single engine's internal
score. Vina/smina/AD4 have no model of coordination chemistry at all, and
even PandaDock's own coordination score used fallback/approximate metal
parameters in the runs seen so far -- this is a second, independent check
computed directly from atomic coordinates.

Method: read the metal ion's position from cofactor_report.json, read the
ligand's top-ranked Vina pose (first MODEL block in the docked PDBQT),
measure distances from the metal to every ligand heavy atom, and count how
many fall within a plausible coordination distance for that metal element.
Flag ligands whose coordination number falls outside the typical range for
that metal.

IDEAL_DISTANCE_RANGES and TYPICAL_COORDINATION_NUMBERS below are
literature-informed approximate defaults, not rigorously calibrated per
target -- treat flags as "worth a closer look", not a certified verdict.

Currently only checks Vina's pose (the PDBQT format already reliably
parsed elsewhere in this pipeline). Extending to smina/AD4 poses is a
natural next step once this is validated against a known case.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import json
import math

# Approximate ideal metal-to-donor-atom distance range, in Angstroms.
# Sources: general inorganic/bioinorganic chemistry references -- these are
# reasonable defaults, not target-specific calibrations.
IDEAL_DISTANCE_RANGES = {
    "Ca": (2.0, 2.8),
    "Zn": (1.9, 2.3),
    "Mg": (1.9, 2.3),
    "Mn": (1.9, 2.4),
    "Fe": (1.9, 2.3),
    "Cu": (1.9, 2.4),
    "Ni": (1.9, 2.3),
    "Co": (1.9, 2.3),
}

# Typical coordination number range for each metal (how many donor atoms
# normally surround it). Ca2+ is notably flexible (6-8 is common).
TYPICAL_COORDINATION_NUMBERS = {
    "Ca": (6, 8),
    "Zn": (4, 6),
    "Mg": (6, 6),
    "Mn": (5, 6),
    "Fe": (5, 6),
    "Cu": (4, 6),
    "Ni": (4, 6),
    "Co": (4, 6),
}

DEFAULT_DISTANCE_RANGE = (1.9, 2.8)
DEFAULT_COORDINATION_RANGE = (4, 8)

# Atom types that don't participate in metal coordination -- skip these
# when scanning the ligand pose.
NON_COORDINATING_TYPES = {"HD", "H", "C", "A"}


def load_metal_position(cofactor_report_path):
    with open(cofactor_report_path) as f:
        report = json.load(f)
    if not report["metals_found"]:
        return None
    metal = report["metals_found"][0]  # first detected metal, matching the rest of this pipeline
    return metal["symbol"], (metal["x"], metal["y"], metal["z"])


def load_top_pose_atoms(pdbqt_path):
    """Read only the FIRST model block (top-ranked pose) from a multi-model PDBQT."""
    atoms = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                atom_type = parts[-1]
                if atom_type in NON_COORDINATING_TYPES:
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_type, (x, y, z)))
    return atoms


def distance(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def evaluate_ligand(metal_symbol, metal_pos, pdbqt_path):
    dist_min, dist_max = IDEAL_DISTANCE_RANGES.get(metal_symbol, DEFAULT_DISTANCE_RANGE)
    coord_min, coord_max = TYPICAL_COORDINATION_NUMBERS.get(metal_symbol, DEFAULT_COORDINATION_RANGE)

    atoms = load_top_pose_atoms(pdbqt_path)
    if not atoms:
        return None

    distances_within_range = []
    closest_distance = min(distance(metal_pos, pos) for _, pos in atoms)

    for atom_type, pos in atoms:
        d = distance(metal_pos, pos)
        if dist_min <= d <= dist_max:
            distances_within_range.append((atom_type, round(d, 2)))

    coordination_number = len(distances_within_range)

    if coordination_number < coord_min:
        flag = "LOW_COORDINATION"
    elif coordination_number > coord_max:
        flag = "HIGH_COORDINATION"
    else:
        flag = "OK"

    return {
        "coordination_number": coordination_number,
        "expected_range": f"{coord_min}-{coord_max}",
        "closest_atom_distance": round(closest_distance, 2),
        "coordinating_atoms": "; ".join(f"{t}@{d}A" for t, d in distances_within_range),
        "flag": flag,
    }


def main():
    cofactor_report_path = str(snakemake.input.report)  # noqa: F821
    pdbqt_paths = list(snakemake.input.poses)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    metal_info = load_metal_position(cofactor_report_path)
    if metal_info is None:
        # Shouldn't happen -- this rule only runs for targets with a detected
        # metal -- but fail gracefully rather than crash the whole pipeline.
        with open(output_path, "w") as out:
            out.write("ligand_id\terror\n")
            out.write("ALL\tNo metal found in cofactor_report.json\n")
        print("[validate_metal_geometry] WARN: no metal found, wrote empty-ish report")
        return

    metal_symbol, metal_pos = metal_info

    rows = []
    for pdbqt_path in pdbqt_paths:
        import os
        ligand_id = os.path.basename(pdbqt_path).replace("_out.pdbqt", "")
        result = evaluate_ligand(metal_symbol, metal_pos, pdbqt_path)
        if result is None:
            rows.append((ligand_id, None))
        else:
            rows.append((ligand_id, result))

    rows.sort(key=lambda r: (r[1] is None, -(r[1]["coordination_number"] if r[1] else 0)))

    with open(output_path, "w") as out:
        out.write(
            "ligand_id\tcoordination_number\texpected_range\tclosest_atom_distance_A\t"
            "coordinating_atoms\tflag\n"
        )
        for ligand_id, result in rows:
            if result is None:
                out.write(f"{ligand_id}\tNA\tNA\tNA\tNA\tPARSE_FAILED\n")
            else:
                out.write(
                    f"{ligand_id}\t{result['coordination_number']}\t{result['expected_range']}\t"
                    f"{result['closest_atom_distance']}\t{result['coordinating_atoms']}\t"
                    f"{result['flag']}\n"
                )

    print(
        f"[validate_metal_geometry] Wrote independent coordination-geometry check "
        f"for {len(rows)} ligands ({metal_symbol} site) to {output_path}"
    )


if __name__ == "__main__":
    main()
