"""
Redocking benchmark: for each real protein-ligand complex, extract the
receptor and native crystal ligand, re-dock the ligand back into its own
pocket with Vina, and measure RMSD between the docked pose and the real
crystal position. This is the standard "redocking accuracy" validation
used throughout the docking literature (including PDBbind-based studies).

IMPORTANT ON DATA SOURCE: true PDBbind requires a license/account to
download. Rather than pretend to use it, this benchmarks against publicly
downloadable RCSB structures that are widely cited in the docking
literature as validation cases -- 1HSG is the exact complex this project
already validated by hand (0.39 A RMSD, see project history). Treat this
as a small proof-of-concept redocking study, not a full CASF-2016-scale
benchmark.

LIGAND AUTO-DETECTION: rather than hardcode each complex's native ligand
residue code from memory (which, if wrong, would silently corrupt the
whole benchmark -- exactly the class of mistake this project has hit
repeatedly with other tools), this script identifies the native ligand
automatically: the largest HETATM group in the structure that isn't a
known solvent/buffer/common-ion artifact. It PRINTS what it detected for
every complex so you can visually sanity-check before trusting the result.

Requires: obabel, vina, meeko (mk_prepare_receptor.py / mk_prepare_ligand.py) --
same tools already used by the main pipeline.

Usage:
    python3 run_pdbbind_benchmark.py
"""

import csv
import math
import os
import subprocess
import urllib.request

# A small, diverse set of real, publicly-downloadable complexes widely used
# as redocking benchmark cases in the docking literature. 1HSG matches this
# project's own hand-validated case (0.39 A RMSD). Ligand identity for the
# rest is AUTO-DETECTED and printed at runtime, not assumed.
BENCHMARK_PDB_IDS = [
    "1HSG",  # HIV-1 protease -- already hand-validated in this project
    "3PTB",  # Beta-trypsin
    "1STP",  # Streptavidin
    "4DFR",  # Dihydrofolate reductase
    "1HVR",  # HIV-1 protease, second inhibitor
    "2CPP",  # Cytochrome P450cam
]

# HETATM residue names to exclude when auto-detecting the native ligand --
# common crystallization solvents/buffers/ions, never the actual ligand of interest.
NON_LIGAND_RESIDUES = {
    "HOH", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "TRS", "PEG", "1PE",
    "MPD", "CIT", "FMT", "IPA", "MES", "BME", "NA", "CL", "MG", "CA", "ZN",
    "MN", "FE", "CU", "NI", "CO", "K", "BR", "IOD", "NH4", "UNK",
}

WORKDIR = "benchmark_work"
RESULTS_TSV = "benchmark_results.tsv"


def run(cmd, **kwargs):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kwargs)


def download_pdb(pdb_id, out_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, out_path)


def detect_native_ligand(pdb_path):
    """
    Returns (resname, chain, resnum, n_atoms) for the largest non-solvent
    HETATM group, or None if nothing qualifies.
    """
    groups = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname in NON_LIGAND_RESIDUES:
                continue
            chain = line[21].strip()
            resnum = line[22:26].strip()
            key = (resname, chain, resnum)
            groups.setdefault(key, 0)
            groups[key] += 1

    if not groups:
        return None

    (resname, chain, resnum), n_atoms = max(groups.items(), key=lambda kv: kv[1])
    return resname, chain, resnum, n_atoms


def extract_receptor_and_ligand(pdb_path, resname, chain, resnum, receptor_out, ligand_out):
    with open(pdb_path) as f, open(receptor_out, "w") as rec, open(ligand_out, "w") as lig:
        for line in f:
            if line.startswith("ATOM"):
                rec.write(line)
            elif line.startswith("HETATM"):
                this_resname = line[17:20].strip()
                this_chain = line[21].strip()
                this_resnum = line[22:26].strip()
                if (this_resname, this_chain, this_resnum) == (resname, chain, resnum):
                    lig.write(line)


def compute_box_from_ligand(ligand_pdb_path, padding=8.0):
    coords = []
    with open(ligand_pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))

    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    xs, ys, zs = [c[0] for c in coords], [c[1] for c in coords], [c[2] for c in coords]
    sx = (max(xs) - min(xs)) + padding
    sy = (max(ys) - min(ys)) + padding
    sz = (max(zs) - min(zs)) + padding
    return (cx, cy, cz), (sx, sy, sz)


def run_one_complex(pdb_id):
    complex_dir = os.path.join(WORKDIR, pdb_id)
    os.makedirs(complex_dir, exist_ok=True)

    raw_pdb = os.path.join(complex_dir, f"{pdb_id}.pdb")
    if not os.path.exists(raw_pdb):
        download_pdb(pdb_id, raw_pdb)

    detected = detect_native_ligand(raw_pdb)
    if detected is None:
        print(f"  [{pdb_id}] SKIP: no suitable native ligand detected")
        return None
    resname, chain, resnum, n_atoms = detected
    print(f"  [{pdb_id}] Detected native ligand: {resname} chain {chain} resnum {resnum} ({n_atoms} atoms)")

    receptor_pdb = os.path.join(complex_dir, "receptor.pdb")
    native_ligand_pdb = os.path.join(complex_dir, "native_ligand.pdb")
    extract_receptor_and_ligand(raw_pdb, resname, chain, resnum, receptor_pdb, native_ligand_pdb)

    center, size = compute_box_from_ligand(native_ligand_pdb)
    print(f"  [{pdb_id}] Box center={tuple(round(c, 2) for c in center)}, "
          f"size={tuple(round(s, 2) for s in size)}")

    # Prep with Meeko (same tools as the main pipeline)
    receptor_pdbqt = os.path.join(complex_dir, "receptor.pdbqt")
    run(f"mk_prepare_receptor.py --read_pdb {receptor_pdb} -p {receptor_pdbqt}")

    native_ligand_sdf = os.path.join(complex_dir, "native_ligand.sdf")
    run(f"obabel {native_ligand_pdb} -O {native_ligand_sdf}")
    ligand_pdbqt = os.path.join(complex_dir, "ligand.pdbqt")
    run(f"mk_prepare_ligand.py -i {native_ligand_sdf} -o {ligand_pdbqt}")

    # Dock
    docked_pdbqt = os.path.join(complex_dir, "docked_out.pdbqt")
    run(
        f"vina --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} "
        f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} "
        f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} "
        f"--exhaustiveness 8 --out {docked_pdbqt}"
    )

    # Extract top pose, compute RMSD to the real crystal ligand
    best_pose_pdb = os.path.join(complex_dir, "best_pose.pdb")
    run(f"obabel {docked_pdbqt} -O {best_pose_pdb} -f 1 -l 1")

    rmsd_result = run(f"obrms {native_ligand_pdb} {best_pose_pdb}")
    rmsd_line = rmsd_result.stdout.strip().splitlines()[-1]
    rmsd = float(rmsd_line.split()[-1])

    print(f"  [{pdb_id}] RMSD to crystal pose: {rmsd:.2f} A")
    return {
        "pdb_id": pdb_id,
        "ligand_resname": resname,
        "ligand_n_atoms": n_atoms,
        "rmsd_angstrom": rmsd,
        "success": rmsd < 2.0,  # standard docking-literature success threshold
    }


def main():
    os.makedirs(WORKDIR, exist_ok=True)
    results = []

    print(f"Running redocking benchmark on {len(BENCHMARK_PDB_IDS)} complexes...\n")
    for pdb_id in BENCHMARK_PDB_IDS:
        print(f"=== {pdb_id} ===")
        try:
            result = run_one_complex(pdb_id)
            if result:
                results.append(result)
        except subprocess.CalledProcessError as e:
            print(f"  [{pdb_id}] FAILED: {e.stderr[:300] if e.stderr else e}")
        except Exception as e:
            print(f"  [{pdb_id}] FAILED: {e}")
        print()

    with open(RESULTS_TSV, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["pdb_id", "ligand_resname", "ligand_n_atoms", "rmsd_angstrom", "success"])
        for r in results:
            writer.writerow([r["pdb_id"], r["ligand_resname"], r["ligand_n_atoms"],
                              f"{r['rmsd_angstrom']:.3f}", r["success"]])

    n_success = sum(1 for r in results if r["success"])
    n_total = len(results)
    print("=" * 50)
    print(f"BENCHMARK COMPLETE: {n_success}/{n_total} complexes redocked within 2.0 A RMSD "
          f"({100 * n_success / n_total:.0f}% success)" if n_total else "No complexes completed successfully.")
    print(f"Full results written to {RESULTS_TSV}")


if __name__ == "__main__":
    main()
