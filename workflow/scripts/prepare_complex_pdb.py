"""
Combine a target's cleaned receptor PDB with one ligand's top-ranked Vina
pose into a single PDB file -- the input format PLIP (and most
protein-ligand interaction tools) actually expect: one file, protein as
ATOM records, ligand as a HETATM block with its own residue name/chain.

Uses the same proven technique as elsewhere in this pipeline: extract
model 1 (top pose) from the multi-model Vina PDBQT via OpenBabel, which
correctly reconstructs bond connectivity from the PDBQT's own torsion-tree
structure (validated earlier in this project against 1HSG/indinavir).

The ligand is force-labeled HETATM / residue "LIG" / chain "Z" regardless
of whatever OpenBabel's default PDB writer chooses, so downstream tools
have a guaranteed-consistent way to identify it.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import subprocess
import tempfile
import os


def extract_top_pose_as_pdb(pdbqt_path, out_pdb_path):
    """Use OpenBabel to pull model 1 (top pose) out as a proper PDB with bonds intact."""
    subprocess.run(
        ["obabel", pdbqt_path, "-O", out_pdb_path, "-f", "1", "-l", "1"],
        check=True, capture_output=True
    )


def relabel_ligand_pdb_lines(pdb_path):
    """Force every atom line to HETATM / residue LIG / chain Z / resnum 1."""
    relabeled = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # Rebuild the line with forced record type, residue name, chain, resnum
                # while preserving atom serial/name/coordinates/etc. columns.
                record = "HETATM" + line[6:17] + "LIG" + " Z" + "   1" + line[26:]
                relabeled.append(record)
            elif line.startswith(("CONECT",)):
                relabeled.append(line)  # keep bond records if OpenBabel wrote any
    return relabeled


def main():
    receptor_path = str(snakemake.input.receptor)  # noqa: F821
    pose_pdbqt_path = str(snakemake.input.pose)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    with tempfile.TemporaryDirectory() as tmpdir:
        ligand_pdb_path = os.path.join(tmpdir, "ligand_pose.pdb")
        extract_top_pose_as_pdb(pose_pdbqt_path, ligand_pdb_path)
        ligand_lines = relabel_ligand_pdb_lines(ligand_pdb_path)

    with open(receptor_path) as f:
        receptor_lines = [l for l in f if l.startswith(("ATOM", "HETATM", "TER"))]

    with open(output_path, "w") as out:
        out.writelines(receptor_lines)
        out.write("TER\n")
        out.writelines(ligand_lines)
        out.write("END\n")

    print(f"[prepare_complex_pdb] Wrote combined receptor+ligand complex to {output_path}")


if __name__ == "__main__":
    main()
