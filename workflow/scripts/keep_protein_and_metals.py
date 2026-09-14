"""
Build a cleaned receptor PDB: keep all protein ATOM lines, plus HETATM lines
for known FUNCTIONAL metal cofactors only (drop waters, buffer ions,
crystallization salts, and any native co-crystallized ligand).

Used by BOTH the standard branch (prep_receptor_clean, feeding Vina/smina/AD4)
and the metal branch (feeding PandaDock) -- every engine should see the same
cleaned structure, differing only in whether it's converted to PDBQT
afterward.

Deliberately excludes alkali metals (Na, K): these show up constantly in
crystal structures as crystallization-buffer artifacts, essentially never as
a genuine, druggable coordination site the way Ca/Zn/Mg/Mn/Fe/Cu are. If your
target genuinely needs one, add it back to METAL_CODES explicitly.

Run only via Snakemake (uses the injected `snakemake` object).
"""

METAL_CODES = {"ZN", "CA", "MG", "FE", "FE2", "MN", "CU", "CU1", "NI", "CO"}


def clean_receptor(in_path, out_path):
    """Keep protein ATOM lines + whitelisted metal HETATM lines only. Returns count kept."""
    kept = 0
    with open(in_path) as f, open(out_path, "w") as out:
        for line in f:
            if line.startswith("ATOM"):
                out.write(line)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in METAL_CODES:
                    out.write(line)
                    kept += 1
            elif line.startswith(("TER", "END")):
                out.write(line)
    return kept


def main():
    in_path = str(snakemake.input.receptor)  # noqa: F821
    out_path = str(snakemake.output.pdb)  # noqa: F821

    kept = clean_receptor(in_path, out_path)

    print(f"[prep_receptor_clean] Wrote {out_path} (kept {kept} metal HETATM line(s))")


if __name__ == "__main__":
    main()
