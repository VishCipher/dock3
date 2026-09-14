"""
Compute basic ADMET/drug-likeness descriptors for every ligand in the
library, using RDKit (already a dependency via Meeko/PandaDock -- no new
install needed).

A strong docking score means little if the compound itself is a poor drug
candidate. This is a cheap, well-established sanity filter (Lipinski's
Rule of Five plus a couple of common extensions), not a replacement for
real ADMET prediction software -- it flags obviously problematic
compounds, it doesn't clear anything for development.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import os

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def compute_properties(sdf_path):
    mol = Chem.MolFromMolFile(sdf_path)
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Lipinski.NumRotatableBonds(mol)

    violations = []
    if mw > 500:
        violations.append("MW>500")
    if logp > 5:
        violations.append("LogP>5")
    if hbd > 5:
        violations.append("HBD>5")
    if hba > 10:
        violations.append("HBA>10")

    return {
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotatable,
        "lipinski_violations": len(violations),
        "violation_flags": ";".join(violations) if violations else "none",
    }


def main():
    ligand_paths = list(snakemake.input.ligands)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    rows = []
    for path in ligand_paths:
        ligand_id = os.path.splitext(os.path.basename(path))[0]
        props = compute_properties(path)
        if props is None:
            print(f"[compute_admet] WARN: RDKit couldn't parse {path}")
            rows.append((ligand_id, None))
        else:
            rows.append((ligand_id, props))

    with open(output_path, "w") as out:
        out.write(
            "ligand_id\tmw\tlogp\thbd\thba\ttpsa\trotatable_bonds\t"
            "lipinski_violations\tviolation_flags\n"
        )
        for ligand_id, props in rows:
            if props is None:
                out.write(f"{ligand_id}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tPARSE_FAILED\n")
            else:
                out.write(
                    f"{ligand_id}\t{props['mw']}\t{props['logp']}\t{props['hbd']}\t"
                    f"{props['hba']}\t{props['tpsa']}\t{props['rotatable_bonds']}\t"
                    f"{props['lipinski_violations']}\t{props['violation_flags']}\n"
                )

    print(f"[compute_admet] Wrote drug-likeness properties for {len(rows)} ligands to {output_path}")


if __name__ == "__main__":
    main()
