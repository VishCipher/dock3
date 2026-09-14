"""
Scan a receptor PDB for known metal-ion HETATM records.

Run only via Snakemake (uses the injected `snakemake` object) — not meant to
be run standalone.

Extend METAL_CODES if your target has a metal we don't already list.
"""

import json

# PDB residue code -> element symbol (as expected by pandadock-metal's
# --metal-type flag). Add more here if you hit an unrecognized cofactor.
# Deliberately excludes alkali metals (Na, K) -- see keep_protein_and_metals.py
# for why: they're near-always crystallization artifacts, not functional sites.
METAL_CODES = {
    "ZN": "Zn",
    "CA": "Ca",
    "MG": "Mg",
    "FE": "Fe",
    "FE2": "Fe",
    "MN": "Mn",
    "CU": "Cu",
    "CU1": "Cu",
    "NI": "Ni",
    "CO": "Co",
}


def scan_receptor(receptor_path):
    metals_found = []
    seen = set()  # dedupe by (chain, resnum, resname) in case of multi-line entries

    with open(receptor_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue

            resname = line[17:20].strip()
            if resname not in METAL_CODES:
                continue

            chain = line[21].strip()
            resnum = line[22:26].strip()
            key = (chain, resnum, resname)
            if key in seen:
                continue
            seen.add(key)

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            metals_found.append(
                {
                    "resname": resname,
                    "symbol": METAL_CODES[resname],
                    "chain": chain,
                    "resnum": resnum,
                    "x": x,
                    "y": y,
                    "z": z,
                }
            )

    return metals_found


def main():
    receptor_path = str(snakemake.input.receptor)  # noqa: F821
    output_path = str(snakemake.output.report)  # noqa: F821

    metals_found = scan_receptor(receptor_path)
    report = {
        "has_metal": len(metals_found) > 0,
        "metals_found": metals_found,
    }

    with open(output_path, "w") as out:
        json.dump(report, out, indent=2)

    if metals_found:
        symbols = ", ".join(f"{m['resname']} ({m['chain']}:{m['resnum']})" for m in metals_found)
        print(f"[detect_cofactors] Found {len(metals_found)} metal ion(s): {symbols}")
    else:
        print("[detect_cofactors] No known metal ions found — standard branch only.")


if __name__ == "__main__":
    main()
