"""
Parse PLIP's XML interaction report for one ligand into a compact summary.

Verified against a real PLIP 3.0.1 run: the report contains one
<bindingsite> block per hetero group in the structure (every calcium ion
AND our docked ligand each get their own block) -- this specifically finds
the block where <hetid> is "LIG" (the residue name prepare_complex_pdb.py
force-assigns to the docked pose) and ignores every other hetero group's
block (e.g. the calcium ions' own binding-site entries).

REAL FINDING from initial testing: the docked EDTA pose's own
<metal_complexes> block came back EMPTY -- PLIP did not detect the docked
Vina pose as directly coordinating the calcium, even though EDTA is a
genuine chelator. This is consistent with validate_metal_geometry.py's
independent LOW_COORDINATION finding for the same pose, AND with
EDTA's hydrogen bonds/salt bridges landing on the exact same residues
(GLY301, PRO407, ASP408) that coordinate the real calcium in its own
binding-site block. Three independently-computed checks agreeing: Vina
places the ligand in the correct general pocket without achieving true
metal coordination -- a real limitation of classical docking on this
ligand class, not a pipeline bug.

Run only via Snakemake (uses the injected `snakemake` object).
"""

import os
import xml.etree.ElementTree as ET


def find_ligand_bindingsite(xml_path, target_hetid="LIG"):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    for bindingsite in root.findall("bindingsite"):
        hetid = bindingsite.findtext("identifiers/hetid")
        if hetid == target_hetid:
            return bindingsite
    return None


def summarize_bindingsite(bindingsite):
    interactions = bindingsite.find("interactions")

    def count(tag):
        block = interactions.find(tag)
        return 0 if block is None else len(list(block))

    hbonds = count("hydrogen_bonds")
    salt_bridges = count("salt_bridges")
    hydrophobic = count("hydrophobic_interactions")
    water_bridges = count("water_bridges")
    pi_stacks = count("pi_stacks")
    pi_cation = count("pi_cation_interactions")
    halogen = count("halogen_bonds")
    metal_complexes = count("metal_complexes")

    total = (hbonds + salt_bridges + hydrophobic + water_bridges
             + pi_stacks + pi_cation + halogen + metal_complexes)

    # Residues involved in H-bonds/salt bridges -- the two interaction
    # types this pipeline's carboxylic-acid-heavy ligand library actually uses.
    contacted_residues = set()
    for hb in interactions.findall("hydrogen_bonds/hydrogen_bond"):
        restype, resnr = hb.findtext("restype"), hb.findtext("resnr")
        if restype and resnr:
            contacted_residues.add(f"{restype}{resnr}")
    for sb in interactions.findall("salt_bridges/salt_bridge"):
        restype, resnr = sb.findtext("restype"), sb.findtext("resnr")
        if restype and resnr:
            contacted_residues.add(f"{restype}{resnr}")

    return {
        "hydrogen_bonds": hbonds,
        "salt_bridges": salt_bridges,
        "hydrophobic_interactions": hydrophobic,
        "direct_metal_coordination": metal_complexes,
        "total_interactions": total,
        "contacted_residues": ";".join(sorted(contacted_residues)),
    }


EMPTY_SUMMARY = {
    "hydrogen_bonds": 0, "salt_bridges": 0, "hydrophobic_interactions": 0,
    "direct_metal_coordination": 0, "total_interactions": 0, "contacted_residues": "",
}


def ligand_id_from_outdir(outdir_path):
    return os.path.basename(outdir_path.rstrip("/"))


def main():
    outdirs = list(snakemake.input.outdirs)  # noqa: F821
    output_path = str(snakemake.output)  # noqa: F821

    rows = []
    for outdir in outdirs:
        ligand_id = ligand_id_from_outdir(outdir)
        xml_path = os.path.join(outdir, "plip_report.xml")

        if not os.path.exists(xml_path):
            print(f"[collect_plip_interactions] WARN: no PLIP XML found for {ligand_id}")
            rows.append((ligand_id, None))
            continue

        bindingsite = find_ligand_bindingsite(xml_path)
        # bindingsite is None when PLIP found literally zero interactions for
        # this ligand at all -- a valid outcome (weak/no binder), not an error.
        summary = summarize_bindingsite(bindingsite) if bindingsite is not None else EMPTY_SUMMARY
        rows.append((ligand_id, summary))

    rows.sort(key=lambda r: (r[1] is None, -(r[1]["total_interactions"] if r[1] else 0)))

    with open(output_path, "w") as out:
        out.write(
            "rank\tligand_id\ttotal_interactions\thydrogen_bonds\tsalt_bridges\t"
            "hydrophobic_interactions\tdirect_metal_coordination\tcontacted_residues\n"
        )
        for rank, (ligand_id, s) in enumerate(rows, start=1):
            if s is None:
                out.write(f"{rank}\t{ligand_id}\tNA\tNA\tNA\tNA\tNA\tPARSE_FAILED\n")
            else:
                out.write(
                    f"{rank}\t{ligand_id}\t{s['total_interactions']}\t{s['hydrogen_bonds']}\t"
                    f"{s['salt_bridges']}\t{s['hydrophobic_interactions']}\t"
                    f"{s['direct_metal_coordination']}\t{s['contacted_residues']}\n"
                )

    print(f"[collect_plip_interactions] Wrote PLIP interaction summary for {len(rows)} ligands to {output_path}")


if __name__ == "__main__":
    main()
