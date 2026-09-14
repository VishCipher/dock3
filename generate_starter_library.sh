#!/bin/bash
# Generates a 12-compound starter library for testing the docking pipeline
# against a calcium-binding site (e.g. 3BH4).
#
# Includes:
#   - Known calcium chelators (positive controls — expect strong scores)
#       EDTA, EGTA, citric acid, oxalic acid, malonic acid, tartaric acid
#   - Unrelated drug-like molecules (negative controls — expect weak scores)
#       aspirin, caffeine, ascorbic acid, glycine, phosphoric acid, salicylic acid
#
# Requires OpenBabel (already installed from earlier steps).
# Run from your project root — it creates/populates a `ligands/` folder.

set -e
mkdir -p ligands

declare -A LIGANDS=(
  ["EDTA"]="OC(=O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O"
  ["EGTA"]="OC(=O)CN(CC(=O)O)CCOCCOCCN(CC(=O)O)CC(=O)O"
  ["citric_acid"]="OC(=O)CC(O)(CC(=O)O)C(=O)O"
  ["oxalic_acid"]="OC(=O)C(=O)O"
  ["malonic_acid"]="OC(=O)CC(=O)O"
  ["tartaric_acid"]="OC(=O)C(O)C(O)C(=O)O"
  ["aspirin"]="CC(=O)OC1=CC=CC=C1C(=O)O"
  ["caffeine"]="CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
  ["ascorbic_acid"]="OC1=C(O)C(=O)OC1C(O)CO"
  ["glycine"]="NCC(=O)O"
  ["phosphoric_acid"]="OP(=O)(O)O"
  ["salicylic_acid"]="OC(=O)C1=CC=CC=C1O"
)

echo "Generating ${#LIGANDS[@]} ligand SDFs into ./ligands/ ..."
for name in "${!LIGANDS[@]}"; do
  smiles="${LIGANDS[$name]}"
  echo "  $name"
  echo "$smiles" | obabel -ismi -O "ligands/${name}.sdf" --gen3d -p 7.4 2> "ligands/${name}.log"
done

echo "Done. Verify with: ls ligands/*.sdf | wc -l   (should print 12)"
