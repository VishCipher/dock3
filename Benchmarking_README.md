# Redocking Benchmark

Validates the pipeline's docking accuracy against real, publicly-available
protein-ligand complexes, using the standard "redocking" methodology from
the docking literature: extract a crystal structure's native ligand, dock
it back into its own pocket from scratch, and measure how close the
predicted pose lands to the real crystallographic position.

## Honest methodology notes

- **Not literal PDBbind.** True PDBbind requires a license/account to
  download. This benchmarks against publicly downloadable RCSB structures
  that are widely cited in the docking literature as validation cases —
  same methodology, openly-sourced data. Treat this as a small
  proof-of-concept study (6 complexes), not a full CASF-2016-scale
  benchmark (285 complexes).
- **Ligand identity is auto-detected, not assumed.** Rather than hardcode
  each complex's native ligand residue code from memory — which, if wrong,
  would silently corrupt the result — the script identifies the largest
  non-solvent HETATM group in each structure and **prints what it found**.
  Sanity-check that output before trusting any complex's result.
- **Success threshold**: RMSD < 2.0 Å between the docked pose and the real
  crystal ligand position, the standard threshold used throughout the
  docking literature for "the method found the right pose."
- **Uses only Vina** (not smina/AD4/PandaDock) for now, to keep the first
  version simple. Extending to all engines is a natural next step —
  particularly interesting to see whether AD4/smina redocking accuracy
  differs meaningfully from Vina's.

## Requirements

Same tools as the main pipeline: `obabel`, `vina`, and Meeko
(`mk_prepare_receptor.py`/`mk_prepare_ligand.py`), plus OpenBabel's `obrms`
for RMSD calculation (bundled with OpenBabel — already installed if the
main pipeline works).

## Running it

```bash
cd benchmark
python3 run_pdbbind_benchmark.py
```

Each complex takes a similar amount of time to a single Vina docking run
in the main pipeline (seconds to a couple minutes). Output:

- `benchmark_work/{PDB_ID}/` — intermediate files per complex (receptor,
  native ligand, docked pose) for manual inspection
- `benchmark_results.tsv` — one row per complex: detected ligand, RMSD,
  pass/fail
- A final summary line printed to the terminal: `N/6 complexes redocked
  within 2.0 A RMSD (X% success)`

## Extending this

- Add more PDB IDs to `BENCHMARK_PDB_IDS` in `run_pdbbind_benchmark.py` —
  no other code changes needed, since ligand detection is automatic.
- Add smina/AD4 as additional engines per complex, reporting a per-engine
  success rate for direct comparison.
- If pursuing this for a paper, cross-reference against literature-reported
  Vina redocking success rates on the same PDB IDs, if available, as an
  extra sanity check that this implementation's baseline behavior is
  reasonable.
