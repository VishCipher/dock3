# Docking Pipeline (Snakemake) - multi-target, multi-engine, cofactor-aware

Screen a ligand library against **one or more protein targets**, using
**three docking engines** (Vina, smina, AutoDock4), with **automatic
detection and metal-aware handling of functional cofactors** (via
PandaDock) wherever a target's receptor contains one, plus independent
geometry validation, drug-likeness filtering, and run provenance tracking.

DiffDock (blind docking) lives in a **separate, standalone pipeline**
(`diffdock-pipeline/`), its dependencies proved too heavy/fragile to
safely integrate here. See that project's own README for status.

## What this pipeline actually does

For every protein listed under `targets:` in `config.yaml`:

1. Cleans the raw receptor PDB (strips waters/crystallization salts, keeps
   protein + any functional metal cofactor).
2. Docks your entire ligand library against it with **Vina**, **smina**,
   and **AutoDock4**; three independent scoring functions, so you can
   look for agreement across engines rather than trusting a single number.
3. **Auto-detects** whether that target's receptor contains a functional
   metal (Ca, Zn, Mg, Mn, Fe, Cu, Ni, Co) and, if so, **automatically**
   also docks against it with **PandaDock**, which models metal
   coordination geometry explicitly.
4. **Independently re-checks** PandaDock's own coordination-geometry call
   by measuring actual atomic distances from the metal to each ligand's
   top Vina pose, a second, from-scratch verification rather than trusting
   any single engine's internal score.
5. Computes basic **drug-likeness (ADMET) descriptors** for every ligand
   (Lipinski's Rule of Five and extensions), independent of docking.
6. Merges every engine's results per target, then stitches every target
   together into one master **protein x ligand x engine** results matrix.
7. Records **run provenance** (timestamp, targets, ligand count, git
   commit if available) for every run.

## Directory layout

```
docking-pipeline-project/
├── Snakefile                          # all pipeline rules
├── config.yaml                        # targets, ligand_dir, box, engine settings
├── data/
│   ├── 1HSG_receptor.pdb              # one raw PDB per target, named to match config.yaml
│   └── 3BH4_receptor.pdb
├── ligands/                           # one .sdf per ligand — screened against ALL targets
│   ├── EDTA.sdf
│   ├── EGTA.sdf
│   └── ...
├── generate_starter_library.sh        # optional: builds a small test ligand set from SMILES
├── workflow/
│   └── scripts/
│       ├── collect_scores.py              # Vina score parsing
│       ├── collect_smina_scores.py        # smina score parsing
│       ├── collect_ad4_scores.py          # AutoDock4 score parsing
│       ├── merge_engine_scores.py         # combines vina/smina/ad4 per target
│       ├── detect_cofactors.py            # checkpoint: scans receptor for metal ions
│       ├── keep_protein_and_metals.py     # cleans receptor (used by all branches)
│       ├── collect_metal_scores.py        # PandaDock score parsing
│       ├── validate_metal_geometry.py     # independent coordination-geometry check
│       ├── compute_admet.py               # RDKit drug-likeness descriptors
│       ├── write_run_provenance.py        # run metadata capture
│       └── build_protein_ligand_matrix.py # final cross-target matrix
├── logs/                               # one log per job: {target}_{rule}_{ligand}.log
└── results/
    ├── ligands_pdbqt/                  # shared ligand PDBQTs (built once, used by every target)
    ├── {target}/
    │   ├── receptor/
    │   │   ├── receptor_clean.pdb      # cleaned PDB (protein + functional metal only)
    │   │   └── receptor.pdbqt          # Meeko-prepared PDBQT (used by Vina/smina/AD4)
    │   ├── ad4_grids/                  # AutoGrid maps, computed ONCE per target
    │   ├── cofactors/cofactor_report.json  # what metal(s), if any, were detected
    │   ├── docking/
    │   │   ├── vina/{ligand}_out.pdbqt + _log.txt
    │   │   ├── smina/{ligand}_out.pdbqt + _log.txt
    │   │   ├── ad4/{ligand}.dlg
    │   │   └── metal/{ligand}/         # PandaDock output (only if metal detected)
    │   └── summary/
    │       ├── all_scores.tsv                  # Vina, ranked
    │       ├── smina_scores.tsv                # smina, ranked
    │       ├── ad4_scores.tsv                  # AutoDock4, ranked
    │       ├── merged_all_engines.tsv          # vina+smina+ad4 side by side
    │       ├── metal_site_scores.tsv           # only if metal detected
    │       └── metal_geometry_validation.tsv   # only if metal detected
    └── summary/
        ├── protein_ligand_matrix.tsv    # THE master output: every target x ligand x engine
        ├── admet_properties.tsv         # drug-likeness, all ligands, target-independent
        └── run_provenance.json          # what this run actually consisted of
```

## Requirements

Core engines (conda/mamba is easiest), all in **one shared environment**:

```bash
conda create -n docking -c conda-forge -c bioconda snakemake-minimal openbabel autodock-vina smina python=3.11
conda activate docking
pip install meeko
conda install -n docking -c bioconda autodock   # provides autodock4 / autogrid4 binaries
```

**Meeko (not OpenBabel) prepares every receptor/ligand PDBQT the pipeline
actually docks with.** It's built by the same lab that maintains
AutoDock/Vina, specifically to produce PDBQT geometry that AutoGrid4 (used
by the AD4 branch) handles correctly. OpenBabel is still useful to have
around (e.g. `generate_starter_library.sh` uses it for SMILES -> 3D SDF),
but it's no longer in the receptor/ligand -> PDBQT path.

**No MGLTools needed anywhere.** AD4's grid/docking parameter files
(GPF/DPF) are generated directly by the `Snakefile` in plain Python.

### Metal/cofactor branch (auto-enabled per target, no setup flag needed)

```bash
git clone https://github.com/pritampanda15/PandaDock.git
cd PandaDock
pip install -e .
```

Whenever a target's receptor contains a recognized metal (see
`METAL_CODES` in `workflow/scripts/detect_cofactors.py`), that target
automatically gets docked with `pandadock-metal` AND independently
geometry-checked. Nothing to configure per-target.

## Configuring `config.yaml`

```yaml
targets:
  1HSG:
    receptor_pdb: "data/1HSG_receptor.pdb"
    box:
      center_x: 13.073
      center_y: 22.467
      center_z: 5.557
      size_x: 23
      size_y: 18
      size_z: 17
  # Add more targets the same way — each needs its own receptor_pdb + box.

ligand_dir: "ligands"
exhaustiveness: 8            # Vina/smina search thoroughness

ad4_ligand_types: "C,A,N,NA,OA,SA,HD,F,Cl,Br,I,S,P"
ad4_ga_run: 10                # AutoDockTools' standard default GA-LS settings
ad4_ga_num_evals: 2500000     # lower both for a faster, less thorough screen
```

**The box coordinates are the single most important — and most
error-prone — setting in this whole file.** Vina/smina/AD4/PandaDock all
search *only* within that box; if the center is wrong (even a placeholder
like `0,0,0` left unfilled), every engine will happily dock into empty
space and hand back plausible-looking but meaningless numbers, often with
no error at all. Compute real coordinates from a known ligand or metal
ion's position in the receptor PDB — never leave a target's box as an
unfilled placeholder.

## Running it

```bash
cd docking-pipeline-project

snakemake -n                          # always dry-run first
snakemake --cores 8 --keep-going      # real run
```

`--keep-going` is worth keeping on by default: one problematic ligand
shouldn't block docking for the rest of your library.

If a run partially fails, check `logs/{target}_{rule}_{ligand}.log` for
that specific job.

## Reading the results

Start with `results/summary/protein_ligand_matrix.tsv`:

```python
import pandas as pd
df = pd.read_csv("results/summary/protein_ligand_matrix.tsv", sep="\t")
df[df.engine == "vina"].pivot(index="ligand_id", columns="target", values="score_kcal_mol")
```

**Vina/smina/AD4 scores are not on directly comparable absolute scales**
(different scoring functions). Use the matrix to look for *agreement in
ranking* across engines, not to average the numbers together.

For metal sites, cross-reference two independent checks:
- `results/{target}/summary/metal_site_scores.tsv` — PandaDock's own
  coordination/binding scores
- `results/{target}/summary/metal_geometry_validation.tsv` - an
  independently computed coordination number and distance check from
  Vina's pose

Agreement between the two is real evidence a hit is chemically sensible;
disagreement tells you exactly where to look closer. Known chelators
(EDTA, EGTA, citric acid, if in your library) should score meaningfully
better than unrelated compounds on both checks, if everything comes back
identical or flagged, re-check that target's box coordinates first.

Check `results/summary/admet_properties.tsv` before getting attached to
any hit — a strong docking score means little if `lipinski_violations` is
high.

## If something looks wrong

See **`TROUBLESHOOTING.md`** — a checklist built directly from real bugs
hit during this project's development (wrong box coordinates, silent
fallback parameters, files landing in the wrong folder, and more). Work
through it before assuming a docking engine itself is broken; in every
case so far, it wasn't.

## Known limitations / honest caveats

- `IDEAL_DISTANCE_RANGES` and `TYPICAL_COORDINATION_NUMBERS` in
  `validate_metal_geometry.py` are literature-informed defaults, not
  calibrated per target, treat flags as "worth a closer look," not a
  certified verdict.
- Geometry validation currently only checks **Vina's** pose. Extending to
  smina/AD4 poses is a natural next step.
- PandaDock's coordination scores were computed using **fallback/generic
  metal parameters** in testing (its own log warns about this) the
  independent geometry check exists specifically because that score alone
  isn't fully trustworthy yet.

## Next steps to extend

- **PLIP interaction fingerprinting**: broader than the metal-geometry
  filter - automatically tabulate H-bonds, salt bridges, hydrophobic
  contacts, π-stacking for any pocket, not just metal sites.
- **Pharmacophore pre-filter**: cheaply screen a large library against a
  pharmacophore built from known actives before docking everything blind.
- **Auto-generated visual reports**: PyMOL/ChimeraX session files or
  rendered images per top hit - right now every result is a TSV, nothing
  shows you what the pose actually looks like.
- **Regression/self-test suite**: bake in validated cases (1HSG/indinavir
  RMSD, a hand-checked calcium site) as an automated check on every change.
- **Multi-site screening**: if a target has multiple distinct metal ions
  (e.g. 3BH4's chain A alone has four separate calcium sites), add each as
  its own `targets:` entry with a different box center.
- **Extend geometry validation to smina/AD4 poses**, not just Vina.
- **MD refinement of top hits**: dock everything cheaply first, then
  short-MD-refine only the top 3-5 per target on a GPU workstation/Colab
  (same two-stage local-prep / external-run pattern as the DiffDock
  pipeline, once GPU access is available again).
