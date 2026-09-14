# Docking Pipeline (Snakemake) — multi-target, multi-engine

Screen a ligand library against **one or more protein targets** using Vina,
smina, and AutoDock4 in parallel, with automatic detection and handling of
metal-ion cofactors (via PandaDock) per target.

## Directory layout

```
docking-pipeline/
├── Snakefile
├── config.yaml            # define your targets + box + ligand_dir here
├── data/
│   ├── 1HSG_receptor.pdb  # one PDB file per target, named to match config.yaml
│   └── 3BH4_receptor.pdb
├── ligands/*.sdf           # one SDF per ligand — screened against ALL targets
└── workflow/scripts/
```

## Requirements

Install these (conda/mamba is easiest):

```bash
conda create -n docking -c conda-forge -c bioconda snakemake-minimal openbabel autodock-vina smina python=3.11
conda activate docking
```

### AutoDock4 setup

No MGLTools needed. The pipeline generates AD4's grid/docking parameter
files (GPF/DPF) itself directly in Python — those are just documented text
config files, and MGLTools' `prepare_gpf4.py`/`prepare_dpf4.py` scripts (and
their finicky bundled Python 2 interpreter) aren't actually required to
produce them. You only need the `autodock4`/`autogrid4` binaries themselves:

```bash
conda install -n docking -c bioconda autodock
```

Two settings in `config.yaml` control the search thoroughness AD4 uses
(`ad4_ga_run`, `ad4_ga_num_evals`) — the defaults match AutoDockTools'
standard settings; lower them for a faster, less thorough screen.

### Metal/cofactor branch (optional, auto-enabled per target)

If a target's receptor contains a metal ion (Ca, Zn, Mg, Fe, Mn, Cu, Ni, Co), that
target automatically also gets routed through PandaDock's `pandadock-metal`
command, since Vina/smina/AD4 don't model metal coordination geometry. Install
separately:

```bash
git clone https://github.com/pritampanda15/PandaDock.git
cd PandaDock
pip install -e .
```

**Before trusting `collect_metal_scores.py`:** run `pandadock-metal --help`
(or dock a validated case by hand once) and check what output file/columns
it actually produces — the parsing script is a template.

## Running it

1. Add each target to `config.yaml` under `targets:` — a name, its receptor
   PDB path, and its docking box. Put the actual PDB file at that path under
   `data/`.
2. Drop ligand SDFs into `ligands/` — screened against every target listed.
3. Fill in `mgltools_pythonsh` / `mgltools_utils_dir`.
4. Dry-run to check the plan:
   ```bash
   snakemake -n
   ```
5. Run for real:
   ```bash
   snakemake --cores 8
   ```
   On a SLURM cluster instead:
   ```bash
   snakemake --executor slurm --jobs 100 --default-resources slurm_partition=<partition>
   ```

Results land in:
- `results/{target}/summary/all_scores.tsv` / `smina_scores.tsv` / `ad4_scores.tsv` —
  per-engine, per-target ranked tables
- `results/{target}/summary/merged_all_engines.tsv` — the three engines side by
  side for that one target
- `results/{target}/summary/metal_site_scores.tsv` — only if that target's
  receptor contains a detected metal ion
- **`results/summary/protein_ligand_matrix.tsv`** — the actual screening
  matrix: every target x every ligand x every engine, in one long-format
  table (`target, ligand_id, engine, score_kcal_mol`). Pivot it in
  pandas/Excel to get a wide ligand-by-protein view per engine.

**Note on comparing engines:** Vina/smina/AD4 scores are not on directly
comparable absolute scales (different scoring functions). Use the matrix to
look for *agreement in ranking* across engines and across targets, not to
average the numbers together.

## Next steps to extend

- Add a cavity-detection step (e.g. fpocket) upstream of `prep_receptor` for
  targets where you don't already know the binding site coordinates.
- Containerize each engine (Docker locally / Apptainer on HPC), so
  dependencies never collide as you add more targets/engines.
- If a target has *more than one* distinct metal site, extend
  `dock_metal_site` to loop over all entries in that target's
  `cofactor_report.json` (it currently only docks against the first metal
  found per target) — add a `metal_site` wildcard alongside `ligand`.
- Add a post-docking geometry-validation step: check each metal-site pose's
  coordination distances/angles against known ideal values for that metal,
  and flag poses that score well but violate real coordination chemistry.
- Bake in a regression test using your already-validated cases (1HSG/indinavir,
  3BH4/EDTA-Ca): re-dock them automatically whenever you change a rule, and
  assert RMSD-to-crystal stays low.
