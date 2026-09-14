# Troubleshooting Guide

This pipeline has, over its development, hit nearly every category of bug
this kind of tool can produce. Every issue below is a **real bug that
actually occurred**, not a hypothetical — this guide exists so the same
mistakes don't get repeated, by you or anyone else who picks this up.

---

## "My scores don't match" / "these numbers look wrong" — start here

This is the single most common thing that will go wrong, and it is almost
never the docking engine's fault. Work through this checklist **in order**
before assuming a tool is broken.

### Step 1: Is the box actually pointing at the right place?

**This caused more wasted debugging time in this project than every other
bug combined.** A wrong or unfilled box center doesn't error — every engine
will happily dock into empty space and hand back plausible-looking,
completely meaningless numbers.

```bash
grep -A10 "TARGET_NAME:" config.yaml
```

Red flags:
- `center_x: 0.0` / `center_y: 0.0` / `center_z: 0.0` (the placeholder default)
- A `# TODO` comment still present next to the box block
- Coordinates that don't match a known ligand or metal position in the raw PDB

**How to verify a box center is real**, using the metal-ion technique from
this project's own history:
```bash
grep "^HETATM" data/YOUR_RECEPTOR.pdb | awk '$4=="CA"'   # or whatever residue code
```
Cross-check the printed x/y/z against `config.yaml`'s `center_x/y/z` for
that target. If a receptor has multiple copies of the same ion (this
project's 3BH4 has eight separate calcium sites across two chains), make
sure you picked the specific one you meant to.

**If you find a wrong/placeholder box**: fix `config.yaml`, then wipe
**everything** downstream for that target, not just the summary tables —
Snakemake tracks file existence/timestamps, not whether a config value
changed, so stale results from the old (wrong) box can silently persist:
```bash
rm -rf results/TARGET_NAME
snakemake --cores 4 --keep-going
```

### Step 2: Is a tool silently using fallback/default parameters?

Some tools degrade gracefully instead of erroring when they're missing a
proper config file — which is worse than an error, because it looks like a
successful run. Always read the **first ~30 lines** of a job's log, not
just the last few:

```bash
head -30 logs/TARGET_dock_metal_LIGAND.log
```

Known example from this project: PandaDock silently falls back to generic
metal parameters if `--metal-params` isn't supplied, producing a flat
`0.000` score for every ligand regardless of actual chemistry — no error,
just meaningless uniform output. Fixed by pointing `--metal-params` at a
real AutoDock4 parameter file (see `metal_params_file` in `config.yaml`).

**General rule: if every ligand gets the identical score, that is itself
a red flag** — real chemistry rarely produces exact ties across a diverse
library. Treat suspiciously uniform output as a sign the tool never really
differentiated your inputs at all.

### Step 3: Did your file edit actually land where you think it did?

This project hit this exact bug at least three separate times:
- `config.yaml` saved into `workflow/scripts/` instead of the project root
- A corrected script downloaded, but the *old* version never actually
  overwritten locally, so a since-fixed bug kept recurring
- A Snakefile edit made, but a stale copy re-run from a different terminal
  tab/directory

**Before debugging code logic, always verify you're looking at the file you
think you are:**
```bash
pwd                                   # confirm working directory
ls -la config.yaml Snakefile          # confirm they exist HERE, not elsewhere
grep "SOME_STRING_YOU_JUST_ADDED" Snakefile   # confirm the edit is actually present
```

### Step 4: Cross-check with an independent method

Never trust a single engine's score in isolation, especially for anything
metal-related. This pipeline deliberately computes results two independent
ways for exactly this reason:

- `results/{target}/summary/metal_site_scores.tsv` — PandaDock's own
  coordination/binding score
- `results/{target}/summary/metal_geometry_validation.tsv` — an
  independently computed coordination number and closest-contact distance,
  measured directly from atomic coordinates (currently from Vina's pose)

**If the two disagree**, don't assume either is "correct" by default —
figure out *why*:
- If PandaDock shows uniform/zero scores but the geometry table shows real,
  differentiated distances → suspect PandaDock's parameter file (Step 2).
- If the geometry table shows everything at 50+ Å → suspect the box
  (Step 1) — no docking engine can legitimately place a pose that far
  outside its own search box.
- If both agree (known chelators show short distances + good PandaDock
  scores; unrelated compounds show long distances + poor scores) → that
  agreement is real evidence the result is trustworthy.

---

## Other known issues and fixes

### Consensus tier assignments flip between runs for borderline ligands
`compute_consensus_ranking.py` gates a ligand into "real metal contact"
tier 1 based on whether its top-ranked Vina pose has `coordination_number
>= 1`. This is a **single pose from a stochastic search** -- Vina doesn't
guarantee the same top pose across separate runs, especially for a
compound sitting right at the coordination-distance boundary (e.g.
`citric_acid` moved from tier 1 to tier 2 between two full pipeline runs
in this project's own testing, with no code or config change in between).

This is not a bug to "fix" so much as an inherent limitation of judging
metal engagement from one pose. Don't treat a borderline tier assignment
(especially rank changes right at the tier boundary) as final without
checking pose-to-pose consistency. A natural refinement: re-dock each
ligand multiple times (varying Vina's random seed) and require the metal
contact to appear in most/all replicates before trusting tier 1
placement, rather than relying on a single pose.

### AutoGrid4 error: "Unknown receptor type: X"
The receptor PDB contains a HETATM the AD4 force field has no parameters
for — almost always a crystallization artifact (sodium, potassium,
chloride) that shouldn't be there in the first place, not a genuine
cofactor. Fix: `keep_protein_and_metals.py`'s `METAL_CODES` whitelist
should exclude it; the receptor-cleaning step (`prep_receptor_clean`)
should strip it before it ever reaches AutoGrid.

### AutoGrid4 error: "no closestH atom was found"
AutoGrid4 is stricter about hydrogen-placement geometry than Vina/smina.
OpenBabel-added hydrogens don't reliably satisfy this. Fixed in this
pipeline by using Meeko (not OpenBabel) for all receptor/ligand -> PDBQT
conversion — see `prep_receptor`/`prep_ligand` rules.

### A ligand fails with `ModuleNotFoundError` partway through a batch pip install
When you `pip install pkgA pkgB pkgC` as one command and any package fails
to build, **none** of them get installed, even ones that would've worked
fine alone. Check what's actually present:
```bash
pip list | grep -iE "package_name_1|package_name_2"
```
If most/all expected packages are missing, don't patch one
`ModuleNotFoundError` at a time — reinstall everything **individually**
(one `pip install` per package) so a single bad package can't sink the rest.

### A rule's Python script errors with an `AttributeError` on `snakemake.params.X`
The `Snakefile` rule stopped passing that param (e.g. after removing a
feature), but the script file wasn't updated to match — a version mismatch
between the two files. Make sure both the `Snakefile` rule and its
corresponding `workflow/scripts/*.py` file are updated together; they're
tightly coupled, not independent.

### A `RuleException`/`CalledProcessError` mid-run
Read the **log file it points you to**, not just the terminal summary:
```bash
cat logs/TARGET_RULE_LIGAND.log
```
The terminal only tells you *that* something failed; the log tells you *why*.

### General debugging workflow, in order
1. `snakemake -n` (dry run) — confirms the DAG builds at all, catches
   missing files/config errors before spending real compute.
2. If a real run fails, find and read the specific job's log file.
3. Check `config.yaml` values are what you think they are (Step 1/3 above).
4. Check the script referenced by the failing rule actually matches the
   rule's current inputs/params.
5. When in doubt, verify a tool's actual CLI/output by running it directly
   and reading `--help`, rather than assuming — every external tool in
   this pipeline (PandaDock, Meeko, AutoGrid4) had at least one detail that
   turned out different from documentation/memory.
