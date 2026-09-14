"""
Multi-target docking pipeline (Snakemake) with cofactor-aware branching
------------------------------------------------------------------------
Every protein target gets its own receptor PDB + docking box, defined under
config["targets"]. The SAME ligand library is screened against ALL targets
in one pipeline invocation, producing a full protein x ligand result matrix.

Per target (wildcard {target}), standard branch (always runs):
  1. prep_receptor  -> convert receptor PDB to PDBQT (Meeko)
  2. dock_vina / dock_smina  -> dock each ligand (scatter over {ligand})
  3. ad4_prep_grids  -> AutoGrid maps, computed ONCE per target (not per ligand)
  4. dock_ad4        -> AutoDock4 per ligand, against that target's shared maps
  5. collect_*_scores + merge_engine_scores -> per-target comparison table

AD4's grid/dock parameter files (GPF/DPF) are generated directly by this
Snakefile in plain Python -- NOT via MGLTools' prepare_gpf4.py/prepare_dpf4.py.
Those scripts just write a documented text format; MGLTools' finicky Python 2
install isn't actually required, only the `autodock4`/`autogrid4` binaries
(installable via conda, see README).

Ligand prep (prep_ligand) is target-independent and runs ONCE regardless of
how many targets you screen against, since the ligand's own 3D structure
doesn't depend on which protein it's being docked into.

Per target, metal/cofactor branch (auto-enabled only if that target's
receptor contains a metal ion):
  A. detect_cofactors     -> checkpoint, scans for known metal ions
  B. dock_metal_site       -> PandaDock's pandadock-metal, per ligand, using
                              the SAME cleaned receptor as the standard branch
  C. collect_metal_scores  -> per-target metal-site score table

Finally, build_protein_ligand_matrix concatenates every target's merged
engine table into one long-format results/summary/protein_ligand_matrix.tsv
(columns: target, ligand_id, engine, score_kcal_mol).

NOTE on pandadock-metal: verify its exact output filename/format for your
installed version before trusting collect_metal_scores.py blindly -- the
parsing there is a template.

NOTE on the hand-written DPF: it uses AutoDockTools' standard default GA-LS
search parameters (ga_run 10, ga_num_evals 2500000). Override via
config["ad4_ga_run"] / config["ad4_ga_num_evals"] if you want a faster,
lower-thoroughness screen.

DiffDock (blind docking) lives in a SEPARATE, standalone pipeline
(../diffdock-pipeline/), not here -- its dependencies proved heavy and
fragile enough that keeping it fully decoupled from this working pipeline
was safer than integrating it as an opt-in branch.
"""


configfile: "config.yaml"

import glob
import json
import os
import shutil
import subprocess

TARGETS = config["targets"]
TARGET_IDS = list(TARGETS.keys())

if not TARGET_IDS:
    raise ValueError("config['targets'] is empty -- add at least one protein target.")

LIGAND_DIR = config["ligand_dir"]
LIGAND_IDS = sorted(
    os.path.splitext(os.path.basename(f))[0]
    for f in glob.glob(f"{LIGAND_DIR}/*.sdf")
)

if not LIGAND_IDS:
    raise ValueError(
        f"No .sdf files found in '{LIGAND_DIR}'. "
        "Add ligand SDFs there or update ligand_dir in config.yaml."
    )

AD4_LIGAND_TYPES = config.get("ad4_ligand_types", "C,A,N,NA,OA,SA,HD,F,Cl,Br,I,S,P").split(",")
AD4_GA_RUN = config.get("ad4_ga_run", 10)
AD4_GA_NUM_EVALS = config.get("ad4_ga_num_evals", 2500000)


def _box(wc, key):
    return TARGETS[wc.target]["box"][key]


def _even_npts(size_angstrom, spacing=0.375):
    """AutoGrid wants an even number of grid points per axis."""
    n = int(size_angstrom / spacing)
    return n if n % 2 == 0 else n + 1


def _parse_pdbqt_atom_types(pdbqt_path):
    """Read the AD4 atom-type column (last whitespace-separated field) off every ATOM/HETATM line."""
    types = set()
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if parts:
                    types.add(parts[-1])
    return sorted(types)


def _parse_ligand_pdbqt_info(pdbqt_path):
    """Return (centroid_x, centroid_y, centroid_z, torsdof) from a ligand PDBQT."""
    coords = []
    torsdof = 0
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            elif line.startswith("TORSDOF"):
                torsdof = int(line.split()[-1])
    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    return cx, cy, cz, torsdof


def _write_gpf(gpf_path, receptor_pdbqt_name, receptor_types, ligand_types, center, npts,
               spacing=0.375, smooth=0.5, dielectric=-0.1465):
    prefix = os.path.splitext(receptor_pdbqt_name)[0]
    lines = [
        f"npts {npts[0]} {npts[1]} {npts[2]}",
        f"gridfld {prefix}.maps.fld",
        f"spacing {spacing}",
        f"receptor_types {' '.join(receptor_types)}",
        f"ligand_types {' '.join(ligand_types)}",
        f"receptor {receptor_pdbqt_name}",
        f"gridcenter {center[0]} {center[1]} {center[2]}",
        f"smooth {smooth}",
    ]
    for lt in ligand_types:
        lines.append(f"map {prefix}.{lt}.map")
    lines.append(f"elecmap {prefix}.e.map")
    lines.append(f"dsolvmap {prefix}.d.map")
    lines.append(f"dielectric {dielectric}")
    with open(gpf_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_dpf(dpf_path, receptor_prefix, ligand_pdbqt_name, ligand_types, torsdof, about,
               ga_run=10, ga_num_evals=2500000):
    lines = [
        "autodock_parameter_version 4.2",
        "outlev 1",
        "intelec",
        "seed pid time",
        f"ligand_types {' '.join(ligand_types)}",
        f"fld {receptor_prefix}.maps.fld",
    ]
    for lt in ligand_types:
        lines.append(f"map {receptor_prefix}.{lt}.map")
    lines.append(f"elecmap {receptor_prefix}.e.map")
    lines.append(f"desolvmap {receptor_prefix}.d.map")
    lines.append(f"move {ligand_pdbqt_name}")
    lines.append(f"about {about[0]} {about[1]} {about[2]}")
    lines += [
        "tran0 random", "quat0 random", "dihe0 random",
        "tstep 2.0", "qstep 50.0", "dstep 50.0",
        f"torsdof {torsdof}",
        "rmstol 2.0", "extnrg 1000.0", "e0max 0.0 10000",
        "ga_pop_size 150", f"ga_num_evals {ga_num_evals}", "ga_num_generations 27000",
        "ga_elitism 1", "ga_mutation_rate 0.02", "ga_crossover_rate 0.8",
        "ga_window_size 10", "ga_cauchy_alpha 0.0", "ga_cauchy_beta 1.0",
        "set_ga",
        "sw_max_its 300", "sw_max_succ 4", "sw_max_fail 4",
        "sw_rho 1.0", "sw_lb_rho 0.01", "ls_search_freq 0.06",
        "set_psw1",
        "unbound_model bound",
        f"ga_run {ga_run}",
        "analysis",
    ]
    with open(dpf_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def get_final_targets(wildcards):
    """
    For every protein target: always want its standard-engine comparison
    table; additionally want its metal-site table if that target's
    detect_cofactors checkpoint found a metal ion. Finally, want the
    cross-target matrix that stitches every target's results together.
    """
    targets_needed = []
    for target in TARGET_IDS:
        targets_needed.append(f"results/{target}/summary/merged_all_engines.tsv")
        targets_needed.append(f"results/{target}/summary/plip_interactions.tsv")
        targets_needed.append(f"results/{target}/summary/consensus_ranking.tsv")

        report_path = checkpoints.detect_cofactors.get(target=target).output.report
        with open(report_path) as f:
            cofactor_report = json.load(f)
        if cofactor_report["has_metal"]:
            targets_needed.append(f"results/{target}/summary/metal_site_scores.tsv")
            targets_needed.append(f"results/{target}/summary/metal_geometry_validation.tsv")

    targets_needed.append("results/summary/protein_ligand_matrix.tsv")
    targets_needed.append("results/summary/admet_properties.tsv")
    targets_needed.append("results/summary/run_provenance.json")
    return targets_needed


rule all:
    input:
        get_final_targets


# =============================================================================
# Ligand prep (target-independent, runs once)
# =============================================================================

rule prep_ligand:
    """
    Convert one ligand SDF -> PDBQT using Meeko (mk_prepare_ligand.py)
    instead of OpenBabel. Meeko is built by the same lab that maintains
    AutoDock/Vina, specifically to produce PDBQT geometry that AutoGrid4
    handles correctly -- OpenBabel's hydrogen placement, while chemically
    valid, doesn't reliably satisfy AutoGrid4's stricter tolerances (see
    the "no closestH atom was found" errors this replaced).
    Shared across all targets, same as before.
    """
    input:
        sdf=os.path.join(LIGAND_DIR, "{ligand}.sdf")
    output:
        pdbqt="results/ligands_pdbqt/{ligand}.pdbqt"
    log:
        "logs/prep_ligand_{ligand}.log"
    shell:
        """
        mkdir -p results/ligands_pdbqt logs
        mk_prepare_ligand.py -i {input.sdf} -o {output.pdbqt} > {log} 2>&1
        """


# =============================================================================
# Per-target receptor prep
# =============================================================================

rule prep_receptor_clean:
    """
    Strip crystallization artifacts (waters, buffer ions, native co-crystallized
    ligands) from the raw receptor PDB, keeping protein atoms plus any
    recognized FUNCTIONAL metal cofactor (Ca, Zn, Mg, Mn, Fe, Cu, ...).

    This feeds BOTH the standard branch (via prep_receptor below) and the
    metal branch (dock_metal_site) -- every engine sees the same cleaned
    structure. Without this step, incidental HETATMs in the raw PDB (e.g. a
    crystallization-buffer Na+ ion) can reach AutoGrid4, which has no force
    field parameters for most non-functional ions and will hard-fail with
    "Unknown receptor type".
    """
    input:
        receptor=lambda wc: TARGETS[wc.target]["receptor_pdb"]
    output:
        pdb="results/{target}/receptor/receptor_clean.pdb"
    script:
        "workflow/scripts/keep_protein_and_metals.py"


rule prep_receptor:
    """
    Convert this target's CLEANED receptor PDB -> PDBQT using Meeko
    (mk_prepare_receptor.py) instead of OpenBabel, for the same
    AutoGrid4-compatibility reason as prep_ligand above.
    """
    input:
        receptor="results/{target}/receptor/receptor_clean.pdb"
    output:
        pdbqt="results/{target}/receptor/receptor.pdbqt"
    log:
        "logs/{target}_prep_receptor.log"
    shell:
        """
        mkdir -p results/{wildcards.target}/receptor logs
        mk_prepare_receptor.py --read_pdb {input.receptor} -p {output.pdbqt} > {log} 2>&1
        """


# =============================================================================
# Standard engines: Vina, smina
# =============================================================================

rule dock_vina:
    """Dock one ligand against one target's receptor with AutoDock Vina."""
    input:
        receptor="results/{target}/receptor/receptor.pdbqt",
        ligand="results/ligands_pdbqt/{ligand}.pdbqt"
    output:
        pose="results/{target}/docking/vina/{ligand}_out.pdbqt",
        log="results/{target}/docking/vina/{ligand}_log.txt"
    params:
        cx=lambda wc: _box(wc, "center_x"),
        cy=lambda wc: _box(wc, "center_y"),
        cz=lambda wc: _box(wc, "center_z"),
        sx=lambda wc: _box(wc, "size_x"),
        sy=lambda wc: _box(wc, "size_y"),
        sz=lambda wc: _box(wc, "size_z"),
        exhaustiveness=config.get("exhaustiveness", 8)
    shell:
        """
        mkdir -p results/{wildcards.target}/docking/vina
        vina --receptor {input.receptor} --ligand {input.ligand} \
             --center_x {params.cx} --center_y {params.cy} --center_z {params.cz} \
             --size_x {params.sx} --size_y {params.sy} --size_z {params.sz} \
             --exhaustiveness {params.exhaustiveness} \
             --out {output.pose} > {output.log} 2>&1
        """


rule collect_scores:
    """Parse this target's Vina logs, rank, write TSV."""
    input:
        logs=lambda wc: expand(
            "results/{target}/docking/vina/{ligand}_log.txt",
            target=wc.target, ligand=LIGAND_IDS
        )
    output:
        "results/{target}/summary/all_scores.tsv"
    script:
        "workflow/scripts/collect_scores.py"


rule dock_smina:
    """Dock one ligand against one target's receptor with smina."""
    input:
        receptor="results/{target}/receptor/receptor.pdbqt",
        ligand="results/ligands_pdbqt/{ligand}.pdbqt"
    output:
        pose="results/{target}/docking/smina/{ligand}_out.pdbqt",
        log="results/{target}/docking/smina/{ligand}_log.txt"
    params:
        cx=lambda wc: _box(wc, "center_x"),
        cy=lambda wc: _box(wc, "center_y"),
        cz=lambda wc: _box(wc, "center_z"),
        sx=lambda wc: _box(wc, "size_x"),
        sy=lambda wc: _box(wc, "size_y"),
        sz=lambda wc: _box(wc, "size_z"),
        exhaustiveness=config.get("exhaustiveness", 8)
    shell:
        """
        mkdir -p results/{wildcards.target}/docking/smina
        smina --receptor {input.receptor} --ligand {input.ligand} \
              --center_x {params.cx} --center_y {params.cy} --center_z {params.cz} \
              --size_x {params.sx} --size_y {params.sy} --size_z {params.sz} \
              --exhaustiveness {params.exhaustiveness} \
              --out {output.pose} > {output.log} 2>&1
        """


rule collect_smina_scores:
    """Parse this target's smina logs, rank, write TSV."""
    input:
        logs=lambda wc: expand(
            "results/{target}/docking/smina/{ligand}_log.txt",
            target=wc.target, ligand=LIGAND_IDS
        )
    output:
        "results/{target}/summary/smina_scores.tsv"
    script:
        "workflow/scripts/collect_smina_scores.py"


# =============================================================================
# AutoDock4 -- grid maps computed once per target, GPF/DPF written directly
# (no MGLTools). Only needs the autodock4/autogrid4 binaries on PATH.
# =============================================================================

rule ad4_prep_grids:
    """
    Compute AutoGrid maps once for this target's receptor+box. Writes the
    GPF ourselves (scanning the receptor PDBQT for atom types actually
    present) instead of calling MGLTools' prepare_gpf4.py.
    """
    input:
        receptor="results/{target}/receptor/receptor.pdbqt"
    output:
        gpf="results/{target}/ad4_grids/receptor.gpf",
        glg="results/{target}/ad4_grids/receptor.glg",
        marker="results/{target}/ad4_grids/.grids_done"
    params:
        cx=lambda wc: _box(wc, "center_x"),
        cy=lambda wc: _box(wc, "center_y"),
        cz=lambda wc: _box(wc, "center_z"),
        npx=lambda wc: _even_npts(_box(wc, "size_x")),
        npy=lambda wc: _even_npts(_box(wc, "size_y")),
        npz=lambda wc: _even_npts(_box(wc, "size_z")),
    run:
        grid_dir = f"results/{wildcards.target}/ad4_grids"
        os.makedirs(grid_dir, exist_ok=True)
        os.makedirs("logs", exist_ok=True)

        local_receptor = os.path.join(grid_dir, "receptor.pdbqt")
        shutil.copyfile(input.receptor, local_receptor)

        receptor_types = _parse_pdbqt_atom_types(local_receptor)

        _write_gpf(
            gpf_path=output.gpf,
            receptor_pdbqt_name="receptor.pdbqt",
            receptor_types=receptor_types,
            ligand_types=AD4_LIGAND_TYPES,
            center=(params.cx, params.cy, params.cz),
            npts=(params.npx, params.npy, params.npz),
        )

        logfile = os.path.abspath(f"logs/{wildcards.target}_ad4_prep_grids.log")
        with open(logfile, "w") as lf:
            subprocess.run(
                ["autogrid4", "-p", "receptor.gpf", "-l", "receptor.glg"],
                cwd=grid_dir, stdout=lf, stderr=subprocess.STDOUT, check=True
            )

        with open(output.marker, "w") as f:
            f.write("done\n")


rule dock_ad4:
    """
    Dock one ligand against this target's precomputed AD4 grid maps. Writes
    the DPF ourselves (reading ligand atom types, TORSDOF, and centroid
    straight out of the ligand PDBQT) instead of calling MGLTools'
    prepare_dpf4.py.
    """
    input:
        ligand="results/ligands_pdbqt/{ligand}.pdbqt",
        grids_marker="results/{target}/ad4_grids/.grids_done"
    output:
        dlg="results/{target}/docking/ad4/{ligand}.dlg"
    params:
        ga_run=AD4_GA_RUN,
        ga_num_evals=AD4_GA_NUM_EVALS
    run:
        grid_dir = f"results/{wildcards.target}/ad4_grids"
        os.makedirs(f"results/{wildcards.target}/docking/ad4", exist_ok=True)
        os.makedirs("logs", exist_ok=True)

        local_ligand_name = f"{wildcards.ligand}.pdbqt"
        local_ligand = os.path.join(grid_dir, local_ligand_name)
        shutil.copyfile(input.ligand, local_ligand)

        ligand_types = _parse_pdbqt_atom_types(local_ligand)
        cx, cy, cz, torsdof = _parse_ligand_pdbqt_info(local_ligand)

        dpf_name = f"{wildcards.ligand}.dpf"
        dlg_name = f"{wildcards.ligand}.dlg"

        _write_dpf(
            dpf_path=os.path.join(grid_dir, dpf_name),
            receptor_prefix="receptor",
            ligand_pdbqt_name=local_ligand_name,
            ligand_types=ligand_types,
            torsdof=torsdof,
            about=(cx, cy, cz),
            ga_run=params.ga_run,
            ga_num_evals=params.ga_num_evals,
        )

        logfile = os.path.abspath(f"logs/{wildcards.target}_dock_ad4_{wildcards.ligand}.log")
        with open(logfile, "w") as lf:
            subprocess.run(
                ["autodock4", "-p", dpf_name, "-l", dlg_name],
                cwd=grid_dir, stdout=lf, stderr=subprocess.STDOUT, check=True
            )

        shutil.copyfile(os.path.join(grid_dir, dlg_name), output.dlg)


rule collect_ad4_scores:
    """Parse this target's AD4 .dlg files, rank, write TSV."""
    input:
        dlgs=lambda wc: expand(
            "results/{target}/docking/ad4/{ligand}.dlg",
            target=wc.target, ligand=LIGAND_IDS
        )
    output:
        "results/{target}/summary/ad4_scores.tsv"
    script:
        "workflow/scripts/collect_ad4_scores.py"


rule merge_engine_scores:
    """Combine this target's Vina/smina/AD4 tables into one long-format table."""
    input:
        vina="results/{target}/summary/all_scores.tsv",
        smina="results/{target}/summary/smina_scores.tsv",
        ad4="results/{target}/summary/ad4_scores.tsv"
    output:
        "results/{target}/summary/merged_all_engines.tsv"
    script:
        "workflow/scripts/merge_engine_scores.py"


# =============================================================================
# Metal / cofactor branch (per target, auto-enabled)
# =============================================================================

checkpoint detect_cofactors:
    """Scan this target's receptor PDB for known metal-ion HETATM records."""
    input:
        receptor=lambda wc: TARGETS[wc.target]["receptor_pdb"]
    output:
        report="results/{target}/cofactors/cofactor_report.json"
    script:
        "workflow/scripts/detect_cofactors.py"


rule dock_metal_site:
    """
    Dock one ligand against this target's metal site with PandaDock.

    Uses the SAME cleaned receptor PDB as the standard branch
    (prep_receptor_clean) -- no separate metal-specific receptor-prep rule
    needed, since both branches want identical cleaning (protein + functional
    metal, artifacts stripped).

    Verified against the installed CLI (`pandadock-metal dock --help`):
    it's a subcommand (`dock`), takes -r/-l for receptor/ligand, --center
    and --box for the grid, and auto-detects the metal directly from the
    receptor PDB -- there is no --metal-type/--metal-residue/
    --coordination-geometry flag, unlike what an earlier, unverified
    version of this rule assumed.

    Passes --metal-params pointing at AD4.1_bound.dat, the real official
    AutoDock4 atomic parameter file (from ccsb-scripps/AutoDock4 on GitHub,
    GPL-licensed, includes verified Ca/Zn/Mg/Mn/Fe entries) -- without it,
    PandaDock silently falls back to generic approximations and its
    avg_metal_binding_score/best_coordination_score come back as a flat
    0.000 for every ligand regardless of chemistry (confirmed: EDTA and
    aspirin scored identically). The independent geometry check in
    validate_metal_geometry.py exists specifically because of this gap --
    worth re-comparing the two once this fix is in.
    """
    input:
        receptor="results/{target}/receptor/receptor_clean.pdb",
        ligand=os.path.join(LIGAND_DIR, "{ligand}.sdf"),
        report="results/{target}/cofactors/cofactor_report.json"
    output:
        outdir=directory("results/{target}/docking/metal/{ligand}")
    params:
        cx=lambda wc: _box(wc, "center_x"),
        cy=lambda wc: _box(wc, "center_y"),
        cz=lambda wc: _box(wc, "center_z"),
        sx=lambda wc: _box(wc, "size_x"),
        sy=lambda wc: _box(wc, "size_y"),
        sz=lambda wc: _box(wc, "size_z"),
        metal_params=config.get("metal_params_file", "AD4.1_bound.dat"),
    log:
        "logs/{target}_dock_metal_{ligand}.log"
    shell:
        """
        mkdir -p {output.outdir} logs
        pandadock-metal dock \
            -r {input.receptor} \
            -l {input.ligand} \
            --center {params.cx} {params.cy} {params.cz} \
            --box {params.sx} {params.sy} {params.sz} \
            --metal-params {params.metal_params} \
            --output-dir {output.outdir} > {log} 2>&1
        """


rule collect_metal_scores:
    """
    Gather this target's PandaDock metal-site results into a ranked TSV.
    Parses the real, differentiating search energy from each job's log
    file (see collect_metal_scores.py docstring) rather than trusting
    metal_docking_summary.json's own binding/coordination score fields,
    which were confirmed to come back as a flat 0.000 regardless of
    chemistry even after supplying a real AutoDock4 parameter file.
    """
    input:
        outdirs=lambda wc: expand(
            "results/{target}/docking/metal/{ligand}",
            target=wc.target, ligand=LIGAND_IDS
        ),
        logs=lambda wc: expand(
            "logs/{target}_dock_metal_{ligand}.log",
            target=wc.target, ligand=LIGAND_IDS
        )
    output:
        "results/{target}/summary/metal_site_scores.tsv"
    script:
        "workflow/scripts/collect_metal_scores.py"


rule validate_metal_geometry:
    """
    Independently verify metal coordination geometry from each ligand's
    top-ranked Vina pose, rather than trusting any single engine's internal
    score alone. See validate_metal_geometry.py's docstring for method and
    the important caveat that this only checks Vina's pose so far.
    """
    input:
        report="results/{target}/cofactors/cofactor_report.json",
        poses=lambda wc: expand(
            "results/{target}/docking/vina/{ligand}_out.pdbqt",
            target=wc.target, ligand=LIGAND_IDS
        )
    output:
        "results/{target}/summary/metal_geometry_validation.tsv"
    script:
        "workflow/scripts/validate_metal_geometry.py"


# =============================================================================
# Cross-target aggregation: the actual protein x ligand matrix
# =============================================================================

rule build_protein_ligand_matrix:
    """Concatenate every target's merged engine table into one master matrix."""
    input:
        merged=expand("results/{target}/summary/merged_all_engines.tsv", target=TARGET_IDS)
    output:
        "results/summary/protein_ligand_matrix.tsv"
    script:
        "workflow/scripts/build_protein_ligand_matrix.py"


# =============================================================================
# ADMET / drug-likeness filter (target-independent, runs once)
# =============================================================================

rule compute_admet:
    """
    Compute basic drug-likeness descriptors (Lipinski's Rule of Five plus a
    couple of common extensions) for every ligand in the library. Pure
    RDKit -- already a dependency via Meeko/PandaDock, no new install.
    Independent of which targets are configured; runs once regardless of
    how many proteins you're screening against.
    """
    input:
        ligands=expand(os.path.join(LIGAND_DIR, "{ligand}.sdf"), ligand=LIGAND_IDS)
    output:
        "results/summary/admet_properties.tsv"
    script:
        "workflow/scripts/compute_admet.py"


# =============================================================================
# Run provenance (cheap reproducibility insurance)
# =============================================================================

rule write_run_provenance:
    """
    Record what this run actually consisted of: timestamp, which targets
    and ligands were involved, which engine branches were active, and the
    git commit this Snakefile was at (if the project is under git) --
    given how many config/file-sync bugs came up during development, being
    able to answer "which exact version of the pipeline produced this
    results folder" is worth having for free.
    """
    input:
        "results/summary/protein_ligand_matrix.tsv"
    output:
        "results/summary/run_provenance.json"
    params:
        targets=TARGET_IDS,
        ligand_ids=LIGAND_IDS
    script:
        "workflow/scripts/write_run_provenance.py"


# =============================================================================
# PLIP interaction fingerprinting (in progress -- see TROUBLESHOOTING.md)
# =============================================================================
# prepare_plip_complex is fully built and safe to run now (pure Python +
# OpenBabel, no external tool schema to guess). The actual `plip` call and
# its output-parsing script are deliberately NOT added yet -- every other
# external tool integrated into this pipeline (PandaDock, Meeko, AutoGrid4)
# had at least one detail that turned out different from documentation, so
# `plip --help` gets verified against the real installed version first.

rule prepare_plip_complex:
    """
    Combine this target's cleaned receptor with one ligand's top-ranked
    Vina pose into a single PDB file -- the input format PLIP (and most
    interaction-analysis tools) expect: one file, ligand as its own HETATM
    residue, not two separate files.
    """
    input:
        receptor="results/{target}/receptor/receptor_clean.pdb",
        pose="results/{target}/docking/vina/{ligand}_out.pdbqt"
    output:
        "results/{target}/plip_complexes/{ligand}_complex.pdb"
    script:
        "workflow/scripts/prepare_complex_pdb.py"


rule run_plip:
    """
    Run PLIP on one target/ligand's combined complex, generating both TXT
    and XML interaction reports. Verified against a real installed run
    (`plip --help`, PLIP 3.0.1): an earlier attempt without -t/-x produced
    no report at all, only PLIP's internal "fixed" intermediate PDB files
    -- -t/-x are what actually trigger report generation.
    """
    input:
        complex_pdb="results/{target}/plip_complexes/{ligand}_complex.pdb"
    output:
        outdir=directory("results/{target}/plip/{ligand}")
    log:
        "logs/{target}_plip_{ligand}.log"
    shell:
        """
        mkdir -p {output.outdir} logs
        plip -f {input.complex_pdb} -o {output.outdir} -x -t --name plip_report > {log} 2>&1
        """


rule collect_plip_interactions:
    """
    Summarize PLIP's interaction fingerprint for every ligand into one
    ranked TSV, sorted by total interaction count. See
    collect_plip_interactions.py's docstring for the real finding that
    direct_metal_coordination came back 0 even for EDTA in initial testing
    -- consistent with the independent geometry validator's
    LOW_COORDINATION result for the same Vina pose.
    """
    input:
        outdirs=lambda wc: expand(
            "results/{target}/plip/{ligand}",
            target=wc.target, ligand=LIGAND_IDS
        )
    output:
        "results/{target}/summary/plip_interactions.tsv"
    script:
        "workflow/scripts/collect_plip_interactions.py"


# =============================================================================
# Consensus ranking across all independently-computed metrics
# =============================================================================

def consensus_inputs(wildcards):
    """
    Always includes the three standard engines + PLIP. Adds the metal-site
    energy and independent geometry check only for targets where
    detect_cofactors actually found a cofactor -- same checkpoint pattern
    used by get_final_targets.
    """
    inputs = {
        "vina": f"results/{wildcards.target}/summary/all_scores.tsv",
        "smina": f"results/{wildcards.target}/summary/smina_scores.tsv",
        "ad4": f"results/{wildcards.target}/summary/ad4_scores.tsv",
        "plip": f"results/{wildcards.target}/summary/plip_interactions.tsv",
    }
    report_path = checkpoints.detect_cofactors.get(target=wildcards.target).output.report
    with open(report_path) as f:
        cofactor_report = json.load(f)
    if cofactor_report["has_metal"]:
        inputs["metal_energy"] = f"results/{wildcards.target}/summary/metal_site_scores.tsv"
        inputs["geometry"] = f"results/{wildcards.target}/summary/metal_geometry_validation.tsv"
    return inputs


rule compute_consensus_ranking:
    """
    Combine every independently-computed metric for this target into one
    consensus ranking (mean rank across available metrics -- see
    compute_consensus_ranking.py's docstring for the rationale, straight
    from this project's own results, for why agreement across independent
    checks matters more than trusting any single score).
    """
    input:
        unpack(consensus_inputs)
    output:
        "results/{target}/summary/consensus_ranking.tsv"
    script:
        "workflow/scripts/compute_consensus_ranking.py"
