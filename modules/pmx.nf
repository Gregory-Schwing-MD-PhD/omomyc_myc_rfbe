nextflow.enable.dsl=2

process downloadInputs {
    cache = true
    publishDir "${params.output_folder}/inputs/", mode: 'copy'

    input:
    val urls

    output:
    path "stateA.tpr"     , emit: stateA_tpr
    path "stateA_1ns.xtc" , emit: stateA_traj
    path "stateB.tpr"     , emit: stateB_tpr
    path "stateB_1ns.xtc" , emit: stateB_traj
    path "dhdlA.zip"      , emit: dhdlA
    path "dhdlB.zip"      , emit: dhdlB
    path "schema.png"     , emit: schema

    script:
    """
    for url in ${urls.join(" ")}; do
        fname=\$(basename \$url)
        [ -s \$fname ] || wget -q \$url -O \$fname
    done
    """
}

process extractSnapshots {
    cache = true
    publishDir "${params.output_folder}/snapshots/", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"

    input:
    path stateA_traj
    path stateA_tpr
    path stateB_traj
    path stateB_tpr

    output:
    path "stateA_frames.zip", emit: stateA_frames
    path "stateB_frames.zip", emit: stateB_frames
    path "frameA*.pdb", emit: stateA_pdbs
    path "frameB*.pdb", emit: stateB_pdbs
    path "frame*.pdb", emit: states_pdbs

    script:
    """
    python3 <<'EOF'
import zipfile
from biobb_analysis.gromacs.gmx_trjconv_str_ens import gmx_trjconv_str_ens

#### State A ####
output_framesA = 'stateA_frames.zip'

prop = {
    'selection'   : "${params.selection}",
    'start'       : int("${params.start}"),
    'end'         : int("${params.end}"),
    'dt'          : int("${params.dt}") if int("${params.dt}") > 0 else max(1, (int("${params.end}") - int("${params.start}")) // int("${params.n_snapshots}")),
    'output_name' : 'frameA',
    'output_type' : "${params.output_type}"
}

gmx_trjconv_str_ens(input_traj_path="${stateA_traj}",
                 input_top_path="${stateA_tpr}",
                 output_str_ens_path=output_framesA,
                 properties=prop)

with zipfile.ZipFile(output_framesA, 'r') as zip_f:
    zip_f.extractall()
    stateA_pdb_list = zip_f.namelist()

#### State B ####
output_framesB = 'stateB_frames.zip'

prop = {
    'selection'   : "${params.selection}",
    'start'       : int("${params.start}"),
    'end'         : int("${params.end}"),
    'dt'          : int("${params.dt}") if int("${params.dt}") > 0 else max(1, (int("${params.end}") - int("${params.start}")) // int("${params.n_snapshots}")),
    'output_name' : 'frameB',
    'output_type' : "${params.output_type}"
}

gmx_trjconv_str_ens(input_traj_path="${stateB_traj}",
                 input_top_path="${stateB_tpr}",
                 output_str_ens_path=output_framesB,
                 properties=prop)

with zipfile.ZipFile(output_framesB, 'r') as zip_f:
    zip_f.extractall()
    stateB_pdb_list = zip_f.namelist()
EOF
    """
}

process pmxMutateFrames {
    cache = true
    publishDir "${params.output_folder}/mutated/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(pdb), val(state), val(index), val(mutation)
    
    output:
    path "*_mut*.pdb", emit: mut_pdbs
    tuple path("*_mut*.pdb"), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_pmx.pmxbiobb.pmxmutate import pmxmutate
import os

# Input files from Nextflow tuple
pdb = '${pdb}'
state = '${state}'
mutation = '${mutation}'

# Create output filename
output_structure_mut = os.path.basename(pdb).replace('.pdb', '_mut' + state + '.pdb')

prop = {
    'force_field' : 'amber99sb-star-ildn-mut',
    'mutation_list': mutation,
    'binary_path' : 'pmx',
    'gmx_lib' : '${params.gmxlib}'
}

pmxmutate(input_structure_path=pdb,
          output_structure_path=output_structure_mut,
          properties=prop)
EOF
    """
}
process pdb2gmxFrames {
    cache = true
    debug = false
    publishDir "${params.output_folder}/topologies/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(input_pdb), val(state), val(index), val(mutation)
    
    output:
    tuple path("*_MUT_top.zip"), path("*_MUT.gro"), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx
import os

# Input
pdb_file = '${input_pdb}'

# Outputs with _MUT appended
basename = os.path.basename(pdb_file).replace('.pdb', '')
output_gro = f"{basename}_MUT.gro"
output_top = f"{basename}_MUT_top.zip"

# Properties
props = {
    'force_field': 'amber99sb-star-ildn-mut',
    'gmx_lib': '${params.gmxlib}'
}

# Run pdb2gmx
pdb2gmx(input_pdb_path=pdb_file,
        output_gro_path=output_gro,
        output_top_zip_path=output_top,
        properties=props)
EOF
    """
}

process pmxGentop {
    cache = true
    debug = false
    publishDir "${params.output_folder}/topologies/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(input_top_zip), path(input_gro), val(state), val(index), val(mutation)
    
    output:
    tuple path(input_top_zip), path("*_MUT_top.zip"), path(input_gro), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_pmx.pmxbiobb.pmxgentop import pmxgentop
import os

# Input
top_zip = '${input_top_zip}'

# Outputs with _MUT appended
basename = os.path.basename(top_zip).replace('.zip', '')
output_top_zip = f"{basename}_MUT_top.zip"
output_log = f"{basename}_MUT_top.log"
print(output_log)

# Properties
prop = {
    'force_field': 'amber99sb-star-ildn-mut',
    'binary_path': 'pmx',
    'gmx_lib': '${params.gmxlib}'
}

# Run pmxgentop
pmxgentop(input_top_zip_path=top_zip,
          output_top_zip_path=output_top_zip,
          output_log_path=output_log,
          properties=prop)
EOF
    """
}

process gmxMakeNdx {
    cache = true
    debug = false
    publishDir "${params.output_folder}/topologies/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), val(state), val(index), val(mutation)
    
    output:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path("*_ndx.ndx"), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_gromacs.gromacs.make_ndx import make_ndx
import os

# Input
gro_file = '${input_gro}'

# Output .ndx file with _ndx appended
basename = os.path.basename(gro_file).replace('.gro', '')
output_ndx = f"{basename}_ndx.ndx"

# Properties
prop = {
    'selection': 'a D*\\n0 & ! 19\\nname 20 FREEZE'
}

# Run make_ndx
make_ndx(input_structure_path=gro_file,
         output_ndx_path=output_ndx,
         properties=prop)
EOF
    """
}
process gmxMinimize {
    cache = true
    debug = false
    publishDir "${params.output_folder}/minimization/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path(input_ndx), val(state), val(index), val(mutation)
    
    output:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path(input_ndx), path("emout.gro"), path("emout.trr"), path("emout.edr"), path("emout.log"), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.mdrun import mdrun
import os

# Input files
gro_file = '${input_gro}'
top_zip = '${processed_top_zip}'
ndx_file = '${input_ndx}'

# Check if FREEZE group exists in the index file
freeze_group_exists = False
try:
    with open(ndx_file, 'r') as f:
        content = f.read()
        if '[ FREEZE ]' in content:
            freeze_group_exists = True
            print("FREEZE group found in index file")
        else:
            print("FREEZE group NOT found in index file")
except:
    print("Could not read index file")

# Grompp: Creating portable binary run file for energy minimization
output_tpr_min = 'em.tpr'

# Base MDP parameters
mdp_params = {
    'integrator': 'steep',
    'emtol': '100',
    'dt': '0.001',
    'nsteps': '10000',
    'nstcomm': '1',
    'nstcalcenergy': '1'
}

# Add freeze parameters only if FREEZE group exists
if freeze_group_exists:
    mdp_params['freezegrps'] = 'FREEZE'
    mdp_params['freezedim'] = 'Y Y Y'
    print("Adding freeze constraints to MDP")
else:
    print("Running without freeze constraints")

prop_grompp = {
    'gmx_lib': '${params.gmxlib}',
    'mdp': mdp_params,
    'simulation_type': 'minimization'
}

# Run grompp
grompp(input_gro_path=gro_file,
       input_top_zip_path=top_zip,
       input_ndx_path=ndx_file,
       output_tpr_path=output_tpr_min,
       properties=prop_grompp)

# Mdrun: Running minimization
output_min_trr = 'emout.trr'
output_min_gro = 'emout.gro'
output_min_edr = 'emout.edr'
output_min_log = 'emout.log'

# Run mdrun
mdrun(input_tpr_path=output_tpr_min,
      output_trr_path=output_min_trr,
      output_gro_path=output_min_gro,
      output_edr_path=output_min_edr,
      output_log_path=output_min_log)
EOF
    """
}

process gmxEnergyAnalysis {
    cache = true
    debug = false
    publishDir "${params.output_folder}/analysis/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path(input_ndx), path(min_gro), path(min_trr), path(min_edr), path(min_log), val(state), val(index), val(mutation)
    
    output:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path(input_ndx), path(min_gro), path(min_trr), path(min_edr), path(min_log), path("*.xvg"), path("*.png"), path("*.dat"), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_analysis.gromacs.gmx_energy import gmx_energy
import matplotlib.pyplot as plt
import numpy as np
import os

# Input files
edr_file = '${min_edr}'

# Create prop dict and inputs/outputs
output_min_ene_xvg = 'min_ene.xvg'
prop = {
    'terms': ["Potential"]
}

# Create and launch bb
gmx_energy(input_energy_path=edr_file, 
          output_xvg_path=output_min_ene_xvg, 
          properties=prop)

# Read data from XVG file and filter energy values higher than 1000 kJ/mol
x_data = []
y_data = []

with open(output_min_ene_xvg, 'r') as energy_file:
    for line in energy_file:
        if not line.startswith(("#", "@")):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    x_val = float(parts[0])
                    y_val = float(parts[1])
                    if y_val < 1000:  # Filter high energy values
                        x_data.append(x_val)
                        y_data.append(y_val)
                except ValueError:
                    continue

# Create matplotlib plot and save as PNG
plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, 'b-', linewidth=1.5, markersize=2)
plt.xlabel('Energy Minimization Step')
plt.ylabel('Potential Energy (kJ/mol)')
plt.title('Energy Minimization')
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save PNG
basename = os.path.basename(edr_file).replace('.edr', '')
png_filename = f"{basename}_energy_plot.png"
plt.savefig(png_filename, dpi=300, bbox_inches='tight')
plt.close()

# Create data file for gnuplot/xmgrace (simple two-column format)
dat_filename = f"{basename}_energy_data.dat"
with open(dat_filename, 'w') as f:
    f.write("# Energy Minimization Data\\n")
    f.write("# Column 1: Step\\n")
    f.write("# Column 2: Potential Energy (kJ/mol)\\n")
    for x, y in zip(x_data, y_data):
        f.write(f"{x:12.6f} {y:12.6f}\\n")

print(f"Generated {png_filename} and {dat_filename}")
print(f"XVG file: {output_min_ene_xvg}")
print(f"Data points: {len(x_data)}")
if len(y_data) > 0:
    print(f"Energy range: {min(y_data):.2f} to {max(y_data):.2f} kJ/mol")
EOF
    """
}
process gmxEquilibrate {
    cache = true
    debug = false
    publishDir "${params.output_folder}/equilibration/state_${state}/frame_${index}", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"
    
    input:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path(input_ndx), path(min_gro), path(min_trr), path(min_edr), path(min_log), path(energy_xvg), path(energy_png), path(energy_dat), val(state), val(index), val(mutation)
    
    output:
    tuple path(input_top_zip), path(processed_top_zip), path(input_gro), path(input_ndx), path(min_gro), path(min_trr), path(min_edr), path(min_log), path(energy_xvg), path(energy_png), path(energy_dat), path("eqout*.gro"), path("eqout*.trr"), path("eqout*.edr"), path("eqout*.log"), val(state), val(index), val(mutation), emit: system
    
    script:
    """
    python3 <<'EOF'
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.mdrun import mdrun
import os

# Input files
min_gro = '${min_gro}'
processed_top_zip = '${processed_top_zip}'
state = '${state}'
gmxlib = '${params.gmxlib}'

print(f"Starting equilibration for state {state}")
print(f"Using structure: {min_gro}")
print(f"Using topology: {processed_top_zip}")

# Grompp: Creating portable binary run file for system equilibration
output_tpr_eq = f'eq{state}_20ps.tpr'

# MDP parameters for equilibration
mdp_params = {
    'nsteps': '10000',      # 10000 steps x 0.001 ps = 10 ps
    'dt': '0.001',          # 1 fs timestep for proper dummy atom equilibration
    'nstcomm': '1',
    'nstcalcenergy': '1'
}

prop_grompp = {
    'gmx_lib': gmxlib,
    'mdp': mdp_params,
    'simulation_type': 'free'
}

# Run grompp (without index file - not needed for basic equilibration)
print(f"Running grompp for state {state}...")
grompp(input_gro_path=min_gro,
       input_top_zip_path=processed_top_zip,
       output_tpr_path=output_tpr_eq,
       properties=prop_grompp)

print(f"Generated TPR file: {output_tpr_eq}")

# Mdrun: Running equilibration
output_eq_trr = f'eqout{state}.trr'
output_eq_gro = f'eqout{state}.gro'
output_eq_edr = f'eqout{state}.edr'
output_eq_log = f'eqout{state}.log'

print(f"Running mdrun for state {state}...")
mdrun(input_tpr_path=output_tpr_eq,
      output_trr_path=output_eq_trr,
      output_gro_path=output_eq_gro,
      output_edr_path=output_eq_edr,
      output_log_path=output_eq_log)

print(f"Equilibration completed for state {state}")
print(f"Generated files:")
print(f"  Structure: {output_eq_gro}")
print(f"  Trajectory: {output_eq_trr}")
print(f"  Energy: {output_eq_edr}")
print(f"  Log: {output_eq_log}")
EOF
    """
}

workflow test {
    main:
    // Step 1: define all URLs in one channel
    urls_ch = Channel.value([
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateA.tpr",
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateA_1ns.xtc",
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateB.tpr",
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateB_1ns.xtc",
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/dhdlA.zip",
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/dhdlB.zip",
        "https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/schema.png"
    ])

    // Step 2: run downloadInputs
    inputs = downloadInputs(urls_ch)

    // Step 3: run extractSnapshots
    snapshots = extractSnapshots(
        inputs.stateA_traj,
        inputs.stateA_tpr,
        inputs.stateB_traj,
        inputs.stateB_tpr
    )


    // Define mutation dictionary
    def mutations = [
        'A': '10Ala',
        'B': '10Ile'
    ]

    stateA_sorted = snapshots.stateA_pdbs
        .collect()
        .map { files ->
            files.sort { f -> (f.name =~ /\d+/)[0].toInteger() }
        }

    stateB_sorted = snapshots.stateB_pdbs
        .collect()
        .map { files ->
            files.sort { f -> (f.name =~ /\d+/)[0].toInteger() }
        }

// Create indexed channels with mutation info and take first n entries
chA_indexed = stateA_sorted
    .flatMap { files ->
        files.withIndex().collect { file, idx -> [file, 'A', idx, mutations['A']] }
    }
    .take(params.n_snapshots_to_run)

chB_indexed = stateB_sorted
    .flatMap { files ->
        files.withIndex().collect { file, idx -> [file, 'B', idx, mutations['B']] }
    }
    .take(params.n_snapshots_to_run)

// Combine both channels for the pipeline
running_snaps = chA_indexed.mix(chB_indexed)

// Step 6: run pmxMutateFrames on all pairs
running_snaps.view()
pmxMutateFrames(running_snaps)
pdb2gmxFrames(pmxMutateFrames.output.system)
pmxGentop(pdb2gmxFrames.output.system)
gmxMakeNdx(pmxGentop.output.system)
gmxMinimize(gmxMakeNdx.output.system)
gmxEnergyAnalysis(gmxMinimize.output.system)
gmxEquilibrate(gmxEnergyAnalysis.output.system)
}

