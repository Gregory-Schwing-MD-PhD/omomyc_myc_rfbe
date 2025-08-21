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
    publishDir "${params.output_folder}/mutated/", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"

    input:
    tuple path(pdbA), path(pdbB)

    output:
    path "*_mutA.pdb", emit: mutA_pdbs
    path "*_mutB.pdb", emit: mutB_pdbs
    tuple path("*_mutA.pdb"), path("*_mutB.pdb"), emit: mut_pdbs

    script:
    """
    python3 <<'EOF'
from biobb_pmx.pmxbiobb.pmxmutate import pmxmutate
import os

# Input files from Nextflow tuple
pdbA = '${pdbA}'
pdbB = '${pdbB}'

# State A: WT -> Mut
output_structure_mutA = os.path.basename(pdbA).replace('.pdb', '_mutA.pdb')
propA = {
    'force_field'  : 'amber99sb-star-ildn-mut',
    'mutation_list': '10Ala',
    'binary_path'  : 'pmx',
    'gmx_lib'      : '${params.gmxlib}'
}
pmxmutate(input_structure_path=pdbA,
          output_structure_path=output_structure_mutA,
          properties=propA)

# State B: Mut -> WT
output_structure_mutB = os.path.basename(pdbB).replace('.pdb', '_mutB.pdb')
propB = {
    'force_field'  : 'amber99sb-star-ildn-mut',
    'mutation_list': '10Ile',
    'binary_path'  : 'pmx',
    'gmx_lib'      : '${params.gmxlib}'
}
pmxmutate(input_structure_path=pdbB,
          output_structure_path=output_structure_mutB,
          properties=propB)
EOF
    """
}
process pdb2gmxFrames {
    cache = true
    debug = false
    publishDir "${params.output_folder}/topologies/", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"

    input:
    path input_pdb

    output:
    path "*_MUT.gro", emit: gro_files
    path "*_MUT_top.zip", emit: top_files

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
    publishDir "${params.output_folder}/topologies/", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"

    input:
    path input_top_zip

    output:
    path "*_MUT_top.zip", emit: top_files
    //path "*_MUT_top.log", emit: log_files

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
    publishDir "${params.output_folder}/topologies/", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"

    input:
    path input_gro  // input structure, e.g., output from pdb2gmx

    output:
    path "*_ndx.ndx", emit: ndx_files

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

process gmxGromppMin {
    cache = true
    debug = true
    publishDir "${params.output_folder}/minimization/", mode: 'copy', overwrite: true
    container "${params.container__biobb_pmx}"

    input:
    path input_gro
    path input_top_zip
    path input_ndx

    output:
    path "*_em.tpr", emit: tpr_files

    script:
    """
    python3 <<'EOF'
from biobb_gromacs.gromacs.grompp import grompp
import os

# Inputs
gro_file = '${input_gro}'
top_zip = '${input_top_zip}'
ndx_file = '${input_ndx}'

# Output filename
basename = os.path.basename(gro_file).replace('.gro', '')
output_tpr = f"{basename}_em.tpr"

# Properties
prop = {
    'gmx_lib': '${params.gmxlib}',
    'mdp': {
        'integrator': 'steep',
        'emtol': '100',
        'dt': '0.001',
        'nsteps': '10000',
        'nstcomm': '1',
        'nstcalcenergy': '1',
        'freezegrps': 'FREEZE',
        'freezedim': 'Y Y Y'
    },
    'simulation_type': 'minimization'
}

# Run grompp
grompp(input_gro_path=gro_file,
       input_top_zip_path=top_zip,
       input_ndx_path=ndx_file,
       output_tpr_path=output_tpr,
       properties=prop)
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

    // Debug / verify outputs
    //snapshots.stateA_pdbs.view { "Extracted StateA snapshot: $it" }
    //snapshots.stateB_pdbs.view { "Extracted StateB snapshot: $it" }

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

    // Debug: let's see what combine actually produces
    stateA_sorted.combine(stateB_sorted).view { "Combined result: $it" }

    // Try a different approach - zip the indexed channels
    chA_indexed = stateA_sorted
        .flatMap { files ->
            files.withIndex().collect { file, idx -> [file, 'A', idx] }
        }

    chB_indexed = stateB_sorted
        .flatMap { files ->
            files.withIndex().collect { file, idx -> [file, 'B', idx] }
        }

    // Group by index and collect pairs
    ch_pbc = chA_indexed
        .mix(chB_indexed)
        .groupTuple(by: 2)  // group by index (3rd element)
        .map { grouped ->
            // grouped will be [file_list, label_list, index]
            def files = grouped[0]
            def labels = grouped[1] 
            def index = grouped[2]
            
            // Create pairs for this index
            (0..<files.size()).collect { i ->
                [files[i], labels[i], index]
            }
        }

    // How many snapshots to run
    running_snaps = ch_pbc.take(params.n_snapshots_to_run).flatMap { it }
    running_snaps.view()
    // Step 6: run pmxMutateFrames on all pairs
    pmxMutateFrames(running_snaps)
    pdb2gmxFrames(pmxMutateFrames.output.mut_pdbs.flatten())
    pmxGentop(pdb2gmxFrames.output.top_files)
    gmxMakeNdx(pdb2gmxFrames.output.gro_files)
}


workflow prod {
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

    // Debug / verify outputs
    //snapshots.stateA_pdbs.view { "Extracted StateA snapshot: $it" }
    //snapshots.stateB_pdbs.view { "Extracted StateB snapshot: $it" }

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

// Debug: let's see what combine actually produces
stateA_sorted.combine(stateB_sorted).view { "Combined result: $it" }

// Try a different approach - zip the indexed channels
chA_indexed = stateA_sorted
    .flatMap { files ->
        files.withIndex().collect { file, idx -> [file, 'A', idx] }
    }

chB_indexed = stateB_sorted
    .flatMap { files ->
        files.withIndex().collect { file, idx -> [file, 'B', idx] }
    }

// Group by index and collect pairs
ch_pbc = chA_indexed
    .mix(chB_indexed)
    .groupTuple(by: 2)  // group by index (3rd element)
    .map { grouped ->
        // grouped will be [file_list, label_list, index]
        def files = grouped[0]
        def labels = grouped[1] 
        def index = grouped[2]
        
        // Create pairs for this index
        (0..<files.size()).collect { i ->
            [files[i], labels[i], index]
        }
    }

ch_pbc.view()

return


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

// Add index and label A to stateA_sorted - correct tuple order: [file, label, index]
chA_indexed = stateA_sorted
    .flatMap { files ->
        files.withIndex().collect { file, idx -> [file, 'A', idx] }
    }

// Add index and label B to stateB_sorted - correct tuple order: [file, label, index]  
chB_indexed = stateB_sorted
    .flatMap { files ->
        files.withIndex().collect { file, idx -> [file, 'B', idx] }
    }

// Merge the channels properly
ch_pbc = chA_indexed.mix(chB_indexed)
ch_pbc.view()

return

    // Step 6: run pmxMutateFrames on all pairs
    pmxMutateFrames(ch_pbc.take(1))
    pdb2gmxFrames(pmxMutateFrames.output.mut_pdbs.flatten())
    pmxGentop(pdb2gmxFrames.output.top_files)
    gmxMakeNdx(pdb2gmxFrames.output.gro_files)
}
