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

    //stateA_sorted.view()

    stateB_sorted = snapshots.stateB_pdbs
        .collect()
        .map { files ->
            files.sort { f -> (f.name =~ /\d+/)[0].toInteger() }
        }
    //stateB_sorted.view()
    chA_indexed = stateA_sorted.flatten()
    chB_indexed = stateB_sorted.flatten()
    chA_indexed.merge(chB_indexed).set { ch_pbc }

    ch_pbc.take(1).view()

    // Step 6: run pmxMutateFrames on all pairs
    pmxMutateFrames(ch_pbc.take(1))
    pdb2gmxFrames(pmxMutateFrames.output.mut_pdbs.flatten())
    pmxGentop(pdb2gmxFrames.output.top_files)

}
