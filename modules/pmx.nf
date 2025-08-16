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
    path "stateA_frame*.pdb", emit: stateA_pdbs
    path "stateB_frame*.pdb", emit: stateB_pdbs

    script:
    """
    python3 <<'EOF'
import zipfile
from biobb_analysis.gromacs.gmx_trjconv_str_ens import gmx_trjconv_str_ens

#### State A ####
output_framesA = 'stateA_frames.zip'
propA = {
    'selection' : 'System',
    'start': 1,
    'end': 1000,
    'dt': 200,
    'output_name': 'stateA_frame',
    'output_type': 'pdb'
}
gmx_trjconv_str_ens(input_traj_path="${stateA_traj}",
                    input_top_path="${stateA_tpr}",
                    output_str_ens_path=output_framesA,
                    properties=propA)

with zipfile.ZipFile(output_framesA, 'r') as zip_f:
    zip_f.extractall()

#### State B ####
output_framesB = 'stateB_frames.zip'
propB = {
    'selection' : 'System',
    'start': 1,
    'end': 1000,
    'dt': 200,
    'output_name': 'stateB_frame',
    'output_type': 'pdb'
}
gmx_trjconv_str_ens(input_traj_path="${stateB_traj}",
                    input_top_path="${stateB_tpr}",
                    output_str_ens_path=output_framesB,
                    properties=propB)

with zipfile.ZipFile(output_framesB, 'r') as zip_f:
    zip_f.extractall()
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
    snapshots.stateA_pdbs.view { "Extracted StateA snapshot: $it" }
    snapshots.stateB_pdbs.view { "Extracted StateB snapshot: $it" }
}
