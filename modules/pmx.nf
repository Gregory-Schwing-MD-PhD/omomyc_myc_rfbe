process downloadInputs {
    publishDir "${params.output_folder}/inputs/", mode: 'copy', overwrite: true

    output:
    path "stateA.tpr",     emit: stateA_tpr
    path "stateA_1ns.xtc", emit: stateA_traj
    path "stateB.tpr",     emit: stateB_tpr
    path "stateB_1ns.xtc", emit: stateB_traj
    path "dhdlA.zip",      emit: dhdlA
    path "dhdlB.zip",      emit: dhdlB
    path "schema.png",     emit: schema

    script:
    """
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateA.tpr     -O stateA.tpr
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateA_1ns.xtc -O stateA_1ns.xtc
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateB.tpr     -O stateB.tpr
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/stateB_1ns.xtc -O stateB_1ns.xtc
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/dhdlA.zip     -O dhdlA.zip
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/dhdlB.zip     -O dhdlB.zip
    wget -q https://github.com/bioexcel/biobb_workflows/raw/main/biobb_wf_pmx_tutorial/docker/pmx_tutorial/schema.png    -O schema.png
    """
}


workflow test {
    main:
    downloadInputs()
}