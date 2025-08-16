#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { test } from './modules/pmx'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run [-profile local/docker/singularity/slurm] .

# Docker
nextflow run . -profile docker

    """.stripIndent()
}


// Main workflow
workflow {

    results = test()

}