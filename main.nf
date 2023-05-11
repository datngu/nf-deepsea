#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-rasqual
========================================================================================
                DeepSea Pipeline with nextflow.
                https://github.com/datngu/nf-deepsea
                Author: Dat T Nguyen
                Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/





/*
 Define the default parameters
*/ 
params.genome          = "$baseDir/data/ref/genome.fa"
params.peaks           = "$baseDir/data/peak/*"

params.outdir          = "results"
params.trace_dir       = "trace_dir"

// running options
params.chrom           = 29 
params.window          = 200
params.seqlen          = 1000 



// pipeline options
params.generate_data  = true
params.train          = false
params.evaluate       = false
params.predict        = false


log.info """\
================================================================
                        nf-deepsea
================================================================
    genome              : $params.genome
    peaks               : $params.peaks
    outdir              : $params.outdir
    trace_dir           : $params.trace_dir


    chrom               : $params.chrom
    window              : $params.window
    seqlen              : $params.seqlen


    generate_data       : $params.generate_data
    train               : $params.train
    evaluate            : $params.evaluate
    predict             : $params.predict


================================================================
"""

nextflow.enable.dsl=2


//include { ATAC_deltaSVM_slipt_bed; ATAC_deltaSVM_gen_null_seqs; ATAC_deltaSVM_train; ATAC_deltaSVM_merge_models; ATAC_deltaSVM_gen_10mers; ATAC_deltaSVM_score_10mers; ATAC_deltaSVM_average_weights; ATAC_deltaSVM_input_generator; ATAC_deltaSVM } from './module/deltaSVM'



workflow {

    INDEX_genome(params.genome)
    BIN_genome(INDEX_genome.out)

    all_peaks = Channel.fromPath($params.peaks)
    all_peaks.view()

}






// preprocessing data



process INDEX_genome {
    container 'ndatth/deepsea:v0.0.0'
    publishDir "${params.outdir}/genome", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path "genome.fa"

    output:
    path("genome.fa*")


    script:
    """
    samtools faidx genome.fa
    """
}


process BIN_genome {
    container 'ndatth/deepsea:v0.0.0'
    publishDir "${params.outdir}/bed_files", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path genome

    output:
    path("genome_window.bed")


    script:
    """
    generate_coordinate_onebed.py --genome genome.fa.fai --out genome_window.bed --window $params.window --chrom $params.chrom
    """
}


process BIN_genome {
    container 'ndatth/deepsea:v0.0.0'
    publishDir "${params.outdir}/bed_files", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path genome

    output:
    path("genome_window.bed")


    script:
    """
    generate_coordinate_onebed.py --genome genome.fa.fai --out genome_window.bed --window $params.window --chrom $params.chrom
    """
}


