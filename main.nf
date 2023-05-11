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
params.chrom           = 1..29 
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

    // channel general processing
    chrom_list_ch = channel.from(params.chrom)



}

// template

process RNA_COMPUTE_rasqual_emperical_pvalues {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_emperical_pvalues", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path merged_results
    path permuation_merged_results

    output:
    path("rasqual_emperical_pvalues.txt")


    script:
    """
    rasqual_emperical_pvalues.R rasqual_emperical_pvalues.txt $merged_results $permuation_merged_results
    """
}



