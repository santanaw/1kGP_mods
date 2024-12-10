#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BASECALL_READS    as BASECALL_READS    } from './modules/local/basecall_reads'
include { SUMMARYZE_READS   as SUMMARYZE_READS   } from './modules/local/summaryze_reads'
include { TRIM_READS        as TRIM_READS        } from './modules/local/trim_reads'



workflow {
    // Create input channel
    Channel.fromPath( params.sample_sheet  )
           .splitCsv( header:true,sep:"\t" )
           .map { 
            row -> 
            def meta = [
                sampleid:   row.sampleid,
                flowcellid: row.flowcellid,
                input_dir:  row.input_dir ]
            def pod5_dir = "${row.input_dir}/pod5"
            [ meta, pod5_dir ] 
           }
           .set{ basecall_reads_ch }

    // Basecall reads
    BASECALL_READS( basecall_reads_ch )

    // Trim basecalled reads
    TRIM_READS( BASECALL_READS.out )

    // Summary of reads
    SUMMARYZE_READS( TRIM_READS.out )

}
