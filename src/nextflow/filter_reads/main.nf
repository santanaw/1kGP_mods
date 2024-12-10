#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FILTER_READS  as FILTER_READS    } from './modules/local/filter_reads'

workflow {
    // Create input channel
    Channel.fromPath( params.sample_sheet  )
           .splitCsv( header:true,sep:"\t" )
           .map { 
            row -> 
            def meta = [
                sampleid:    row.sampleid,
                flowcellid:  row.flowcellid]
            [ meta, row.input_file ] 
           }
           .set{ filter_reads_in }

    // Basecall reads
    FILTER_READS( filter_reads_in )

}
