#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MAP_READS          as MAP_READS          }  from './modules/local/map_reads'
include { FILTER_ALIGNMENT   as FILTER_ALIGNMENT   }  from './modules/local/filter_alignment'


workflow {
    // Create input channel
    Channel.fromPath( params.sample_sheet  )
           .splitCsv( header:true,sep:"\t" )
           .map { 
            row -> 
            def meta = [
                sampleid:    row.sampleid,
                flowcellid:  row.flowcellid]
            [ meta, row.input_file, params.genome ] 
           }
           .set{ map_reads_in }

    // Map reads
    MAP_READS( map_reads_in )

    // Filter alignment reads
    FILTER_ALIGNMENT( MAP_READS.out )

}
