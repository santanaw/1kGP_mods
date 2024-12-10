#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TAG_ALIGNMENT             as TAG_ALIGNMENT              } from './modules/local/tag_alignment'
include { MERGE_FLOWCELL_BAM        as MERGE_FLOWCELL_BAM         } from './modules/local/merge_flowcell_bam'
include { SUBSET_ALIGNMENT          as SUBSET_ALIGNMENT           } from './modules/local/subset_alignment'
include { FILTER_HAPLOTYPE_READS    as FILTER_HAPLOTYPE1_READS    } from './modules/local/filter_haplotype_reads'
include { FILTER_HAPLOTYPE_READS    as FILTER_HAPLOTYPE2_READS    } from './modules/local/filter_haplotype_reads'
include { PILEUP_MODS_HAPLOTYPE     as PILEUP_MODS_HAPLOTYPE1     } from './modules/local/pileup_mods_haplotype'
include { PILEUP_MODS_HAPLOTYPE     as PILEUP_MODS_HAPLOTYPE2     } from './modules/local/pileup_mods_haplotype'
include { ANNOTATE_PILEUP_HAPLOTYPE as ANNOTATE_PILEUP_HAPLOTYPE1 } from './modules/local/annotate_pileup_haplotype'
include { ANNOTATE_PILEUP_HAPLOTYPE as ANNOTATE_PILEUP_HAPLOTYPE2 } from './modules/local/annotate_pileup_haplotype'


workflow {

    // Create main input channel
    Channel.fromPath( params.sample_sheet  )
           .splitCsv( header:true,sep:"\t" )
           .map { 
            row -> 
            def meta = [
                sampleid:    row.sampleid,
                flowcellid:  row.flowcellid
	    ]
            [ meta, row.input_file ] 
           }
           .set{ tag_alignment_in }

    // Create read haplotype channel
    Channel.fromPath( params.haplotype_read_sheet  )
           .splitCsv( header:true,sep:"\t" )
           .map { 
            row -> 
            def meta = [
                sampleid:          row.sampleid,
                haplotype1_reads:  row.haplotype1_reads,
                haplotype2_reads:  row.haplotype2_reads
	    ]
            [ row.sampleid, row.haplotype1_reads, row.haplotype2_reads ] 
           }
           .set{ haplotype_reads }
    // Count number of flowcells per sample
    tag_alignment_in.map{
                        meta, input_file ->
                        [meta.sampleid, meta.flowcellid]
                     }
                     .groupTuple()
                     .map{ 
                        sampleid, flowcells -> 
                        [(sampleid): flowcells.size()]
                     }
                     .set{ flowcell_counts }

    // Tag alignment
    TAG_ALIGNMENT( tag_alignment_in )

    // VIEW
    TAG_ALIGNMENT.out.view()
    flowcell_counts.view()

    // Group input alignments to merge
    TAG_ALIGNMENT.out
                 .map{
                        meta, input_file -> 
                        [ meta.sampleid, input_file ]
                 }
                 .groupTuple()
                 .set{ merge_flowcell_bam_in }

    // VIEW
    merge_flowcell_bam_in.view()

    // Merge all alignments for each sample
    MERGE_FLOWCELL_BAM( merge_flowcell_bam_in )

    // Add bed region file to channels
    MERGE_FLOWCELL_BAM.out
                      .map{
                        sampleid, alignment, alignment_index ->
                        [sampleid, alignment, alignment_index, params.gene_regions]
                      }
                      .set{ subset_alignment_in }

    // Subset alignment by gene regions
    SUBSET_ALIGNMENT( subset_alignment_in )

    // Add read ids for haplotypes 1 and 2
    SUBSET_ALIGNMENT.out
                    .join(haplotype_reads)
                    .map{
                        sampleid, alignment, alignment_index, 
                        haplotype1_reads, haplotype2_reads ->
                        [sampleid, "hap1", alignment, alignment_index,haplotype1_reads]
                    }
                    .set{ filter_haploype1_reads_in }

    SUBSET_ALIGNMENT.out
                    .join(haplotype_reads)
                    .map{
                        sampleid, alignment, alignment_index, 
                        haplotype1_reads, haplotype2_reads ->
                        [sampleid, "hap2", alignment, alignment_index,haplotype2_reads]
                    }
                    .set{ filter_haploype2_reads_in }

    // Filter by haplotype 1
    FILTER_HAPLOTYPE1_READS( filter_haploype1_reads_in )

    // Filter by haplotype 2
    FILTER_HAPLOTYPE2_READS( filter_haploype2_reads_in )

    // Add Genome to pileup
    FILTER_HAPLOTYPE1_READS.out
                           .map{ 
                                sampleid, hap_n, alignment, alignment_index -> 
                                [sampleid, hap_n, alignment, alignment_index, params.genome ] 
                           }
                           .set{ pileup_mods_haplotype1_in }
    FILTER_HAPLOTYPE2_READS.out
                           .map{ 
                                sampleid, hap_n, alignment, alignment_index -> 
                                [sampleid, hap_n, alignment, alignment_index, params.genome ] 
                           }
                           .set{ pileup_mods_haplotype2_in }

    // Pileup mods for haplotype 1
    PILEUP_MODS_HAPLOTYPE1( pileup_mods_haplotype1_in )
    // Pileup mods for haplotype 2
    PILEUP_MODS_HAPLOTYPE2( pileup_mods_haplotype2_in )

    // Add bed region file to channels
    PILEUP_MODS_HAPLOTYPE1.out
                           .map{ 
                                sampleid, hap_n, pileup_log, pileup  -> 
                                [sampleid, hap_n, pileup, params.gene_regions ] 
                           }
                           .set{ annotate_pileup_haplotype1_in }
    PILEUP_MODS_HAPLOTYPE2.out
                           .map{ 
                                sampleid, hap_n, pileup_log, pileup  -> 
                                [sampleid, hap_n, pileup, params.gene_regions ] 
                           }
                           .set{ annotate_pileup_haplotype2_in }

    // Annotate pileups of haplotype 1 with genes
    ANNOTATE_PILEUP_HAPLOTYPE1( annotate_pileup_haplotype1_in )
    // Annotate pileups of haplotype 2 with genes
    ANNOTATE_PILEUP_HAPLOTYPE2( annotate_pileup_haplotype2_in )


}
