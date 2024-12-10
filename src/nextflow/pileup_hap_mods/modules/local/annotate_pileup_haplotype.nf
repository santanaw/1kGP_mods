// Annotate pileups per haplotype
process ANNOTATE_PILEUP_HAPLOTYPE {
        tag "$sampleid"
        memory '20 GB'
        cpus '1'
        time '23h'

        input:
        tuple(
        val(sampleid), 
        val(hap_n),
        path(pileup), 
        path(gene_regions) 
        )

        output:
        tuple( 
        val(sampleid), 
        val(hap_n),
	path("${sampleid}.${hap_n}.genes.mod.bed")
        )

        script:
        def annotated_bed      = "${sampleid}.${hap_n}.genes.mod.bed"

        """
        bedtools intersect \\
        -a ${pileup} \\
        -b ${gene_regions} -wb > \\
        ${annotated_bed}
        """
        
        stub:
        """
        touch ${sampleid}.${hap_n}.genes.mod.bed
        """
}
