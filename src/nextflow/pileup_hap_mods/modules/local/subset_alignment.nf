// Subset alignments for each bed region
process SUBSET_ALIGNMENT {
        tag "$sampleid"
        memory '10 GB'
        cpus '5'
        time '23h'

        input:
        tuple(
        val(sampleid), 
        path(alignment), 
        path(alignment_index), 
        path(bed_regions)
        )

        output:
        tuple( 
        val(sampleid), 
        path("${sampleid}.all.genes.sorted.mod.bam"),
        path("${sampleid}.all.genes.sorted.mod.bam.csi")        
        )

        script:
        def gene_alignment         = "${sampleid}.genes.mod.bam"
        def alignment_index        = "${sampleid}.genes.mod.bam.csi"
        def sorted_gene_alignment  = "${sampleid}.all.genes.sorted.mod.bam"
        def alignment_sorted_index = "${sampleid}.all.genes.sorted.mod.bam.csi"

        """
        samtools view -@ ${task.cpus} -bh \\
        --target-file ${bed_regions} \\
        ${alignment} > \\
        ${gene_alignment}; \\
        samtools sort -@ ${task.cpus} \\
        -o ${sorted_gene_alignment} \\
        ${gene_alignment}; \\
        samtools index -@ ${task.cpus} \\
        -o ${alignment_sorted_index} \\
        ${sorted_gene_alignment}
        """
        
        stub:
        """
        touch ${sampleid}.genes.mod.bam
        touch ${sampleid}.genes.mod.bam.csi
        touch ${sampleid}.all.genes.sorted.mod.bam
        touch ${sampleid}.all.genes.sorted.mod.bam.csi
        """
}
