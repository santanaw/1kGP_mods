// Merge bams of different flowcells for same sample
process FILTER_HAPLOTYPE_READS {
        tag "$sampleid"
        memory '20 GB'
        cpus '10'
        time '23h'

        input:
        tuple(
        val(sampleid), 
        val(hap_n),
        path(alignment), 
        path(alignment_index), 
        path(haplotype_reads)
        )

        output:
        tuple( 
        val(sampleid), 
        val(hap_n),
        path("${sampleid}.${hap_n}.genes.sorted.mod.bam"),
        path("${sampleid}.${hap_n}.genes.sorted.mod.bam.csi")        
        )

        script:
        def hap_alignment                = "${sampleid}.${hap_n}.genes.mod.bam"
        def sorted_hap_alignment         = "${sampleid}.${hap_n}.genes.sorted.mod.bam"
        def sorted_hap_alignment_index   = "${sampleid}.${hap_n}.genes.sorted.mod.bam.csi"

        """
        samtools view -bh -@ ${task.cpus} \\
        -N ${haplotype_reads} \\
        ${alignment} > \\
        ${hap_alignment}; \\
        samtools sort -@ ${task.cpus} \\
        -o ${sorted_hap_alignment} \\
        ${hap_alignment}; \\
        samtools index -@ ${task.cpus} \\
        -o ${sorted_hap_alignment_index} \\
        ${sorted_hap_alignment}
        """
        
        stub:
        """
        touch ${sampleid}.${hap_n}.genes.mod.bam
        touch ${sampleid}.${hap_n}.genes.sorted.mod.bam
        touch ${sampleid}.${hap_n}.genes.sorted.mod.bam.csi
        """
}
