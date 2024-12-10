
// Merge bams of different flowcells for same sample
process MERGE_FLOWCELL_BAM {
        tag "$sampleid"
        memory '20 GB'
        cpus '10'
        time '23h'

        input:
        tuple(val(sampleid), path(alignment))

        output:
        tuple( 
        val(sampleid), 
        path("${sampleid}.all.mod.bam"),
        path("${sampleid}.all.mod.bam.csi")        
        )

        script:
        def sampleid_alignment  = "${sampleid}.all.mod.bam"
        def sampleid_index      = "${sampleid}.all.mod.bam.csi"

        """
        samtools merge -@ ${task.cpus} \\
        ${alignment} \\
        -o ${sampleid_alignment};
        samtools index -@ ${task.cpus} \\
        -o ${sampleid_index} \\
        ${sampleid_alignment}
        """
        
        stub:
        """
        touch ${sampleid}.all.mod.bam
        touch ${sampleid}.all.mod.bam.csi
        """
}
