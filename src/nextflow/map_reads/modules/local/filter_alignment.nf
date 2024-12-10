// Filter read alignment
process FILTER_ALIGNMENT {
        tag "$meta.flowcellid"
        memory '24 GB'
        cpus '12'
        time '23h'

        input:
        tuple(val(meta), path(alignment), path(alignment_idx))

        output:
        tuple(
        val(meta), 
        path("${meta.sampleid}.${meta.flowcellid}.final.sorted.mod.bam"),
        path("${meta.sampleid}.${meta.flowcellid}.final.sorted.mod.bam.bai")
        )

        script:
        def filtered_aln  = "${meta.sampleid}.${meta.flowcellid}.final.mod.bam"
        def sorted_aln    = "${meta.sampleid}.${meta.flowcellid}.final.sorted.mod.bam"
        def aln_index     = "${meta.sampleid}.${meta.flowcellid}.final.sorted.mod.bam.bai"       

        """
        samtools view \\
        -@ ${task.cpus} -F 0x004 -bh \\
        ${alignment} > \\
        ${filtered_aln};\\
        samtools sort -@ ${task.cpus} \\
        -o ${sorted_aln} \\
        ${filtered_aln};\\
        samtools index -@ ${task.cpus} \\
        -o ${aln_index} \\
        ${sorted_aln}
        """
        
        stub:

        """
        touch ${meta.sampleid}.${meta.flowcellid}.final.sorted.mod.bam
        touch ${meta.sampleid}.${meta.flowcellid}.final.sorted.mod.bam.bai
        """
}