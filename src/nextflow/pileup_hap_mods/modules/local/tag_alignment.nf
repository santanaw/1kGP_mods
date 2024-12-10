// Tag alignment merge bams
process TAG_ALIGNMENT {
        tag "$meta.flowcellid"
        memory '20 GB'
        cpus '10'
        time '23h'

        input:
        tuple(val(meta), path(alignment))

        output:
        tuple( 
        val(meta), 
        path("${meta.sampleid}.${meta.flowcellid}.tag.sorted.mod.bam")
        )

        script:
        def tag_alignment  = "${meta.sampleid}.${meta.flowcellid}.tag.sorted.mod.bam"

        """
        samtools addreplacerg -@ ${task.cpus} --write-index \\
        -r "@RG\tID:${meta.flowcellid}\tSM:${meta.sampleid}\tPL:ONT" \\
        -o ${tag_alignment} \\
        ${alignment}
        """
        
        stub:

        """
        touch ${meta.sampleid}.${meta.flowcellid}.tag.sorted.mod.bam
        """
}