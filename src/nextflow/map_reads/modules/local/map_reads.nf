// Map reads to genome
process MAP_READS {
        tag "$meta.flowcellid"
        memory '48 GB'
        cpus '24'
        time '23h'

        input:
        tuple(val(meta), path(reads), path(genome))

        output:
        tuple( 
        val(meta), 
        path("${meta.sampleid}.${meta.flowcellid}.sorted.mod.bam"),
        path("${meta.sampleid}.${meta.flowcellid}.sorted.mod.bam.bai")
        )

        script:
        def alignment  = "${meta.sampleid}.${meta.flowcellid}.mod.bam"
        def sorted_aln = "${meta.sampleid}.${meta.flowcellid}.sorted.mod.bam"
        def aln_index  = "${meta.sampleid}.${meta.flowcellid}.sorted.mod.bam.bai"

        """
        samtools fastq -T '*' \\
        ${reads} | \\
        minimap2 -a -x map-ont --rmq=yes --MD --cs -L -y -t ${task.cpus} ${genome} - | \\
        samtools view -bh > ${alignment};\\
        samtools sort -@ ${task.cpus} \\
        -o ${sorted_aln} \\
        ${alignment};\\
        samtools index -@ ${task.cpus} \\
        -o ${aln_index} \\
        ${sorted_aln}
        """
        
        stub:

        """
        touch ${meta.sampleid}.${meta.flowcellid}.sorted.mod.bam
        touch ${meta.sampleid}.${meta.flowcellid}.sorted.mod.bam.bai
        """
}