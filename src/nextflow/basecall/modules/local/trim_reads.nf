// Trim reads
process TRIM_READS {
        tag "$meta.flowcellid"
        memory '20 GB'
        cpus '20'
        time '23h'

        input:
        tuple(val(meta), path(reads))

        output:
        tuple(val(meta), path("${meta.sampleid}.${meta.flowcellid}.trim.mod.bam"))

        script:
        
        def args = task.ext.args   ?: ""  
        """
        dorado trim \\
        ${args} \\
        ${meta.sampleid}.${meta.flowcellid}.raw.mod.bam > \\
        ${meta.sampleid}.${meta.flowcellid}.trim.mod.bam
        """
        
        stub:
        """
        touch ${meta.sampleid}.${meta.flowcellid}.trim.mod.bam
        """
}