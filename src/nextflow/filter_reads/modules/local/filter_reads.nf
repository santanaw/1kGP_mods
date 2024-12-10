// Filter reads in mod bam
process FILTER_READS {
        tag "$meta.flowcellid"
        memory '8 GB'
        cpus '1'
        time '23h'

        input:
        tuple(val(meta), path(reads))

        output:
        tuple(val(meta), path("${meta.sampleid}.${meta.flowcellid}.filtered.mod.pass.bam"))

        script:

        def args           = task.ext.args   ?: ""  
        def filtered_reads = "${meta.sampleid}.${meta.flowcellid}.filtered.mod.bam"
        
        """
        ont_tools.py filter_bam \\
        ${args} \\
        -i ${reads} \\
        -o ${filtered_reads}
        """
        
        stub:
        def filtered_reads = "${meta.sampleid}.${meta.flowcellid}.filtered.mod.pass.bam"

        """
        touch ${filtered_reads}
        """
}