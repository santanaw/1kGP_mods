// Summary reads
process SUMMARYZE_READS {
        tag "$meta.flowcellid"
        memory '5 GB'
        cpus '1'
        time '23h'

        input:
        tuple(val(meta), path(reads))

        output:
        tuple(val(meta), path("summary.${meta.sampleid}.${meta.flowcellid}.txt.gz"))

        script:

        """
        dorado summary \\
        ${reads} | gzip -  > summary.${meta.sampleid}.${meta.flowcellid}.txt.gz
        """
        
        stub:
        """
        touch summary.${meta.sampleid}.${meta.flowcellid}.txt.gz
        """
}