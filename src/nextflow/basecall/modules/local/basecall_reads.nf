// Call modified nucleotides from pod5 files
process BASECALL_READS {
        tag "$meta.flowcellid"
        memory '240 GB'
        cpus '24'
        time '167h'

        input:
        tuple(val(meta), path(pod5_dir))

        output:
        tuple(val(meta), path("${meta.sampleid}.${meta.flowcellid}.raw.mod.bam"))

        script:

        def args = task.ext.args   ?: ""  

        """
        dorado basecaller \\
        --emit-moves --device "cuda:all" \\
        --no-trim \\
        ${args} \\
        ${pod5_dir}/ \\
        --modified-bases '5mCG_5hmCG' > \\
        ${meta.sampleid}.${meta.flowcellid}.raw.mod.bam
        """
        
        stub:
        """
        touch ${meta.sampleid}.${meta.flowcellid}.raw.mod.bam
        """
}