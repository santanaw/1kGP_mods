// Pileup modifications per haplotype
process PILEUP_MODS_HAPLOTYPE {
        tag "$sampleid"
        memory '20 GB'
        cpus '5'
        time '23h'

        input:
        tuple(
        val(sampleid), 
        val(hap_n),
        path(alignment), 
        path(alignment_index), 
        path(genome)
        )

        output:
        tuple( 
        val(sampleid), 
        val(hap_n),
        path("${sampleid}.${hap_n}.genes.pileup.log"),
        path("${sampleid}.${hap_n}.genes.sorted.mod.bedmethyl")        
        )

        script:
        def pileup_log         = "${sampleid}.${hap_n}.genes.pileup.log"
        def pileup_bedmethyl   = "${sampleid}.${hap_n}.genes.sorted.mod.bedmethyl"

        """
        modkit pileup \\
        --sampling-frac 1 --cpg --threads ${task.cpus} --seed 1234 \\
        --ref ${genome} \\
        --log-filepath ${pileup_log} \\
        ${alignment} \\
        ${pileup_bedmethyl}
        """
        
        stub:
        """
        touch ${sampleid}.${hap_n}.genes.pileup.log
        touch ${sampleid}.${hap_n}.genes.sorted.mod.bedmethyl
        """
}
