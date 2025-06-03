process TRIM_GALORE {
    publishDir "TRIMMED", mode:'copy'
    
    input:
        tuple val(sampleid), path(reads)
    
    output:
        path "*"
        path "*trimmed*.fq.gz", emit:trimmed
    
    script:
        """
        trim_galore --paired -q 20 --gzip --basename ${sampleid}_trimmed ${reads}
        """
}


workflow {
    ref_fasta=Channel.fromPath(params.ref_fasta)
    ref_gtf=Channel.fromPath(params.ref_gtf)
    fastq_ch=Channel.fromFilePairs(params.reads)
    strand=Channel.of(params.strand)

    TRIM_GALORE(fastq_ch).set{trimmed}
}