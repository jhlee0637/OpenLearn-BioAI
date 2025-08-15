workflow {
    ref_fasta=Channel.fromPath(params.ref_fasta)
    ref_gtf=Channel.fromPath(params.ref_gtf)
    fastq_ch=Channel.fromFilePairs(params.reads)
    strand=Channel.of(params.strand)

    ref_fasta.view()
    ref_gtf.view()
    strand.view()
    fastq_ch.view()
}