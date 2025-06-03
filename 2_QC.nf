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

process QC {
    publishDir "QC_REPORT", mode:'copy'

    input:
        path(reads)
    
    output:
        path "*"
    
    script:
        """
        fastqc ${reads} 
        multiqc *fastqc*

        mkdir FASTQC
        mv *fastqc* FASTQC
        """
}

workflow {
    ref_fasta=Channel.fromPath(params.ref_fasta)
    ref_gtf=Channel.fromPath(params.ref_gtf)
    fastq_ch=Channel.fromFilePairs(params.reads)
    strand=Channel.of(params.strand)

    TRIM_GALORE(fastq_ch).set{trimmed}

    raw_fastq=fastq_ch.map{items -> items[1]}.flatten().collect()
    trimmed_fastq=trimmed.trimmed.flatten().collect()
    raw_fastq.mix(trimmed_fastq).collect() | QC
}