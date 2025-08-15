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

process STAR_INDEX{
    publishDir "INDEX", mode: 'copy'

    input:
        path(fasta)
        path(gtf)

    output:
        path "*", emit:index
    
    script:
        """
        STAR --runThreadN 8 \\
             --runMode genomeGenerate \\
             --genomeDir index \\
             --genomeFastaFiles ${fasta} \\
             --sjdbGTFfile ${gtf} \\
             --genomeSAindexNbases 12
        """
}

process STAR_MAPPING{
    publishDir "MAPPING", mode:'copy'
    cpus 6  //my pc: 8 cores

    input:
        tuple val(sampleid), path(read1), path(read2), path(index)
    
    output:
        path "*"
        path "*.bam", emit:bams

    script:
        """
        STAR --runThreadN 4 --genomeDir ${index} \\
             --readFilesIn ${read1} ${read2} \\
             --outSAMtype BAM SortedByCoordinate \\
             --outFileNamePrefix ${sampleid} \\
             --readFilesCommand zcat
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

    STAR_INDEX(ref_fasta, ref_gtf).set{star_index}
    trimmed.trimmed.map{read1, read2 -> tuple("${read1.getFileName()}".split("_trimmed")[0], read1, read2)} | combine(star_index.index) | STAR_MAPPING | set{bams}
}