RNA-seq Data Analysis Tutorial Using Nextflow    
[한국어](./README_KR.md)

### Reference
1) https://youtu.be/78YDqvg2ADA?si=7s-QkpVHaBQ_h-BD
2) https://github.com/vappiah

## Analysis Overview
1. FASTQ data download
2. Trimming FASTQ
3. Comparing Trimmed FASTQ <-> Raw FASTQ 
4. Mapping Trimmed FASTQ 
5. Reading Mapped FASTQ 

## My Execution Environment
#### Container Environment
|name|version|memo|
|:---:|:---:|:---|
|OS|Rocky linux 9.3|
|conda|25.3.1|miniconda<div>[conda channel priority modification](https://bioconda.github.io/#usage)

#### Conda env Environment
|name|version|how did I install|memo|
|:---:|:---:|:---|:---|
|python|3.9.21|
|java|17.0.10|[nextflow guide](https://www.nextflow.io/docs/latest/install.html#requirements)| unzip, zip installation required
|nextflow|25.4.2|`pip install nextflow`
|nf-core|3.2.1|`pip install nf-core`|git installation required
|fastqc|0.12.1|`conda install bioconda::fastqc`|If errors occur during QC process<div> consider `yum install freetype freetype-devel`
|multiqc|1.19|`conda install bioconda::multiqc`|
|trim_galore|0.6.10|`conda install bioconda::trim-galore` | libxcrypt installation required|
|STAR|2.7.10b|`conda install bioconda::star`|
|subread|2.1.1|`conda install bioconda::subread`|

## FASTQ Data Download
 * Download approximately 2GB file

```bash
#!/bin/bash
#WARNING: This script is for education purposes. You modify this script at your own risk

#we are interested in the raw-data(samples) and reference genome(fasta and gff) and phenotype information. 
#So we will put them in a directory called resources.All other files will be removed

wget -c ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
tar xvfz chrX_data.tar.gz


mv -v chrX_data/samples ./fastq 
mkdir ref
mv chrX_data/genes/chrX.gtf ref
mv chrX_data/genome/chrX.fa ref


#remove unwanted files
rm -fr chrX_data

echo "Data has been downloaded. Check your current directory " 
```
[Code Reference (@vappiah)](https://github.com/vappiah/bioinformatics-tutorials/blob/main/workflow-managers/nextflow/rna-seq/episode-1-star-feature-count/get_data.sh)  
[FASTQ Data Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC5032908/#:~:text=Downloading%20and%20organizing%20required%20data)  

## Nextflow Pipeline
#### Memo
* In Nextflow, 'planning' is important.
* You need to know in advance which files will be modified for what purpose and applied to which tools.
* The reason for deciding to perform trimming on FASTQ files in this exercise is that after examining the reads of the FASTQ files, we determined that read quality could be improved through trimming. This consideration process is omitted in this exercise.

#### Code Execution

```bash
$ nextflow run 4_FEATURE_COUNT.nf --ref_fasta ref/chrX.fa --ref_gtf ref/chrX.gtf --strand 0 --reads 'fastq/*{1, 2}*fastq.gz'
```
>:bulb: The codes below show the development process from `0_channel_check.nf` to `4_FEATURE_COUNT.nf`.    
To see results quickly, you can run only `4_FEATURE_COUNT.nf`.    
If you want to see results for each code, just replace the `4_FEATURE_COUNT.nf` part in the above command line.

### [0_channel_check.nf](./0_channel_check.nf)
* Check channels to configure the pipeline
### [1_TRIM_GALORE.nf](./1_TRIM_GALORE.nf)
* Trim raw FASTQ files
### [2_QC.nf](./2_QC.nf)
* Perform FASTQC, MultiQC on trimmed FASTQ and raw FASTQ
* Compare FASTQC results
    * QC_REPORT/multiqc_report.html
    * QC_REPORT/FASTQC/{sample_name}.html
### [3_STAR.nf](./3_STAR.nf)
* Perform indexing and mapping using STAR
### [4_FEATURE_COUNT.nf](./4_FEATURE_COUNT.nf)
* Generate exon count data using featureCounts from the subRead package