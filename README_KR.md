Nextflow를 이용한 rna-seq 데이터 분석 따라하기
### Reference
1) https://youtu.be/78YDqvg2ADA?si=7s-QkpVHaBQ_h-BD
2) https://github.com/vappiah

## 분석의 개요
1. FASTQ 데이터 다운로드
2. Trimming FASTQ
3. Trimmed FASTQ <-> Raw FASTQ 비교
4. Mapping Trimmed FASTQ 
5. Reading Mapped FASTQ 

## 내 실행환경
#### 컨테이너 환경    
|name|version|memo|
|:---:|:---:|:---|
|OS|Rocky linux 9.3|
|conda|25.3.1|miniconda<div>[conda 채널 우선순위 변경](https://bioconda.github.io/#usage)

#### Conda env 환경
|name|version|how did I install|memo|
|:---:|:---:|:---|:---|
|python|3.9.21|
|java|17.0.10|[nextflow 가이드](https://www.nextflow.io/docs/latest/install.html#requirements)| unzip, zip 설치 필요
|nextflow|25.4.2|`pip install nextflow`
|nf-core|3.2.1|`pip install nf-core`|git 설치 필요
|fastqc|0.12.1|`conda install bioconda::fastqc`|QC 과정에서 에러나는 경우<div> `yum install freetype freetype-devel` 고려하기
|multiqc|1.19|`conda install bioconda::multiqc`|
|trim_galore|0.6.10|`conda install bioconda::trim-galore` | libxcrypt 설치 필요|
|STAR|2.7.10b|`conda install bioconda::star`|
|subread|2.1.1|`conda install bioconda::subread`|

## FASTQ 데이터 다운로드
 * 약 2Gb 파일 다운로드

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

## Nextflow 파이프라인
#### 메모
* Nextflow에서는 '계획'을 세우는 것이 중요합니다.
* 어떤 파일을 어떤 목적으로 수정하고 툴에 적용할지 미리 알아야 합니다.
* 이 연습에서 FASTQ 파일에 대한 trimming 과정을 수행하기로 결정한 이유는, 해당 FASTQ 파일의 reads를 살펴본 결과 trimming을 통해서 read quality가 향상될 수 있다고 판단했기 때문입니다. 이 연습에서는 그런 고민 과정이 생략되어있습니다.

#### 코드의 실행

```bash
$ nextflow run 4_FEATURE_COUNT.nf --ref_fasta ref/chrX.fa --ref_gtf ref/chrX.gtf --strand 0 --reads 'fastq/*{1, 2}*fastq.gz'
```
>:bulb: 아래의 코드들은 `0_channel_check.nf`부터 시작하여 `4_FEATURE_COUNT.nf` 코드까지 코드가 발전하는 과정을 보여줍니다.    
빠르게 결과를 보기 위해서는 `4_FEATURE_COUNT.nf`만 실행하여도 괜찮습니다.    
각 코드 별 결과를 보고싶다면 위의 명령줄에서 `4_FEATURE_COUNT.nf` 부분만 교체하여 수행하면 됩니다.

### [0_channel_check.nf](./0_channel_check.nf)
* 파이프라인을 구성할 channel을 확인하기
### [1_TRIM_GALORE.nf](./1_TRIM_GALORE.nf)
* raw FASTQ 파일을 trimming하기
### [2_QC.nf](./2_QC.nf)
* trimmed FASTQ, raw FASTQ에 대하여 FASTQC, MultiQC 수행하기
* FASTQC 결과를 비교하기
    * QC_REPORT/multiqc_report.html
    * QC_REPORT/FASTQC/{샘플명}.html
### [3_STAR.nf](./3_STAR.nf)
* STAR를 이용하여 indexing, mapping하기
### [4_FEATURE_COUNT.nf](./4_FEATURE_COUNT.nf)
* subRead 패키지의 featureCounts 사용하여 exon count 데이터 생산하기