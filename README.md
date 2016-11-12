# SpliceFisher
Multiple Fisher's exact tests for differential alternative splicing detection using RNA-seq data


Method
------

![SpliceFisher.method.png](SpliceFisher.method.png)

1. Counting reads from BAM files
  - Exon skipping
    - *a*, *b* exon-junction reads
    - *c* exon-skipping reads
    - *d* exon-mapping reads
    - *e* gene-mapping reads
  - Intron retention
    - *a*, *b* exon-intron reads
    - *c* exon-exon reads
    - *d* intron-mapping reads
    - *e* gene-mapping reads
2. Fisher's exact test
  - (control *a* / control *c*) / (test *a* / test *c*)
  - (control *b* / control *c*) / (test *b* / test *c*)
  - (control *d* / control *e*) / (test *d* / test *e*)
3. Adjustment of p-values by false discovery rate (FDR) method


Requirements
------------

1. Perl - https://www.perl.org
2. Perl module "Bio::DB::Sam" - http://search.cpan.org/~lds/Bio-SamTools-1.43/lib/Bio/DB/Sam.pm
3. R (Rscript) - https://www.r-project.org
4. Common linux commands: bash, cut, sort, uniq, ...


Install
-------

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.
```
git clone https://github.com/jiwoongbio/SpliceFisher.git
```


Usages
------

1. Preparation of exon/intron regions
  ```
  ./prepare.sh <gene.gtf>
  ```
  - Example: Before running the following command, download "Homo_sapiens_UCSC_hg19.tar.gz" from http://support.illumina.com/sequencing/sequencing_software/igenome.html and unzip the file.
  ```
  ./prepare.sh Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  ```
  - The repository already contains prebuilt exon/intron regions in hg19 (exon.unique.txt, intron.unique.txt).

2. Detection of differential alternative splicing
  ```
  ./SpliceFisher.sh <outputPrefix> <alpha> <control.bam> <test.bam>
  ```
  - Example: It uses SRA data of https://www.ncbi.nlm.nih.gov/pubmed/24840202. The reads mapped to 11p13 were extracted to generate the input BAM files.
  ```
  ./SpliceFisher.sh hnRNPM_KD 0.05 Control_rep1.bam,Control_rep2.bam hnRNPM_KD_rep1.bam,hnRNPM_KD_rep2.bam
  ```
  - Example result: exon skipping and intron retention of CD44

![hnRNPM_KD.CD44.png](hnRNPM_KD.CD44.png)
