# SpliceFisher
Multiple Fisher's exact tests for differential alternative splicing detection using RNA-seq data


Method
------

![method](SpliceFisher.method.png)

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
  - (control *a* x test *c*) / (control *c* / test *a*)
  - (control *b* x test *c*) / (control *c* / test *b*)
  - (control *d* x test *e*) / (control *e* / test *d*)
3. Adjustment of p-values by false discovery rate (FDR) method


Requirements
------------

1. Perl - https://www.perl.org
2. Perl module "Bio::DB::Sam" - http://search.cpan.org/~lds/Bio-SamTools-1.43/lib/Bio/DB/Sam.pm
3. R - https://www.r-project.org
