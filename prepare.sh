# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

codeDir=`dirname $0`
fastaFile=$1
gtfFile=$2

if [ -z "$gtfFile" ]; then
	echo 'Usage: ./prepare.sh <genome.fasta> <gene.gtf>' 1>&2
	exit 1
fi

type=exon
perl $codeDir/prepare_$type.pl $gtfFile | perl $codeDir/sort_by_reference.pl - $fastaFile 0 1 2 > $type.txt
cut -f1,2,3,4,5 $type.txt | sort -u | perl $codeDir/sort_by_reference.pl - $fastaFile 0 1 2 > $type.unique.txt

type=intron
perl $codeDir/prepare_$type.pl $gtfFile | perl $codeDir/sort_by_reference.pl - $fastaFile 0 1 2 > $type.txt
cut -f1,2,3,4,5 $type.txt | sort -u | perl $codeDir/sort_by_reference.pl - $fastaFile 0 1 2 > $type.unique.txt
