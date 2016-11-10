# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

codeDir=`dirname $0`
outputPrefix=$1
alpha=$2
bamFiles=${@:3}

if [ -z "$bamFiles" ]; then
	echo 'Usage: ./SpliceFisher.sh <outputPrefix> <alpha> <control.bam> <test.bam>' 1>&2
	exit 1
fi
perl -MBio::DB::Sam -e '' || exit 1

perl $codeDir/SpliceFisher_gene.pl $codeDir/exon.unique.txt $bamFiles > $outputPrefix.gene.count.txt

type=exon
perl $codeDir/SpliceFisher_$type.pl $codeDir/$type.unique.txt $outputPrefix.gene.count.txt $bamFiles > $outputPrefix.$type.count.txt
Rscript $codeDir/SpliceFisher.R $outputPrefix.$type.count.txt $outputPrefix.$type.txt
perl $codeDir/SpliceFisher_filter.pl $outputPrefix.$type.txt $alpha > $outputPrefix.$type.filtered.txt

type=intron
perl $codeDir/SpliceFisher_$type.pl $codeDir/$type.unique.txt $outputPrefix.gene.count.txt $bamFiles > $outputPrefix.$type.count.txt
Rscript $codeDir/SpliceFisher.R $outputPrefix.$type.count.txt $outputPrefix.$type.txt
perl $codeDir/SpliceFisher_filter.pl $outputPrefix.$type.txt $alpha > $outputPrefix.$type.filtered.txt
