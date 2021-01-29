# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

codeDir=`dirname $0`
gtfFile=$1

if [ -z "$gtfFile" ]; then
	echo 'Usage: ./prepare.sh <gene.gtf>' 1>&2
	exit 1
fi

type=exon
perl $codeDir/prepare_$type.pl $gtfFile | sort --field-separator=$'\t' -k1,1 -k2,2n -k3,3n > $type.txt
cut -f1,2,3,4,5 $type.txt | sort --field-separator=$'\t' -k1,1 -k2,2n -k3,3n | uniq > $type.unique.txt

type=intron
perl $codeDir/prepare_$type.pl $gtfFile | sort --field-separator=$'\t' -k1,1 -k2,2n -k3,3n > $type.txt
cut -f1,2,3,4,5 $type.txt | sort --field-separator=$'\t' -k1,1 -k2,2n -k3,3n | uniq > $type.unique.txt

type=exon_pair
perl $codeDir/prepare_$type.pl $gtfFile | sort --field-separator=$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5n > $type.txt
