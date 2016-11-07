# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;

chomp(my $directory = `dirname $0`);
my ($gtfFile) = @ARGV;
my ($geneTranscriptChromosomeStrand, @exonStartEndList) = ('');
my @columnList = ('chromosome', 'start', 'end', 'strand', 'gene_id', 'transcript_id');
open(my $reader, "perl $directory/gtf_extract.pl -f feature=exon $gtfFile @columnList | sort --field-separator='\t' -k1,1 -k4 -k2,2n -k3,3n |");
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{@columnList} = split(/\t/, $line);
	if($geneTranscriptChromosomeStrand ne (my $currentGeneTranscriptChromosomeStrand = join("\t", @tokenHash{'gene_id', 'transcript_id', 'chromosome', 'strand'}))) {
		printExon(split(/\t/, $geneTranscriptChromosomeStrand), @exonStartEndList) if($geneTranscriptChromosomeStrand ne '');
		($geneTranscriptChromosomeStrand, @exonStartEndList) = ($currentGeneTranscriptChromosomeStrand);
	}
	push(@exonStartEndList, [@tokenHash{'start', 'end'}]);
}
close($reader);
printExon(split(/\t/, $geneTranscriptChromosomeStrand), @exonStartEndList) if($geneTranscriptChromosomeStrand ne '');

sub printExon {
	my ($gene, $transcript, $chromosome, $strand, @exonStartEndList) = @_;
	@exonStartEndList = reverse(@exonStartEndList) if($strand eq '-');
	print join("\t", $chromosome, @{$exonStartEndList[$_]}, $strand, $gene, $transcript, $_ + 1), "\n" foreach(0 .. $#exonStartEndList);
}
