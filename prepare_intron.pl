# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

if(scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl prepare_intron.pl gene.gtf

EOF
}
chomp(my $directory = `dirname $0`);
my ($gtfFile) = @ARGV;
my ($geneTranscriptChromosomeStrand, @exonStartEndList) = ('');
my @columnList = ('chromosome', 'start', 'end', 'strand', 'gene_id', 'transcript_id');
open(my $reader, "perl $directory/gtf_extract.pl -f feature=exon $gtfFile @columnList | sort --field-separator='\t' -k1,1 -k4,4 -k5,5 -k6,6 -k2,2n -k3,3n |");
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{@columnList} = split(/\t/, $line);
	if($geneTranscriptChromosomeStrand ne (my $currentGeneTranscriptChromosomeStrand = join("\t", @tokenHash{'gene_id', 'transcript_id', 'chromosome', 'strand'}))) {
		printIntron(split(/\t/, $geneTranscriptChromosomeStrand), @exonStartEndList) if($geneTranscriptChromosomeStrand ne '');
		($geneTranscriptChromosomeStrand, @exonStartEndList) = ($currentGeneTranscriptChromosomeStrand);
	}
	push(@exonStartEndList, [@tokenHash{'start', 'end'}]);
}
close($reader);
printIntron(split(/\t/, $geneTranscriptChromosomeStrand), @exonStartEndList) if($geneTranscriptChromosomeStrand ne '');

sub printIntron {
	my ($gene, $transcript, $chromosome, $strand, @exonStartEndList) = @_;
	if(join(',', map {join('-', @$_)} @exonStartEndList) =~ /^[0-9]+-(.*)-[0-9]+$/) {
		my @intronStartEndList = map {[$_->[0] + 1, $_->[1] - 1]} map {[split(/,/)]} split(/-/, $1);
		@intronStartEndList = reverse(@intronStartEndList) if($strand eq '-');
		print join("\t", $chromosome, @{$intronStartEndList[$_]}, $strand, $gene, $transcript, $_ + 1), "\n" foreach(0 .. $#intronStartEndList);
	}
}
