# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

if(scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl prepare_exon_pair.pl gene.gtf

EOF
}
chomp(my $directory = `dirname $0`);
my ($gtfFile) = @ARGV;
my ($geneChromosomeStrand, %transcriptExonStartEndListHash) = ('');
my @columnList = ('chromosome', 'start', 'end', 'strand', 'gene_id', 'transcript_id');
open(my $reader, "perl $directory/gtf_extract.pl -f feature=exon $gtfFile @columnList | sort --field-separator='\t' -k1,1 -k4,4 -k5,5 -k6,6 -k2,2n -k3,3n |");
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{@columnList} = split(/\t/, $line);
	if($geneChromosomeStrand ne (my $currentGeneChromosomeStrand = join("\t", @tokenHash{'gene_id', 'chromosome', 'strand'}))) {
		printExonPair(split(/\t/, $geneChromosomeStrand), %transcriptExonStartEndListHash) if($geneChromosomeStrand ne '');
		($geneChromosomeStrand, %transcriptExonStartEndListHash) = ($currentGeneChromosomeStrand);
	}
	push(@{$transcriptExonStartEndListHash{$tokenHash{'transcript_id'}}}, [@tokenHash{'start', 'end'}]);
}
close($reader);
printExonPair(split(/\t/, $geneChromosomeStrand), %transcriptExonStartEndListHash) if($geneChromosomeStrand ne '');

sub printExonPair {
	my ($gene, $chromosome, $strand, %transcriptExonStartEndListHash) = @_;
	my %exonTranscriptListHash = ();
	my %end2startExonListHash = ();
	my %start2endExonListHash = ();
	foreach my $transcript (sort keys %transcriptExonStartEndListHash) {
		my @exonStartEndList = map {[join('-', @$_), @$_]} @{$transcriptExonStartEndListHash{$transcript}};
		push(@{$exonTranscriptListHash{$_->[0]}}, $transcript) foreach(@exonStartEndList);
		foreach my $index (1 .. $#exonStartEndList) {
			push(@{$end2startExonListHash{$exonStartEndList[$index - 1]->[2]}->{$exonStartEndList[$index]->[1]}}, $exonStartEndList[$index]->[0]);
			push(@{$start2endExonListHash{$exonStartEndList[$index]->[1]}->{$exonStartEndList[$index - 1]->[2]}}, $exonStartEndList[$index - 1]->[0]);
		}
	}
	my %exonPairTypeHash = ();
	foreach(map {[@$_{sort {$a <=> $b} keys %$_}]} @end2startExonListHash{sort {$a <=> $b} keys %end2startExonListHash}) {
		my @exonListList = @$_;
		foreach my $index1 (0 .. $#exonListList) {
			foreach my $index2 ($index1 + 1 .. $#exonListList) {
				foreach my $exon1 (@{$exonListList[$index1]}) {
					foreach my $exon2 (@{$exonListList[$index2]}) {
						my $exonPair = join(',', $exon1, $exon2);
						$exonPairTypeHash{$exonPair}->{'end2start'} = 1;
					}
				}
			}
		}
	}
	foreach(map {[@$_{sort {$a <=> $b} keys %$_}]} @start2endExonListHash{sort {$a <=> $b} keys %start2endExonListHash}) {
		my @exonListList = @$_;
		foreach my $index1 (0 .. $#exonListList) {
			foreach my $index2 ($index1 + 1 .. $#exonListList) {
				foreach my $exon1 (@{$exonListList[$index1]}) {
					foreach my $exon2 (@{$exonListList[$index2]}) {
						my $exonPair = join(',', $exon1, $exon2);
						$exonPairTypeHash{$exonPair}->{'start2end'} = 1;
					}
				}
			}
		}
	}
	foreach my $exonPair (keys %exonPairTypeHash) {
		my ($exon1, $exon2) = split(/,/, $exonPair);
		my ($start1, $end1) = split(/-/, $exon1);
		my ($start2, $end2) = split(/-/, $exon2);
		if($exonPairTypeHash{$exonPair}->{'end2start'} && $exonPairTypeHash{$exonPair}->{'start2end'}) {
			unless($start1 <= $end2 && $start2 <= $end1) {
				my $transcript1 = join(',', @{$exonTranscriptListHash{$exon1}});
				my $transcript2 = join(',', @{$exonTranscriptListHash{$exon2}});
				my $pairType = 'MXE';
				print join("\t", $chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, $gene, $transcript1, $transcript2), "\n";
			}
		} elsif($exonPairTypeHash{$exonPair}->{'end2start'}) {
			if($start1 <= $end2 && $start2 <= $end1) {
				my $transcript1 = join(',', @{$exonTranscriptListHash{$exon1}});
				my $transcript2 = join(',', @{$exonTranscriptListHash{$exon2}});
				my $pairType = '';
				$pairType = 'A3SS' if($strand eq '+');
				$pairType = 'A5SS' if($strand eq '-');
				print join("\t", $chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, $gene, $transcript1, $transcript2), "\n";
			}
		} elsif($exonPairTypeHash{$exonPair}->{'start2end'}) {
			if($start1 <= $end2 && $start2 <= $end1) {
				my $transcript1 = join(',', @{$exonTranscriptListHash{$exon1}});
				my $transcript2 = join(',', @{$exonTranscriptListHash{$exon2}});
				my $pairType = '';
				$pairType = 'A5SS' if($strand eq '+');
				$pairType = 'A3SS' if($strand eq '-');
				print join("\t", $chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, $gene, $transcript1, $transcript2), "\n";
			}
		}
	}
}
