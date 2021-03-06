# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Sam;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'q=i' => \(my $minimumMappingQuality = 0),
	's=s' => \(my $stranded = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl SpliceFisher_gene.pl [options] exon.txt group1.bam,[...] group2.bam,[...] > gene.count.txt

Options: -h       display this help message
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -s       stranded, "f" or "r"

EOF
}
my ($exonFile, @bamFilesList) = @ARGV;
my @samListList = ();
foreach my $bamFiles (@bamFilesList) {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} split(/,/, $bamFiles);
	push(@samListList, \@samList);
}
my %chromosomeHash = ();
$chromosomeHash{$_} = 1 foreach(map {$_->seq_ids} map {@$_} @samListList);
{
	my ($previousGene, @chromosomeStartEndStrandList) = ('');
	open(my $reader, "sort --field-separator='\t' -k5 $exonFile |");
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $start, $end, $strand, $gene) = @tokenList;
		next unless($chromosomeHash{$chromosome});
		if($gene ne $previousGene) {
			print join("\t", $previousGene, getReadCountsList(@chromosomeStartEndStrandList)), "\n" if($previousGene ne '');
			($previousGene, @chromosomeStartEndStrandList) = ($gene);
		}
		push(@chromosomeStartEndStrandList, [$chromosome, $start, $end, $strand]);
	}
	close($reader);
	print join("\t", $previousGene, getReadCountsList(@chromosomeStartEndStrandList)), "\n" if($previousGene ne '');
}

sub getReadCountsList {
	my @chromosomeStartEndStrandList = @_;
	my @readCountListList = ();
	foreach(@samListList) {
		my @samList = @$_;
		my @readCountList = ();
		foreach my $index (0 .. $#samList) {
			my %readCountHash = ();
			foreach(@chromosomeStartEndStrandList) {
				my ($chromosome, $start, $end, $strand) = @$_;
				foreach my $alignment ($samList[$index]->get_features_by_location(-seq_id => $chromosome, -start => $start, -end => $end)) {
					next if($alignment->qual < $minimumMappingQuality);
					if($stranded eq 'f') {
						if($alignment->flag & 1) {
							next if($strand eq '+' && scalar(grep {$_ == $alignment->flag} (99, 147)) == 0);
							next if($strand eq '-' && scalar(grep {$_ == $alignment->flag} (83, 163)) == 0);
						} else {
							next if($strand eq '+' && $alignment->flag != 0);
							next if($strand eq '-' && $alignment->flag != 16);
						}
					}
					if($stranded eq 'r') {
						if($alignment->flag & 1) {
							next if($strand eq '+' && scalar(grep {$_ == $alignment->flag} (83, 163)) == 0);
							next if($strand eq '-' && scalar(grep {$_ == $alignment->flag} (99, 147)) == 0);
						} else {
							next if($strand eq '+' && $alignment->flag != 16);
							next if($strand eq '-' && $alignment->flag != 0);
						}
					}
					my @junctionStartEndList = getJunctionStartEndList($alignment->start, $alignment->cigar_str);
					$readCountHash{$alignment->qname} += 1 if(scalar(grep {$_->[0] <= $start && $end <= $_->[1]} @junctionStartEndList) == 0);
				}
			}
			$readCountList[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
		}
		push(@readCountListList, \@readCountList);
	}
	return map {join(',', @$_)} @readCountListList;
}

sub getJunctionStartEndList {
	my ($position, $cigar) = @_;
	my @junctionStartEndList = ();
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		$position += $length if($operation eq 'M');
		$position += $length if($operation eq 'D');
		if($operation eq 'N') {
			push(@junctionStartEndList, [$position, $position + $length - 1]);
			$position += $length;
		}
	}
	return @junctionStartEndList;
}
