# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use List::Util qw(sum);
#use Bio::DB::Sam;
use Getopt::Long;

GetOptions('q=i' => \(my $minimumMappingQuality = 0), 's=s' => \(my $stranded = ''));
my ($regionFile, $geneFile, @bamFilesList) = @ARGV;
my @samListList = ();
foreach my $bamFiles (@bamFilesList) {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} split(/,/, $bamFiles);
	push(@samListList, \@samList);
}

my %geneReadCountListListHash = ();
{
	open(my $reader, $geneFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my ($gene, @readCountsList) = @tokenList;
		$geneReadCountListListHash{$gene} = [map {[split(/,/, $_)]} @readCountsList];
	}
	close($reader);
}
my @columnList = ('chromosome', 'start', 'end', 'strand', 'gene');
push(@columnList, map {"count$_->[1]$_->[2]_$_->[0]"} getCombinationList([map {$_ + 1} 0 .. $#samListList], ['Head', 'Tail', 'Body'], [1, 2]));
print join("\t", @columnList), "\n";
open(my $reader, $regionFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my ($chromosome, $start, $end, $strand, $gene) = @tokenList;

	my @readCountListListList = map {[getExonReadCountListList($chromosome, $start, $end, $strand, @$_)]} @samListList;
	my @geneReadCountListList = @{$geneReadCountListListHash{$gene}};
	foreach my $index (0 .. $#samListList) {
		$readCountListListList[$index]->[5]->[$_] = $geneReadCountListList[$index]->[$_] - $readCountListListList[$index]->[4]->[$_] foreach(0 .. $#{$samListList[$index]});
	}

	my @geneReadCountList = map {sum(@$_)} @geneReadCountListList;
	my @regionReadCountList = map {sum(@{$_->[4]})} @readCountListListList;
	my $regionLength = $end - $start + 1;
	print join("\t", $chromosome, $start, $end, $strand, $gene, (map {join(',', @$_)} map {@$_} @readCountListListList)), "\n";
}
close($reader);

sub getExonReadCountListList {
	my ($chromosome, $start, $end, $strand, @samList) = @_;
	my @readCountListList = ();
	foreach my $index (0 .. $#samList) {
		$readCountListList[$_]->[$index] = 0 foreach(0 .. 3);
		my %readCountHash = ();
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
			my @junctionStartEndList = ($alignment->cigar_str =~ /[0-9]+N/) ? getJunctionStartEndList($alignment->start, $alignment->cigar_str) : ();
			if(grep {$_->[0] <= $start && $end <= $_->[1]} @junctionStartEndList) {
				$readCountListList[1]->[$index] += 1;
				$readCountListList[3]->[$index] += 1;
			} else {
				$readCountListList[0]->[$index] += 1 if(grep {$_->[1] == $start - 1} @junctionStartEndList);
				$readCountListList[2]->[$index] += 1 if(grep {$_->[0] == $end + 1} @junctionStartEndList);
			}
			$readCountHash{$alignment->qname} += 1 unless(grep {$_->[0] <= $start && $end <= $_->[1]} @junctionStartEndList);
		}
		$readCountListList[4]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
	}
	return @readCountListList;
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

sub getCombinationList {
	return if(scalar(@_) == 0);
	return map {[$_]} @{$_[0]} if(scalar(@_) == 1);
	my @combinationList = ();
	foreach my $element (@{$_[0]}) {
		push(@combinationList, map {[$element, @$_]} getCombinationList(@_[1 .. $#_]));
	}
	return @combinationList;
}
