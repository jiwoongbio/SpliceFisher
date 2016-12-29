# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use List::Util qw(sum);
use Bio::DB::Sam;
use Getopt::Long;

GetOptions('q=i' => \(my $minimumMappingQuality = 0), 's=s' => \(my $stranded = ''));
my ($intronFile, $geneReadCountFile, @bamFilesList) = @ARGV;
my @samListList = ();
foreach my $bamFiles (@bamFilesList) {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} split(/,/, $bamFiles);
	push(@samListList, \@samList);
}

my %geneReadCountListListHash = ();
{
	open(my $reader, $geneReadCountFile);
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
open(my $reader, $intronFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my ($chromosome, $start, $end, $strand, $gene) = @tokenList;

	my @readCountListListList = map {[getIntronReadCountListList($chromosome, $start, $end, $strand, @$_)]} @samListList;
	my @geneReadCountListList = @{$geneReadCountListListHash{$gene}};
	foreach my $index (0 .. $#samListList) {
		$readCountListListList[$index]->[5]->[$_] = $geneReadCountListList[$index]->[$_] foreach(0 .. $#{$samListList[$index]});
	}
	print join("\t", $chromosome, $start, $end, $strand, $gene, (map {join(',', @$_)} map {@$_} @readCountListListList)), "\n";
}
close($reader);

sub getIntronReadCountListList {
	my ($chromosome, $start, $end, $strand, @samList) = @_;
	my @readCountListList = ();
	foreach my $index (0 .. $#samList) {
		$readCountListList[$_]->[$index] = 0 foreach(0 .. 3);
		my %readCountHash = ();
		foreach my $alignment ($samList[$index]->get_features_by_location(-seq_id => $chromosome, -start => $start, -end => $end)) {
			next if($alignment->qual < $minimumMappingQuality);
			next if($stranded ne '' && getReadStrand($alignment->flag) ne $strand);
			my @junctionStartEndList = ($alignment->cigar_str =~ /[0-9]+N/) ? getJunctionStartEndList($alignment->start, $alignment->cigar_str) : ();
			if(@junctionStartEndList) {
				$readCountListList[1]->[$index] += 1 if(grep {$_->[0] == $start} @junctionStartEndList);
				$readCountListList[3]->[$index] += 1 if(grep {$_->[1] == $end} @junctionStartEndList);
			} else {
				$readCountListList[0]->[$index] += 1 if($alignment->start < $start);
				$readCountListList[2]->[$index] += 1 if($alignment->end > $end);
			}
			$readCountHash{$alignment->qname} += 1 unless(grep {$_->[0] <= $start && $end <= $_->[1]} @junctionStartEndList);
		}
		$readCountListList[4]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
	}
	@readCountListList[0, 1, 2, 3] = @readCountListList[2, 3, 0, 1] if($strand eq '-');
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

sub getReadStrand {
	my ($flag) = @_;
	if($stranded eq 'f') {
		if($flag & 1) {
			return '+' if(grep {$_ == $flag} (99, 147));
			return '-' if(grep {$_ == $flag} (83, 163));
		} else {
			return '+' if($flag == 0);
			return '-' if($flag == 16);
		}
	}
	if($stranded eq 'r') {
		if($flag & 1) {
			return '+' if(grep {$_ == $flag} (83, 163));
			return '-' if(grep {$_ == $flag} (99, 147));
		} else {
			return '+' if($flag == 16);
			return '-' if($flag == 0);
		}
	}
	return '';
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
