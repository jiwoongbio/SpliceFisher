# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use List::Util qw(sum);
use Bio::DB::Sam;
use Getopt::Long;

GetOptions('q=i' => \(my $minimumMappingQuality = 0), 's=s' => \(my $stranded = ''));
my ($exonPairFile, @bamFilesList) = @ARGV;
my @samListList = ();
foreach my $bamFiles (@bamFilesList) {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} split(/,/, $bamFiles);
	push(@samListList, \@samList);
}

my @columnList = ('chromosome', 'start1', 'end1', 'start2', 'end2', 'strand', 'pairType', 'gene');
push(@columnList, map {"count$_->[1]$_->[2]_$_->[0]"} getCombinationList([map {$_ + 1} 0 .. $#samListList], ['Head', 'Tail', 'Body'], [1, 2]));
print join("\t", @columnList), "\n";
open(my $reader, $exonPairFile);
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line, -1);
	my ($chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, $gene) = @tokenList;

	my @readCountListListList = map {[getExonPairReadCountListList($chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, @$_)]} @samListList;
	print join("\t", $chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, $gene, (map {join(',', @$_)} map {@$_} @readCountListListList)), "\n";
}
close($reader);

sub getExonPairReadCountListList {
	my ($chromosome, $start1, $end1, $start2, $end2, $strand, $pairType, @samList) = @_;
	if($pairType eq 'MXE') { # alternative exons
		($start1, $end1, $start2, $end2) = ($start2, $end2, $start1, $end1) if($strand eq '-');
	} elsif(($pairType eq 'A5SS' && $strand eq '+') || ($pairType eq 'A3SS' && $strand eq '-')) { # alternative ends
		$start2 = $end1 + 1;
	} elsif(($pairType eq 'A5SS' && $strand eq '-') || ($pairType eq 'A3SS' && $strand eq '+')) { # alternative starts
		$end1 = $start2 - 1;
	}
	my @readCountListList = ();
	foreach my $index (0 .. $#samList) {
		$readCountListList[$_]->[$index] = 0 foreach(0 .. 3);
		{
			my %readCountHash = ();
			foreach my $alignment ($samList[$index]->get_features_by_location(-seq_id => $chromosome, -start => $start1, -end => $end1)) {
				next if($alignment->qual < $minimumMappingQuality);
				next if($stranded ne '' && getReadStrand($alignment->flag) ne $strand);
				my @junctionStartEndList = ($alignment->cigar_str =~ /[0-9]+N/) ? getJunctionStartEndList($alignment->start, $alignment->cigar_str) : ();
				if($pairType eq 'MXE') { # alternative exons
					$readCountListList[0]->[$index] += 1 if(grep {$_->[1] == $start1 - 1} @junctionStartEndList); # upstream junction
					$readCountListList[2]->[$index] += 1 if(grep {$_->[0] == $end1 + 1} @junctionStartEndList); # downstream junction
				} elsif(($pairType eq 'A5SS' && $strand eq '+') || ($pairType eq 'A3SS' && $strand eq '-')) { # alternative ends
					$readCountListList[1]->[$index] += 1 if(grep {$_->[0] == $end1 + 1} @junctionStartEndList);
					$readCountListList[2]->[$index] += 1 if(scalar(@junctionStartEndList) == 0 && $alignment->end > $end1);
				} elsif(($pairType eq 'A5SS' && $strand eq '-') || ($pairType eq 'A3SS' && $strand eq '+')) { # alternative starts
					$readCountListList[0]->[$index] += 1 if(grep {$_->[1] == $start1 - 1} @junctionStartEndList);
				}
				$readCountHash{$alignment->qname} += 1 unless(grep {$_->[0] <= $start1 && $end1 <= $_->[1]} @junctionStartEndList);
			}
			if($pairType eq 'MXE') { # alternative exons
				$readCountListList[4]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
			} elsif(($pairType eq 'A5SS' && $strand eq '+') || ($pairType eq 'A3SS' && $strand eq '-')) { # alternative ends
				$readCountListList[3]->[$index] = $readCountListList[1]->[$index];
				$readCountListList[5]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
			} elsif(($pairType eq 'A5SS' && $strand eq '-') || ($pairType eq 'A3SS' && $strand eq '+')) { # alternative starts
				$readCountListList[4]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
			}
		}
		{
			my %readCountHash = ();
			foreach my $alignment ($samList[$index]->get_features_by_location(-seq_id => $chromosome, -start => $start2, -end => $end2)) {
				next if($alignment->qual < $minimumMappingQuality);
				next if($stranded ne '' && getReadStrand($alignment->flag) ne $strand);
				my @junctionStartEndList = ($alignment->cigar_str =~ /[0-9]+N/) ? getJunctionStartEndList($alignment->start, $alignment->cigar_str) : ();
				if($pairType eq 'MXE') { # alternative exons
					$readCountListList[1]->[$index] += 1 if(grep {$_->[1] == $start2 - 1} @junctionStartEndList); # upstream junction
					$readCountListList[3]->[$index] += 1 if(grep {$_->[0] == $end2 + 1} @junctionStartEndList); # downstream junction
				} elsif(($pairType eq 'A5SS' && $strand eq '+') || ($pairType eq 'A3SS' && $strand eq '-')) { # alternative ends
					$readCountListList[0]->[$index] += 1 if(grep {$_->[0] == $end2 + 1} @junctionStartEndList);
				} elsif(($pairType eq 'A5SS' && $strand eq '-') || ($pairType eq 'A3SS' && $strand eq '+')) { # alternative starts
					$readCountListList[1]->[$index] += 1 if(grep {$_->[1] == $start2 - 1} @junctionStartEndList);
					$readCountListList[2]->[$index] += 1 if(scalar(@junctionStartEndList) == 0 && $alignment->start < $start2);
				}
				$readCountHash{$alignment->qname} += 1 unless(grep {$_->[0] <= $start2 && $end2 <= $_->[1]} @junctionStartEndList);
			}
			if($pairType eq 'MXE') { # alternative exons
				$readCountListList[5]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
			} elsif(($pairType eq 'A5SS' && $strand eq '+') || ($pairType eq 'A3SS' && $strand eq '-')) { # alternative ends
				$readCountListList[4]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
			} elsif(($pairType eq 'A5SS' && $strand eq '-') || ($pairType eq 'A3SS' && $strand eq '+')) { # alternative starts
				$readCountListList[3]->[$index] = $readCountListList[1]->[$index];
				$readCountListList[5]->[$index] = scalar(grep {$readCountHash{$_} > 0} keys %readCountHash);
			}
		}
	}
	if($pairType eq 'MXE') { # alternative exons
		@readCountListList[0, 1, 2, 3] = @readCountListList[2, 3, 0, 1] if($strand eq '-');
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
