# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($exonCountFile) = @ARGV;
open(my $reader, $exonCountFile);
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line, -1);
my %indexCountHash = ();
foreach my $column (@columnList) {
	$indexCountHash{$1} += 1 if($column =~ /_([0-9]+)$/);
}
my @indexList = sort {$a <=> $b} keys %indexCountHash;
my @printColumnList = ('chromosome', 'start', 'end', 'strand', 'gene', map {("inclusion_reads_$_", "exclusion_reads_$_", "percent_spliced_in_$_")} @indexList);
print join("\t", @printColumnList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{@columnList} = split(/\t/, $line, -1);
	foreach my $index (@indexList) {
		$tokenHash{"inclusion_reads_$index"} = $tokenHash{"countBody1_$index"};
		$tokenHash{"exclusion_reads_$index"} = $tokenHash{"countTail2_$index"};
		my $divisor = $tokenHash{"inclusion_reads_$index"} + $tokenHash{"exclusion_reads_$index"};
		$tokenHash{"percent_spliced_in_$index"} = $divisor > 0 ? $tokenHash{"inclusion_reads_$index"} / $divisor : 0;
	}
	print join("\t", @tokenHash{@printColumnList}), "\n";
}
close($reader);
