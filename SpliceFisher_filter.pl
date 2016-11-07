use strict;
use warnings;
use Getopt::Long;

GetOptions('noadjust' => \(my $noadjust = 0),
	'head' => \(my $checkHead = 0),
	'body' => \(my $checkBody = 0),
	'tail' => \(my $checkTail = 0),
	'kruskal' => \(my $kruskal = 0));
my $prefix = $noadjust ? 'pvalue' : 'padjust';
$prefix = "kruskal_$prefix" if($kruskal);
my @targetList = ();
push(@targetList, "${prefix}Head") if($checkHead);
push(@targetList, "${prefix}Body") if($checkBody);
push(@targetList, "${prefix}Tail") if($checkTail);
if($checkHead + $checkBody + $checkTail == 0) {
	push(@targetList, "${prefix}Head");
	push(@targetList, "${prefix}Body");
	push(@targetList, "${prefix}Tail");
}
my ($inputFile, $value) = @ARGV;
open(my $reader, $inputFile);
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line);
print join("\t", @columnList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{@columnList} = split(/\t/, $line);
	if(scalar(grep {$_ ne "NA" && $_ < $value} @tokenHash{@targetList}) == scalar(@targetList)) {
		print join("\t", @tokenHash{@columnList}), "\n";
	}
}
close($reader);
