# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(all);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'r' => \(my $raw = ''),
	'k' => \(my $kruskal = ''),
	'H' => \(my $headOnly = ''),
	'B' => \(my $bodyOnly = ''),
	'T' => \(my $tailOnly = ''),
);
my $prefix = $raw ? 'pvalue' : 'padjust';
$prefix = "kruskal_$prefix" if($kruskal);
my @targetList = ();
push(@targetList, "${prefix}Head") if($headOnly);
push(@targetList, "${prefix}Body") if($bodyOnly);
push(@targetList, "${prefix}Tail") if($tailOnly);
if(all {$_ eq ''} ($headOnly, $bodyOnly, $tailOnly)) {
	push(@targetList, "${prefix}Head");
	push(@targetList, "${prefix}Body");
	push(@targetList, "${prefix}Tail");
}
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl SpliceFisher_filter.pl [options] output.txt alpha > output.filtered.txt

Options: -h       display this help message
         -r       use raw p-values instead of FDR-adjusted p-values
         -k       use Kruskal-Wallis rank sum test p-values
         -H       use head p-values only
         -B       use body p-values only
         -T       use tail p-values only

EOF
}
my ($inputFile, $alpha) = @ARGV;

open(my $reader, $inputFile);
chomp(my $line = <$reader>);
my @columnList = split(/\t/, $line);
print join("\t", @columnList), "\n";
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{@columnList} = split(/\t/, $line);
	if(all {$_ ne "NA" && $_ < $alpha} @tokenHash{@targetList}) {
		print join("\t", @tokenHash{@columnList}), "\n";
	}
}
close($reader);
