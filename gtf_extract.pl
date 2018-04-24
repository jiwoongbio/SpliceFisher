# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(all);

my @filterList = ();
use Getopt::Long qw(:config no_ignore_case);
GetOptions(
	'h' => \(my $help = ''),
	'f=s' => \@filterList,
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl gtf_extract.pl [options] gene.gtf [column [...]] > extract.txt

Options: -h       display this help message
         -f STR   filter e.g. feature=exon

EOF
}
@filterList = map {[split(/=/, $_, 2)]} @filterList;
my ($gtfFile, @columnList) = @ARGV;

open(my $reader, $gtfFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/);
	my %tokenHash = ();
	@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
	$tokenHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);
	print join("\t", @tokenHash{@columnList}), "\n" if(all {$tokenHash{$_->[0]} eq $_->[1]} @filterList);
}
close($reader);
