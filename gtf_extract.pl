# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Getopt::Long;

my @filterList = ();
GetOptions('f=s' => \@filterList);
@filterList = map {[split(/=/, $_)]} @filterList;

my ($gtfFile, @columnList) = @ARGV;
open(my $reader, $gtfFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/);
	my %tokenHash = ();
	@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
	$tokenHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);

	my %filterHash = ();
	$filterHash{$_->[0]} = $filterHash{$_->[0]} || $tokenHash{$_->[0]} eq $_->[1] foreach(@filterList);
	print join("\t", @tokenHash{@columnList}), "\n" unless(grep {!$_} values %filterHash);
}
close($reader);
