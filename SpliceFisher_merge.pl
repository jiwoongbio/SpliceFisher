# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(all);

if(scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl merge.pl group1.count.txt,[...] group2.count.txt,[...] > count.txt

EOF
}
my (@filesList) = @ARGV;
my @readerList = ();
my @indexListList = ();
my $index = 0;
foreach my $files (@filesList) {
	my @fileList = split(/,/, $files);
	my @indexList = ();
	foreach my $file (@fileList) {
		open(my $reader, $file);
		push(@readerList, $reader);
		push(@indexList, $index);
		$index += 1;
	}
	push(@indexListList, \@indexList);
}
my $line = '';
my @columnList = ();
{
	my @lineList = map {$line = <$_>} @readerList;
	chomp($_) foreach(@lineList);
	my %lineHash = ();
	$lineHash{$_} = 1 foreach(@lineList);
	@columnList = split(/\t/, $line) if(scalar(($line) = keys %lineHash) == 1);
	s/_[0-9]+$// foreach(@columnList);
}
my @eventColumnList = grep {$_ !~ /^count/} @columnList;
my @countColumnList = grep {$_ =~ /^count/} @columnList;
{
	my @columnList = @eventColumnList;
	foreach my $group (1 .. scalar(@indexListList)) {
		push(@columnList, map {"$_\_$group"} @countColumnList);
	}
	print join("\t", @columnList), "\n";
}
while(all {defined($_)} (my @lineList = map {$line = <$_>} @readerList)) {
	chomp($_) foreach(@lineList);
	my @tokenHashList = ();
	foreach my $line (@lineList) {
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		push(@tokenHashList, \%tokenHash);
	}
	my %eventHash = ();
	$eventHash{$_} = 1 foreach(map {join("\t", @$_{@eventColumnList})} @tokenHashList);
	if(scalar((my $event) = keys %eventHash) == 1) {
		my @countsList = ();
		foreach(@indexListList) {
			my @countListList = map {[@$_{@countColumnList}]} @tokenHashList[@$_];
			foreach my $index (0 .. $#countColumnList) {
				push(@countsList, join(',', map {$_->[$index]} @countListList));
			}
		}
		print join("\t", $event, @countsList), "\n";
	}
}
