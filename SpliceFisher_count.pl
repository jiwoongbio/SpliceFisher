# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

if(scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl SpliceFisher_count.pl [name=output.filtered.txt [...]]

EOF
}
my (@fileList) = @ARGV;
foreach(map {defined($_->[1]) ? $_ : [$_->[0], $_->[0]]} map {[split(/=/, $_, 2)]} @fileList) {
	my ($name, $file) = @$_;
	my %changeCountHash = ();
	open(my $reader, $file);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		$changeCountHash{$tokenHash{'change'}} += 1;
	}
	close($reader);
	print join("\t", $name, map {defined($_) ? $_ : 0} @changeCountHash{-1, 0, 1}), "\n";
}
