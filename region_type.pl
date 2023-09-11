#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'S=s' => \(my $stranded = ''),
);

my ($gtfFile) = @ARGV;
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2n -k3,3n");
{
	my @columnList = ('chromosome', 'start', 'end', 'strand', 'gene_id', 'transcript_id');
	chomp(my $directory = `dirname $0`);
	my ($geneTranscriptChromosomeStrand, @exonStartEndList) = ('');
	open(my $reader, "perl $directory/gtf_extract.pl -f feature=exon $gtfFile @columnList | sort --field-separator='\t' -k1,1 -k4,4 -k5,5 -k6,6 -k2,2n -k3,3n |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		next if(($stranded eq 'f' || $stranded eq 'forward') && $tokenHash{'strand'} eq '-');
		next if(($stranded eq 'r' || $stranded eq 'reverse') && $tokenHash{'strand'} eq '+');
		if($geneTranscriptChromosomeStrand ne (my $currentGeneTranscriptChromosomeStrand = join("\t", @tokenHash{'gene_id', 'transcript_id', 'chromosome', 'strand'}))) {
			printExonIntron(split(/\t/, $geneTranscriptChromosomeStrand), @exonStartEndList) if($geneTranscriptChromosomeStrand ne '');
			($geneTranscriptChromosomeStrand, @exonStartEndList) = ($currentGeneTranscriptChromosomeStrand);
		}
		push(@exonStartEndList, [@tokenHash{'start', 'end'}]);
	}
	close($reader);
	printExonIntron(split(/\t/, $geneTranscriptChromosomeStrand), @exonStartEndList) if($geneTranscriptChromosomeStrand ne '');

	sub printExonIntron {
		my ($gene, $transcript, $chromosome, $strand, @exonStartEndList) = @_;
		my @intronStartEndList = map {[$exonStartEndList[$_ - 1]->[1] + 1, $exonStartEndList[$_]->[0] - 1]} 1 .. $#exonStartEndList;
		if($strand eq '-') {
			@exonStartEndList = reverse(@exonStartEndList);
			@intronStartEndList = reverse(@intronStartEndList);
		}
		print $writer join("\t", $chromosome, @{$exonStartEndList[$_]}, $strand, $gene, $transcript, 'exon', $_ + 1), "\n" foreach(0 .. $#exonStartEndList);
		print $writer join("\t", $chromosome, @{$intronStartEndList[$_]}, $strand, $gene, $transcript, 'intron', $_ + 1), "\n" foreach(0 .. $#intronStartEndList);
	}
}
close($writer);
{
	my @columnList = ('chromosome', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'region', 'number');
	my @tokenHashList = ();
	my ($chromosome, $start, $end) = ('');
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		while(@tokenHashList && ($tokenHash{'chromosome'} ne $chromosome || $tokenHashList[0]->{'end'} < $tokenHash{'start'})) {
			$end = $tokenHashList[0]->{'end'};
			print join("\t", $chromosome, $start, $end, getRegionType()), "\n" if($start <= $end);
			$start = $end + 1;
			shift(@tokenHashList);
		}
		if($tokenHash{'chromosome'} eq $chromosome) {
			$end = $tokenHash{'start'} - 1;
			print join("\t", $chromosome, $start, $end, getRegionType()), "\n" if($start <= $end);
			$start = $end + 1;
		} else {
			($chromosome, $start) = @tokenHash{'chromosome', 'start'};
		}
		push(@tokenHashList, \%tokenHash);
		@tokenHashList = sort {$a->{'end'} <=> $b->{'end'}} @tokenHashList;
	}
	while(@tokenHashList) {
		$end = $tokenHashList[0]->{'end'};
		print join("\t", $chromosome, $start, $end, getRegionType()), "\n" if($start <= $end);
		$start = $end + 1;
		shift(@tokenHashList);
	}

	sub getRegionType {
		if(@tokenHashList) {
			foreach my $region (map {$_->{'region'}} @tokenHashList) {
				return $region if($region eq 'exon');
			}
			return 'intron';
		} else {
			return 'intergenic';
		}
	}
}
close($reader);
waitpid($pid, 0);
