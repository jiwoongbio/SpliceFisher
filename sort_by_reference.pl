# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Getopt::Long;

my ($append, $comment, $header) = (0, 0, 0);
GetOptions('a' => \$append, 'c' => \$comment, 'h' => \$header);
my ($tableFile, $referenceFastaFile, $chromosomeIndex, @positionIndexList) = @ARGV;
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
chomp(my $hostname = `hostname`);
system(sprintf('rm -rf %s', getTemporaryFile('*')));
my %chromosomeWriterHash = ();
open(my $reader, $tableFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/ && (!$comment || (print "$line\n")));
	next if($header && !($header = 0) && (print "$line\n"));
	my @tokenList = split("\t", $line, -1);
	my $chromosome = $tokenList[$chromosomeIndex];
	$chromosome =~ s/[ |>;()\$].*$//;
	if(defined($chromosomeWriterHash{$chromosome})) {
		print {$chromosomeWriterHash{$chromosome}} "$line\n";
	} else {
		open(my $writer, ($append ? '>>' : '>'), getTemporaryFile($chromosome));
		print $writer "$line\n";
		if($append) {
			close($writer);
		} else {
			$chromosomeWriterHash{$chromosome} = $writer;
		}
	}
}
close($reader);
foreach my $chromosome (keys %chromosomeWriterHash) {
	close($chromosomeWriterHash{$chromosome});
}
chomp(my @chromosomeList = `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null | xargs cat | cut -f1`);
chomp(@chromosomeList = `grep '^>' $referenceFastaFile | sed 's/^>//'`) unless(@chromosomeList);
s/[ |>;()\$].*$// foreach(@chromosomeList);
my @sortOptionList = map {sprintf('-k%d,%dn', $_ + 1, $_ + 1)} @positionIndexList;
@sortOptionList = (sprintf('-k%d,%d', $chromosomeIndex + 1, $chromosomeIndex + 1), @sortOptionList);
foreach my $chromosome (@chromosomeList) {
	system("sort --field-separator='\t' @sortOptionList " . getTemporaryFile($chromosome)) if(-e getTemporaryFile($chromosome));
}
system(sprintf('rm -rf %s', getTemporaryFile('*')));

sub getTemporaryFile {
	my ($chromosome) = @_;
	return "$temporaryDirectory/sort_by_reference.$hostname.$$.$chromosome";
}
