#!/usr/bin/perl -w
#Estimating parental ancestry proportion in mapping data (20% tails + 25% parental differences)
use strict;
use List::Util qw( min max sum );
#my @c= ('X','2L','2R','3L','3R');
my $i=0;
my $matrix = '';
my $names = '';
my @AllFiles = ();
my $cmd = '';
my $vcf=shift;

$cmd = 'module load vcftools; vcftools --vcf ' . $vcf . ' --plink --remove-indels --out plink';
system($cmd);
$cmd = 'module load plink/1.07; plink --file plink --distance-matrix --out distance --cluster';
system($cmd);
$cmd = 'cut -f 1 -d \' \' distance.cluster2 > samples.txt';
system($cmd);

#open directory (SWITCH CHROMOSOME ARM HERE)
opendir DIR, "." or die "couldn't open directory\n";
@AllFiles = readdir(DIR);
closedir DIR;

#find correct input files
for ($i = 0; $i < @AllFiles; $i++) {
  if ($AllFiles[$i] eq 'distance.mdist') {
    $matrix = $AllFiles[$i];
  }
  if ($AllFiles[$i] eq 'samples.txt') {
    $names = $AllFiles[$i];
  }
}

my $outfile = ">distance.matrix";
open(O,$outfile);

my @samples = ();

open (H, $names);
while (<H>) {
	chomp;
	push @samples, $_;
}
close H;

my @distances = ();

open (C, $matrix);
while (<C>) {
	chomp;
	push @distances, $_;
}
close C;

my @letters = ();

print O scalar@samples . "\n";

for ($a = 0; $a < @samples; $a++) {
	@letters = (split '', $samples[$a]);
	my $chars = 0;
	for ($b = 0; $b < @letters; $b++) {
		print O $letters[$b];
		$chars++;
		last if ($chars == 10);
	}
	while ($chars < 11) {
		print O ' ';
		$chars++;
	}
	print O "$distances[$a]\n";
}
close O;

$cmd = 'module load phylip; neighbor << EOF' . "\n" . 'distance.matrix' . "\n" . "N\nY\nEOF";
system($cmd);
$cmd = 'module load phylip; drawtree << EOF' . "\nouttree\n" . '/data/CCBR/local/lib/fontfile' . "\nP\nW\n3000\n3000\nC\n" . '0.75' . "\nL\nR\nY\nEOF";
system($cmd);
$cmd = 'mv plotfile sample_network.bmp';
system($cmd);