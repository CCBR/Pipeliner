#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my @germ = ();
my @tumor = ();
my @line = ();
my $sample = $ARGV[0];
my @chrom = ();
my $knowns = '';
my @samples = ();
my @fusion=();
my $genea='';
my $geneb='';

my $input = '';

my $infile= 'fusioncatcher/fusioninspector/' . $sample . '/' . $sample . "_fusionInspector.input";
open I, ">$infile";
$input = 'fusioncatcher/' . $sample . '/final-list_candidate-fusion-genes.txt';
	
open (J, $input);
while (<J>) {
	chomp;
 		last if m/^$/;
	@line = split;
	next if ($line[0] =~ m'Gene_1_symbol');
	$genea = ucfirst lc $line[0];
	$geneb = ucfirst lc $line[1];
	push @fusion, ($genea . '--' . $geneb);
}
	
my $d = 0;
	
for ($d=0; $d<@fusion; $d++) {
	if ($d == 0) {
		print I "$fusion[$d]\n";	
	}
	else {
		if ($fusion[$d] ne $fusion[$d-1]) {
			print I "$fusion[$d]\n";	
		}
	}
}