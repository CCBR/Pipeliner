#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my @germ = ();
my @tumor = ();
my @line = ();
my $sample = $ARGV[0];
my $type = $ARGV[1];
my @chrom = ();
my $knowns = '';
my @samples = ();
my @fusion=();
my $chr = '';
my $pos = '';
my $input = '';

my $outfile= $type . '/oncofuse/' . $sample . '/' . $sample . ".oncofuse.input";
open I, ">$outfile";
$input = $type . '/' . 'fusioninspector/' . $sample . '/' . $sample . '.fusion_predictions.final';
@fusion = ();
	
open (J, $input);
while (<J>) {
	chomp;
 	last if m/^$/;
	@line = split;
	next if ($line[0] =~ m'#FusionName');
	$chr = ((split ':', $line[5])[0]);
	$pos = ((split ':', $line[5])[1]);
	print I "$chr\t$pos\t";
	$chr = ((split ':', $line[8])[0]);
	$pos = ((split ':', $line[8])[1]);
	print I "$chr\t$pos\tAVG\n";
}
close I;
close J;