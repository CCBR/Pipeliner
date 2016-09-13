#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $notrimSAM = '';
my $trimSAM= '';

my $outfile = 'exome_targets.bed'; #to fix...
open C, ">$outfile";
#print C "Chromosome\tRead_1_Start\tRead_2_Start\tHighCounts\tLowCounts\tHomozygous_informative\tHeterozygous_informative\tComps\tHighProp\tLowProp\tWinner\n";

my $infile = $ARGV[0]; #to fix...
my @line = ();

open U, "<$infile";
while (<U>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	if ($line[0] ne '#'){
#		$chrom = ((split /hr/, $line[0])[1]);
		print C "$line[0]\t$line[1]\t$line[2]\t$line[3]\0\t.\n";
	}
}