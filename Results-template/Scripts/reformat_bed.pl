#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $notrimSAM = '';
my $trimSAM= '';

my $outfile = 'exome_targets.bed';
my $freecout = 'freec_targets.bed';
open C, ">$outfile";
open D, ">$freecout";

my $infile = $ARGV[0];
my @line = ();

open U, "<$infile";
while (<U>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	if ($line[0] ne '#'){
#		$chrom = ((split /hr/, $line[0])[1]);
		print C "$line[0]\t$line[1]\t$line[2]\t$line[3]\t0\t.\n";
		if ($ARGV[1] eq 'hg19') {
			print D "chr" . "$line[0]\t$line[1]\t$line[2]\n";
		}
		else {
			print D "$line[0]\t$line[1]\t$line[2]\n";
		}
	}
}