#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT
my $headers = $ARGV[0];
my $vcf= $ARGV[1];
my $annotations='full_annot.txt';
my @variants = ();
my $a = 0;
my @line=();

open C, ">variants.database";

open B, "<$headers";
while (<B>){
	chomp;
 	@line = split;
 	last if ($a>0);
	if ($line[0] !~ m/##/) {
	 	print C "$line[0]\t$line[1]\t$line[3]\t$line[4]";	
		for ($a = 9; $a < @line; $a++) {
			print C "\t$line[$a].genotypes";
		}
		print C "\t";
		$a++;
	}
}

open D, "<$vcf";
while (<D>){
	chomp;
 	last if m/^$/;
 	@line = split;
	push @variants, $_;
}

my $b=0;
my $c=0;

open E, "<$annotations";
while (<E>){
	chomp;
 	last if m/^$/;
 	if ($b == 0) {
 		print C "$_\n";
 		$b++;
 	}
 	else {
 		print C "$variants[$c]\t$_\n";
 		$c++;
 	}
 }