#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $key = $ARGV[0];
open C, ">$ARGV[1]";

my $pedfile = $ARGV[2]; #to fix...
my @line = ();
my @sample=();
my @pop=();

open H, "<$key";
while (<H>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne '#') {
		push @sample, $line[0];
		push @pop, $line[1];
	}
}

open G, "<$pedfile";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne '#') {
		my $a=0;
		my $status=0;
		for ($a=0; $a<@sample; $a++) {
			if ($line[0] eq $sample[$a]) {
				print C "$pop[$a]\n";
				$status++;
				last;
			}
		}
		if ($status == 0) {
			print C '-' . "\n";
		}
	}
}