#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT
my $infile=shift;
my $position = '';
my $indel = ''; #amount of time (seconds) to wait between pings to AVIA
my $stop = '';
my @line = ();

open C, ">variants.bed";

open D, "<$infile";
while (<D>){
	chomp;
 	last if m/^$/;
 	@line = split;
  	next if (($line[0] =~ m/#/) || ($line[0] =~ m/GL00/));
	print C "$line[0] ";
  	if ((length$line[3] == 1) && (length$line[4] == 1)) {
		print C "$line[1] $line[1] $line[3] $line[4]\n";
	}
	elsif (length$line[3] > 1) {
		$indel = substr $line[3], 1;
		$position = ($line[1] + 1);
		$stop = (($position + length$indel) - 1);
		print C "$position $stop $indel -\n";
	}
	elsif (length$line[4] > 1) {
		if ($line[4] !~ m',') {
			$indel = substr $line[4], 1;
			$position = ($line[1] + 1);
			print C "$position $position - $indel\n";
		}
		else {
			print C "$line[1] $line[1] $line[3] $line[4]\n";
		}
	}
}