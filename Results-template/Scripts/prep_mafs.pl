#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $finalmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_merged.maf'; #to fix...
open C, ">$finalmaf";
#print C "Chromosome\tRead_1_Start\tRead_2_Start\tHighCounts\tLowCounts\tHomozygous_informative\tHeterozygous_informative\tComps\tHighProp\tLowProp\tWinner\n";

my $maffile = $ARGV[0]; #to fix...
my @line = ();
#my $path = $ARGV[4];
my $header='null';
my @variants=();
my @freq='';
my $calc='';

open G, "<$maffile";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split "\t", $_;
  	next if ($line[0] =~ m'#');
	if ($line[0] =~ m'Hugo') {
		if ($header eq 'null') {
			$header=$_;
		}
	}
	else {
		push @variants, $_;
		$calc=($line[41]/$line[39]);
		push @freq, $calc;
	}
}

print C "$header\ttumor_freq\n";

my $a = 0;
for ($a = 0; $a < @variants; $a++) {
	print C "$variants[$a]\t$freq[$a]\n";
}
close C;
close G;