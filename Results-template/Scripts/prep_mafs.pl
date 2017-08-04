#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $finalmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_merged.maf'; #to fix...
open C, ">$finalmaf";
#print C "Chromosome\tRead_1_Start\tRead_2_Start\tHighCounts\tLowCounts\tHomozygous_informative\tHeterozygous_informative\tComps\tHighProp\tLowProp\tWinner\n";

my $maffile = $ARGV[0]; #to fix...
my $fixedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_variants_fixed.maf';
my @line = ();
#my $path = $ARGV[4];
my $header='null';
my @variants=();
my @freq='';
my $calc='';

my $cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $maffile . ' > ' . $fixedmaf;
system($cmd);

open G, "<$fixedmaf";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split "\t", $_;
	if ($line[0] =~ m'Hugo') {
		if ($header eq 'null') {
			$header=$_;
		}
	}
	elsif ($line[0] !~ m'#') {
		push @variants, $_;
		if ($ARGV[1] eq 'mutect2'){
		$calc=($line[41]/$line[39]);
		push @freq, $calc;
		}
	}
}

if ($ARGV[1] eq 'mutect2'){

print C "$header\ttumor_freq\n";

my $a = 0;
for ($a = 0; $a < @variants; $a++) {
	print C "$variants[$a]\t$freq[$a+1]\n";
}
}

else {

print C "$header\n";

my $a = 0;
for ($a = 0; $a < @variants; $a++) {
	print C "$variants[$a]\n";
}
}
close C;
close G;