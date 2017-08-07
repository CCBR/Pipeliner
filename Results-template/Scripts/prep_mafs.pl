#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $finalmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_merged.maf'; #to fix...
open C, ">$finalmaf";
#print C "Chromosome\tRead_1_Start\tRead_2_Start\tHighCounts\tLowCounts\tHomozygous_informative\tHeterozygous_informative\tComps\tHighProp\tLowProp\tWinner\n";

my $maffile = $ARGV[0]; #to fix...
my $fixedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_variants_fixed.maf';
my $sortedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_variants_fixed_geneSorted.maf';
my @line = ();
my $header='null';
my @variants=();
my @freq='';
my $calc='';
my $gene='null';
my $muts=0;
my @hugoid=();
my @count=();
my $sample='null';

my $cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $maffile . ' > ' . $fixedmaf;
system($cmd);
$cmd = 'sort -k1,1 ' . $fixedmaf . ' > ' . $sortedmaf;
system($cmd);

open H, "<$sortedmaf";
while (<H>){
	chomp;
  	last if m/^$/;
  	@line = split "\t", $_;
  	next if ($line[0] =~ m'#');
  	next if ($line[8] =~ m/RNA|Silent|Intron|IGR|Flank|UTR/);
	if ($line[0] !~ m/Hugo/) {
		if (($line[0] eq $gene) || ($gene eq 'null')) {
			if (($line[15] ne $sample) || ($sample eq 'null')) {
				$muts++;
				$gene = $line[0];
				$sample = $line[15];
			}
		}
		else {
			push @hugoid, $gene;
			push @count, $muts;
			$gene = $line[0];
			$sample = $line[15];
			$muts=1;
		}
	}
}
push @hugoid, $gene;
push @count, $muts;


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
#		if ($ARGV[1] eq 'mutect2'){
		$calc=($line[41]/$line[39]);
		push @freq, $calc;
#		}
	}
}

#if ($ARGV[1] eq 'mutect2'){

print C "$header\ttumor_freq\tNonsilent_gene_mutation_count\n";

my $a = 0;
for ($a = 0; $a < @variants; $a++) {
	my $d=0;
	my $genecount=0;
	for ($d=0; $d<@hugoid; $d++) {
		@line=split "\t", $variants[$a];
		if ($hugoid[$d] eq $line[0]) {
			$genecount=$count[$d];
			last;
		}
	}
	print C "$variants[$a]\t$freq[$a+1]\t$genecount\n";
}
#}

#elsif ($ARGV[1] eq 'mutect'){

#print C "$header\ttumor_freq\tNonsilent_gene_mutation_count\n";

#my $a = 0;
#for ($a = 0; $a < @variants; $a++) {
#	my $d=0;
#	my $genecount=0;
#	for ($d=0; $d<@hugoid; $d++) {
#		@line=split "\t", $variants[$a];
#		if ($hugoid[$d] eq $line[0]) {
#			$genecount=$count[$d];
#			last;
#		}
#	}
#	print C "$variants[$a]\t$freq[$a+1]\t$genecount\n";
#}
#}

#else {

#print C "$header\tNonsilent_gene_mutation_count\n";

#my $a = 0;
#for ($a = 0; $a < @variants; $a++) {
#	my $d=0;
#	my $genecount=0;
#	for ($d=0; $d<@hugoid; $d++) {
#		@line=split "\t", $variants[$a];
#		if ($hugoid[$d] eq $line[0]) {
#			$genecount=$count[$d];
#			last;
#		}
#	}
#	print C "$variants[$a]\t$genecount\n";
#}
#}
close C;
close G;