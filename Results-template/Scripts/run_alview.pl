#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $outfile = $ARGV[3] . '_variants_alview.input'; #to fix...
open C, ">$outfile";
#print C "Chromosome\tRead_1_Start\tRead_2_Start\tHighCounts\tLowCounts\tHomozygous_informative\tHeterozygous_informative\tComps\tHighProp\tLowProp\tWinner\n";

my $maffile = $ARGV[0]; #to fix...
my @line = ();
my $path = $ARGV[4];
my $sortmaf = $ARGV[2] . $ARGV[3] . '_variants_sorted.maf';
my $fixedmaf = $ARGV[2] . $ARGV[3] . '_variants_fixed.maf';
my $sample = 'null';
my $muts=0;
my $name = 'null';
my $gene = 'null';
my @type = ();
my @full = ();
my $a = 0;
my @chrom=();
my @position=();
my @tumorsample=();
my @germsample=();

my $cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $maffile . ' > ' . $fixedmaf;
system($cmd);
$cmd = 'sort -k1,1 -k17,17 ' . $fixedmaf . ' > ' . $sortmaf;
system($cmd);
#$cmd = 'sort -k1,1 -k17,17 ' . $maffile . ' | awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' > ' . $sortmaf;

open G, "<$sortmaf";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	next if (($line[0] =~ m'#') || ($line[0] =~ m"Hugo_Symbol") || ($line[8] =~ m"Silent") || ($line[8] =~ m"IGR") || ($line[8] =~ m"Intron"));
	if ($line[0] !~ m'#') {
		if (($line[0] eq $gene) || ($gene eq 'null')) {
			if (($line[16] ne $sample) || ($sample eq 'null')) {
				$gene = $line[0];
				$sample = $line[16];
				push @chrom, $line[4];
				push @position, $line[5];
				push @tumorsample, $line[15];
				push @germsample, $line[16];
				push @type, $line[8];
				push @full, $_;
				$muts++;
			}
			else {
				push @chrom, $line[4];
				push @position, $line[5];
				push @tumorsample, $line[15];
				push @germsample, $line[16];
				push @type, $line[8];
				push @full, $_;
			}
		}
		else {
			if ($muts > 1) {
				$a = 0;
				for ($a = 0; $a < @type; $a++) {
					print C "$chrom[$a]\t$position[$a]\t" . $path . '/' . "$tumorsample[$a]\t$gene\t$muts\t$type[$a]\n";
					print C "$chrom[$a]\t$position[$a]\t" . $path . '/' . "$germsample[$a]\t$gene\t$muts\t$type[$a]\n";
				}
			}
			@type = ();
			@full = ();
			@chrom = ();
			@position=();
			@tumorsample=();
			@germsample=();
			$gene = $line[0];
			$sample = $line[16];
			push @chrom, $line[4];
			push @position, $line[5];
			push @tumorsample, $line[15];
			push @germsample, $line[16];				
			push @type, $line[8];
			push @full, $_;
			$muts=1;
		}
	}
}
$a = 0;
for ($a = 0; $a < @type; $a++) {
	print C "$chrom[$a]\t$position[$a]\t" . $path . '/' . "$tumorsample[$a]\t$gene\t$muts\t$type[$a]\n";
	print C "$chrom[$a]\t$position[$a]\t" . $path . '/' . "$germsample[$a]\t$gene\t$muts\t$type[$a]\n";
}
close C;
close G;