#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

#my $mergedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_merged.maf'; #to fix...
#open C, ">$mergedmaf";

my $finalmaf = $ARGV[1] . '_out/oncotator_out/final_filtered.maf'; #to fix...
open C, ">$finalmaf";

my $maffile = $ARGV[0]; #to fix...
my $fixedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_variants_fixed.maf';
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

open G, "<$fixedmaf";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split "\t", $_;
	if ($line[0] =~ m'Hugo') {
		if ($header eq 'null') {
			$header=$_;
			print C "$header\ttumor_freq\n";
		}
	}
	elsif ($line[0] !~ m'#') {
		next if ($line[0] =~ /TTN|MUC16|OBSCN|AHNAK2|SYNE1|FLG|MUC5B|DNAH17|PLEC|DST|SYNE2|NEB|HSPG2|LAMA5|AHNAK|HMCN1|USH2A|DNAH11|MACF1|MUC17|DNAH5|GPR98|FAT1|PKD1|MDN1|RNF213|RYR1|DNAH2|DNAH3|DNAH8|DNAH1|DNAH9|ABCA13|SRRM2|CUBN|SPTBN5|PKHD1|LRP2|FBN3|CDH23|DNAH10|FAT4|RYR3|PKHD1L1|FAT2|CSMD1|PCNT|COL6A3|FRAS1|FCGBP|RYR2|HYDIN|XIRP2|LAMA1/);
		if (($line[44] eq '-') || ($line[44] < 2)) {
			if (($line[41] > 2) && ($line[39] > 9)) {
				if ((($line[123] eq '-') || ($line[123] < 0.001)) && (($line[76] eq '-') || ($line[76] < 0.01)) && (($line[99] eq '-') || ($line[99] < 0.001)) && (($line[41]/$line[39]) > 0.05) {
					print C "$_\t";
					if ($line[39] != 0){
						$calc=($line[41]/$line[39]);
						print C "$calc\n";
					}
					else {
						print "0\n";		
					}
				}
			}
		}
	}
}