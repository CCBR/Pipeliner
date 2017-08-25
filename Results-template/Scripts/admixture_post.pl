#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $key = $ARGV[0];
open C, ">$ARGV[1]";

my $ancestry = $ARGV[2];
my $ref = $ARGV[3];
my $ped = $ARGV[4];
my @line = ();
my @sample=();
my @pop=();
my @pops=();
my @name=();

if (($ref eq 'hg19') || ($ref eq 'hg38')) {
	@pops = ('EUR','EAS','AMR','SAS','AFR');
}
elsif ($ref eq 'mm10') {
	@pops = ('129P2_OlaHsd','129S1_SvImJ','129S5SvEvBrd','AKR_J','A_J','BALB_cJ','BTBR_T__Itpr3tf_J','BUB_BnJ','C3H_HeH','C3H_HeJ','C57BL_10J','C57BL_6NJ','C57BR_cdJ','C57L_J','C58_J','CAST_EiJ','CBA_J','DBA_1J','DBA_2J','FVB_NJ','I_LnJ','JF1_MsJ','KK_HiJ','LEWES_EiJ','LG_J','LP_J','MOLF_EiJ','NOD_ShiLtJ','NZB_B1NJ','NZO_HlLtJ','NZW_LacJ','PWK_PhJ','RF_J','SEA_GnJ','SJL_J','SPRET_EiJ','ST_bJ','WSB_EiJ','ZALENDE_EiJ');
}

print C "Sample";

my $d=0;

for ($d=0; $d<@pops; $d++) {
	print C "\t$pops[$d]";
}
print C "\n";

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

open I, "$ped";
while (<I>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne '#') {
		push @name, $line[0];
	}
}

my $step=0;

open G, "<$ancestry";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne '#') {
		my $a=0;
		my $status=0;
		for ($a=0; $a<@sample; $a++) {
			if ($name[$step] eq $sample[$a]) {
				$status++;
				last;
			}
		}
		if ($status == 0) {
			my $e=0;
			print C "$name[$step]";
			for ($e=0; $e<@line; $e++) {
				print C "\t$line[$e]";
			}
			print C "\n";
		}
	}
	$step++;
}