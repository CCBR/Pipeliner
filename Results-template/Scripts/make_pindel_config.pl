#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $i = 0;
my $a=0;
my $cmd = "";
my $chr = 10;

my @AllFiles = ();
my @FasFiles = ();
my @line = ();
my $out = 'pindel_out/pindel_config';

#open directory (SWITCH CHROMOSOME ARM HERE)
opendir DIR, "." or die "couldn't open directory\n";
@AllFiles = readdir(DIR);
closedir DIR;

#find correct input files
for ($i = 0; $i < @AllFiles; $i++){
  if ($AllFiles[$i] =~ m'.recal.bam'){
  	if ($AllFiles[$i] !~ m'.recal.bam.bai'){
    	push (@FasFiles, (split /.recal/, $AllFiles[$i])[0]);
	}
  }
}

open G, ">$out";

for ($a=0; $a<@FasFiles; $a++) {
	my $file='QC/' . $FasFiles[$a] . '.qualimapReport/genome_results.txt';
	my $size = '';
	open (H, $file);
	while (<H>) {
		chomp;
#  		last if m/^$/;
  		if ($_ =~ m'median insert size') {
#  			print $_;
  			$size=((split ' = ', $_)[1]);
  		}
  	}
  	close H;
  	print G "$FasFiles[$a]" . '.recal.bam' . "\t$size\t$FasFiles[$a]\n";
}
close G