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
my $name=$ARGV[0];
my $out = $name . '_readlength.txt';

open G, ">$out";

my $size = '';
my $temp='';

$cmd= 'unzip ' . $name . '.R1_fastqc.zip';
system($cmd);
my $file= $name . '.R1_fastqc/fastqc_data.txt';

open (H, $file);
while (<H>) {
	chomp;
#  	last if m/^$/;
  	if ($_ =~ m'Sequence length') {
#  		print $_;
  		$temp=((split ' ', $_)[2]);
  		if ($temp =~ m'-') {
			$size=((split '-', $temp)[1]);
	  	}
	  	else {
	  		$size=$temp;
	  	}
  	}
}
close H;
print G "$size\n";
close G;