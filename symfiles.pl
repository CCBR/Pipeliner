#!/usr/bin/perl
##############################################################
##  Script     : symfiles
##  Author     : chenx3
##  Date       : 09/30/2016
##  Last Edited: 10/01/2016, chenx3
##  Description: Create symlinks with sample information for 
##				CCBR_Pipeliner
###############################################################
## Purpose:
## - create a directory called symdata in the raw_data folder, 
##	use provided sample information to create symlink with 
##	changed name. CCBR_Pipeliner then can use these files
##	to the downstream analysis.
## Requirements:
## - tab delimited file, old_file prefix<tab>new_file prefix
## - path to the raw_data directory
## - old file name could contain _001, or other junk
## - but it has to be fastq.gz file, and has R1/R2 information
## - example: To change Julio_1_S19_R1_100.fastq.gz 
##   It will be:	Julio_1_S19<tab>L6_IgG2_1 
##
## Syntax: symfiles [path_to_raw_data] [tab_file] [path_to_work_dir]
##
## Example: perl symfiles.pl /data/ccbr741/raw_data filetab.txt
###############################################################
use strict;

my $dir = $ARGV[0];
my $filetab = $ARGV[1];
my $outd = $ARGV[2];
 
opendir(DH, $dir) || die "Can not open $dir: $!";
my @files=readdir DH;
close DH;

## mkdir "$dir\/symdata", 0755;

open IN, $filetab or die "no input file";
my ($old, $new, %repl);
while(<IN>){
chomp;
($old, $new) = split /\t/;
$repl{$old} = $new;
}
close IN;

for my $oldfile (@files){
	for my $oold (keys %repl){
		if (($oldfile =~ /$oold[_|\.]/) && (my ($a,$r) = $oldfile =~ /(\.|_)(R\d)/) && ($oldfile =~ /fastq\.gz/)){
		my $newfile = $repl{$oold} . "." . $r . ".fastq.gz";
		symlink ("$dir\/$oldfile", "$outd\/$newfile");
		last;
		}
	}
}

