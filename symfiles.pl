#!/usr/bin/perl
#############################################################
##
## Syntax: symfiles [path_to_raw_data] [tab_file] [path_to_work_dir]
##
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

#for my $oldfile (@files){
#	for my $oold (keys %repl){
#		if (($oldfile =~ /$oold[_|\.]/) && (my ($a,$r) = $oldfile =~ /(\.|_)(R\d)/) && ($oldfile =~ /fastq\.gz/)){
#		my $newfile = $repl{$oold} . "." . $r . ".fastq.gz";
#		symlink ("$dir\/$oldfile", "$outd\/$newfile");
#		last;
#		}
#	}
#}

for my $oldfile (@files){
        for my $oold (keys %repl){
                if ($oldfile eq $oold) {
                my $newfile = $repl{$oold};
                symlink ("$dir\/$oldfile", "$outd\/$newfile");
                last;
                }
         }
}       
