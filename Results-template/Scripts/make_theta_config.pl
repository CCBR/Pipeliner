#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $germ = $ARGV[0];
my $tumor = $ARGV[1];
my $sites = $ARGV[2];
my $path = $ARGV[3];

my $tumorconfig = $tumor . '.config';
my $germconfig = $germ . '.config';
open F, ">$path" . '/' . "$tumorconfig";
print F 'SNP_FILE=' . $sites . "\n";
print F 'OUTPUT_PREFIX=' . $tumor . "\n";
print F 'BAM_FILE=' . $tumor . '.recal.bam';
close F;
open H, ">$path" . '/' . "$germconfig";
print H 'SNP_FILE=' . $sites . "\n";
print H 'OUTPUT_PREFIX=' . $germ . "\n";
print H 'BAM_FILE=' . $germ . '.recal.bam';
close H;