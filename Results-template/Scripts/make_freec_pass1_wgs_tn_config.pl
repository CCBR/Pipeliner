#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

#my $mergedmaf = $ARGV[1] . '_out/oncotator_out/' . $ARGV[1] . '_merged.maf'; #to fix...
#open C, ">$mergedmaf";

my $outfile = $ARGV[0] . '/freec_wgs_config.txt';
my $chrLenFile = $ARGV[1];
my $chrFiles = $ARGV[2];
my $tumormateFile = $ARGV[3];
my $controlmateFile = $ARGV[4];
my $makePileup = $ARGV[5];
my $fastaFile = $ARGV[6];
my $SNPfile = $ARGV[7];

open C, ">$outfile";

print C '[general]' . "\n\n";
 
print C "chrLenFile = $chrLenFile\n"
print C "ploidy = 2,3,4,5,6\nbreakPointThreshold = 0.8\nwindow = 1000\n";
print C "chrFiles = $chrFiles\n";
print C "minimalSubclonePresence = 20\ncontaminationAdjustment = TRUE\nmaxThreads = 24\nnumberOfProcesses = 24\n";
print C "outputDir = $ARGV[0]\n";
 
print C '[sample]' . "\n\n";
 
print C "mateFile = $tumormateFile\n";
print C "inputFormat = BAM\nmateOrientation = FR\n\n";
 
print C '[control]' . "\n\n";
 
print C "mateFile = $controlmateFile\n";
print C "inputFormat = BAM\nmateOrientation = FR\n";
 
print C '[target]' . "\n\n";

print C '[BAF]' . "\n\n";

print C "makePileup = $makePileup\n";
print C "fastaFile = $fastaFile\n";
print C "minimalCoveragePerPosition = 20\nminimalQualityPerPosition = 20\n";
print C "SNPfile = $SNPfile";
