#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my @format = ();
my $files = 'calls';
my $newformat = '';
my $first = '';
my $depth = '';
my $ref = '';
my $alt = '';
my $knowns = '';
my $newsample = '';
my @line = ();
my @file=();
my @name = ();
my $loc = '';

my $cmd = '';
$cmd = 'ls *calls.cns > calls';
system($cmd);

open (H, $files);
while (<H>) {
	chomp;
  	last if m/^$/;
	@line = split;
	push (@file, (split '_calls', $line[0])[0]);
}
my $c=0;

#All Homozygous Deletions

my $outfile = 'HomozygousDeletions.txt';
open G, ">$outfile";

for ($c = 0; $c < @file; $c++) {
my $vcf = $file[$c] . '_calls.cns';
my $fixedvcf = $file[$c] . '_fixed.calls.cns';

$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $vcf . ' > ' . $fixedvcf;
system($cmd);

my @genes = ();

open (H, $fixedvcf);
while (<H>) {
	chomp;
  	last if m/^$/;
	@line = split;
	next if ($line[0] =~ m'chromosome');
	next if ($line[0] eq 'X');
	next if ($line[0] eq 'Y');
	next if ($line[6] > 0);
	@genes = (split ',', $line[3]);
	my $a = 0;
	for ($a = 0; $a < @genes; $a++) {
		$loc = $line[0] . ':' . $line[1] . '-' . $line[2];
		print G "$file[$c]\t$genes[$a]\t$line[6]\t$loc\n";
	}
}

close H;
}
close G;

#All Amplifications >=5

$outfile = 'HighLevelAmplifications.txt';
open G, ">$outfile";

for ($c = 0; $c < @file; $c++) {
my $vcf = $file[$c] . '_calls.cns';
my $fixedvcf = $file[$c] . '_fixed.calls.cns';

$cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $vcf . ' > ' . $fixedvcf;
system($cmd);

my @genes = ();

open (H, $fixedvcf);
while (<H>) {
	chomp;
  	last if m/^$/;
	@line = split;
	next if ($line[0] =~ m'chromosome');
	next if ($line[0] eq 'X');
	next if ($line[0] eq 'Y');
	next if ($line[6] < 5);
	@genes = (split ',', $line[3]);
	my $a = 0;
	for ($a = 0; $a < @genes; $a++) {
		$loc = $line[0] . ':' . $line[1] . '-' . $line[2];
		print G "$file[$c]\t$genes[$a]\t$line[6]\t$loc\n";
	}
}

close H;
}
close G;

#All LOH Events

$outfile = 'LOHevents.txt';
open G, ">$outfile";

for ($c = 0; $c < @file; $c++) {
my $vcf = $file[$c] . '_calls.cns';
my $fixedvcf = $file[$c] . '_fixed.calls.cns';

$cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $vcf . ' > ' . $fixedvcf;
system($cmd);

my @genes = ();

open (H, $fixedvcf);
while (<H>) {
	chomp;
  	last if m/^$/;
	@line = split;
	next if ($line[0] =~ m'chromosome');
	next if ($line[0] eq 'X');
	next if ($line[0] eq 'Y');
	next if (($line[7] eq '-') || ($line[8] eq '-'));
	next if ((($line[7] > 0) && ($line[8] > 0)) || ($line[6] == 0));
	@genes = (split ',', $line[3]);
	my $a = 0;
	for ($a = 0; $a < @genes; $a++) {
		$loc = $line[0] . ':' . $line[1] . '-' . $line[2];
		print G "$file[$c]\t$genes[$a]\t$line[6]\t$loc\n";
	}
}

close H;
}
close G;

#AllCNVevents
$outfile = 'AllNonDiploidEvents.txt';
open G, ">$outfile";

for ($c = 0; $c < @file; $c++) {
my $vcf = $file[$c] . '_calls.cns';
my $fixedvcf = $file[$c] . '_fixed.calls.cns';

$cmd = '';
$cmd = 'awk \'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i ="-" }; 1\' ' . $vcf . ' > ' . $fixedvcf;
system($cmd);

my @genes = ();

open (H, $fixedvcf);
while (<H>) {
	chomp;
  	last if m/^$/;
	@line = split;
	next if ($line[0] =~ m'chromosome');
	next if ($line[0] eq 'X');
	next if ($line[0] eq 'Y');
	next if ($line[6] == 2);
	@genes = (split ',', $line[3]);
	my $a = 0;
	for ($a = 0; $a < @genes; $a++) {
		$loc = $line[0] . ':' . $line[1] . '-' . $line[2];
		print G "$file[$c]\t$genes[$a]\t$line[6]\t$loc\n";
	}
}

close H;
}
close G;

#SummarizeGeneLevelCounts

$cmd='';

$cmd = 'sort -k2,2 -k1,1 HighLevelAmplifications.txt > SortedHighLevelAmplifications.txt';
system($cmd);
$cmd = 'sort -k2,2 -k1,1 HomozygousDeletions.txt > SortedHomozygousDeletions.txt';
system($cmd);
$cmd = 'sort -k2,2 -k1,1 AllNonDiploidEvents.txt > SortedAllNonDiploidEvents.txt';
system($cmd);
$cmd = 'sort -k2,2 -k1,1 LOHevents.txt > SortedLOHevents.txt';
system($cmd);

#SummarizeAmplificationGeneLevelCounts

my $sample = 'null';
my $gene = 'null';
my $muts = 0;
my @tumors=();
my @copynumber=();
my @locs=();

$outfile = 'AmplifiedGenes.txt'; #to fix...
open C, ">$outfile";
print C "Gene\tSampleCount\tSamples\tCopyNumber\tSegmentLocations\n";

open G, "<SortedHighLevelAmplifications.txt";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	next if ($line[1] eq m/SampleName/);
	if (($line[1] eq $gene) || ($gene eq 'null')) {
		if (($line[0] ne $sample) || ($sample eq 'null') || ($gene eq '-')) {
			$muts++;
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
		}
	}
	else {
		print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
			@tumors = ();
			@copynumber = ();
			@locs=();
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
			$muts=1;
	}
}

print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
close C;
close G;

#SummarizeDeletedGeneLevelCounts

$outfile = 'DeletedGenes.txt'; #to fix...
open C, ">$outfile";
print C "Gene\tSampleCount\tSamples\tCopyNumber\tSegmentLocations\n";

$sample = 'null';
$gene = 'null';
$muts = 0;
@tumors=();
@copynumber=();
@locs=();

open G, "<SortedHomozygousDeletions.txt";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	next if ($line[1] eq m/SampleName/);
	if (($line[1] eq $gene) || ($gene eq 'null')) {
		if (($line[0] ne $sample) || ($sample eq 'null') || ($gene eq '-')) {
			$muts++;
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
		}
	}
	else {
		print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
			@tumors = ();
			@copynumber = ();
			@locs=();
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
			$muts=1;
	}
}

print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
close C;
close G;

#SummarizeLOHGeneLevelCounts

$outfile = 'LOHGenes.txt'; #to fix...
open C, ">$outfile";
print C "Gene\tSampleCount\tSamples\tCopyNumber\tSegmentLocations\n";

$sample = 'null';
$gene = 'null';
$muts = 0;
@tumors=();
@copynumber=();
@locs=();

open G, "<SortedLOHevents.txt";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	next if ($line[1] eq m/SampleName/);
	if (($line[1] eq $gene) || ($gene eq 'null')) {
		if (($line[0] ne $sample) || ($sample eq 'null') || ($gene eq '-')) {
			$muts++;
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
		}
	}
	else {
		print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
			@tumors = ();
			@copynumber = ();
			@locs=();
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
			$muts=1;
	}
}

print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
close C;
close G;

#SummarizeAllCNVsGeneLevelCounts

$outfile = 'AllCNValteredGenes.txt'; #to fix...
open C, ">$outfile";
print C "Gene\tSampleCount\tSamples\tCopyNumber\tSegmentLocations\n";

$sample = 'null';
$gene = 'null';
$muts = 0;
@tumors=();
@copynumber=();
@locs=();

open G, "<SortedAllNonDiploidEvents.txt";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
  	next if ($line[1] eq m/SampleName/);
	if (($line[1] eq $gene) || ($gene eq 'null')) {
		if (($line[0] ne $sample) || ($sample eq 'null') || ($gene eq '-')) {
			$muts++;
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
		}
	}
	else {
		print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
			@tumors = ();
			@copynumber = ();
			@locs=();
			$gene = $line[1];
			$sample = $line[0];
			push @tumors, $line[0];
			push @copynumber, $line[2];
			push @locs, $line[3];
			$muts=1;
	}
}

print C "$gene\t$muts\t@tumors\t@copynumber\t@locs\n";
close C;
close G;

$cmd = 'rm *_fixed.calls.cns SortedLOHevents.txt SortedHomozygousDeletions.txt SortedHighLevelAmplifications.txt SortedAllNonDiploidEvents.txt AllNonDiploidEvents.txt LOHevents.txt HighLevelAmplifications.txt HomozygousDeletions.txt calls';
system($cmd);