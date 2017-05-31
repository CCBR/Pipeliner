#!/bin/perl
use List::Util qw / max min /;
local $/=undef;
my $gtffile   = $ARGV[0];
#open FILE, "/fdb/GENCODE/Gencode_human/release_19/gencode.v19.annotation.gtf" or die "Couldn't open file: $!";
open FILE, $gtffile or die "Couldn't open file: $!";
$sample = <FILE>;
close FILE;

my @sample_lines = split /\n/, $sample;
my $path;

my %coords;
my %chr;
my %strand;
foreach my $line (@sample_lines) {
  @t=split/\t/,$line;
  if ($t[2]=="gene") {
  ($genename) = $line =~ m/gene_name \"([^s]+)\";/g;
  $genename =~ s/\s+//g;
  next if length $genename == 0;
  $chr{$genename}=$t[0];
  $strand{$genename}=$t[6];
  $coords{$genename}=[] unless exists $coords{$genename};
  push(@{$coords{$genename}},int($t[3]));
  push(@{$coords{$genename}},int($t[4]));
  }
}
foreach my $genename (keys %coords) { 
my $start=min @{$coords{$genename}};
my $end=max @{$coords{$genename}};
print "$chr{$genename}\t$start\t$end\t$strand{$genename}\t$genename\n";
}
