#!/usr/local/bin/perl

($infile, $out1file) = @ARGV;

open IN, "<$infile" or die "Could not open for read :$!\n";
open OUT, ">$out1file" or die "Failed to open $out1file file\n";

while (<IN>){
  chomp;
  next if /^all/;
@array = split (/\t/, $_);
($chr,$start,$stop,$gene,$n_reads,$n_bases,$exon_width)=@array;

$exon_coord = "$gene\t$chr\t$start\t$stop";

if (!exists($exon_count{$exon_coord})){
$total_exon_length{$gene} += $exon_width;
$exon_count{$exon_coord}++;
}

$genes{$gene}++;

if ($n_reads >= 1 ) {
    $coverage_1{$gene} += $n_bases;}

if ($n_reads >= 5 ) {
    $coverage_5{$gene} += $n_bases;}

if ($n_reads >= 10 ) {
    $coverage_10{$gene} += $n_bases;}

if ($n_reads >= 20 ) {
    $coverage_20{$gene} += $n_bases;}

if ($n_reads >= 50 ) {
    $coverage_50{$gene} += $n_bases;}

if ($n_reads >= 100) {
    $coverage_100{$gene} += $n_bases;}

  }

close IN;

print OUT "gene\tlength\tpct1x\tpct5x\tpct10x\tpct20x\tpct50x\tpct100x\n";

foreach $gene (sort keys %genes){

  $pct_cov1{$gene} = sprintf("%d",($coverage_1{$gene}/$total_exon_length{$gene})*100);
  $pct_cov5{$gene} = sprintf("%d",($coverage_5{$gene}/$total_exon_length{$gene})*100);
  $pct_cov10{$gene} = sprintf("%d",($coverage_10{$gene}/$total_exon_length{$gene})*100);
  $pct_cov20{$gene} = sprintf("%d",($coverage_20{$gene}/$total_exon_length{$gene})*100);
  $pct_cov50{$gene} = sprintf("%d",($coverage_50{$gene}/$total_exon_length{$gene})*100);
  $pct_cov100{$gene} = sprintf("%d",($coverage_100{$gene}/$total_exon_length{$gene})*100);

print OUT "$gene\t$total_exon_length{$gene}\t$pct_cov1{$gene}\t$pct_cov5{$gene}\t$pct_cov10{$gene}\t$pct_cov20{$gene}\t$pct_cov50{$gene}\t$pct_cov100{$gene}\n";
}

close OUT;
