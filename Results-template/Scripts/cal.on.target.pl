#!/bin/env perl  

# Script for calcualting exome-seq on-target rate by extracting reads that are
# mapped in the target region including reads that overhang the target BED specified feature region.
# The on_target_ratio = total_reads_on_target/total_mapped_reads_on_genome 
#
# Written by Yongmei Zhao
###############################

#use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IO::File;
use File::Basename;
use File::Spec;

$o=`module load bamtools`;
$o=`module load bedtools`;

my $BAMTOOLS = "bamtools";
#my $BAMTOOLS = "/opt/nasapps/production/bamtools/2.2.3/bin/bamtools";
my $intersectBed = "intersectBed";
#my $intersectBed = "/opt/nasapps/production/bedtools/2.17.0/bin/intersectBed";

my $USAGE = "Usage: cal_on_target_rate.pl [options] bam_file \n";
my $HELPTEXT = qq~
                Program for calculating on_target percentage for exome-seq data.
                $USAGE
                options: --target_bed file   target bed file, which is required in bed3 format

        Example:   cal_on_target_rate.pl --target_bed bed_file input_bam_file
~;
my ($help,  $target_bed);
my $result = GetOptions( 'target_bed|b=s' => \$target_bed,
			 'help|h'    => \$help,
    ) || &print_help();
my $bam = $ARGV[0];
&print_help if $help;

unless (defined $bam && defined $target_bed ){
    die "must specify bam file and --target_bed file as input\n";
}

$bam = File::Spec->rel2abs($bam);
my ($file_name, $path) = fileparse($bam);
my $out1 = "$file_name\.bam_stats";

#Calculate original BAM file statistics:
my $cmd = "$BAMTOOLS stats -in $bam -insert >$out1 2>$out1\.err";
my $ret = system ($cmd);
if ($ret){
    die "Error run cmd $cmd: $! \n";
}

#generate intersection of target
my $target_bam = $file_name . ".onTarget.bam";;
$cmd = "$intersectBed -abam $bam -b $target_bed 1>$target_bam 2>$target_bam\.err" ;
$ret = system ($cmd);
if ($ret){
    die "Error run cmd $cmd: $! \n";
}

# find the stats of the onTarget BAM file
my $out2 = $file_name .".onTarget.bam_stats";;
$cmd = "$BAMTOOLS stats -in $target_bam -insert > $out2 2>$out2\.err";
$ret = system ($cmd);
if ($ret){
    die "Error run cmd $cmd: $! \n";
}

my $target_stats_ref = &parse_bam_stats($out2);
my $original_stats_ref = &parse_bam_stats($out1);
my $on_target_rate = sprintf '%.2f',  $target_stats_ref->{mapped_reads}*100/$original_stats_ref->{mapped_reads};
print STDOUT "PERCENT_ON_TARGET:\t$file_name\t$on_target_rate\n";

exit 0;

sub parse_bam_stats ($){
    my $file = shift;
    my %info = ();
    
    open (IN, $file) or die "Cannot open file $file: $!\n";

    while (my $line =<IN>){
	chomp $line;
	my @cols = split (/:\s+/, $line);
	next unless defined $cols[0];
	if ($cols[0] =~ /Total\sreads/){
	    $cols[1] =~ s/\s//g;
	    $info{total_reads} =$cols[1] ;
	}elsif ($cols[0]=~ /Mapped\sreads/){
	    $cols[1]=~s/\)//g;
	    my @counts = split(/\s\(/, $cols[1]);
	    $info{mapped_reads} = $counts[0];
	    $info{mapped_perc} =$counts[1];
	}
    }
    close IN;
    return \%info;
}

sub print_help {
    print "$HELPTEXT\n";
    exit 0;
}
