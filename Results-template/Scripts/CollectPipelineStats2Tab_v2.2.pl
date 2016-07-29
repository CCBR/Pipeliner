#!/usr/bin/perl

use strict;
use Cwd;
use File::Spec;

# Gathers all pipeline metrics from a project and outputs two tab-delimited report files, ready for conversion to excel format.


# ChangeLog:
#    5.11.2015: -added options for print (demultiplex stats and/or summary stats) and pipeline (ChIPseq or RNAseq) mode.
#    5.13.2015: -changed output file format to include project name, cellflow name, report type, and an optional user designated suffix.
#               -added output field '% Uniquely aligned' to the chipseq pipeline summary report.
#    5.15.2015: -new/improved option parsing routine.  Doesn't break with hyphenated file names.
#               -print mode 3 (analysis stats) is now independent (i.e. you can print analysis stats without a corresponding demultiplex directory).
#                However, yield and Q30 will be unavailable in this mode.
#    5/19/2015: -added option to collect Hiseq RNAseq statistics.
#    5/26/2015: -added paired-end/single-end option.  However, default behavior of paired-end RNAseq and
#                hiseq ChIPseq, and single-end nextseq ChIPseq are the only currently tested modes.
#    6/18/2015: -bebugged/added ability to correctly parse hiseq samples on multiple lanes
#    7/17/2015: -mode 3 added to parse exome-seq pipeline stats.
#    9/18/2015: -fixed bug that messed up output when there was no sample barcode
#    10/2/2015: -v1.1: added gDNA pipeline (-r 4) for BWA.  Uses Exomeseq pipeline without on-target stats (i.e. run_bwa_pipeline.pl)
#               -fixed a major html parsing bug related to empty fields in the sample sheet causing cascading hash errors
#    2/1/2016:  -v2.0: code re-write for bcl2fastq 2.17 and up.  Hiseq, Miseq, and Nextseq demultiplex
#                formats are now the same.
#    4/15/16:   -2.1: added sample summary output to demultiplex only reports.  This has added benifit of not breaking Tab2Excel.pl
#    5/18/16:   -added additional sammary output field 'Median Insert Size' to gDNA and Exome pipelines.
#               -added additional error detection.
#    6/14/2016: -added Star RNAseq pipeline option
#               -added back trim reads, total mapped reads, and uniquely mapped read stats (they were removed in last edit)

#known issues:
#               -print mode 3 (i.e. Summary output only) is broken as of v2.0.


# defaults, initializations, and constants
my $help = "\nCollectPipelineStats2Tab v2.1\nDescription:\n  Gathers all pipeline metrics from a project and outputs two tab-delimited report files.\n".
          "\te.g. from ChipSeq run directory:\n\t :~\$ CollectPipelineStats2Tab.pl -p Hager_Kim_18289_10Chip\n\t     -u Unaligned/Hager_Kim__18289_10_Chip_11_30_15\n\t     -d analysis/project_Hager_Kim__18289_10_Chip/ -r 1".
          "\nOptions/Switches:\n".
          "\t-o  Option: specify output file name suffix.  Optional,\n\t    format='[ProjectName]_[FlowcellName]_{Suffix}_[ReportID].txt'.\n".
          "\t-u  Option: specify unaligned read directory.  Optional,\n\t    default='Unaligned/[Project Name]/.\n".
          "\t-p  Option: specify project name.  Required.\n".
          "\t-f  Option: specify flowcell.  Optional, default=derived from\n\t    run directory.\n".
          "\t-d  Option: specify project analysis directory path.  Optional,\n\t    default='analysis/'.\n".
          "\t-m  Option: output mode:\n\t\t1: collect all stats.\n\t\t2: collect demultiplex stats only.\n\t\t3: collect pipeline stats only.\n\t    Optional, default is m=1.\n".
          "\t-r  Option: pipeline mode:\n\t\t1: collect ChIPseq stats.\n\t\t2: collect Tophat RNAseq stats.\n\t\t3: collect ExomeSeq stats.\n\t\t4: collect gDNAseq stats.\n\t\t5: collect Star RNAseq stats.\n\t    Optional, default is r=1\n".
          "\t-e  Option: read type:\n\t\t1: single-end\n\t\t2: paired-end\n\t    Optional, default behavior is paired-end for all RNAseq\n\t    and Hiseq ChIPseq, single-end for Nextseq ChIPseq.\n".
          "\n************\nAuthor: J. Cristobal Vera, email: cris.vera\@nih.gov\n";
my $usage = "\nUsage:\nCollectPipelineStats2Tab.pl -o [Output File Name Suffix] -p [Project Name] -u [Unaligned Path] -d [Project Path] -f [Flowcell] -m [Output Mode] -r [Pipeline Mode] -e [Read Type]\n";
my $dirbase = getcwd;
my $cellname = $dirbase;
my ($projectname,$outfile,$demultend);
my $projectpath = './';  #default path to project
my $demultpath = 'laneBarcode.html'; ####default partial path to nextseq demultiplex html report.
my $printmode = 1;     #output mode default=both
my $pipeline = 3; #pipeline mode default=ExomeSeq
my $paired = 2; #single or paired-end, see below, default = single-end
my $z = 0;
my $x = 0;
my $y = 0;

#process command line custom script options/switches
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
$outfile = $cmds{'o'} if ($cmds{'o'});
$outfile = '_'.$outfile if ($outfile);
$projectname = $cmds{'p'} if ($cmds{'p'});
$projectpath = $cmds{'d'} if ($cmds{'d'});
$demultpath = $cmds{'u'} if ($cmds{'u'});
if ($demultpath =~ s/\/([^\/]+)\/*$//) {
  $demultend = $1;
}
else{
  $demultend = $projectname;
}
$pipeline = $cmds{'r'} if ($cmds{'r'});
$printmode = $cmds{'m'} if ($cmds{'m'});
$paired = 2 if ($pipeline == 2 or $pipeline == 3 or $pipeline == 4); #paired end if any RNAseq or if hiseq chipseq
$paired = $cmds{'e'} if ($cmds{'e'});
$paired = 0 if ($paired == 1);
#get flowcell name from working dir if available
if ($cellname =~ m/^.+_[AB]/) {
  $cellname =~ s/^.+_[AB]//;
  $cellname =~ s/\/$//;
}
elsif ($cellname =~ m/_0+-(A[A-Z0-9]+)\/*$/){
  $cellname = $1;
}
$cellname = $cmds{'f'} if ($cmds{'f'});

print STDERR "\nRead Type:  Paired-end\n" if ($paired);
print STDERR "\nRead Type:  Single-end\n" if (!$paired);

die "\nERROR: must specify the project name.\n" if (!$projectname);

#get absolute paths
$demultpath = File::Spec->rel2abs($demultpath);
$projectpath = File::Spec->rel2abs($projectpath);

#die "\nError: demultiplex directory not found: $demultpath\n\n" if (-d $demultpath);
#die "\nError: project directory not found: $projectpath\n\n" if (-d $projectpath);

#gather html project demultiplex stats
my $start = 0;
my $fieldcount = 0;
my $samplecount = 0;
my $tReads1 = my $tReads2 = my $tYield = 0;
my $cReads = my $cYield = my $cQ30 = 0;
my ($sample,$samplename,$lane,$fields,$pf);
#demultiplex titles
my @fields = ('Lane','Sample','Barcode Sequence','PF Clusters','Percent of Lane','Percent Perfect Barcode',
              'Percent One Mismatch Barcode','Yield (Mbases)','Percent PF Clusters','Percent Bases >= Q30','Mean Quality Score',
              );
my %demultstats;
if ($printmode == 1 or $printmode == 2) {
  $demultpath = 'laneBarcode.html';
  open(DMT,"<$demultpath") or die "Cannot open $demultpath: $!\n";
  while (my $line = <DMT>) {
    chomp $line;
    next if ($line eq '');
    if ($line =~ m/^\<h2\>Flowcell Summary/) {
      $start = 1;
    }
    elsif ($line =~ m/^\<h2\>Lane Summary/){
      $start = 2;
      $fieldcount = 0;
    }
    elsif ($start == 1 and $line =~ m/^\<td\>(.+)\<\/td\>/){
      $fieldcount += 1;
      $tReads1 = $1 if ($fieldcount == 1);
      $tReads2 = $1 if ($fieldcount == 2);
      $tYield = $1 if ($fieldcount == 3);
    }
    elsif ($start == 2 and $line =~ m/^\<\/tr\>/){
      $fieldcount = 0;
    }
    elsif ($start == 2 and $line =~ m/^\<td.*\>(.+)\<\/td\>/){
      $fieldcount += 1;
      my $tmp = $1;
      $tmp =~ s/,//g;
      $lane = $tmp if ($fieldcount == 1);
      $sample = $tmp if ($fieldcount == 2);
      $demultstats{$sample}{'Lane'} .= "\t$lane" if ($fieldcount == 2 and exists $demultstats{$sample}{'Lane'});
      $demultstats{$sample}{'Lane'} = $lane if ($fieldcount == 2 and !exists $demultstats{$sample}{'Lane'});
      $demultstats{$sample}{'Barcode Sequence'} .= "\t$tmp" if ($fieldcount == 3 and exists $demultstats{$sample}{'Barcode Sequence'});
      $demultstats{$sample}{'Barcode Sequence'} = $tmp if ($fieldcount == 3 and !exists $demultstats{$sample}{'Barcode Sequence'});
      $demultstats{$sample}{'PF Clusters'} .= "\t$tmp" if ($fieldcount == 4 and exists $demultstats{$sample}{'PF Clusters'});
      $demultstats{$sample}{'PF Clusters'} = $tmp if ($fieldcount == 4 and !exists $demultstats{$sample}{'PF Clusters'});
      $cReads += $tmp if ($fieldcount == 4);
      $demultstats{$sample}{'Percent of Lane'} .= "\t$tmp" if ($fieldcount == 5 and exists $demultstats{$sample}{'Percent of Lane'});
      $demultstats{$sample}{'Percent of Lane'} = $tmp if ($fieldcount == 5 and !exists $demultstats{$sample}{'Percent of Lane'});
      $demultstats{$sample}{'Percent Perfect Barcode'} .= "\t$tmp" if ($fieldcount == 6 and exists $demultstats{$sample}{'Percent Perfect Barcode'});
      $demultstats{$sample}{'Percent Perfect Barcode'} = $tmp if ($fieldcount == 6 and !exists $demultstats{$sample}{'Percent Perfect Barcode'});
      $demultstats{$sample}{'Percent One Mismatch Barcode'} .= "\t$tmp" if ($fieldcount == 7 and exists $demultstats{$sample}{'Percent One Mismatch Barcode'});
      $demultstats{$sample}{'Percent One Mismatch Barcode'} = $tmp if ($fieldcount == 7 and !exists $demultstats{$sample}{'Percent One Mismatch Barcode'});
      $demultstats{$sample}{'Yield (Mbases)'} .= "\t$tmp" if ($fieldcount == 8 and exists $demultstats{$sample}{'Yield (Mbases)'});
      $demultstats{$sample}{'Yield (Mbases)'} = $tmp if ($fieldcount == 8 and !exists $demultstats{$sample}{'Yield (Mbases)'});
      $cYield += $tmp if ($fieldcount == 8);
      $demultstats{$sample}{'Percent PF Clusters'} .= "\t$tmp" if ($fieldcount == 9 and exists $demultstats{$sample}{'Percent PF Clusters'});
      $demultstats{$sample}{'Percent PF Clusters'} = $tmp if ($fieldcount == 9 and !exists $demultstats{$sample}{'Percent PF Clusters'});
      $demultstats{$sample}{'Percent Bases >= Q30'} .= "\t$tmp" if ($fieldcount == 10 and exists $demultstats{$sample}{'Percent Bases >= Q30'});
      $demultstats{$sample}{'Percent Bases >= Q30'} = $tmp if ($fieldcount == 10 and !exists $demultstats{$sample}{'Percent Bases >= Q30'});
      $cQ30 += $tmp if ($fieldcount == 10);
      $demultstats{$sample}{'Mean Quality Score'} .= "\t$tmp" if ($fieldcount == 11 and exists $demultstats{$sample}{'Mean Quality Score'});
      $demultstats{$sample}{'Mean Quality Score'} = $tmp if ($fieldcount == 11 and !exists $demultstats{$sample}{'Mean Quality Score'});
    }
  }
  close(DMT);
}
elsif ($printmode == 3){
  opendir(DIR,$demultpath) or die "Can't open directory: $projectpath: $!\n";
  my @dirs = grep {-d "$projectpath/$_"} readdir(DIR);
  close(DIR);
  foreach my $dir (@dirs){
    if ($dir eq 'QC' or $dir =~ m/^\./) {}    
    elsif ($dir) {
        $demultstats{$dir} = 1;
    }
  }
}

#print  multiplex and summary tab-delimited text file headers
open (OUT1, ">$projectname\_$cellname$outfile\_Demultiplex.txt") or die "Cannot create $projectname\_$cellname$outfile\_Demultiplex.txt: $!\n" if ($printmode == 1 or $printmode == 2);
open (OUT2, ">$projectname\_$cellname$outfile\_Summary.txt") or die "Cannot create $projectname\_$cellname$outfile\_Summary.txt: $!\n";
open (OUT3, ">$projectname\_$cellname$outfile\_RnaStats.txt") or die "Cannot create $projectname\_$cellname$outfile\_RnaStats.txt: $!\n" if (($pipeline == 2 or $pipeline == 5) and ($printmode == 1 or $printmode == 3));

#print multiplex titles and totals
if ($printmode == 1 or $printmode == 2){
  print OUT1 "Flowcell Summary: $cellname\n\n";
  print OUT1 "Cluster (Raw)\tClusters (PF)\tYield (MBases)\n";
  print OUT1 "$tReads1\t$tReads2\t$tYield\n\n";
  print OUT1 "Lane Summary\n\n";
  print OUT1 "Lane\tSample\tBarcode Sequence\tPF Clusters\tPercent of Lane\tPercent Perfect Barcode\tPercent One Mismatch Barcode\tYield (Mbases)\tPercent PF Clusters\tPercent Bases >= Q30\tMean Quality Score\n";
}

#print summary and RNA stats titles
if ($printmode == 1 or $printmode == 3) {
  print OUT2 "Sample ID\tSample Yield (Mbases)\tPercent of (PF) Bases >= Q30\tTotal Reads (PF)\tTotal Reads After Trimming\tPercent Total Reads after Trimming\tTotal Mapped Reads (Trimmed)\tPercent Total Mapped Reads (Trimmed)\tUniquely Mapped Reads (Trimmed)\tPercent Uniquely Mapped Reads (Trimmed)\tPercent Non-duplicated Reads (Trimmed)" if ($pipeline == 1 or $pipeline == 2 or $pipeline == 5);
  print OUT3 "Sample ID\tPF_BASES\tPF_ALIGNED_BASES\tRIBOSOMAL_BASES\tCODING_BASES\tUTR_BASES\tINTRONIC_BASES\tINTERGENIC_BASES\tCORRECT_STRAND_READS\tINCORRECT_STRAND_READS\tPCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tPCT_MRNA_BASES\tPCT_USABLE_BASES\tPCT_CORRECT_STRAND_READS\tMEDIAN_CV_COVERAGE\tMEDIAN_5PRIME_BIAS\tMEDIAN_3PRIME_BIAS\tMEDIAN_5PRIME_TO_3PRIME_BIAS\n" if ($pipeline == 2 or $pipeline == 5);
  print OUT2 "\tPCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tPCT_MRNA_BASES\tPCT_CORRECT_STRAND_READS\tMEDIAN_5PRIME_TO_3PRIME_BIAS" if ($pipeline == 2 or $pipeline == 5);
  print OUT2 "Sample ID\tSample Yield (Mbases)\tPercent of (PF) Bases >= Q30\tTotal Reads (PF)\tTotal Reads After Trimming\tPercent Total Reads after Trimming\tTotal Mapped Reads (Trimmed)\tPercent Total Mapped Reads (Trimmed)\tPercent Non-duplicated Reads (Trimmed)\tPercent Reads Mapped On Target\tMean Coverage Inside Target (X)\tPercent Coverage >= 30x\tPercent GC\tMedian Insert Size" if ($pipeline == 3);
  print OUT2 "Sample ID\tSample Yield (Mbases)\tPercent of (PF) Bases >= Q30\tTotal Reads (PF)\tTotal Reads After Trimming\tPercent Total Reads after Trimming\tTotal Mapped Reads (Trimmed)\tPercent Total Mapped Reads (Trimmed)\tPercent Non-duplicated Reads (Trimmed)\tMean Coverage Depth\tMean Coverage Depth SD\tPercent of Coverage >= 5X\tPercent GC\tMedian Insert Size" if ($pipeline == 4);
  print OUT2 "\n";
}
print OUT2 "Sample ID\tSample Yield (Mbases)\tPercent of (PF) Bases >= Q30\tTotal Reads (PF)\n" if ($printmode == 2);

#get summary stats and print all
my ($trimreads,$percenttrimmed,$covSD,$mappedreads,$meancov,$percentcov30,$targetreads,$gc,$percenttarget,$percentdups,$rnastats,
    $minrnastats,$percentnondups,$percentmapped,$mult,$unimappedreads,$percentunimapped,$insert);
foreach my $sample (sort {lc($a) cmp lc($b)} keys %demultstats){
  $z += 1;
  #finish printing to demutliplex file
  my $subYield = my $subQ30 = my $subReads = 0;
  if ($printmode == 1 or $printmode == 2) {
    my $n = scalar(split /\t/,$demultstats{$sample}{'Lane'});
    for (my $i = 1; $i <= $n; $i++){
      foreach my $field (@fields){
        my @string = split /\t/,$demultstats{$sample}{$field};
        print OUT1 $string[($i-1)] if ($field eq 'Lane');
        print OUT1 "\t$sample" if ($field eq 'Sample');
        print OUT1 "\t$string[($i-1)]" if ($field ne 'Sample' and $field ne 'Lane');
        $subYield += $string[($i-1)] if ($field eq 'Yield (Mbases)');
        $subQ30 += $string[($i-1)] if ($field =~ m/Q30/);
        $subReads += $string[($i-1)] if ($field eq 'PF Clusters');
      }
      print OUT1 "\n";
    }
    $subQ30 = Round1000th($subQ30/$n);
    #print sample summaries
    print OUT2 "$sample\t$subYield\t$subQ30\t$subReads\n" if ($printmode == 2);
  }
  $subYield = $subQ30 = $subReads = 0 if ($printmode == 3);
  
  #find alignment stats in analysis sample folders
  my $samplefile = $sample; 
  $samplefile = 'Sample_'.$samplefile if (($pipeline == 1 or $pipeline == 2) and $samplefile !~ m/^Sample_/); #should be adding 'Sample_'. to all sample names for future
  my $markdfile = $samplefile.'_MARKEDUPmetrics.txt'; # for ChIPseq and Star RNAseq duplication stats
  $markdfile = $samplefile.'_MKDUPmetrics.txt' if ($pipeline == 2);
  $markdfile = $samplefile.'.sorted.txt' if ($pipeline == 3 or $pipeline == 4);
  my $bamfile = $samplefile.'_bowtie2.err';  # for ChIPseq stats
  my $trimfile = 'QC/'.$samplefile.'_run_trimmomatic.err';  # for RNAseq and exome trimming stats
  my $picardfile = $samplefile.'_RnaSeqMetrics.txt'; # for RNAseq RNA stats
  my $tophatfile = 'align_summary.txt';  # for tophat RNSseq alignment stats
  my $starfile = $samplefile.'_Log.final.out'; #for star RNAseq alignment stats
  my $exomefile1 = $samplefile.'.dedup.bam.bam_stats';  ###exome alignment and markdup stats
  my $exomefile2 = $samplefile.'.dedup.bam.onTarget.bam_stats'; ###exome on-target stats
  my $exomecovfile = 'QC/'.$samplefile.'.qualimapReport/genome_results.txt';
  my $samplefolder = "$projectpath";
  my $sampleexfolder = "$projectpath";
  if ($printmode == 1 or $printmode == 3) {
    if ($pipeline == 1) {
      open (BAM, "<$samplefolder".$bamfile) or die "Cannot open $samplefolder$bamfile:$!\n";
      chomp (my @bams = <BAM>);
      close BAM;
      if ($printmode == 3) {  ###get read trimming stats from trimmomatic if print mode 3
        open (TRM, "<$samplefolder".$trimfile) or die "Cannot open $samplefolder$trimfile:$!\n";
        chomp(my @trims = <TRM>);
        close TRM;
        foreach my $trim (@trims){
          if ($trim =~ m/^Input Reads: +(\d+)/) {
            $subReads = $1;
            #$paired = 0;  ##maybe make this obligatory?
          }
          elsif($trim =~ m/^Input Read Pairs: +(\d+)/){
            $subReads = $1; #multiply by two to get total reads when hiseq chipseq
            #$paired = 1;  ##maybe make this obligatory?
          }
        }
      }
      elsif ($paired){
        $subReads = $subReads*2;
      }
      my ($mp2,$nomap);
      $unimappedreads = 0;
      foreach my $bam (@bams){
        if ($bam =~ m/^([0-9]+) reads; of these:/) {
          $trimreads = $1;
          $trimreads *= 2 if ($paired);
        }
        elsif (!$paired and $bam =~ m/^ +([0-9]+).+ aligned exactly 1 time/){
          $unimappedreads = $1;
        }
        elsif (!$paired and $bam =~ m/ +([0-9]+).+ aligned \>1 times/){
          $mp2 = $1;
        }
        elsif ($paired and $bam =~ m/ +([0-9]+).+ aligned concordantly exactly 1 time/){
          $unimappedreads += $1*2;
        }
        elsif ($paired and $bam =~ m/ +([0-9]+).+ aligned discordantly 1 time/){
          $unimappedreads += $1*2;
        }
        elsif ($paired and $bam =~ m/ +([0-9]+).+ aligned exactly 1 time/){
          $unimappedreads += $1;
        }
        elsif ($paired and $bam =~ m/^ +([0-9]+) .+ aligned 0 times$/){
          $nomap = $1;  ##$nomap is reads, not read pairs
        }
      }
      $mappedreads = $unimappedreads + $mp2 if (!$paired);
      $mappedreads = $trimreads - $nomap if ($paired);
    }
    elsif ($pipeline == 2){
      open (TRM, "<$samplefolder".$trimfile) or die "Cannot open $samplefolder$trimfile:$!\n";
      chomp(my @trims = <TRM>);
      close TRM;
      open (TOP, "<$samplefolder".'tophat_alignment/'.$tophatfile) or die "Cannot open $samplefolder"."tophat_alignment/$tophatfile:$!\n";
      chomp(my @tops = <TOP>);
      close TOP;
      open (PIC, "<$samplefolder".'tophat_alignment/'.$picardfile) or die "Cannot open $samplefolder"."tophat_alignment/$picardfile:$!\n";
      chomp(my @pics = <PIC>);
      close PIC;
      #RNAseq is always assumed to be paired-end, must add single-end if this becomes untrue
      foreach my $trim (@trims){
        if ($trim =~ m/ Both Surviving: +(\d+)/) {
          $trimreads = $1;
          $trimreads *= 2 if ($paired); #convert paired-end to total reads
        }
        if ($trim =~ m/^Input Read.{0,7}: +(\d+)/) {
            $subReads = $1;
            $subReads *= 2 if ($paired);
          }
      }
      my $next = 0;
      $mappedreads = 0;
      $mult = 0;
      foreach my $top (@tops){
        if ($next < 2 and $top =~ m/^ *Mapped *: +(\d+)/) {
          $mappedreads += $1;
        }
        elsif ($next < 2 and $top =~ m/^ *of these: +(\d+)/){
          $mult += $1;
          $next += 1;
        }
      }
      $next = 0;
      foreach my $pic (@pics){
        if ($pic =~ m/^PF_BASES/) {
          $next = 1;
        }
        elsif ($next){
          my @rnastats = my @minrnastats = split /\t/,$pic;
          splice @rnastats,7,1;
          splice @minrnastats,0,10;
          splice @minrnastats,6,1;
          splice @minrnastats,7,3;
          $rnastats = join "\t",@rnastats;
          $minrnastats = join "\t",@minrnastats;
          $next = 0;
        }
      }
    }
    elsif ($pipeline == 5){
      #open files, get lines
      open (TRM, "<$samplefolder".'logs/'.$trimfile) or die "Cannot open $samplefolder"."logs/$trimfile:$!\n";
      chomp(my @trims = <TRM>);
      close TRM;
      open (STAR, "<$samplefolder".$starfile) or die "Cannot open $samplefolder$starfile:$!\n";
      chomp(my @stars = <STAR>);
      close STAR;
      open (PIC, "<$samplefolder".$picardfile) or die "Cannot open $samplefolder$picardfile:$!\n";
      chomp(my @pics = <PIC>);
      close PIC;
      #parse trim file lines
      foreach my $trim (@trims){
        if ($trim =~ m/ Both Surviving: +(\d+)/) {
          $trimreads = $1;
          $trimreads *= 2 if ($paired); #convert paired-end to total reads
        }
        if ($trim =~ m/^Input Read.{0,7}: +(\d+)/) {
          $subReads = $1;
          $subReads *= 2 if ($paired);
        }
      }
      $unimappedreads = 0;
      $mult = $mappedreads = 0;
      foreach my $star (@stars){
        if ($star =~ m/^ *Uniquely mapped reads number \|\t([0-9]+)/) {
          $unimappedreads = $1;
          $unimappedreads *= 2 if ($paired);
        }
        elsif ($star =~ m/^ *Number of reads mapped to multiple loci \|\t([0-9]+)/){
          $mult = $1;
          $mult *= 2 if ($paired);
        }
      }
      $mappedreads = $unimappedreads + $mult;
      my $next = 0;
      foreach my $pic (@pics){
        if ($pic =~ m/^PF_BASES/) {
          $next = 1;
        }
        elsif ($next){
          my @rnastats = my @minrnastats = split /\t/,$pic;
          splice @rnastats,7,1;
          splice @minrnastats,0,10;
          splice @minrnastats,6,1;
          splice @minrnastats,7,3;
          $rnastats = join "\t",@rnastats;
          $minrnastats = join "\t",@minrnastats;
          $next = 0;
        }
      }
    }
    elsif ($pipeline == 3 or $pipeline == 4){
      open (EX1, "<".$exomefile1) or die "Cannot open $exomefile1: $!\n";
      chomp(my @exomes = <EX1>);
      close EX1;
      my @targets;
      if ($pipeline == 3) {
        open (EX2, "<".$exomefile2) or die "Cannot open $exomefile2: $!\n";
        chomp(@targets = <EX2>);
        close EX2;
      }
      open (TRM, "<".$trimfile) or die "Cannot open $projectpath/$trimfile: $!\n";
      chomp(my @trims = <TRM>);
      close TRM;
      open (COV, "<$sampleexfolder/".$exomecovfile) or die "Cannot open $sampleexfolder/$exomecovfile: $!\n";
      chomp(my @coverages = <COV>);
      close COV;
      #Exomeseq and gDNA 
      foreach my $trim (@trims){
        if ($trim =~ m/ Both Surviving: +(\d+)/) {
          $trimreads = $1;
          $trimreads *= 2 if ($paired); #convert paired-end to total reads
        }
        if ($trim =~ m/^Input Read.{0,7}: +(\d+)/) {
          $subReads = $1;
          $subReads *= 2 if ($paired);
        }
      }
      foreach my $exome (@exomes){
        if ($exome =~ m/Mapped reads: +([0-9]+)/) {
          $mappedreads = $1;
        }
      }
      if ($pipeline == 3) {
        foreach my $target (@targets){
          if ($target =~ m/Mapped reads: +([0-9]+)/) {
            $targetreads = $1;
          }
        }
      }
      foreach my $coverage (@coverages){
        if ($coverage =~ m/mean coverageData = ([0-9]+\.[0-9]+)X/) {
          $meancov = $1;
        }
        if ($pipeline == 3 and $coverage =~ m/There is a ([0-9]+\.[0-9]+)%.*\>\= 30X/) {
          $percentcov30 = $1;
        }
        if ($pipeline == 4 and $coverage =~ m/There is a ([0-9]+\.[0-9]+)%.*\>\= 5X/) {
          $percentcov30 = $1;
        }
        if ($pipeline == 4 and $coverage =~ m/std coverageData = ([0-9]+\.[0-9]+)X/) {
          $covSD = $1;
        }
        if ($coverage =~ m/GC percentage \= ([0-9]+\.[0-9]+)\%/) {
          $gc = $1;
        }
        if ($coverage =~ m/median insert size \= ([0-9\.]+)/) {
          $insert = $1;
        }
      }
    }
    
    #get markdup stats
    open (MRK, "<".$markdfile) or die "Cannot open $samplefolder$markdfile: $!\n" if ($pipeline == 1 or $pipeline == 2 or $pipeline == 5);
    open (MRK, "<".$markdfile) or die "Cannot open $projectpath$markdfile: $!\n" if ($pipeline == 3 or $pipeline == 4);
    chomp (my @marks = <MRK>);
    close MRK;
    my $mtitle = 0;
    foreach my $mark (@marks){
      my @line = split /\t/,$mark;
      if ($mtitle) {
        $percentdups = $line[$mtitle];
        $mtitle = 0;
      }
      elsif ($mark =~ m/\t*PERCENT_DUPLICATION\t*/){
        my $y = 0;
        foreach my $field (@line){
          $mtitle = $y if ($field eq 'PERCENT_DUPLICATION');
          $y += 1;
        }
      }
    }
    $percenttrimmed = Round1000th(($trimreads/$subReads)*100);
    $percentnondups = Round1000th((1 - $percentdups)*100);
    $percentmapped = Round1000th(($mappedreads/$trimreads)*100);
    $unimappedreads = ($mappedreads - $mult) if ($pipeline == 2);
    $percentunimapped = Round1000th(($unimappedreads/$trimreads)*100);
    $percenttarget = Round1000th(($targetreads/$mappedreads)*100);
  }
  
  #print to summary file
  $subYield = 'n/a' if ($printmode == 3);
  $subQ30 = 'n/a' if ($printmode == 3);
  print OUT2 "$sample\t$subYield\t$subQ30\t$subReads\t$trimreads\t$percenttrimmed\t$mappedreads\t$percentmapped\t$unimappedreads\t$percentunimapped\t$percentnondups\n" if ($pipeline == 1 and ($printmode == 1 or $printmode == 3));
  print OUT2 "$sample\t$subYield\t$subQ30\t$subReads\t$trimreads\t$percenttrimmed\t$mappedreads\t$percentmapped\t$unimappedreads\t$percentunimapped\t$percentnondups\t$minrnastats\n" if (($pipeline == 2 or $pipeline == 5) and ($printmode == 1 or $printmode == 3));
  print OUT2 "$sample\t$subYield\t$subQ30\t$subReads\t$trimreads\t$percenttrimmed\t$mappedreads\t$percentmapped\t$percentnondups\t$percenttarget\t$meancov\t$percentcov30\t$gc\t$insert\n" if ($pipeline == 3 and ($printmode == 1 or $printmode == 3));
  print OUT2 "$sample\t$subYield\t$subQ30\t$subReads\t$trimreads\t$percenttrimmed\t$mappedreads\t$percentmapped\t$percentnondups\t$meancov\t$covSD\t$percentcov30\t$gc\t$insert\n" if ($pipeline == 4 and ($printmode == 1 or $printmode == 3));
  print OUT3 "$sample\t$rnastats\n" if (($pipeline == 2 or $pipeline == 5) and ($printmode == 1 or $printmode == 3))
}

print STDERR "\nTotal samples parsed $z\n";

#round to the nearest 1000th decimal place
sub Round1000th{
  my ($number) = @_;
  if ($number =~ s/^([0-9]+\.[0-9][0-9][0-9])([0-9])[0-9]*/$1/) {
    if ($2 >= 5) {
      $number += 0.001;
    }
    return $number;
  }
  else{
    return $number;
  }
}

sub ReturnCmds{
  my (@cmds) = @_;
  my ($opt);
  my %cmds;
  foreach my $cmd (@cmds){
    if (!$opt and $cmd =~ m/^-([a-zA-Z])/) {
      $opt = $1;
    }
    elsif ($opt and $cmd =~ m/^-([a-zA-Z])/){
      $cmds{$opt} = 1;
      $opt = $1;
    }
    elsif ($opt){
      $cmds{$opt} = $cmd;
      $opt = '';
    }
  }
  $cmds{$opt} = 1 if ($opt);
  return %cmds;
}