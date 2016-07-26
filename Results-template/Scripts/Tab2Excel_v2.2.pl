#!/usr/bin/perl

use lib "/is2/projects/CCR-SF/scratch/illumina/Processing/ANALYSIS/DATA/verajc/perl/lib/share/perl/5.14.2";

use Cwd;
use Excel::Writer::XLSX;

#version 2 is for bcl2fastq 2.17 and up
#version 2.2: added star RNAseq pipeline and additional fields

# defaults, initializations, and constants
my $help = "Tab2Excel v2.1\nDescription:\nConverts CollectPipelineStats2Excel output files (tab-delimied text files) to".
            "excel file (.xlsx) and creates the apporapriate figures.\n".
            "\nOptions/Switches:\n".
            "\t-i  Option: specify input project/file name ([Project]_[Flowcell]).\n\t    Required.\n".
            "\t-o  Option: specify output file name. Optional,\n\t    default is to use input name for output.\n".
            "\t-r  Option: pipeline mode:\n\t\t1: collect ChIPseq stats.\n\t\t2: collect Tophat RNAseq stats.\n\t\t3: collect Exome-seq stats.\n\t\t4: collect BWA gDNA stats.\n\t\t5: collect STAR RNAseq stats.\n\t    Optional, default is r=1\n".
            "\n************\nAuthor: J. Cristobal Vera, email: cris.vera\@nih.gov\n";      
my $usage = "\nUsage:\nTab2Excel.pl -i [Input Project/File Name] -o [Output File Name] -r [Pipeline Mode]\n";
my $dirbase = getcwd;
my ($infile1,$infile2,$infile3,$outfile);
my (@dems,@sums,@rna);
my $pipeline = 1; #default to Chipseq

#process command line custom script options/switches
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
$pipeline = $cmds{'r'} if ($cmds{'r'});
if ($cmds{'i'}) {
  $cmds{'i'} =~ s/_Demultiplex//;
  $cmds{'i'} =~ s/_Summary//;
  $cmds{'i'} =~ s/_RnaStats//;
  $infile1 = $infile2 = $infile3 = $cmds{'i'};
  if ($cmds{'i'} =~ m/\./) {
    $infile1 =~ s/\.(.+)$/_Demultiplex.$1/;
    $infile2 =~ s/\.(.+)$/_Summary.$1/;
    $infile3 =~ s/\.(.+)$/_RnaStats.$1/;
  }
  else{
    $infile1 .= "_Demultiplex.txt";
    $infile2 .= "_Summary.txt";
    $infile3 .= "_RnaStats.txt";
  }
  open ('IN', "<$infile1") or die "Cannot open $infile1: $!\n";
  chomp(@dems = <IN>);
  close('IN');
  open ('IN2', "<$infile2") or die "Cannot open $infile2: $!\n";
  chomp(@sums = <IN2>);
  close('IN2');
  if ($pipeline == 2 or $pipeline == 5) {
    open ('IN3', "<$infile3") or die "Cannot open $infile3: $!\n";
    chomp(@rna = <IN3>);
    close('IN3');
  }
}
else{
    die "ERROR: must specify input file/project name.\n"
}
if ($cmds{'o'}) {
  $outfile = $cmds{'o'};
}
else{
    $outfile = $cmds{'i'};
}
$outfile =~ s/\..+$//;




#create workbook and two or three sheets
my $report = Excel::Writer::XLSX->new("$outfile.xlsx");
my $demulti = $report->add_worksheet('Demultiplex');
my $summary = $report->add_worksheet('Summary');
my $rna = $report->add_worksheet('RNA Statistics') if ($pipeline == 2 or $pipeline == 5);

#define title formats
my $titleformat = $report->add_format( bold => 1, size => 18, font => 'Calibri' );
my $subtitleformat = $report->add_format( bold => 1, size => 12, font => 'Calibri', border => 1, align => 'center', valign => 'vcenter', text_wrap => 1 );
#my $subtitleformat2 = $report->add_format( bold => 1, size => 12, font => 'Calibri', border => 1, align => 'center' );
#row formats
my $charformat = $report->add_format( bold => 0, size => 11, font => 'Calibri', border => 1 );
my $numformat = $report->add_format( bold => 0, size => 11, font => 'Calibri', border => 1, num_format => 3);
my $performat = $report->add_format( bold => 0, size => 11, font => 'Calibri', border => 1, num_format => 2 );
my $decformat = $report->add_format( bold => 0, size => 11, font => 'Calibri', border => 1, num_format => '0.0000' );
my $charformatb = $report->add_format( bold => 1, size => 11, font => 'Calibri' );
my $numformatb = $report->add_format( bold => 1, size => 11, font => 'Calibri', num_format => 3 );
my $performatb = $report->add_format( bold => 1, size => 11, font => 'Calibri', num_format => 2 );

#write demliplex worksheet
my $r = 0;
my $start = 0;
foreach my $row (@dems){
  $r += 1;
  my @row = split /\t/,$row;
  my $c = 0;
  foreach my $cell (@row){
    $c += 1;
    if ($r == 1 or $r == 6){
      $demulti->write(($r-1),($c-1),$cell,$titleformat);
    }
    elsif ($r == 3 or $r == 8){
      $demulti->write(($r-1),($c-1),$cell,$subtitleformat);
      $start = 1;
    }
    elsif ($r == 4){
      $demulti->write(($r-1),($c-1),$cell,$numformat);
    }
    elsif ($start){
      $demulti->write(($r-1),($c-1),$cell,$numformat) if ($c == 1 or $c == 4 or $c == 8);
      $demulti->write(($r-1),($c-1),$cell,$performat) if (($c >= 5 and $c <= 7) or ($c >= 9 and $c <= 11));
      $demulti->write(($r-1),($c-1),$cell,$charformat) if ($c == 2 or $c <= 3);
    }
  }
}

#write RNA Statistics worksheet for RNAseq
if ($pipeline == 2 or $pipeline == 5) {
    $r = 0;
  foreach my $row (@rna){
    $r += 1;
    my @row = split /\t/,$row;
    my $c = 0;
    foreach my $cell (@row){
      $c += 1;
      $cell *= 100 if ($r > 1 and ($c >= 11 and $c <= 18));
      $rna->write(($r-1),($c-1),$cell,$subtitleformat) if ($r == 1);
      $rna->write(($r-1),($c-1),$cell,$charformat) if ($r > 1 and $c == 1);
      $rna->write(($r-1),($c-1),$cell,$numformat) if ($r > 1 and ($c >= 2 and $c <= 10));
      $rna->write(($r-1),($c-1),$cell,$performat) if ($r > 1 and ($c >= 11 and $c <= 18));
      $rna->write(($r-1),($c-1),$cell,$decformat) if ($r > 1 and ($c >= 19));
    }
  }
}

#write summary worksheet for all pipelines
$r = 0;
foreach my $row (@sums){
  $r += 1;
  my @row = split /\t/,$row;
  my $c = 0;
  foreach my $cell (@row){
    $c += 1;
    $summary->write(($r-1),($c-1),$cell,$subtitleformat) if ($r == 1);
    if ($pipeline == 1) {   #ChIPseq
      $summary->write(($r-1),($c-1),$cell,$charformat) if ($r > 1 and $c == 1);
      $summary->write(($r-1),($c-1),$cell,$numformat) if ($r > 1 and ($c == 2 or $c == 4 or $c == 5 or $c == 7 or $c == 9));
      $summary->write(($r-1),($c-1),$cell,$performat) if ($r > 1 and ($c == 3 or $c == 6 or $c == 8 or $c == 10 or $c == 11));
    }
    elsif ($pipeline == 2 or $pipeline == 5){ #Tophat and Star RNAseq
      $cell *= 100 if ($r > 1 and ($c >= 12 and $c <= 18));
      $summary->write(($r-1),($c-1),$cell,$charformat) if ($r > 1 and $c == 1);
      $summary->write(($r-1),($c-1),$cell,$numformat) if ($r > 1 and ($c == 2 or $c == 4 or $c == 5 or $c == 7 or $c == 9));
      $summary->write(($r-1),($c-1),$cell,$performat) if ($r > 1 and ($c == 3 or $c == 6 or $c == 8 or ($c >= 10 and $c <= 18)));
      $summary->write(($r-1),($c-1),$cell,$decformat) if ($r > 1 and ($c == 19));
    }
   elsif ($pipeline == 3 or $pipeline == 4){  #Exome-seq and gDNA
      $summary->write(($r-1),($c-1),$cell,$charformat) if ($r > 1 and $c == 1);
      $summary->write(($r-1),($c-1),$cell,$numformat) if ($r > 1 and ($c == 2 or $c == 4 or $c == 5 or $c == 7 or $c == 14));
      $summary->write(($r-1),($c-1),$cell,$performat) if ($r > 1 and ($c == 3 or $c == 6 or ($c >= 8 and $c <= 13)));
    }
  }
}
$r -= 1;

#insert charts
#Chipseq, RNAseq, and Exome-seq/gDNA chart 1.  combo column and line graph
my $chart1a = $report->add_chart( type => 'column', embedded => 1 );
$chart1a->add_series(
    name    => ['Summary',0,1],
    values  => ['Summary',1,$r,1,1],
    categories => ['Summary',1,$r,0,0],
    data_labels => {value => 1, position => 'inside_base' },
);
my $chart1b = $report->add_chart( type => 'line', embedded => 1 );
$chart1b->add_series(
    name    => ['Summary',0,2],
    values  => ['Summary',1,$r,2,2],
    categories => ['Summary',1,$r,0,0],
    marker => {type => 'diamond'},
    data_labels => {value => 1, position => 'above', num_format => '#00' },
    y2_axis => 1,
);
$chart1b->set_y2_axis( min => 70, max => 100, major_unit => 5 );
$chart1a->combine($chart1b);
$chart1a->set_legend( position => 'top' );
$chart1a->set_title( name => 'Sample Yield & Percent of Bases >= Q30' );
$summary->insert_chart( ($r+2),1,$chart1a,0,0,1.25,1.25 );

#Chipseq chart 2. column chart with 2 series
if ($pipeline == 1) {
  my $chart2 = $report->add_chart( type => 'column', embedded => 1 );
  $chart2->add_series(
      name    => ['Summary',0,5],
      values  => ['Summary',1,$r,5,5],
      categories => ['Summary',1,$r,0,0],
      data_labels => {value => 1, position => 'inside_end', num_format => '#00' },
  );
  $chart2->add_series(
      name    => ['Summary',0,7],
      values  => ['Summary',1,$r,7,7],
      categories => ['Summary',1,$r,0,0],
      data_labels => {value => 1, position => 'inside_end', num_format => '#00' },
  );
  $chart2->add_series(
      name    => ['Summary',0,9],
      values  => ['Summary',1,$r,9,9],
      categories => ['Summary',1,$r,0,0],
      data_labels => {value => 1, position => 'inside_end', num_format => '#00' },
  );
  $chart2->add_series(
      name    => ['Summary',0,10],
      values  => ['Summary',1,$r,10,10],
      categories => ['Summary',1,$r,0,0],
      data_labels => {value => 1, position => 'inside_end', num_format => '#00' },
  );
  $chart2->set_legend( position => 'top' );
  $chart2->set_title( name => 'Percent Total Reads After Trimming, Percent Total Mapped Reads, Percent Uniquely Mapped Reads, & Percent Non-duplicated Reads' );
  $chart2->set_y_axis( max => 100 );
  $summary->insert_chart( ($r+2),11,$chart2,0,0,1.25,1.25 );
}
elsif ($pipeline == 2 or $pipeline == 5){
  ##RNAseq line chart with 2 series
  my $chart2 = $report->add_chart( type => 'line', embedded => 1 );
  $chart2->add_series(
      name    => ['Summary',0,7],
      values  => ['Summary',1,$r,7,7],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'diamond'},
  );
  $chart2->add_series(
      name    => ['Summary',0,9],
      values  => ['Summary',1,$r,9,9],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'square'},
  );
  $chart2->add_series(
      name    => ['Summary',0,10],
      values  => ['Summary',1,$r,10,10],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'triangle'},
  );
  $chart2->set_legend( position => 'top' );
  $chart2->set_y_axis( max => 100 );
  $chart2->set_title( name => 'Percent Total Mapped Reads, Percent Uniquely Mapped Reads, & Percent Non-duplicated Reads' );
  $summary->insert_chart( ($r+2),11,$chart2,0,0,1.25,1.25 );
  
  #RNAseq Radar chart
  #pct ribosomal, 
  my $chart3 = $report->add_chart( type => 'radar', embedded => 1 );
  $chart3->add_series(
      name    => ['Summary',0,11],
      values  => ['Summary',1,$r,11,11],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'diamond'},
  );
  $chart3->add_series(
      name    => ['Summary',0,12],
      values  => ['Summary',1,$r,12,12],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'square'},
  );
  $chart3->add_series(
      name    => ['Summary',0,13],
      values  => ['Summary',1,$r,13,13],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'triangle'},
  );
  $chart3->add_series(
      name    => ['Summary',0,14],
      values  => ['Summary',1,$r,14,14],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'x'},
  );
  $chart3->add_series(
      name    => ['Summary',0,15],
      values  => ['Summary',1,$r,15,15],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'star'},
  );
  $chart3->add_series(
      name    => ['Summary',0,16],
      values  => ['Summary',1,$r,16,16],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'circle'},
  );
  $chart3->set_legend( position => 'top' );
  $chart3->set_title( name => 'RNA Statistics' );
  $summary->insert_chart( ($r+2),21,$chart3,0,0,1.50,1.50 );
}
elsif ($pipeline == 3){
  ###Exome-seq second chart
  my $chart2 = $report->add_chart( type => 'line', embedded => 1 );
  $chart2->add_series(
      name    => ['Summary',0,5],
      values  => ['Summary',1,$r,5,5],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'diamond'},
  );
  $chart2->add_series(
      name    => ['Summary',0,7],
      values  => ['Summary',1,$r,7,7],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'square'},
  );
  $chart2->add_series(
      name    => ['Summary',0,8],
      values  => ['Summary',1,$r,8,8],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'triangle'},
  );
  $chart2->add_series(
      name    => ['Summary',0,11],
      values  => ['Summary',1,$r,11,11],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'circle'},
  );
  $chart2->set_legend( position => 'top' );
  $chart2->set_y_axis( max => 100 );
  $chart2->set_title( name => 'Percent Total Reads After Trimmed, Percent Total Mapped, Percent Mapped On-Target, Percent Non-duplicated, & Percent Coverage >= 30X' );
  $summary->insert_chart( ($r+2),11,$chart2,0,0,1.25,1.25 );
}
elsif ($pipeline == 4){
  ###gDNA second chart
  my $chart2 = $report->add_chart( type => 'line', embedded => 1 );
  $chart2->add_series(
      name    => ['Summary',0,5],
      values  => ['Summary',1,$r,5,5],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'diamond'},
  );
  $chart2->add_series(
      name    => ['Summary',0,7],
      values  => ['Summary',1,$r,7,7],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'square'},
  );
  $chart2->add_series(
      name    => ['Summary',0,8],
      values  => ['Summary',1,$r,8,8],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'triangle'},
  );
  $chart2->add_series(
      name    => ['Summary',0,11],
      values  => ['Summary',1,$r,11,11],
      categories => ['Summary',1,$r,0,0],
      marker => {type => 'circle'},
  );
  $chart2->set_legend( position => 'top' );
  $chart2->set_y_axis( max => 100 );
  $chart2->set_title( name => 'Percent Total Reads After Trimming, Percent Total Mapped, Percent Non-duplicated, Percent Coverage >= 5X, & , Percent GC' );
  $summary->insert_chart( ($r+2),11,$chart2,0,0,1.25,1.25 );
}
print STDERR "\nSamples reported: $r\n";

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