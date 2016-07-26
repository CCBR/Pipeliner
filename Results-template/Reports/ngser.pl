#!/bin/env perl

($indir,$outdir)=@ARGV;
$home=`pwd`;
print "Looking for reports in $indir\tPlacing summary in $outdir\n";

$files=`find $indir -name 'output*html'`;
$header="<head><link rel='stylesheet' type='text/css' href='stats.css'><style>div{overflow-x:auto}.autoResizeImage {max-width: 50%;height: auto;width: auto;vertical-align: top;}table.graphics {width: 1200px;table-layout:fixed;text-align:left;vertical-align:top;}</style></head>";
#$header="";
$contents="<H1 id='contents'>Table of Contents</H1><hr><ul>";
print "Found these reports:\n $files\n";


{
    local $/=undef;

for $F (split(/\n/,$files)){
    
    @f=split(/\//,$F);
    $title=pop @f;     
    @tables=();@final_tables=();
    print "Processing report $F\n";
    $contents.="<li><a href=#$title>$title</a></li>";
open(IN,"<$F");
    $R=<IN>;
#    print "Report looks like this: $R\n\n";
    close IN;
    if($header eq ""){$R=~/<head>.+<\/head>/s;$header=$&;}
    while($R=~/<table.+?<\/table>/gs){push @tables, $&;}

    $ntables=scalar @tables ;
    print "Incorporating $ntables tables.\n";

#    $page.="<H1>Report for $F</H1><br>".join("<br><br>",@tables);
    for($i=0;$i<$ntables;$i++){$t[$i]=transpose_table($tables[$i]);}

#    for ($i=0;$i<1;$i++){push @final_tables,$tables[$i];}
    for ($i=2;$i<4;$i++){push @final_tables,$t[$i];}
    @graph_tables=();
    for ($i=4;$i<$ntables;$i++){
	$tables[$i]=~s/<[^i].+?>//gs;
	
	push @graph_tables,$tables[$i];
    }

    $page.="<H1 id=$title>Report for $title</H1><br><br><a href='$F'><h5>Full Report for $title</a></h5>&nbsp&nbsp<a href='#contents'><h5>Back to Table of Contents</a></h5>".join("<br><br>",@final_tables);
    $graphs="<table class='graphics'><tr><td>".join("<td>",@graph_tables)."</tr></table>";
    $page.="<hr>".$graphs;

	}
}

$page="<html>$header<body>".$contents."</ul><hr>".$page."</body></html>";

$page=~s/<img src="/<img class='autoResizeImage' src="$indir\//g;
$page=~s/align=\".+?\"//sg;
$page=~s/width=\".+?\"//sg;
$page=~s/class=\"cnt\"//sg;

open(OUT,">$outdir/ngser.html");
print OUT $page;
close OUT;


sub transpose_table {

    my $table=$_[0];
    my $matrix;
    my $r=0;my $c=0;my $rw; my $co;
    my @R; my $thisRow; my $transposed;

    while($table=~/<tr(.*?)>(.+?)<\/tr>/gs){
	$row=$2;$r=0;$rowstyle=$1;
        while($row=~/<td(.*?)>(.+?)<\/td>/gs){
	    $item=$2;$colstyle=$1;
	    $matrix[$r][$c]=$item;
	    $r++;
	}
	$c++;
    }


    for($rw=0;$rw<$r;$rw++){
        @R=();
	for($co=0;$co<$c;$co++){push @R,$matrix[$rw][$co];}
	$thisRow=join("<td $rowstyle>",@R);
	$transposed.="<tr $colstyle><td>$thisRow</tr>";
    }

    $transposed="<table>$transposed</table>";
    return $transposed;

}
