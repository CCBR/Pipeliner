#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';
use String::Random;
$r=new String::Random;

$id=$r->randpattern("ccccnnnn");

#INPUT

$vcf=shift;

my $status = 0;
my $cmd = '';
my $time = 300; #amount of time (seconds) to wait between pings to AVIA



$cmd="/usr/local/bin/curl  https://avia-abcc.ncifcrf.gov/apps/site/upload_viz -X POST -F user_file=@$vcf -F user.ver=hg19 -F user_inputformat=vcf4 -F user_api=cFMdtEdwm34iVzXOZ6 -F "user_email=justin.lack|nih.gov" --insecure -F user_id=$id";

    system($cmd);

#sleep($time);
open(LOG,">$id.avia.log");
my $pings=0;			
while ($status == 0) {
	$cmd = "curl https://avia-abcc.ncifcrf.gov/apps/site/status/?id=$id > $id.avia_status.txt";
	$pings++;
	print LOG "$pings\n";
	print STDERR "$pings\n";
	system($cmd);
	
	my @line = ();

	open G, "<$id.avia_status.txt";
	while (<G>){
		chomp;
  		last if m/^$/;
  		next if ($_ !~ m/INFO/);
	  	if ($_ =~ m/is ready for download/) {
			$status++;
		}
		else {
			sleep($time);
		}
	}
close G;
}

close OUT;

$cmd = "wget -O $id.full_annot.txt.zip https://avia-abcc.ncifcrf.gov/apps/site/download/?file=$id/$id.full_annot.txt.zip";
system($cmd);
