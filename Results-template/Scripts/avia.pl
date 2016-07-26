#!/bin/env perl
use List::Util 'shuffle';
use String::Random;
$r=new String::Random;

$id="viz-ccbr-".$r->randpattern("ccccnnnn");

#INPUT

$vcf=shift;

$status = 0;
$cmd = '';
$time = 300; #amount of time (seconds) to wait between pings to AVIA

$cmd="/usr/local/bin/curl  https://avia-abcc.ncifcrf.gov/apps/site/upload_viz -X POST -F user_file=\@$vcf -F user.ver=hg19 -F user_inputformat=bed -F user_api=cFMdtEdwm34iVzXOZ6 -F 'user_email=justin.lack|nih.gov' --insecure -F user_id=$id";

print STDERR "Executing command: $cmd\n";

    system($cmd);

print STDERR "Begin waiting for the result.\n";

sleep($time);
open(LOG,">$id.avia.log");
$pings=0;			
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
  		if (($_ =~ m/INFO/) || ($_ =~ m/ERROR/)) {
	  		if ($_ =~ m/is ready for download/) {
				$status++;
			}
			elsif ($_ =~ m/is not accessible/) {
				die;
			}
			else {
				sleep($time);
			}
	}
}
close G;
}

close OUT;

$cmd = "wget -O full_annot.txt.zip https://avia-abcc.ncifcrf.gov/apps/site/download/?file=$id/$id.full_annot.txt.zip";
system($cmd);
