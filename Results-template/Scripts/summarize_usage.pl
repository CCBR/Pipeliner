#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my @line=();
my $count='';
my @jobs=();
my $cmd = '';
my $temp='';

$cmd = 'ls slurmfiles/slurm*out > slurmfiles.txt';
system($cmd);

open H, "<slurmfiles.txt";
while (<H>){
	chomp;
	last if m/^$/;
	@line = split;
	if ($line[0] !~ m/Hugo/) {
		$temp=((split 'urm-', $line[0])[1]);
		push (@jobs, (split '.ou', $temp)[0]);
	}
}
close H;
my $a=0;

for ($a=0; $a<@jobs; $a++) {
	$cmd = 'jobhist ' . $jobs[$a] . ' >> job_usage.txt';
	system($cmd);
}

my $b=0;
my $c=0;
my $name='';
my $on=0;

open C, ">HPC_usage_table.txt";
print C "JobName\tJobid\tPartition\tState\tNodes\tCPUs\tWalltime\tRuntime\tMemReq\tMemUsed\tNodelist\n";

open G, "<job_usage.txt";
while (<G>) {
	chomp;
#	last if m/^$/;
	@line=split /\s+/, $_;
	if ($on==1) {
		$on++;
	}
	if ($_ =~ m/job-name/) {
		$b=0;
		for ($b=0; $b<@line; $b++) {
			if ($line[$b] =~ m/job-name/) {
				$name=((split '=', $line[$b])[1]);
			}
		}
	}
	elsif ($_ =~ m/Nodelist/) {
		$on++;
	}
	if ($on==2) {
		$on=0;
		$c=0;
		print C "$name";
		for ($c=0; $c<@line; $c++) {
			print C "\t$line[$c]";
		}
	print C "\n";
	}
}
close C;
close G;
close H;
$cmd = 'rm slurmfiles.txt job_usage.txt';
system($cmd);