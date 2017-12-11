#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';
use File::Temp;
use File::Copy;
my $file = tmpnam();
my @t = split/\//,$file;
$file = "/scratch/CCBRPipeliner/usage/" . $t[-1] . ".usage.txt";
my $file2 = "/scratch/CCBRPipeliner/usage/" . $t[-1] . ".run.json";
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

#@jobs=('54856894');

my $a=0;

for ($a=0; $a<@jobs; $a++) {
	$cmd = 'jobhist ' . $jobs[$a] . ' >> job_usage.txt';
	system($cmd);
}

my $b=0;
my $c=0;
my $name='';
my $jobid='';
my $on=0;

open C, ">HPC_usage_table.txt";
print C "JobName\tJobid\tPartition\tState\tNodes\tCPUs\tWalltime\tRuntime\tMemReq\tMemUsed\tNodelist\tMaxCPUUsed\tQueuetime\tCPUHours\tAccount\tUsername\n";

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
	elsif ($_ =~ m/^JobId/) {
		$b=0;
		$jobid=$line[2];
	}
	elsif ($_ =~ m/Nodelist/) {
		$on++;
	}
	if ($on==2) {
		$on=0;
		$c=0;
		$cmd="jobdata " . $jobid . " >jobdata.txt";
		system($cmd);
		my $maxcpu='';
		my $qtime='';
		my $cpuhrs='';
		my $account='';
		my $username='';
		open H, "<jobdata.txt";
		while(<H>) {
		chomp;
		my @line2=split /\t/, $_;
		if ($line2[0] =~ m/^max_cpu_used\b/) {
			$maxcpu=$line2[1];
		}
		elsif ($line2[0] =~ m/^queued\b/) {
			$qtime=$line2[1];
		}
		elsif ($line2[0] =~ m/^cpu_hours\b/) {
			$cpuhrs=$line2[1];
		}
		elsif ($line2[0] =~ m/^account\b/) {
			$account=$line2[1];
		}
		elsif ($line2[0] =~ m/^username\b/) {
			$username=$line2[1];
		}
		}
		close H;
		system("rm -f jobdata.txt");
		print C "$name";
		for ($c=0; $c<@line; $c++) {
			print C "\t$line[$c]";
		}
		print C "\t$maxcpu";
		print C "\t$qtime";
		print C "\t$cpuhrs";
		print C "\t$account";
		print C "\t$username";
		print C "\n";
	}
}
close C;
close G;
$cmd = 'rm slurmfiles.txt job_usage.txt';
system($cmd);
copy("HPC_usage_table.txt",$file);
copy("run.json",$file2);


