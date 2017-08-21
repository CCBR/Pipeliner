#!/bin/bash
fastqgz=$1
taxatxt=$2
fastqgz_nreads=`zcat $fastqgz |wc -l`
fastqgz_nreads=$((fastqgz_nreads/4))
#fastqgz_nreads=56853786
samplename=${fastqgz/.fastq.gz/}
samplename=`echo $samplename|awk -F "/" '{print $NF}'`
total=`cat $taxatxt|awk -v n=$fastqgz_nreads '{sum=sum+$1}END{print sum}'`
echo $total|awk -v f=$samplename -v n=$fastqgz_nreads '{printf("%s\t%d\t%.4f\t#total\n",f,$1,100*$1/n)}'
cat $taxatxt|sed "s/|/ /g"|awk -v f=$samplename -v n=$fastqgz_nreads '{p=100*$1/n; if (p>0.005) {printf("%s\t%d\t%.4f\t%s\n",f,$1,100*$1/n,$NF)}}'
