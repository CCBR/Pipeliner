#!/bin/bash
module load ngsplot
bam=$1
genome=$2
name=${bam/.bam/}
for region in tss tes genebody ; do
for go in max km; do
ngs.plot.r -G $genome -GO $go -R $region -C $bam -O ${name}.${region}.${go} -T $bam -L 3000 -FL 225
done
done