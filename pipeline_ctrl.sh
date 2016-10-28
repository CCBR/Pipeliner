#!/bin/bash

module load python/3.5

#PBS -N PipelineTest

cd $SLURM_SUBMIT_DIR
#cd $PBS_O_WORKDIR
D=/data/CCBR/dev/Pipeline/Pipeliner
R=/scratch/dwheeler/test2

snakemake  -s $R/Snakefile -d $R --printshellcmds --cluster-config $D/cluster.json --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem}" -j 50 --rerun-incomplete --stats $R/Reports/initialqc.stats -T 2>&1|tee -a $R/Reports/snakemake.log



cd $R/Reports && python summarize_fastqc_reports.py *fastqc.zip ../QC ./
cd $R/Reports && $R/Reports/ngser.pl ../QC ./
cd $R/Reports && $R/Reports/stats2html.py

