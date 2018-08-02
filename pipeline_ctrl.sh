#!/bin/bash

module load python/3.5

#PBS -N PipelineTest

cd $SLURM_SUBMIT_DIR
#cd $PBS_O_WORKDIR
D=/data/CCBR/dev/Pipeline/Pipeliner
R=/scratch/dwheeler/test2

snakemake  -s $R/Snakefile -d $R --printshellcmds --cluster-config $D/cluster.json --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem}" -j 50 --rerun-incomplete --keep-going --stats $R/Reports/initialqc.stats -T 2>&1|tee -a $R/Reports/snakemake.log

mv slurm*out slurmfiles
if [ -f $R/HPC_usage_table.txt ]; then
	modtime1=`stat -c %y $R/HPC_usage_table.txt|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
	mv $R/HPC_usage_table.txt $R/HPC_usage_table.txt.${modtime1}
fi

module load perl/5.24.3

perl Scripts/summarize_usage.pl
python Scripts/filter_usage_summary.py > $R/HPCusagetable.txt.tmp
mv $R/HPCusagetable.txt.tmp $R/HPC_usage_table.txt

cd $R/Reports && python summarize_fastqc_reports.py *fastqc.zip ../QC ./
cd $R/Reports && $R/Reports/ngser.pl ../QC ./
cd $R/Reports && $R/Reports/stats2html.py

