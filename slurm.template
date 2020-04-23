#!/bin/bash
set -e

module load python/3.5
module load snakemake/5.1.3

#PBS -N PipelineTest

D=/home/fake
R=/home/fake

cd $D
module load git/2.26.2
git log |head -n3 > $R/Reports/ccbr_pipeliner_git_commit_version.txt

cd $SLURM_SUBMIT_DIR
#cd $PBS_O_WORKDIR


snakemake -s $R/Snakefile -d $R --printshellcmds --cluster-config $R/cluster.json --keep-going --restart-times 1 --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname}" -j 500 --rerun-incomplete --stats $R/Reports/initialqc.stats -T 2>&1|tee -a $R/Reports/snakemake.log

# snakemake -s $R/Snakefile -d $R --printshellcmds --cluster-config $D/cluster.json --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem}" -j 100 --rerun-incomplete --stats $R/Reports/initialqc.stats -T 2>&1|tee -a $R/Reports/snakemake.log

cd $R
mv $R/slurm-*.out $R/slurmfiles/
if [ -f $R/HPC_usage_table.txt ]; then
	modtime1=`stat -c %y $R/HPC_usage_table.txt|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
	mv $R/HPC_usage_table.txt $R/HPC_usage_table.txt.${modtime1}
fi
perl Scripts/summarize_usage.pl
python Scripts/filter_usage_summary.py > $R/HPCusagetable.txt.tmp
mv $R/HPCusagetable.txt.tmp $R/HPC_usage_table.txt

#cd $R/Reports && python summarize_fastqc_reports.py *fastqc.zip ../QC ./
cd $R/Reports && $R/Reports/ngser.pl ../QC ./
#cd $R/Reports && $R/Reports/stats2html.py
rm -rf $R/bamstats
rm -rf $R/*bam.cnt
