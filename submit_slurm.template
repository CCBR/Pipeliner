#!/bin/bash


cd $SLURM_SUBMIT_DIR

PrjName=$1

if [[ -z "${param// }" ]];
then
	Slrm_J="-J $PrjName"
	echo "Setting an slurm job submission name option: $Slrm_J"
else 
	Slrm_J=""
	echo "No slurm job submission name set! PrjName=$PrjName"
fi


#. /usr/local/Modules/default/init/bash

module load python/3.5
module load snakemake/5.1.3

D=/home/fake
R=/home/fake

modtime1=`stat -c %y $R/Reports/snakemake.log|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
modtime2=`stat -c %y $R/Reports/makeasnake.log|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
mv $R/Reports/snakemake.log $R/Reports/snakemake.log.$modtime1
mv $R/Reports/makeasnake.log $R/Reports/makeasnake.log.$modtime2

touch $R/Reports/makeasnake.log
touch $R/Reports/snakemake.log
BUY_IN_NODE=$(ACCOUNT_SPONSOR=$(sacctmgr -rn list user | awk '{print $2}') && scontrol show partitions | grep -i $ACCOUNT_SPONSOR -B1 | grep '^PartitionName' | cut -d '=' -f2 | grep -iv 'gpu'| tr '\n' ',' | sed 's/.$//')

if [[ $BUY_IN_NODE =~ "ccr" ]]; then
	BUY_IN_NODE="ccr,norm"
	sed -i 's/\"norm\"/\"ccr,norm\"/g' $R/cluster.json
else
	BUY_IN_NODE="norm"
fi

sbatch $Slrm_J --partition=$BUY_IN_NODE --gres=lscratch:200 --time=120:00:00 --mail-type=BEGIN,END,FAIL $R/pipeline_ctrl.sh
#sbatch $Slrm_J --partition=centos7 --gres=lscratch:200 --time=24:00:00 --mail-type=BEGIN,END,FAIL $R/pipeline_ctrl.sh
#sbatch $Slrm_J --partition=norm --gres=lscratch:200 --time=10-00:00:00 --mail-type=BEGIN,END,FAIL $R/pipeline_ctrl.sh
#sbatch --partition=ccr,norm --gres=lscratch:200 --time=120:00:00  $R/pipeline_ctrl.sh
