#!/bin/bash

#. /usr/local/Modules/default/init/bash

module load python/3.5
module load snakemake/5.1.3
module load graphviz

export pipehome=$(dirname $0)
export PYTHONPATH=$PYTHONPATH:$pipehome/local/lib/

cd  $pipehome;
export pipehome=`pwd`;
python pipeliner2.py &
echo $pipehome;
