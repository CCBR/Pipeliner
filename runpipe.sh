#!/bin/bash

#. /usr/local/Modules/default/init/bash

module load python/3.4.0
module load graphviz

pipehome=$(dirname $0)

cd $pipehome && ./pipeliner.py &
