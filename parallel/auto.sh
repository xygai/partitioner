#!/bin/bash

#SBATCH -n 16        # 8 cores = 1 node on lonsdale
#SBATCH -p compute
#SBATCH -t 6:00:00  # 6 hours
#SBATCH -J xgai

# source the module commands: note the new module location
#source /home/support/spack/paddy/modules.sh

# load the modules used to build the xhpl binary
module load cports
module load gcc/9.2.0-gnu openmpi/3.1.3-gnu

# run it
mpirun ./build/parallel auto.graph auto.part 16

