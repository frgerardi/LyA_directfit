#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --constraint=haswell
#SBATCH --mail-user=francesca.gerardi.19@ucl.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output /global/cfs/projectdirs/desi/users/fgerardi/mpi_sampling/logs/test.log


echo $shell
source activate COBAYA_mpi2
cd /global/cfs/projectdirs/desi/users/fgerardi/mpi_sampling

srun --tasks-per-node=1 cobaya-run info.yaml
