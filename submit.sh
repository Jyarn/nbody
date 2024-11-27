#!/bin/sh
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:10:00
#SBATCH --job-name=nbody_sim
#SBATCH --output=scaling_res_%j.txt

export OMP_NUM_THREADS=16
module load TeachEnv/2022a gcc openmpi

make build
time make run

echo "scaling complete!"
