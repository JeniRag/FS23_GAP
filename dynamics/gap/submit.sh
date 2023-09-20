#!/bin/bash

#SBATCH --job-name=NVT-gap
#SBATCH --get-user-env
#SBATCH --mem-per-cpu=550MB
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --cpus-per-task=1

module load anaconda3
source activate gap

#change to your system
IPI=/home/sragav/data/conda/envs/gap/bin/i-pi
IPI_DRIVER=/home/sragav/data/conda/envs/gap/bin/i-pi-py_driver

DYNAMICS_PATH=/data/sragav/our_pbe_set/dynamics
SCRATCH_PATH=/scratch/sragav/our_pbe_set/dynamics
model=${DYNAMICS_PATH}/models/demo_model/pbe_model_grad_0.9.json

Gap=gap
if [ -e /tmp/ipi_$Gap ]; then
	rm /tmp/ipi_$Gap
	echo "/tmp/ipi_$Gap was removed"
else
	echo "/tmp/ipi_$Gap didn't exist"
fi
 
echo "Model: $model"
echo "IPI NODE ADDRESS IS $HOSTNAME"

sed -e "s:<address>.*</address>:<address>$Gap</address>:" ./input.nvt.xml > initial.xml #change input file in input.nvt.xml
# make sure to adapt total_time. Should be slightly smaller than the computation time you requested, in order to stop early and produce a RESTART file.

time $IPI ${SCRATCH_PATH}/$Gap/initial.xml > ipi.log&

sleep 15

$IPI_DRIVER -m rascal -u -a $Gap -o ${model},${DYNAMICS_PATH}/inputs/idx_0.xyz #idx_0.xyz. Doesn't have to be your input file. It should just have the same atoms and number of atoms as in you input files.
echo "simulation done"
conda deactivate
