#!/bin/bash

#SBATCH --partition=standard
#SBATCH --qos=user_qos_denard
#SBATCH --ntasks=200
#SBATCH --mem-per-cpu=2GB
#SBATCH --job-name=tennesen_no_gnu_parallel
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=denard
#SBATCH --mail-user=jesus.murga@uab.cat

DPATH="/xdisk/denard/mig2020/rsgrps/denard"
module load slim

for i in {1..1000}
do
        srun -n1 --exclusive slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=8.22e-15 -d pposL=0.011536 -d pposH=7.94794e-5 -d p=\'${DPATH}/jmurga/simulations/tennesen/tennesen_0.4_0.3_0.999\' -d nF=$i ${DPATH}/jmurga/src/tennesen.slim &
done
wait

# srun slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.14e-8 -d pposL=0.011536 -d pposH=7.94794e-5 -d p=\'${DPATH}/jmurga/simulations/tennesen/tennesen_0.4_0.3_0.4\' -d nF=$SLURM_ARRAY_TASK_ID ${DPATH}/jmurga/src/tennesen.slim


#seq 24000 25000 | parallel -u -j100 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.14e-8 -d pposL=0.011536 -d pposH=7.94794e-5 -d p=\'${DPATH}/jmurga/simulations/tennesen/tennesen_0.4_0.3_0.4\' -d nF={} ${DPATH}/jmurga/src/tennesen.slim"
