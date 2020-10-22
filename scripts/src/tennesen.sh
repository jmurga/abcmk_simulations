#!/bin/bash

#SBATCH --array=1-10
#SBATCH --partition=standard
#SBATCH --qos=user_qos_denard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=1GB
#SBATCH --job-name=tennesen_model_jmurga
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=denard
#SBATCH --mail-user=jesus.murga@uab.cat

PATH="/xdisk/denard/mig2020/rsgrps/denard/"

slim -d popSize=7310 -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.6e-11 -d pposL=0.01153603001 -d pposH=7.9479354e-05 -d p=\'/home/jmurga/mkt/202004/rawData/simulations/tennesen/tennesen_0.4_0.3_0.999\' -d nF=1 /home/jmurga/mkt/202004/scripts/slimRecipes/tennesen.slim