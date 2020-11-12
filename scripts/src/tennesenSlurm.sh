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

DPATH="/xdisk/denard/mig2020/rsgrps/denard/"

for i in {9000..9500}
srun slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.14e-8 -d pposL=0.011536 -d pposH=7.94794e-5 -d p=\'${DPATH}/jmurga/simulations/tennesen/tennesen_0.4_0.3_0.4\' -d nF=${i} ${DPATH}/jmurga/src/tennesen.slim &
done


#SBATCH --partition=standard
#SBATCH --qos=user_qos_denard
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=60
#SBATCH --ntasks=240
#SBATCH --mem-per-cpu=2GB
#SBATCH --job-name=tennesen
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=denard
#SBATCH --mail-user=jesus.murga@uab.cat

DPATH="/xdisk/denard/mig2020/rsgrps/denard"

module load slim
module load parallel


#seq 9500 10000 | parallel -u -j240 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.14e-8 -d pposL=0.011536 -d pposH=7.94794e-5 -d p=\'${DPATH}/jmurga/simulations/tennesen/tennesen_0.4_0.3_0.4\' -d nF={} ${DPATH}/jmurga/src/tennesen.slim"

~                                        
