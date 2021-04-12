#!/bin/bash

#SBATCH --partition=standard
#SBATCH --job-name=noDemog
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesus.murga@uab.cat
#SBACTH --array=1-50000%100

DPATH="/xdisk/denard/mig2020/rsgrps/denard"

module load slim
module load parallel

# Each command takes about 12hours to run using this configuration

##########BGS 0.2
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.1_0.2\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.2_0.2\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.3_0.2\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"

##########BGS 0.4
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.1_0.4\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.2_0.4\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.3_0.4\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"

##########BGS 0.8
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.1_0.8\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.2_0.8\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.3_0.8\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"

##########BGS 0.999
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.1_0.999\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.2_0.999\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
seq 1 50000 | parallel -u -j200 "slim -d sampleSize=661 -d codingLength=2000  -d weaklyStrength=10 -d strongStrength=500 -d bgsMutationRate=5.1471e-8 -d pposL=0.01153603001 -d pposH=7.9479354e-5 -d p=\'${DPATH}/jmurga/simulations/isolation/isolation_0.4_0.3_0.999\' -d nF={} ${DPATH}/jmurga/src/bgs.slim"
