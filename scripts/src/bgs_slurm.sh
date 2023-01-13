#!/bin/bash

#SBATCH --partition=standard
#SBATCH --job-name=no_demog
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesus.murga@uab.cat
#SBACTH --array=1-50000

DPATH="/xdisk/denard/mig2020/rsgrps/denard"

module load slim
module load parallel

# Each command takes about 12hours to run using this configuration

##########BGS 0.2

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=6.608754e-07 -d pposL=0.003877 -d pposH=0.000311 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.1_0.2\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=6.608754e-07 -d pposL=0.007739 -d pposH=0.000207 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.2_0.2\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=6.608754e-07 -d pposL=0.011586 -d pposH=0.000103 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.3_0.2\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

##########BGS 0.4
parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=3.762519e-07 -d pposL=0.003877 -d pposH=0.000311 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.1_0.4\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=3.762519e-07 -d pposL=0.007739 -d pposH=0.000207 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.2_0.4\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=3.762519e-07 -d pposL=0.011586 -d pposH=0.000103 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.3_0.4\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

##########BGS 0.8
parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=9.162832e-08 -d pposL=0.003877 -d pposH=0.000311 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.1_0.8\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=9.162832e-08 -d pposL=0.007739 -d pposH=0.000207 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.2_0.8\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=9.162832e-08 -d pposL=0.011586 -d pposH=0.000103 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.3_0.8\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

##########BGS 0.999
parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=4.108304e-10 -d pposL=0.003877 -d pposH=0.000311 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.1_0.999\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=4.108304e-10 -d pposL=0.007739 -d pposH=0.000207 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.2_0.999\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`

parallel -j200 -u "slim -d N=500 -d sampleSize=500 -d wS=10 -d sS=500 -d muBgs=4.108304e-10 -d pposL=0.011586 -d pposH=0.000103 -d nF={} -d output=\'${DPATH}/jmurga/simulations/no_demog/no_demog_0.4_0.3_0.999\' ${DPATH}/jmurga/src/bgs.slim" ::: `seq 1 50000`


