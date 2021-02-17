DOFE="${HOME}/mkt/202004/rawData/dofe"
SRC="${HOME}/mkt/202004/scripts/src"
POLYDFE="${HOME}/.conda/envs/abcmk/bin/polydfe"
PARALLEL="${HOME}/.conda/envs/abcmk/bin/parallel"
OUTPUT="${HOME}/mkt/202004/results/polydfe"
MODEL="noDemog"
NTHREADS=20

for CASE in `ls ${DOFE}/${MODEL} | fgrep 0.999`;
do 

    ANALYSIS="${DOFE}/${MODEL}/${CASE}"
    RESULTS="${OUTPUT}/${MODEL}/${CASE}"
    mkdir -p ${RESULTS}

    echo "############${CASE}"
    echo "####ALPHA_DIV"
    time ${PARALLEL} -u -j${NTHREADS} "${POLYDFE} -d ${ANALYSIS}/${CASE}_polydfe_{}_10.tsv -i ${SRC}/init_model_BandC_10.txt 7 -m B > ${RESULTS}/${CASE}_alphadiv_{}.txt" ::: $(seq 1 100)
    echo "####ALPHA_DFE"
    ${PARALLEL} "${POLYDFE} -d ${ANALYSIS}/${CASE}_polydfe_{}_10.tsv -i ${SRC}/init_model_BandC_10.txt 6 -m B -w > ${RESULTS}/${CASE}_alphadfe_{}.txt" ::: $(seq 1 100)
done


MODEL="isolation"

# echo "#################### isolation fixed"
# /home/jmurga/.conda/envs/abcmk/bin/parallel --link /home/jmurga/.conda/envs/abcmk/bin/polydfe -d ${DOFE}/inputs/${MODEL}/{1}_polydfe_20.tsv -i /home/jmurga/mkt/202004/scripts/src/init_model_BandC_fixed.txt {2} -m B ">" ${DOFE}/outputs/${MODEL}/fixed/{1}.polydfe ::: ${MODEL}_0.4_0.1_0.2  ${MODEL}_0.4_0.1_0.8 ${MODEL}_0.4_0.1_0.999 ${MODEL}_0.4_0.2_0.2 ${MODEL}_0.4_0.2_0.8 ${MODEL}_0.4_0.2_0.999 ${MODEL}_0.4_0.3_0.2 ${MODEL}_0.4_0.3_0.8 ${MODEL}_0.4_0.3_0.999 ::: 1 1 1 2 2 2 3 3 3

# echo "#################### isolation r"
# /home/jmurga/.conda/envs/abcmk/bin/parallel --link /home/jmurga/.conda/envs/abcmk/bin/polydfe -d ${DOFE}/inputs/${MODEL}/{1}_polydfe_20.tsv -i /home/jmurga/mkt/202004/scripts/src/init_model_BandC_r.txt {2} -m B ">" ${DOFE}/outputs/${MODEL}/r_i/{1}.polydfe ::: ${MODEL}_0.4_0.1_0.2  ${MODEL}_0.4_0.1_0.8 ${MODEL}_0.4_0.1_0.999 ${MODEL}_0.4_0.2_0.2 ${MODEL}_0.4_0.2_0.8 ${MODEL}_0.4_0.2_0.999 ${MODEL}_0.4_0.3_0.2 ${MODEL}_0.4_0.3_0.8 ${MODEL}_0.4_0.3_0.999 ::: 1 1 1 2 2 2 3 3 3

echo "#################### isolation nonFixed"
/home/jmurga/.conda/envs/abcmk/bin/parallel --link /home/jmurga/.conda/envs/abcmk/bin/polydfe -d ${DOFE}/inputs/${MODEL}/{1}_polydfe_20.tsv -i /home/jmurga/mkt/202004/scripts/src/init_model_BandC.txt {2} -m B ">" ${DOFE}/outputs/${MODEL}/nonFixed/{1}.polydfe ::: ${MODEL}_0.4_0.1_0.2  ${MODEL}_0.4_0.1_0.8 ${MODEL}_0.4_0.1_0.999 ${MODEL}_0.4_0.2_0.2 ${MODEL}_0.4_0.2_0.8 ${MODEL}_0.4_0.2_0.999 ${MODEL}_0.4_0.3_0.2 ${MODEL}_0.4_0.3_0.8 ${MODEL}_0.4_0.3_0.999 ::: 1 1 1 2 2 2 3 3 3

