#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root

MIGRATION_FILE=$1
TRUE_HIST=efficiency_numerator
WARPED_FILE=$2
RECO_HIST=selected_signal_reco
TARGET=$3

PDG=$4

#declare -a arr=("BjorkenX")
#declare -a arr=("pTmu" "pZmu" "BjorkenX" "BjorkenY" "Emu" "Enu" "pTmu_pZmu")
declare -a arr=("pTmu" "pZmu" "BjorkenX" "BjorkenY" "Emu" "Enu" "Erecoil")
#declare -a arr=("pTmu_pZmu")
## loop through above array (quotes are important if your elements may contain spaces)
for i in "${arr[@]}"
do
    num_dim=1
    VARIABLE=$i
    if [[ "$VARIABLE" == "pTmu_pZmu" ]]; then
        num_dim=2
        TRUE_HIST=truth
        RECO_HIST=reco
    fi
    echo $TARGET
    OUTFILE_NAME=$(basename $1)_$(basename $2)_${TARGET}_${VARIABLE}.root
    echo $OUTFILE_NAME
    args=(--output_file Warping_$OUTFILE_NAME --data $RECO_HIST --data_file $WARPED_FILE --data_truth $TRUE_HIST --data_truth_file $WARPED_FILE --migration migration --migration_file $MIGRATION_FILE --reco $RECO_HIST --reco_file $MIGRATION_FILE --truth $TRUE_HIST --truth_file $MIGRATION_FILE --num_iter 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 --num_uni 100 -C 2 --variable $VARIABLE --target_string $TARGET --num_dim $num_dim --PDG $PDG)
    #[[ "$VARIABLE" == "BjorkenX" ]] && args+=(--exclude_bins 7,8)
    TransWarpExtractionMod "${args[@]}"
done

