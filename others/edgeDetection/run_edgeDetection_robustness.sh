#! /bin/bash

#! /bin/bash

declare -a DATASETS=("Birds")
declare -a NOISE_BLUR=("5" "9" "13" "17")
declare -a NOISE_SP=("0.04" "0.08" "0.12" "0.16")

METHOD_NAME=EdgeMaps

for DS_NAME in ${DATASETS[@]}; do

    # AVG BLUR
    for NOISE in ${NOISE_BLUR[@]}; do
    
    	DS_PATH="../../datasets/${DS_NAME}/avgblur/${NOISE}"
    	RESULTS_PATH="../../${METHOD_NAME}/avgblur/${NOISE}/${DS_NAME}"
	mkdir -p ${RESULTS_PATH}
    	matlab -nodisplay -r "cd('./'); edgesMap('${DS_PATH}', '${RESULTS_PATH}'); exit;"
    	
    done
    
    # SALT AND PEPPER
    for NOISE in ${NOISE_SP[@]}; do
    
    	DS_PATH="../../datasets/${DS_NAME}/salt_pepper/${NOISE}"
    	RESULTS_PATH="../../${METHOD_NAME}/salt_pepper/${NOISE}/${DS_NAME}"
	mkdir -p ${RESULTS_PATH}
    	matlab -nodisplay -r "cd('./'); edgesMap('${DS_PATH}', '${RESULTS_PATH}'); exit;"
    	
    done

done


