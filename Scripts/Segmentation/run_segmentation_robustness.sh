#! /bin/bash

declare -a DATASETS=("Sky")
declare -a SUPERPIXELS=("200")
declare -a METHODS=("CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "ETPS" "GMMSP" "GRID" "IBIS" "ISF" "LNSNet" "LSC" "ODISF" "RSS" "SCALP" "SICLE" "SEEDS" "SH" "SLIC" "SNIC")

ROBUSTNESS_FOLDER="Robustness"
declare -a NOISE_BLUR=("5" "9" "13" "17")
BLUR_FOLDER="avgblur"
declare -a NOISE_SP=("0.04" "0.08" "0.12" "0.16")
SP_FOLDER="salt_pepper"

SCRIPTS_PATH="Scripts/Segmentation/Methods"
DATASET_DIR="../../datasets/"
RESULTS_DIR="../../RESULTS/"
SALIENCYMAP_DIR="../../others/SaliencyMaps/"
EDGEMAP_DIR="../../others/EdgeMaps/"

SEGMENTATION_DIR="${RESULTS_DIR}/Segm/${ROBUSTNESS_FOLDER}"
EVALUATION_DIR="${RESULTS_DIR}/Eval"

################################################################

cd ../../methods/

runMethod () {
	if [ "$METHOD_NAME" == 'CRS' ] || [ "$METHOD_NAME" = 'ERGC' ] || [ "$METHOD_NAME" = 'SEEDS' ]; then
		cd stutz_bench
	else
		cd $METHOD_NAME
	fi
	
	PARAMS="${DATASET_DIR}/${DS_NAME}/${NOISE_FOLDER}/${NOISE} "
	PARAMS+="${superpixels} "
	
	if [ "$METHOD_NAME" = 'ODISF' ] || [ "$METHOD_NAME" = 'SICLE' ]; then
		PARAMS+="${SALIENCYMAP_DIR}/${NOISE_FOLDER}/${NOISE}/${DS_NAME}/u2net_results/ "
	elif [ "$METHOD_NAME" = 'SCALP' ]; then
		PARAMS+="${EDGEMAP_DIR}/${NOISE_FOLDER}/${NOISE}/${DS_NAME}/ "
	fi
	
	RESULTS_PATH="${SEGMENTATION_DIR}/${NOISE_FOLDER}/${NOISE}/${METHOD_NAME}/${DS_NAME}/${superpixels}"
	TIME_PATH="${EVALUATION_DIR}/${METHOD_NAME}/${DS_NAME}/${ROBUSTNESS_FOLDER}/${NOISE_FOLDER}/${NOISE}/time"

	mkdir -p "${RESULTS_PATH}"
	mkdir -p "${TIME_PATH}"
	
	PARAMS+="${RESULTS_PATH} "
	PARAMS+="${TIME_PATH}/${METHOD_NAME}-${superpixels}.txt "
	PARAMS+="${TIME_PATH}/${METHOD_NAME}-${DS_NAME}-time.txt "
	
	bash ../../${SCRIPTS_PATH}/run_${METHOD_NAME}.sh ${PARAMS}
	cd ..
}

for DS_NAME in ${DATASETS[@]}; do
	for superpixels in ${SUPERPIXELS[@]}; do

		for NOISE in ${NOISE_BLUR[@]}; do
			NOISE_FOLDER="${BLUR_FOLDER}"
		    	
		    	for METHOD_NAME in ${METHODS[@]}; do
		
				echo  "############### ${METHOD_NAME} ${superpixels} ${DS_NAME} ###############"
				runMethod
			done
		done

		for NOISE in ${NOISE_SP[@]}; do
			NOISE_FOLDER="${SP_FOLDER}"
		    	
		    	for METHOD_NAME in ${METHODS[@]}; do
				echo  "############### ${METHOD_NAME} ${superpixels} ${DS_NAME} ###############"
				runMethod
			done
		done
	done
done
cd ..
