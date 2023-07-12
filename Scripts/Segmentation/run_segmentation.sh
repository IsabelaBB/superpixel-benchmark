#! /bin/bash

declare -a DATASETS=("Sky")
declare -a SUPERPIXELS=("200")
declare -a METHODS=("CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "ETPS" "GMMSP" "GRID" "IBIS" "ISF" "LNSNet" "LSC" "ODISF" "RSS" "SCALP" "SICLE" "SEEDS" "SH" "SLIC" "SNIC")

SCRIPTS_PATH="Scripts/Segmentation/Methods"
DATASET_DIR="../../datasets/"
RESULTS_DIR="../../RESULTS/"
SALIENCYMAP_DIR="../../others/SaliencyMaps/"
EDGEMAP_DIR="../../others/EdgeMaps/"

SEGMENTATION_DIR="${RESULTS_DIR}/Segm"
EVALUATION_DIR="${RESULTS_DIR}/Eval"

################################################################

cd ../../methods/

runMethod () {
	if [ "$METHOD_NAME" == 'CRS' ] || [ "$METHOD_NAME" = 'ERGC' ] || [ "$METHOD_NAME" = 'SEEDS' ]; then
		cd stutz_bench
	else
		cd $METHOD_NAME
	fi
	
	PARAMS="${DATASET_DIR}/${DS_NAME}/images "
	PARAMS+="${superpixels} "
	
	if [ "$METHOD_NAME" = 'ODISF' ] || [ "$METHOD_NAME" = 'SICLE' ]; then
		PARAMS+="${SALIENCYMAP_DIR}/${DS_NAME}/u2net_results/ "
	elif [ "$METHOD_NAME" = 'SCALP' ]; then
		PARAMS+="${EDGEMAP_DIR}/${DS_NAME}/ "
	fi
	
	RESULTS_PATH="${SEGMENTATION_DIR}/${METHOD_NAME}/${DS_NAME}/${superpixels}"
	TIME_PATH="${EVALUATION_DIR}/${METHOD_NAME}/${DS_NAME}/time"

	DTIME_FILE="${TIME_PATH}/${METHOD_NAME}-${superpixels}.txt"
	TIME_FILE="${TIME_PATH}/${METHOD_NAME}-${DS_NAME}-time.txt"

	mkdir -p "${RESULTS_PATH}"
	mkdir -p "${TIME_PATH}"
	
	PARAMS+="${RESULTS_PATH} "
	PARAMS+="${DTIME_FILE} "
	PARAMS+="${TIME_FILE} "
	
	bash ../../${SCRIPTS_PATH}/run_${METHOD_NAME}.sh ${PARAMS}
	cd ..
}

for DS_NAME in ${DATASETS[@]}; do
	for superpixels in ${SUPERPIXELS[@]}; do
		for METHOD_NAME in ${METHODS[@]}; do
			echo  "############### ${METHOD_NAME} ${superpixels} ${DS_NAME} ###############"
			runMethod
		done
	done
done
cd ..
