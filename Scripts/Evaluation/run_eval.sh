#! /bin/bash

####################################
####### MANDATORY PARAMETERS

declare -a DATASETS=("Sky")
declare -a SUPERPIXELS=("200")
declare -a METHODS=("CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "ETPS" "GMMSP" "GRID" "IBIS" "ISF" "LNSNet" "LSC" "ODISF" "RSS" "SCALP" "SICLE" "SEEDS" "SH" "SLIC" "SNIC")
declare -a MEASURES=("BR" "UE" "SIRS" "EV" "CO" "Connectivity" "SuperpixelNumber")

SCRIPTS_PATH="Scripts/Evaluation"
DATASET_DIR="../../datasets/"
RESULTS_DIR="../../RESULTS/"

IMG_FOLDER="images"
GT_FOLDER="gt"

SEGM_RESULTS_DIR="${RESULTS_DIR}/Segm"
EVAL_RESULTS_DIR="${RESULTS_DIR}/Eval"

####################################
####### OPTIONAL PARAMETERS
SAVE_SCORES=0				# txt file with evaluated results (for each image)
SAVE_OVERALL_SCORES=0			# txt file with overall results
SAVE_LABELS=0 				# Save image superpixels after enforce coonectivity/minimum number of superpixels. Used in metrics 6 and 7.

####################################
####### OPTIONAL SIRS/EV PARAMETERS
SAVE_SCORE_IMAGE=1			# Save labeled images whose color map indicates SIRS/EV scores
DRAW_SCORES=1				# Boolean option when using SAVE_SCORE_IMAGE to show score values
SAVE_RECONSTRUCTED_IMAGE=0 		# Save the reconstructed image

SCORES_IMAGE_PATH="${RESULTS_DIR}/ScoreImages"
RECONSTRUCTED_IMAGE_PATH="${RESULTS_DIR}/ReconstructedImages"

####################################

runMetric () {
	PARAMS+="--label ${SEGM_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${superpixels} "
	
	if [ "$SAVE_SCORES" -eq 1 ] || [ "$SAVE_OVERALL_SCORES" -eq 1 ]; then
	    SCORES_PATH="${EVAL_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}"
	    mkdir -p ${SCORES_PATH}
	fi
	
	if [ "$SAVE_SCORES" -eq 1 ]; then
	    PARAMS+="--dlog ${SCORES_PATH}/${METHOD_NAME}-${superpixels}.txt "
	fi
	
	if [ "$SAVE_OVERALL_SCORES" -eq 1 ]; then
	    PARAMS+="--log ${SCORES_PATH}/${METHOD_NAME}-${DS_NAME}-${MEASURE_NAME}.txt "
	fi
    
    	echo "../../evaluation/bin/main ${PARAMS}"
    	../../evaluation/bin/main ${PARAMS}
}

# GET THE LABELED IMAGE EXTENSION
getExt () {
	declare -a PGM_METHODS=("CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "GMMSP" "GRID" "IBIS" "ISF" "LSC" "ODISF" "RSS" "SEEDS" "SH" "SLIC" "SICLE")
	declare -a PNG_METHODS=("ETPS" "LNSNet" "SCALP" "SNIC")
	
	if [[ " ${PGM_METHODS[@]} " =~ " ${METHOD_NAME} " ]]; then
		EXT="pgm"
	elif [[ " ${PNG_METHODS[@]} " =~ " ${METHOD_NAME} " ]]; then
		EXT="png"
	else
		echo "Method ${METHOD_NAME} not founded. Set default extension pgm. \n"
		EXT="pgm"
	fi
}

# SET PARAMETERS ACCORDING TO EACH MEASURE
getParams() {
	if [ "$MEASURE_NAME" = "SIRS" ]; then
	    MEASURE_NUM=1
	fi
	
	if [ "$MEASURE_NAME" = "EV" ]; then
	    MEASURE_NUM=2
	fi
	
	if [ "$MEASURE_NAME" = "BR" ]; then
	    MEASURE_NUM=3
	fi
	
	if [ "$MEASURE_NAME" = "UE" ]; then
	    MEASURE_NUM=4
	fi
	
	if [ "$MEASURE_NAME" = "CO" ]; then
	    MEASURE_NUM=5
	fi
	
	if [ "$MEASURE_NAME" = "Connectivity" ]; then
	    MEASURE_NUM=6
	fi
	
	if [ "$MEASURE_NAME" = "SuperpixelNumber" ]; then
	    MEASURE_NUM=7
	    PARAMS+="--k ${superpixels} "
	fi
	
	PARAMS+="--metric ${MEASURE_NUM} "
	
	if [ "$MEASURE_NAME" = "SIRS" ] || [ "$MEASURE_NAME" = "EV" ] || [ "$MEASURE_NAME" = "SuperpixelNumber" ]; then
		PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${IMG_FOLDER} "
	else
		PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${GT_FOLDER} "
	fi
	
	if ([ "$MEASURE_NAME" = "Connectivity" ] || [ "$MEASURE_NAME" = "SuperpixelNumber" ]) && [ ${SAVE_RECONSTRUCTED_IMAGE} -eq 1 ]; then
		PARAMS+="--save ${RESULTS_DIR}/${MEASURE_NAME}_Images/${DS_NAME}/${superpixels} "
		mkdir -p "${RESULTS_DIR}/${MEASURE_NAME}_Images/${DS_NAME}/${superpixels}"
	fi
	
	if [ ${SAVE_SCORE_IMAGE} -eq 1 ] && ( [ "$MEASURE_NAME" = "SIRS" ] || [ "$MEASURE_NAME" = "EV" ] ); then
	    PARAMS+="--imgScores ${SCORES_IMAGE_PATH}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels} "
	    mkdir -p ${SCORES_IMAGE_PATH}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels}
	    PARAMS+="--drawScores ${DRAW_SCORES} "
	fi 
	    
	if [ ${SAVE_RECONSTRUCTED_IMAGE} -eq 1 ] && ( [ "$MEASURE_NAME" = "SIRS" ] || [ "$MEASURE_NAME" = "EV" ] );then
		PARAMS+="--recon ${RECONSTRUCTED_IMAGE_PATH}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels} "
		mkdir -p ${RECONSTRUCTED_IMAGE_PATH}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels}
	fi
}


for DS_NAME in ${DATASETS[@]}; do
	for superpixels in ${SUPERPIXELS[@]}; do
		for METHOD_NAME in ${METHODS[@]}; do
			getExt
			for MEASURE_NAME in ${MEASURES[@]}; do
				PARAMS="--ext ${EXT} "
				getParams
				runMetric
			done
		done
	done
done

