#! /bin/bash

####################################
####### MANDATORY PARAMETERS

declare -a DATASETS=("Sky")
declare -a SUPERPIXELS=("200")
declare -a MEASURES=("BR" "UE" "SIRS" "EV" "Connectivity")
declare -a METHODS=("CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "ETPS" "GMMSP" "GRID" "IBIS" "ISF" "LNSNet" "LSC" "ODISF" "RSS" "SCALP" "SICLE" "SEEDS" "SH" "SLIC" "SNIC")

ROBUSTNESS_FOLDER="Robustness"
declare -a NOISE_BLUR=("5" "9" "13" "17")
BLUR_FOLDER="avgblur"
declare -a NOISE_SP=("0.04" "0.08" "0.12" "0.16")
SP_FOLDER="salt_pepper"

SCRIPTS_PATH="Scripts/Evaluation"
DATASET_DIR="../../datasets/"
RESULTS_DIR="../../RESULTS/"

GT_FOLDER="gt"

SEGM_RESULTS_DIR="${RESULTS_DIR}/Segm/${ROBUSTNESS_FOLDER}"
EVAL_RESULTS_DIR="${RESULTS_DIR}/Eval"

####################################
####### OPTIONAL PARAMETERS
SAVE_SCORES=1				# txt file with evaluated results (for each image)
SAVE_OVERALL_SCORES=1			# txt file with overall results
SAVE_LABELS=0 				# Save image superpixels after enforce coonectivity/minimum number of superpixels. Used in metrics 6 and 7.

####################################
####### OPTIONAL SIRS/EV PARAMETERS
SAVE_SCORE_IMAGE=0			# Save labeled images whose color map indicates SIRS/EV scores
DRAW_SCORES=0				# Boolean option when using SAVE_SCORE_IMAGE to show score values
SAVE_RECONSTRUCTED_IMAGE=0 		# Save the reconstructed image

SCORES_IMAGE_PATH="${RESULTS_DIR}/ScoreImages/${ROBUSTNESS_FOLDER}"
RECONSTRUCTED_IMAGE_PATH="${RESULTS_DIR}/ReconstructedImages/${ROBUSTNESS_FOLDER}"

####################################

runMetric () {
	PARAMS+="--label ${SEGM_RESULTS_DIR}/${NOISE_FOLDER}/${NOISE}/${METHOD_NAME}/${DS_NAME}/${superpixels} "
	
	if [ "$SAVE_SCORES" -eq 1 ] || [ "$SAVE_OVERALL_SCORES" -eq 1 ]; then
	    SCORES_PATH="${EVAL_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${ROBUSTNESS_FOLDER}/${NOISE_FOLDER}/${NOISE}/${MEASURE_NAME}"
	    mkdir -p ${SCORES_PATH}
	fi
	
	if [ "$SAVE_SCORES" -eq 1 ]; then
	    PARAMS+="--dlog ${SCORES_PATH}/${METHOD_NAME}-${superpixels}.txt "
	fi
	
	if [ "$SAVE_OVERALL_SCORES" -eq 1 ]; then
	    PARAMS+="--log ${SCORES_PATH}/${METHOD_NAME}-${DS_NAME}-${MEASURE_NAME}.txt "
	fi
    
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
		PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${NOISE_FOLDER}/${NOISE} "
	else
		PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${GT_FOLDER} "
	fi
	
	if ([ "$MEASURE_NAME" = "Connectivity" ] || [ "$MEASURE_NAME" = "SuperpixelNumber" ]) && [ ${SAVE_RECONSTRUCTED_IMAGE} -eq 1 ]; then
		PARAMS+="--save ${RESULTS_DIR}/${NOISE_FOLDER}/${NOISE}/${MEASURE_NAME}_Images/${DS_NAME}/${superpixels} "
		mkdir -p "${RESULTS_DIR}/${NOISE_FOLDER}/${NOISE}/${MEASURE_NAME}_Images/${DS_NAME}/${superpixels}"
	fi
	
	if [ ${SAVE_SCORE_IMAGE} -eq 1 ] && ( [ "$MEASURE_NAME" = "SIRS" ] || [ "$MEASURE_NAME" = "EV" ] ); then
	    PARAMS+="--imgScores ${SCORES_IMAGE_PATH}/${NOISE_FOLDER}/${NOISE}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels} "
	    mkdir -p ${SCORES_IMAGE_PATH}/${NOISE_FOLDER}/${NOISE}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels}
	    PARAMS+="--drawScores ${DRAW_SCORES} "
	fi 
	    
	if [ ${SAVE_RECONSTRUCTED_IMAGE} -eq 1 ] && ( [ "$MEASURE_NAME" = "SIRS" ] || [ "$MEASURE_NAME" = "EV" ] );then
		PARAMS+="--recon ${RECONSTRUCTED_IMAGE_PATH}/${NOISE_FOLDER}/${NOISE}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels} "
		mkdir -p ${RECONSTRUCTED_IMAGE_PATH}/${NOISE_FOLDER}/${NOISE}/${METHOD_NAME}/${DS_NAME}/${MEASURE_NAME}/${superpixels}
	fi
}


for DS_NAME in ${DATASETS[@]}; do
	for superpixels in ${SUPERPIXELS[@]}; do

		for NOISE in ${NOISE_BLUR[@]}; do
			NOISE_FOLDER="${BLUR_FOLDER}"
		    	
		    	for METHOD_NAME in ${METHODS[@]}; do
				getExt
				for MEASURE_NAME in ${MEASURES[@]}; do
					PARAMS="--ext ${EXT} "
					getParams
					runMetric
				done
			done
		done
		
		for NOISE in ${NOISE_SP[@]}; do
			NOISE_FOLDER="${SP_FOLDER}"
		    	
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
done

