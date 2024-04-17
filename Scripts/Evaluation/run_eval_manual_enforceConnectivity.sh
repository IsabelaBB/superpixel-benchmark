#! /bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
else

####################################
####### MANDATORY PARAMETERS

# Eval the connectivity in a superpixel segmentation.
# Optionally, enforce the connectivity of superpixels and merge them to maintain the number of labels.

METHOD_NAME=$1 # Name of the method to evaluate
DS_NAME=$2 # Name of the dataset to evaluate

# number of superpixels
declare -a SUPERPIXELS=("25" "50" "75" "100" "200" "300" "400" "500" "600" "700" "800" "900" "1000")

SCRIPTS_PATH="Scripts/Evaluation" # Path to this file
DATASET_DIR="../../datasets/" # Path to the datasets
RESULTS_DIR="../../RESULTS/" # Path to the results

IMG_FOLDER="images" # folder name in DATASET_DIR to the original images
GT_FOLDER="gt"			# folder name in DATASET_DIR to the ground truth images

SEGM_RESULTS_DIR="${RESULTS_DIR}/Segm" # folder of the superpixel segmentation images
EVAL_RESULTS_DIR="${RESULTS_DIR}/Eval" # folder of the output file(s)

####################################
####### OPTIONAL PARAMETERS
SAVE_SCORES=1							# txt file with evaluated results (for each image)
SAVE_OVERALL_SCORES=1			# txt file with overall results

SAVE_LABELS=0 						# Save image superpixels after enforce coonectivity/minimum number of superpixels. 
if [ "$#" -gt 3 ]; then
    SAVE_LABELS=$4
fi

####################################
####### OPTIONAL SIRS/EV PARAMETERS
SAVE_SCORE_IMAGE=0						# Save labeled images whose color map indicates SIRS/EV scores
DRAW_SCORES=0									# Boolean option when using SAVE_SCORE_IMAGE to show score values
SAVE_RECONSTRUCTED_IMAGE=0		# Save the reconstructed image

SCORES_IMAGE_PATH="${RESULTS_DIR}/ScoreImages"
RECONSTRUCTED_IMAGE_PATH="${RESULTS_DIR}/ReconstructedImages"

####################################

MEASURE_NAME="SuperpixelNumber"

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
    
    	../../evaluation/bin/main ${PARAMS}
}

# GET THE LABELED IMAGE EXTENSION
getExt () {
	declare -a PGM_METHODS=("CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "GMMSP" "GRID" "IBIS" "ISF" "LSC" "ODISF" "RSS" "SEEDS" "SH" "SLIC" "SICLE")
	declare -a PNG_METHODS=("ETPS" "LNSNet" "SCALP" "SNIC" "AINET" "SIN" "SSFCN")
	
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
	
	PARAMS+="--eval ${MEASURE_NUM} "
	
	if [ "$MEASURE_NAME" = "SIRS" ] || [ "$MEASURE_NAME" = "EV" ] || [ "$MEASURE_NAME" = "SuperpixelNumber" ]; then
		PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${IMG_FOLDER} "
	else
		PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${GT_FOLDER} "
	fi
	
	if ([ "$MEASURE_NAME" = "Connectivity" ] || [ "$MEASURE_NAME" = "SuperpixelNumber" ]) && [ ${SAVE_LABELS} -eq 1 ]; then
		PARAMS+="--save ${SEGM_RESULTS_DIR}/${MEASURE_NAME}_Images/${METHOD_NAME}/${DS_NAME}/${superpixels} "
		mkdir -p "${SEGM_RESULTS_DIR}/${MEASURE_NAME}_Images/${METHOD_NAME}/${DS_NAME}/${superpixels}"
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
    	

getExt
echo "Dataset:${DS_NAME} Method:${METHOD_NAME} Measure: ${MEASURE_NAME} "
for superpixels in ${SUPERPIXELS[@]}; do		
	PARAMS="--ext ${EXT} "
	getParams
	runMetric
done
fi
