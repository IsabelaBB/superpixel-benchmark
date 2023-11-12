#! /bin/bash

if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
else

####################################
####### MANDATORY PARAMETERS

# Merge two superpixel segmentations of the same image, drawing a diagonal line between them. 
# Also colorize the boundaries of superpixels with a specific RGB color.

METHOD_NAME=$1 # Name of the method to evaluate
DS_NAME=$2 # Name of the dataset to evaluate
SPX1=$3 # Number of superpixels of the top left image
SPX2=$ # Number of superpixels of the bottom right image
DIST=$5 # x,y distance to draw the divisor line. Use -1 for maximum distance

# RGB color in range [0,1]
if [ "$#" -lt 6 ]; then
COLOR="0.7,0,0"
else
COLOR=$6
fi

SCRIPTS_PATH="Scripts/Evaluation" # Path to this file
DATASET_DIR="../../datasets/" # Path to the datasets
RESULTS_DIR="../../RESULTS/" # Path to the results

IMG_FOLDER="images"		# folder name in DATASET_DIR to the original images

SEGM_RESULTS_DIR="${RESULTS_DIR}/Segm" # folder of the superpixel segmentation images
EVAL_RESULTS_DIR="${RESULTS_DIR}/Qualitative" # folder of the output image(s)

####################################

runMetric () {
	PARAMS+="--eval 8 "
	
	PARAMS+="--img ${DATASET_DIR}/${DS_NAME}/${IMG_FOLDER} "
	PARAMS+="--save ${EVAL_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${SPX1}-${SPX2}/ "
	mkdir -p "${EVAL_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${SPX1}-${SPX2}/"
	
	PARAMS+="--label ${SEGM_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${SPX1} "
	PARAMS+="--label2 ${SEGM_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${SPX2} "
	PARAMS+="--distances ${DIST} "
	PARAMS+="--rgb ${COLOR} "
	
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

    	
getExt
echo "Dataset:${DS_NAME} Method:${METHOD_NAME} Measure: Qualitative "
PARAMS="--ext ${EXT} "
runMetric

for path in "${EVAL_RESULTS_DIR}/${METHOD_NAME}/${DS_NAME}/${SPX1}-${SPX2}/"*; do 
	tmp=${path%'.'*}; 
	fileName=${tmp##*/}
	mv "${path}" "${EVAL_RESULTS_DIR}/${fileName}-${METHOD_NAME}-${DS_NAME}-${SPX1}-${SPX2}.png"
done

fi
