#! /bin/bash

declare -a DATASETS=("Sky")

DATASET_DIR="../../datasets/"
EDGEMAP_CODE="../../others/edgeDetection/"
EDGEMAP_RESULT="../../others/EdgeMaps/"
IMG_FOLDER="images"

################################################################

for DS_NAME in ${DATASETS[@]}; do

DS_PATH="${DATASET_DIR}/${DS_NAME}/${IMG_FOLDER}"
EDGEMAP_RESULT_PATH="${EDGEMAP_RESULT}/${DS_NAME}"
    
mkdir -p ${EDGEMAP_RESULT_PATH}

matlab -nodisplay -r "cd('${EDGEMAP_CODE}'); edgesMap('${DS_PATH}', '${EDGEMAP_RESULT_PATH}'); exit;"

done


