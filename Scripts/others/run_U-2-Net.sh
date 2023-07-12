#! /bin/bash

declare -a DATASETS=("Sky")

DATASET_DIR="../../datasets/"
SALIENCYMAP_CODE="../../others/U-2-Net/"
SALIENCYMAP_RESULT="../../others/SaliencyMaps/"

MODEL="${SALIENCYMAP_CODE}/saved_models/u2net/u2net"
IMG_FOLDER="images"
DEVICE="cuda"

################################################################

for DS_NAME in ${DATASETS[@]}; do

DS_PATH="${DATASET_DIR}/${DS_NAME}/${IMG_FOLDER}"
SALIENCYMAP_RESULT_PATH="${SALIENCYMAP_RESULT}/${DS_NAME}"
    
mkdir -p ${SALIENCYMAP_RESULT_PATH}

python3 ${SALIENCYMAP_CODE}/u2net_test.py --img ${DS_PATH} --out ${SALIENCYMAP_RESULT_PATH} --model ${MODEL} --device ${DEVICE}

done
