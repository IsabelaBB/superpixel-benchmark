#! /bin/bash

#! /bin/bash

declare -a DATASETS=("Birds" "Sky" "ECSSD" "Insects")

METHOD_NAME=EdgeMaps

for DS_NAME in ${DATASETS[@]}; do

DS_PATH="../../datasets/${DS_NAME}/images"

RESULTS_PATH="../../${METHOD_NAME}/${DS_NAME}"
    
mkdir -p ${RESULTS_PATH}

echo "${DS_NAME}"

#for FILE_FULLNAME in "${DS_PATH}/"*; do

#IMAGE_NAME="$(basename -- ${FILE_FULLNAME} | cut -f 1 -d '.')" # image name only

matlab -nodisplay -r "cd('./'); edgesMap('${DS_PATH}', '${RESULTS_PATH}'); exit;"

done


