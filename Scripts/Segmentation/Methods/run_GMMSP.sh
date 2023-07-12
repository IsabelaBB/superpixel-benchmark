#! /bin/bash

if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	RESULTS_PATH=$3
	DTIME_FILE=$4
	TIME_FILE=$5
	
	matlab -nodisplay -r "cd('./'); GMMSP('${DS_PATH}', ${superpixels}, '${RESULTS_PATH}', '${DTIME_FILE}', '${TIME_FILE}'); exit;"
fi
