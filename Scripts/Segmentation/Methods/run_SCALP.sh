#! /bin/bash

if [ "$#" -lt 6 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	EDGEMAP_PATH=$3
	RESULTS_PATH=$4
	DTIME_FILE=$5
	TIME_FILE=$6
	
	PARAMS="--img ${DS_PATH} "
	PARAMS+="--k ${superpixels} "
	PARAMS+="--contour ${EDGEMAP_PATH} "
	
	if [[ "$RESULTS_PATH" != "-1" ]]; then
	    PARAMS+="--label ${RESULTS_PATH} "
	fi
	
	if [[ "$DTIME_FILE" != "-1" ]]; then
	    PARAMS+="--dtime ${DTIME_FILE} "
	fi
	
	if [[ "$TIME_FILE" != "-1" ]]; then
	    PARAMS+="--time ${TIME_FILE} "
	fi
	
	./SCALP ${PARAMS}
fi

