#! /bin/bash

if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	RESULTS_PATH=$3
	DTIME_FILE=$4
	TIME_FILE=$5
	COMPACTNESS=10
	
	PARAMS="--img ${DS_PATH} "
	PARAMS+="--k ${superpixels} "
	PARAMS+="--cmp ${COMPACTNESS} "
	
	if [[ "$RESULTS_PATH" != "-1" ]]; then
	    PARAMS+="--label ${RESULTS_PATH} "
	fi
	
	if [[ "$DTIME_FILE" != "-1" ]]; then
	    PARAMS+="--dtime ${DTIME_FILE} "
	fi
	
	if [[ "$TIME_FILE" != "-1" ]]; then
	    PARAMS+="--time ${TIME_FILE} "
	fi
	
	python3 ./snic.py ${PARAMS}
fi

