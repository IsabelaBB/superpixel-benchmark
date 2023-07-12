#! /bin/bash

if [ "$#" -lt 6 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	SALIENCY_PATH=$3
	RESULTS_PATH=$4
	DTIME_FILE=$5
	TIME_FILE=$6
	
	PARAMS="--img ${DS_PATH} "
	PARAMS+="--nf ${superpixels} "
	PARAMS+="--objsm ${SALIENCY_PATH} "
	PARAMS+="--boost 1 "
	
	if [[ "$RESULTS_PATH" != "-1" ]]; then
	    PARAMS+="--out ${RESULTS_PATH} "
	fi
	
	if [[ "$DTIME_FILE" != "-1" ]]; then
	    PARAMS+="--dtime ${DTIME_FILE} "
	fi
	
	if [[ "$TIME_FILE" != "-1" ]]; then
	    PARAMS+="--time ${TIME_FILE} "
	fi
	
	./bin/RunSICLE ${PARAMS}
	
fi

