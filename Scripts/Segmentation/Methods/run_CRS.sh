#! /bin/bash

if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	RESULTS_PATH=$3
	DTIME_FILE=$4
	TIME_FILE=$5
	
	COMPACTNESS=0.045
	CLIQUE=0.3
	ITERATIONS=3
	COLOR=0
	
	PARAMS="-i ${DS_PATH} "
	PARAMS+="-s ${superpixels} "
	
	if [[ "$RESULTS_PATH" != "-1" ]]; then
	    PARAMS+="-o ${RESULTS_PATH} "
	fi
	
	if [[ "$DTIME_FILE" != "-1" ]]; then
	    PARAMS+="--dtime ${DTIME_FILE} "
	fi
	
	if [[ "$TIME_FILE" != "-1" ]]; then
	    PARAMS+="--time ${TIME_FILE} "
	fi
	
	./bin/crs_cli ${PARAMS}
fi


