#! /bin/bash

if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	RESULTS_PATH=$3
	DTIME_FILE=$4
	TIME_FILE=$5
	
	BINS=5
	CONFIDENCE=0.1
	PRIOR=1
	MEANS=1
	ITERATIONS=2
	COLOR=1
	
	PARAMS="-i ${DS_PATH} "
	PARAMS+="-s ${superpixels} "
	PARAMS+="-b ${BINS} "
	PARAMS+="-c ${CONFIDENCE} "
	PARAMS+="-p ${PRIOR} "
	PARAMS+="-m ${MEANS} "
	PARAMS+="-t ${ITERATIONS} "
	PARAMS+="-r ${COLOR} "
	
	if [[ "$RESULTS_PATH" != "-1" ]]; then
	    PARAMS+="-o ${RESULTS_PATH} "
	fi
	
	if [[ "$DTIME_FILE" != "-1" ]]; then
	    PARAMS+="--dtime ${DTIME_FILE} "
	fi
	
	if [[ "$TIME_FILE" != "-1" ]]; then
	    PARAMS+="--time ${TIME_FILE} "
	fi
	
	./bin/seeds_cli ${PARAMS}
fi

