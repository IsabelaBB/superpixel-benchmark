#! /bin/bash

if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
else
	DS_PATH=$1
	superpixels=$2
	RESULTS_PATH=$3
	DTIME_FILE=$4
	TIME_FILE=$5
	
	octave --path ./ --eval "LSC('${DS_PATH}', '${RESULTS_PATH}', ${superpixels}, '${DTIME_FILE}', '${TIME_FILE}');"
	
fi

