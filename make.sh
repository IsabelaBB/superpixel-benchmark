#! /bin/bash

declare -a METHODS_FOLDER=("SEEDS" "CRS" "ERGC" "DAL-HERS" "DISF" "DRW" "ERS" "ETPS" "GMMSP" "GRID" "IBIS" "ISF" "LSC" "ODISF" "RSS" "SCALP" "SH" "SICLE" "SLIC")

cd methods
for METHOD_NAME in ${METHODS_FOLDER[@]}; do
	echo "### ${METHOD_NAME} ###"
	cd ${METHOD_NAME}; make clean; #make; 
	cd ..
done
cd ..

cd evaluation; make clean; #make; 
cd ..

