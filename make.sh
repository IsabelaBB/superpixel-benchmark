#! /bin/bash

declare -a METHODS_FOLDER=("AINET" "CRS" "DAL-HERS" "DISF" "DRW" "ERGC" "ERS" "ETPS" "GMMSP" "GRID" "IBIS" "ISF" "LNSNet" "LSC" "ODISF" "RSS" "SCALP" "SEEDS" "SH" "SICLE" "SIN" "SLIC" "SNIC" "SSFCN")

cd methods
for METHOD_NAME in ${METHODS_FOLDER[@]}; do
	echo "### ${METHOD_NAME} ###"
	cd ${METHOD_NAME}; make clean; 
	make; 
	cd ..
done
cd ..

cd evaluation; make clean; make; 
cd ..

