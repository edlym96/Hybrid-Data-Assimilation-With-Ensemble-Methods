#!/bin/bash

steps=${2:-"10"}
data_filepath=${2:-../data/small3DLSBU/}
xB_filepath=${3:-../data/converted_data/background_state.npz}
y_filepath=${4:-../data/converted_data/observations.npz}
pos_filepath=${5:-../data/converted_data/positions.npz}

./convert_vtu.py -fp $data_filepath

COUNTER=10

while [	$COUNTER -le $(($steps * 10)) ];
do
	echo Calculating the ensemble for ens size of $COUNTER ...
	./optimal_covariance.py -ens --ens_size $COUNTER -fp $xB_filepath
	let COUNTER=COUNTER+10
done

COUNTER=10
while [ $COUNTER  -le $(($steps * 10)) ];
do
	filepath=${5:-../data/matrix_prec_494/matrixVensembleSplit$COUNTER.npz}
	echo Building solution for ensemble size of $COUNTER
	./VarDA_3Dtracers_Covariance.py -Vp $filepath -posp $pos_filepath -xBp $xB_filepath -yp $y_filepath
	let COUNTER=COUNTER+10
done

echo All done
