#!/bin/bash
COUNTER=10
while [ $COUNTER  -le 100 ];
do
	filepath=../data/matrix_prec_494/matrixVensembleSplit"$COUNTER"state.npz
	echo Building solution for ensemble size of $COUNTER
	./VarDA_3Dtracers_Covariance.py -Vp $filepath -xBp ../data/converted_data/background_state.npz -yp ../data/converted_data/observations.npz
	let COUNTER=COUNTER+10
done
echo All done!
