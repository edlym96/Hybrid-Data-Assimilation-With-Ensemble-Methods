#!/bin/bash
COUNTER=10
while [ $COUNTER  -le 100 ];
do
	filepath="../data/matrix_prec_494/matrixVensembleSplit$COUNTER.npz"
	echo Building solution for ensemble size of $COUNTER
	./VarDA-3Dtracers+Covariance.py -Vp $filepath
	let COUNTER=COUNTER+10
done
echo All done!
