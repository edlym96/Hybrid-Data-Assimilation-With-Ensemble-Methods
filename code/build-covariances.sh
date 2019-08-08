#!/bin/bash

COUNTER=10

while [	$COUNTER -le 100 ];
do
	echo Calculating for ens size of $COUNTER ...
	./optimal_covariance.py -ens --ens_size $COUNTER
	let COUNTER=COUNTER+10
done

echo All done
