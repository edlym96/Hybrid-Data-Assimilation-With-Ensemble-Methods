#!/bin/bash

steps=${1:-"10"}
data_filepath=${2:-../data/small3DLSBU/}
pos_filepath=${3:-../data/converted_data/positions.npz}

./convert_vtu.py -fp $data_filepath
./localisation.py -posp $pos_filepath
./optimal_covariance.py -fp ../data/converted_data/background_state.npz -tsvd
./optimal_covariance.py -fp ../data/converted_data/background_velocity.npz -tsvd

./VarDA_3Dtracers_Covariance.py -Vp ../data/matrix_prec_494/matrixVprec145state.npz -xBp ./data/converted_data/background_state.npz -yp ../data/converted_data/observations.npz
./VarDA_3Dtracers_Covariance.py -Vp ../data/matrix_prec_494/matrixVprec145velocity.npz -xBp ./data/converted_data/background_velocity.npz -yp ../data/converted_data/obs_velocity.npz

COUNTER=10
while [	$COUNTER -le $(($steps * 10)) ];
do

    filepath=../data/matrix_prec_494/matrixVensembleSplit"$COUNTER"state.npz
	velocity_filepath=../data/matrix_prec_494/matrixVensembleSplit"$COUNTER"velocity.npz

    background_fp=../data/converted_data/background_state.npz
	background_velocity_fp=../data/converted_data/background_velocity.npz

	obs_fp=../data/converted_data/observations.npz
	obs_velocity_fp=../data/converted_data/obs_velocity.npz

	echo Calculating the ensemble for ens size of $COUNTER ...
	./optimal_covariance.py -ens --ens_size $COUNTER -fp $background_fp
	./optimal_covariance.py -ens --ens_size $COUNTER -fp $background_velocity_fp

	echo Building solution for ensemble size of $COUNTER
	./VarDA_3Dtracers_Covariance.py -Vp $filepath -posp $pos_filepath -xBp $background_fp -yp $obs_fp
	./VarDA_3Dtracers_Covariance.py -Vp $velocity_filepath -posp $pos_filepath -xBp $background_velocity_fp -yp $obs_velocity_fp

	let COUNTER=COUNTER+10
done

    $hlocal_fp = ../data/converted_data/reduced_localisation_h3_200.npz
    $vlocal_fp = ../data/converted_data/reduced_localisation_v3.npz

    ./VarDA_3Dtracers_Covariance.py -Vp ../data/matrix_prec_494/matrixVensembleSplit40state.npz -xBp $background_fp -yp $obs_fp -hlocal $hlocal_fp
	./VarDA_3Dtracers_Covariance.py -Vp ../data/matrix_prec_494/matrixVensembleSplit40state.npz -xBp $background_fp -yp $obs_fp -hlocal $hlocal_fp -vlocal $vlocal_fp

echo All done
