#!/bin/bash

steps=${1:-"10"}
data_filepath=${2:-../data/small3DLSBU/}
pos_filepath=${3:-../data/converted_data/positions.npz}

#./convert_vtu.py -fp $data_filepath

#COUNTER=10

#while [	$COUNTER -le $(($steps * 10)) ];
#do
#	echo Calculating the ensemble for ens size of $COUNTER ...
#	./optimal_covariance.py -ens --ens_size $COUNTER -fp ../data/converted_data/background_state.npz
#	./optimal_covariance.py -ens --ens_size $COUNTER -fp ../data/converted_data/background_pressure.npz
#	./optimal_covariance.py -ens --ens_size $COUNTER -fp ../data/converted_data/background_velocity.npz
#	let COUNTER=COUNTER+10
#done

COUNTER=10
while [ $COUNTER  -le $(($steps * 10)) ];
do
	filepath=../data/matrix_prec_494/matrixVensembleSplit"$COUNTER"state.npz
	pressure_filepath=../data/matrix_prec_494/matrixVensembleSplit"$COUNTER"pressure.npz
	velocity_filepath=../data/matrix_prec_494/matrixVensembleSplit"$COUNTER"velocity.npz

    background_fp=../data/converted_data/background_state.npz
	background_pressure_fp=../data/converted_data/background_pressure.npz
	background_velocity_fp=../data/converted_data/background_velocity.npz

	obs_fp=../data/converted_data/observations.npz
	obs_pressure_fp=../data/converted_data/obs_pressure.npz
	obs_velocity_fp=../data/converted_data/obs_velocity.npz

	echo Building solution for ensemble size of $COUNTER
	./VarDA_3Dtracers_Covariance.py -Vp $filepath -posp $pos_filepath -xBp $background_fp -yp $obs_fp -local
	./VarDA_3Dtracers_Covariance.py -Vp $pressure_filepath -posp $pos_filepath -xBp $background_pressure_fp -yp $obs_pressure_fp
	./VarDA_3Dtracers_Covariance.py -Vp $pressure_filepath -posp $pos_filepath -xBp $background_pressure_fp -yp $obs_pressure_fp -local
	./VarDA_3Dtracers_Covariance.py -Vp $velocity_filepath -posp $pos_filepath -xBp $background_velocity_fp -yp $obs_velocity_fp
	./VarDA_3Dtracers_Covariance.py -Vp $velocity_filepath -posp $pos_filepath -xBp $background_velocity_fp -yp $obs_velocity_fp -local
	let COUNTER=COUNTER+10
done

echo All done
