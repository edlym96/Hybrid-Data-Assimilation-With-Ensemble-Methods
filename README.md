# MSc-Individual-Project---Hybrid-Data-Assimilation-With-Ensemble-Methods

The code is separated into 6 files,

-   `convert_vtu.py` - Extracts the relevant fields and converts data
    from .vtu format to .npz for easier access

-   `optimal_covariance.py` - Builds the background error covariance
    CVT, $\textbf{V}$

-   `localisation.py` - Builds the localisation matrices,
    $\textbf{C}_{h}$ and $\textbf{C}_{v}$

-   `VarDA_3Dtracers_Covariance.py` - Carries out the Data Assimilation

-   `evaluate_DA_solution.py` - Auxiliary file for evaluating output of
    DA

-   `vtktools.py` - Dependency for reading `.vtu` files
\
The files should be ran in this sequence with the exception of
`vtktools.py`. The instructions for running each file and their
arguments are as follows.\
-	`convert_vtu.py`
	-   `filepath`: provide the filepath for the folder where the data is
		stored. The data files must follow the format, `LSBU_i.vtu` where
		$i$ is the timestep index

	-   `ntime`: Number of timesteps for the data (default = $988$)\

	-   **Example Usage**:

	-   `./convert_vtu.py --filepath $FILEPATH`
\
-   `optimal_covariance.py`

    -   `filepath`: provide the filepath where the background state data
        is stored (`.npz` format only)

    -   `tsvd`: Flag to perform TSVD method for building error
        covariance matrix

    -   `ens`: Flag to perform ensemble method for building error
        covariance matrix

    -   `trnc`: Provide truncation parameter for TSVD (default = $145$)

    -   `ens_size`: Provide ensemble size for ensemble (default = $50$)

    -   `ntime`: Number of timesteps of DA (default = $494$)\

    -   **Example Usage**:

    -   `./optimal_covariance.py --filepath $FILEPATH -tsvd -ens`
\
-   `localisation.py`

    -   `posp`: provide the filepath where the position data is stored
        (`.npz` format only)

    -   `rh`: Provide the number of dominant EOFs to keep for
        eigendecomposition for horizontal localisation (default = $3$)

    -   `rv`: Provide the number of dominant EOFs to keep for
        eigendecomposition for vertical localisation (default = $3$)

    -   `h_cutoff`: Provide the cutoff value for horizontal localisation
        (default = $200$)\

    -   **Example Usage**:

    -   `./localisation -posp $POS_FILEPATH`
\
-   `VarDA_3Dtracers_Covariance.py`

    -   `xBp`: provide the filepath where the background state data is
        stored (`.npz` format only)

    -   `yp`: provide the filepath where the observation data is stored
        (`.npz` format only)

    -   `Vp`: provide the filepath where the error covariance is stored
        (`.npz` format only)

    -   `hlocal`: provide the filepath where the horizontal localisation
        matrix is stored (default = $None$)

    -   `vlocal`: provide the filepath where the vertical localisation
        matrix is stored (default = $None$)

    -   `ntime`: Number of timesteps of DA (default = $494$)\

    -   **Example Usage**:

    -   `./VarDA_3Dtracers_Covariance.py -xBp $XBP_FILEPATH -yp $YP_FILEPATH`

    -   `-Vp $VP_FILEPATH -hlocal $HLOCAL_FILEPATH -vlocal $VLOCAL_FILEPATH`
\
-   `evaluate_DA_solution.py`

    -   `results`: provide the filepath where the results of DA are
        stored (`.npz` or `.vtu` format only)

    -   `xBp`: provide the filepath where the background state data is
        stored (`.npz` format only)

    -   `yp`: provide the filepath where the observation data is stored
        (`.npz` format only)\

    -   **Example Usage**:

    -   `./evaluate_DA_solution.py -results $RESULTS_FILEPATH -xBp`

    -   `$XBP_FILEPATH -yp $YP_FILEPATH`
	
Additionally, a bash script `run-DA.sh` is available that runs all the
files in the correct order to achieve the same results as presented in
this paper. The default filepath for the folder storing the data in this
setup is `../data/small3DLSBU/` where all the `.vtu` files within this
folder follow the naming convention `LSBU_i.vtu` where $i$ is the time
index.\
\
Dependency list:

-   `numpy`

-   `scipy`

-   `vtk`