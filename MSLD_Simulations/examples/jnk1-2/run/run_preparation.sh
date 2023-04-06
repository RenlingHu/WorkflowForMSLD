#!/bin/sh

# prepare/overwrite system.inp which 
# defines the perturbation-related parameters and procedures for MSLD simulations
# parameter_complex is set to False for ligand system
cd ../prep/charmminp
python inp_for_BLaDE.py --sysname=jnk1-2 --complex=True 

# recompile WHAM in current path
bash run_compile.sh

# make sure these scripts have executable permissions
cd ../ALF
chmod 777 *

# initialize starting values of the biases, and create the file variables1.inp
cd ../analysis0
python InitVars.py