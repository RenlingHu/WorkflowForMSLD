# MSLD Simulations (Adaptive Landscape Flattening & Production Simulations)

### This folder contains scripts and input files required to automately run MSLD Simulations

05/04/2023

Author: Renling Hu

E-mails: 22260237@zju.edu.cn

#### Note
 - origin from: https://github.com/RyanLeeHayes/ALF/tree/master/blade/template
 - these MLSD Simulations use ****BLaDE** instead CHARMM
 - mainly optimized the interface with the previous stage
 - get BLaDE and more details for using BLaDE: https://github.com/RyanLeeHayes/BLaDE

---
#### INPUT
**Note** : this section contains what needs to add or replace these files (or just change some parameters) to the corresponding path according to systems before MSLD Simulations; the remaining part of preparation had already been automated by scripts, which will be covered in the next section 'USAGE'
 - `./prep`
    - crd, pdb and psf file for complex or ligand system in `./prep`,generated in the "Simulation System Setup" stage
    - MSLD multiple topology in `./prep/toppar`, generated in the "Multiple Topology Creation" stage
    - CHARMM inp in `./prep/charmminp`, generated in the "Multiple Topology Creation" stage
 - `./`
    - `name` `nblocks` `ncentral` `nnodes` `nreps` `nsubs`, generated in the "Multiple Topology Creation" stage
    - set pmegridsize in `nbond.str`
    - change the number of runs in phase1-3 & production in `runset2.sh` `runset3.sh` `runset3_2.sh` `runset4.sh` if needed
       - **ini** sets the earliest run that can be used for this purpose
       - **ifi** is the final step for flattening and controls how many runs will be used in the phase
       - **iri** is generally set to "ini" minus 5 in phase2 and phase3(set to 1 in phase1)
 

#### USAGE
1. Preparation
`cd ./run`
`bash preparation.sh`

2. Simulations
`cd ./run`
`bash run_All`
or step-by-step:
`cd ./run`
`bash run_flat.sh`
`bash run_prod.sh`
`bash run_postprocess.sh`
**Note** : step-by-step is recommended, since after the flattening, one could check the landscape and choose the best parameters

#### OUTPUT
 - `./analysis#/` : ./ALF/PlotFreeEnergy5.m could be copied into an analysis directory and run in matlab to visualize the alchemical free energy profiles
 - `./analysis131/Results.txt` : calculated ddG

---
#### File information
 - see more details in https://github.com/RyanLeeHayes/ALF/tree/master/bladelib/template
 - some additional scripts can be seen in the script for introduction and function