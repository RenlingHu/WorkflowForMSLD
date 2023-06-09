# Multiple Topology Creation

### This folder contains Python scripts to automate Multiple Topology Model Creation and MSLD Simulation Input File Preparation

05/04/2023

Author: Renling Hu

E-mails: 22260237@zju.edu.cn

#### Note
 - origin from: https://github.com/Vilseck-Lab/msld-py-prep
 - mainly modified the msld_mcs module which combined with the LIC_align algorithm to identify the MCS more accurately
 - secondly modified the msld_wrt module to connect to the next 'Simulation System Setup' stage more flawlessly and robustly
 
---
#### INPUT
 - structure file(mol2) for each ligand, **make sure** the atom names within each molecule are unique (atom names across ligands may be reused)
 - toppar/param files(str or rtf&prm) for each ligand, **make sure** to name these files the same as the corresponding .mol2 file. 
 - "mol_list.txt" with mol2 filenames (without suffix)

#### USAGE
`./run`

#### OUTPUT
1. **After the first calling**

`./MCS_for_MSLD.txt`
 - each molecule is divided up into core(anchor) + fragment
 - one could manually change the selection of the core and the setting of the substituent site to explore perturbation patterns of interest
    - _NSUBS_ contains information about how many different substituents are at each site.
    - _REFLIG_ includes the name of the ligand that is used as a reference to match every other ligand’s atoms to. This is not to be confused with the reference ligand used for the subsequent free energy calculations, although they could be the same.
    - _CORE_ section includes a list of all of the ligands and the corresponding atom names that are in the core. The matched atoms are spatially aligned across ligands, meaning that each column corresponds to the atom names of all the ligands whose atoms are a match (or the ‘same’).
    - _ANCHOR_ atom section corresponds to the atoms that connect the core to the different sites. Therefore, the first column after the ligand name corresponds to the anchor atoms for each ligand in the first site, the second column to those in the second site, and so on. If more than two atoms serve as anchors to a specific site, the variable _DUM_ will be used as a placeholder (meaning dummy atom). The anchor atoms are the first atoms that connect to each fragment in a particular site and will become part of the fragment after the charge renormalization step. These atoms will appear in the _CORE_ section as well.   
    - _SITE N FRAGMENTS_ section (where N corresponds to the site number) includes the unique fragments that are identified at each site across all of the queried ligands.

2. **After the second calling**

`./build.sysname`
 - contains Input Files for the next 'Simulation System Setup' and 'MSLD Simulations'
 - see more details in 'Simulation System Setup' and 'MSLD Simulation'

`./translation.txt`
 - New atom name assigned to multiple ligands that contain one common core and different fragments

---

#### File information
`./run`
 - scripts for submitting tasks containing (1) mol_list write (2) Multiple Topology Model Creation and MSLD Simulation Input File Preparation 
 - *manually calling* msld_py_prep.py twice is recommended
 - after first calling, please check the identified MCS either through MCS_for_MSLD.txt or vis_check.py

`./msld_py_prep.py`
 - python, numpy and pandas are required
 - contains msld_chk msld_mcs msld_crn msld_prm msld_wrt module

`./vis_check.py`
 - Usage: open Pymol, then run this script in the command line `run vis_check.py`

`./Lg_Solvate.sh`
 - convpdb.pl from the MMTSB toolset (https://github.com/mmtsb/toolset) is required
 - up-to-date CHARMM toppar files in `./tappar/` are required
 - Usage: `./Lg_Solvate.sh [cutoff_distance]`
 - Not used in workflow for now

`./rename_atoms.py`
 - make the atom names within each molecule unique
