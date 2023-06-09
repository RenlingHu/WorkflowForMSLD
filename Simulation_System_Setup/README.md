# Simulation System Setup

### This folder contains CHARMM scripts and relative input files to automate Simulation System Setup using CHARMM

05/04/2023

Author: Renling Hu

E-mails: 22260237@zju.edu.cn

#### Note
 - this folder is itself an example, so there is no ./example

---
#### INPUT
1. **CHARMM scripts**
 - `./syaname.inp` and `./sysname-lig.inp`
 - generated in the "Multiple Topology Creation" stage, right in the build.sysname folder
   
 **Note**:
 - modify the "Reads coordinates and pdb file for the protein" section according to the starting and ending residues of each subunit of the protein
 - modify the method and steps in the "minimization" section according to the requirements
 - modify or add any other parameters if needed

2. **relative input files**
 - `./prep` contains: 
    (1) protein.pdb, solv.pdb, ions.pdb generated in the "Protein Preparation" stage
    (2) files that define the core (anchor) and alchemy perturbation R groups generated in the "Multiple Topoloy Moder Creation" stage, right in the build.sysname folder

#### Usage
`charmm -i sysname.inp -o complex.out`

#### OUTPUT
 - complex_pbc.psf crd pdb file
 - complex_mini.psf crd pdb file
 - ligand_pbc.psf crd pdb file
 - ligand_mini.psf crd pdb file
