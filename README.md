# Workflow For MSLD

### Brief Intro.
An open-source workflow for MSLD Simulations based on BLaDE, from "Protein and Ligand Preparation" to "MSLD Simulations"
This workflow not only **incorporated** all open source scripts currently available, but also **added** workflows and scripts for protein and ligand preparation(e.g. MCS align) and **optimized** MCSS algorithm and the interface between stages

![Workflow](https://github.com/RenlingHu/WorkflowForMSLD/blob/main/Workflow.jpg)

![Main modifications(red part)](https://github.com/RenlingHu/WorkflowForMSLD/blob/main/Main_Modifications.jpg)

### Overview
 - `Protein_and_Ligand_Preparation/` : (1) workflow for protein and ligand preparation; (2) Python/Jupter Notebook scripts for aligning ligands according to the MCS
 - `Multiple_Topology_Creation/` : modified-Python scripts for creating multiple topology models of ligands and preparing input files needed in next stages
 - `Simulation_System_Setup/` : CHARMM scripts and input files for setting up simulation systems
 - `MSLD_Simulations/` : (1) Python and Shell scripts for running MSLD Simulations; (2) workflow for adaptive landscape flattening and production runs

### Requirements
 - python 3.9.16
 - jupyternotebook 6.5.3
 - rdkit 2022.9.5
 - scipy 1.10.1
 - numpy 1.24.2
 - pandas 1.5.3
 - Schrödinger 2022-4
 - CHARMM c46b2
 - cuda 11.2
 - gcc 9.4.0
 - matlab R2022b
 - BLaDE 3.1

### Usage
 - There are more details about how to use them in the corresponding README.md under each folder

### Acknowledgement
 - Thanks to RyanLeeHayes and Vilseck-Lab for open-source scripts, corresponds to "Multiple Topology Creation" and "MSLD_Simulations" respectively
 - This workflow is based on the scripts of the two, and the contributions of their lab to MSLD
