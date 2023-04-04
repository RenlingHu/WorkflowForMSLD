#! /usr/bin/env python3.9
# numpy and pandas are required

# Executable script 
# to create 'Multiple Topology Model' for ligands and Input Files for the next 'Simulation System Setup' and 'MSLD Simulations'

# Execute this script twice for：
# (1) Maximum Common Substructure Search(MCSS)
# (2) Charge Renormalization & Input Files generation

import msld_chk
import msld_mcs
import msld_crn
import msld_prm
import msld_wrt
import glob
import argparse

parser = argparse.ArgumentParser(description='main parameters')
parser.add_argument('--sysname',type=str,default=None)
parser.add_argument('--cgenff',type=str,default=True)
parser.add_argument('--boxsize',type=int,default=None)
parser.add_argument('--pmegrid',type=int,default=None)
args = parser.parse_args()

#####################################################################
## (1) Define System and File Variables

sysname = args.sysname                    # name of future output files
molfile = "mol_list.txt"                  # list of mol2 file names

mcsout = 'MCS_for_MSLD.txt'               # MCS output filename
outdir = 'build.'+sysname                 # MSLD output directory

cgenff = args.cgenff                      # Are CGenFF/ParamChem parameters being used?
boxsize = args.boxsize
pmegrid = args.pmegrid

inFrag=[[],[]]  # reflig core atoms to include in each fragment at each site (list of nsub lists)
AnCore=[[],[]]  # anchor atoms at each site to include in the core (list of nsub lists)

#####################################################################
if len(glob.glob(mcsout)) == 0:
    ## (2) Check molfile and toppar files before getting started
    msld_chk.MsldCHK(molfile)
    print("chk finished")
    
    #####################################################################
    ## (3) Maximum Common SubStruct Search with bonded-environments
    ## "mcsout" = results of the search; edit this file to manual edit ligand splicing
    ## cutoff = RMSD & distance cutoff to differentiate different atoms
    ## change debug to True to get more stdout printed 
    
    reflig = msld_mcs.MsldMCS(molfile,mcsout,cutoff=0.8,debug=False)
    print("MCS results printed to "+mcsout)
    print("Reference Ligand is "+reflig)
    quit()

#####################################################################
## (4) Perform Charge-Renormalization 

msld_crn.MsldCRN(mcsout,outdir,inFrag,AnCore,ChkQChange=True,verbose=True,debug=False)

#####################################################################
## (5) Write Ligand Parameters & the Charmm ALF input scripts
msld_prm.MsldPRM(outdir,cgenff,verbose=False,debug=False)
msld_wrt.writeALF_Files(sysname,outdir,cgenff,boxsize,pmegrid)

## Final Notes to the user
print("default TOPPAR parameters copied into build."+sysname+". Check to make sure these work for your system!")

## FINISHED