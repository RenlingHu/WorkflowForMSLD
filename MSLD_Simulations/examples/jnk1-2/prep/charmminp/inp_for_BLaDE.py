#! /usr/bin/env python

##
## Using CHARMM inp to prepare BLaDE inp needed for MSLD Simulations
##

import numpy as np
import pandas as pd
import glob,os
import argparse

parser = argparse.ArgumentParser(description='main parameters')
parser.add_argument('--sysname',type=str,default=None)
parser.add_argument('--complex',type=str,default=True)
args = parser.parse_args()

sysname = args.sysname
complex = args.complex

if complex:
    fp = open(sysname+'.inp','r')
    name = 'COMPLEX'
else:
    fp = open(sysname+'-lig.inp','r')
    name = 'LIGAND'
for line in fp.readlines():
    if line[0:7] == 'set box':
        boxsize=line[10:12]

fp2=open("../"+sysname+".inp","w")
fp2.seek(0)
fp2.truncate()

fp2.write("variables set box "+boxsize+"\n\n")

fp2.write("variables set ligseg LIG\n")
fp2.write("variables set resnum 1\n\n")

fp2.write("parameters file prep/toppar/par_all36m_prot.prm \n")
fp2.write("parameters file prep/toppar/par_all36_cgenff.prm\n")
fp2.write("parameters file prep/toppar/full_ligand.prm\n")
fp2.write("parameters file prep/toppar/par_water_ions.prm\n\n")

fp2.write("structure file psf prep/" +name+ "_MINI.psf\n")
fp2.write("coordinates file crd prep/" +name+ "_MINI.crd\n")
fp2.write("coordinates box cubi {box} 0 0   0 {box} 0   0 0 {box}\n")
fp2.write("coordinates velocity {temp}\n\n")

fp2.write("selection limit 200\n\n")

fp.seek(0)
for line in fp.readlines():
    if line[0:11] == 'define site':
        fp2.write("selection define site"+line[11]+"sub"+line[16]+" atomnames")
    if line[3:23] == 'atom @ligseg @resnum':
        fp2.write(" "+line[24:28])
    if line[3:7] == 'none':
        fp2.write("\n")

fp2.write("\nmsld nblocks {nblocks}\n\n")
fp2.write("variables set prevblock 0\n")
fp2.write("variables set ii 1\n")
fp2.write("while <= {ii} {nsites}\n")
fp2.write("   variables set jj 1\n")
fp2.write("   while <= {jj} {nsubs{ii}}\n")
fp2.write("      variables calculate jp0 int + {jj} {prevblock}\n\n")
fp2.write("      msld call {jp0} site{ii}sub{jj}\n\n")
fp2.write("      variables calculate jj int + {jj} 1\n")
fp2.write("   endwhile\n")
fp2.write("   variables calculate prevblock int + {prevblock} {nsubs{ii}}\n")
fp2.write("   variables calculate ii int + {ii} 1\n")
fp2.write("endwhile\n\n")

fp2.write("!               B S   t0  tv  tm  fb   q\n")
fp2.write("   msld initialize 0 0 0 0 5 0 0\n\n")

fp2.write("variables set blockassign 0\n")
fp2.write("variables set prevblock 0\n")
fp2.write("variables set ii 1\n")
fp2.write("while <= {ii} {nsites}\n")
fp2.write("   variables set jj 1\n")
fp2.write("   variables set theta0 1\n")
fp2.write("   while <= {jj} {nsubs{ii}}\n")
fp2.write("      variables calculate jp0 int + {jj} {prevblock}\n\n")

fp2.write("!               B S   t0  tv  tm  fb   q\n")
fp2.write("      msld initialize {jp0} {ii} {theta0} 0 5 {lams{ii}s{jj}} 0\n\n")

fp2.write("      ! variables set theta0 -1 ! start all equal\n")
fp2.write("      variables calculate jj int + {jj} 1\n")
fp2.write("   endwhile\n")
fp2.write("   variables calculate prevblock int + {prevblock} {nsubs{ii}}\n")
fp2.write("   variables calculate ii int + {ii} 1\n")
fp2.write("endwhile\n")
fp2.write("   msld removescaling bond urey angle dihe impr \n\n")
fp2.write("   msld softcore on\n")
fp2.write("   msld softcore14 off\n")
fp2.write("   ! pmel ex ! msld ewaldtype 2 ! default\n\n")

fp2.write("variables set prevblock 0\n")
fp2.write("variables set ii 1\n")
fp2.write("while <= {ii} {nsites}\n")
fp2.write("   variables set jj 1\n")
fp2.write("   while <= {jj} {nsubs{ii}}\n")
fp2.write("      variables calculate jp0 int + {jj} {prevblock}\n")
fp2.write("      variables calculate kk int + {jj} 1\n")
fp2.write("      while <= {kk} {nsubs{ii}}\n")
fp2.write("         variables calculate kp0 int + {kk} {prevblock}\n\n")

fp2.write("         msld bias {jp0} {kp0} 6 0.0 {cs{ii}s{jj}s{ii}s{kk}} 0\n")
fp2.write("         msld bias {jp0} {kp0} 10 -5.56 {xs{ii}s{jj}s{ii}s{kk}} 0\n")
fp2.write("         msld bias {jp0} {kp0} 8 0.017 {ss{ii}s{jj}s{ii}s{kk}} 0\n")
fp2.write("         msld bias {kp0} {jp0} 10 -5.56 {xs{ii}s{kk}s{ii}s{jj}} 0\n")
fp2.write("         msld bias {kp0} {jp0} 8 0.017 {ss{ii}s{kk}s{ii}s{jj}} 0\n")
fp2.write("         variables calculate kk int + {kk} 1\n")
fp2.write("      endwhile\n")
fp2.write("      variables calculate jj int + {jj} 1\n")
fp2.write("   endwhile\n")
fp2.write("   variables calculate prevblock int + {prevblock} {nsubs{ii}}\n")
fp2.write("   variables calculate ii int + {ii} 1\n")
fp2.write("endwhile\n")
fp2.write("\n")

fp2.close()
fp.close()