#! /usr/bin/env python

# Maximum Common Substructure search

import numpy as np
from copy import deepcopy
    
# Use bond-connectivity, atomtype definitions, and distance metrics to identify the core and fragment atoms for a set of supplied molecules
# Requirments for successful use:
#     ** Molecules MUST be spacially aligned
#     ** Atom names (within each mol2 file) must be unique for each atom(atom names can be similar between different mol2 files)
#     ** Atom names cannot have the "+" symbol in their name

def MsldMCS(molfile,mcsout,cutoff=0.8,debug=False):
    ms=[]    # ms is a list of mol2 files
    atoms=[] # atoms is a list of lists of all atom names
    heavy=[] # heavy is a list of lists of HEAVY atom names only
    xyzs=[]  # xyzs is a list of dictionaries for coordinates
    bonds=[] # bonds is a list of np.arrays of bonds (0/1 = no/yes bond)
    types=[] # types is a list of dictionaries of atom names -> atom types
    
    ## (1) Read "mol_list"
    fp=open(molfile,'r')
    for line in fp:
        ms.append(line.rstrip())
    fp.close()
    
    ## (2) Read each mol file to load in the coordinates for each atom and its bonds
    # loop over each file
    for m in range(len(ms)):
        atoms.append([])
        heavy.append([])
        xyzs.append({})
        fp=open(ms[m]+'.mol2','r')
        line=fp.readline()
        while line:
            if line.rstrip() == '@<TRIPOS>ATOM':
                # extract atom names and coordinates
                line=fp.readline()
                while line.rstrip() != '@<TRIPOS>BOND':
                    tmp=line.split()
                    atoms[m].append(tmp[1])
                    xyzs[m][tmp[1]]={'X':float(tmp[2]),'Y':float(tmp[3]),'Z':float(tmp[4]),'Elem':tmp[5].split('.')[0]}
                    if tmp[5].split('.')[0] != 'H':
                        heavy[m].append(tmp[1])
                    line=fp.readline()
            if line.rstrip() == '@<TRIPOS>BOND':
                bonds.append(np.zeros((len(atoms[m]),len(atoms[m])),dtype=int))
                line=fp.readline()
                # create a bond matrix
                while line:
                    if line == '\n':
                        # empty line
                        break
                    if line[0] == '@':
                        # signifies end of "BOND" section of mol2
                        break
                    tmp=line.split()
                    bonds[m][int(tmp[1])-1][int(tmp[2])-1]=1
                    bonds[m][int(tmp[2])-1][int(tmp[1])-1]=1
                    line=fp.readline()
            line=fp.readline()
        fp.close()
    
    if debug:
        for m in range(len(ms)):
            print("Molecule",ms[m])
            print(atoms[m])
            print(xyzs[m])
            print(bonds[m])
            for l in range(bonds[m].shape[0]):
                print(bonds[m][l])
            print("")
            print("")

    ## (3) Read the rtf file to get the atom types
    for m in range(len(ms)):
        types.append({})
        numLPs=0
        fp=open(ms[m]+'.rtf','r')
        line=fp.readline()
        while line:
            if line[0:4] == 'ATOM':
                tmp=line.split()
                # check for LonePair site
                if tmp[1][0:2] == 'LP':
                    numLPs+=1
                else:
                    types[m][tmp[1]] = tmp[2]
            line=fp.readline()
        # chk for same number of atoms
        if len(types[m]) != len(atoms[m]):
            print("ERROR: Inconsistent # of atoms between types[m] and atoms[m] after reading RTF file",ms[m]+'.str')
            print("len(types[m])=",len(types[m]),"; len(atoms[m])=",len(atoms[m]))
            quit()
    
    if debug:
        for m in range(len(types)):
            print(types[m])
            print("")
    
    ## (4) Start finding the core (heavy atom only search)
    # function to create an atom's bond pattern 
    def getBonded(atomnum,molnum):
        # first pull out the atom indices that "atomnum" is bonded to
        b1=[]
        for at in range(bonds[molnum].shape[0]):
            if bonds[molnum][atomnum][at] == 1:
                b1.append(at)
        # now convert those indices to atom types
        b2=[]
        for at in range(len(atoms[molnum])):
            if at in b1:
                b2.append(types[molnum][atoms[molnum][at]])
        bonded=['1','2']
        bonded[0]=types[molnum][atoms[molnum][atomnum]]  # bonded[0] = current atom's type
        bonded[1]=b2                                     # bonded[1] = neighbor types
    
        return bonded
    
    if debug:
        for m in range(len(ms)):
            patt = getBonded(16,m)
            print(patt)
    
    ## (i) Find all atom type matches
    matches=[]
    for m1 in range(len(ms)):
        matches.append({})
        for at1 in range(len(atoms[m1])):
            # compare only heavy atoms
            if xyzs[m1][atoms[m1][at1]]['Elem'] == 'H':
                continue
            matches[-1][atoms[m1][at1]]=[]
            patt1=getBonded(at1,m1)
            for m2 in range(len(ms)):
                matches[-1][atoms[m1][at1]].append([])
                if m2 == m1:
                    continue
                for at2 in range(len(atoms[m2])):
                    if xyzs[m2][atoms[m2][at2]]['Elem'] == 'H':
                        continue
                    patt2=getBonded(at2,m2)
                    # compare patterns 1 and 2 - store if they match
                    # if at1 and at2 types match and they both have the same number of neighbors
                    if (patt1[0] == patt2[0]) and (len(patt1[1]) == len(patt2[1])):
                        atsum=0
                        for at in patt1[1]:
                            if at in patt2[1]:
                                atsum+=1
                                # remove at(atoms) from patt2[1] to prevent duplicate matches
                                for i in range(len(patt2[1])):
                                    if at == patt2[1][i]:
                                        patt2[1].pop(i)
                                        break
                        if atsum == len(patt1[1]):
                            # store possible match
                            matches[m1][atoms[m1][at1]][m2].append(atoms[m2][at2])
    
    # debug functions
    def printListDict(struct):
        for m in range(len(ms)):
            print("Molecule",m)
            keys = struct[m].keys()
            for k in keys:
                print(k,struct[m][k])
            print("")
    def printListList(struct):
        for m in range(len(ms)):
            print("Molecule",m)
            for l in struct[m]:
                print(l)
            print("")
    
    if debug:
        print("FIRST MATCH based on atom_type and connections")
        printListDict(matches)
    
                
    ## (ii) start sorting through the matches
    ## (iiA) remove entries with more than one empty list 
    for m1 in range(len(ms)):
        droplist=[]
        ibuff=-1
        for at1 in range(len(atoms[m1])):
            empty=0
            if xyzs[m1][atoms[m1][at1]]['Elem'] == 'H':
                continue
            ibuff+=1
            for m2 in range(len(ms)):
                if matches[m1][atoms[m1][at1]][m2] == []:
                    empty+=1
            if empty > 1:
                droplist.append(atoms[m1][at1])
        # pop off items from matches dict
        for i in droplist:
            matches[m1].pop(i)
    
    if debug:
        printListDict(matches)
                
    ## (iiB) remove entries with exactly one match and add them to cores
    cores=[] # core atoms and their matches
    for m1 in range(len(ms)):
        droplist=[]
        cores.append({})
        for at1 in range(len(atoms[m1])):
            if atoms[m1][at1] in matches[m1].keys():
                chk=0
                for m2 in range(len(ms)):
                    if m1 != m2:
                        if len(matches[m1][atoms[m1][at1]][m2]) == 1:
                            chk+=1
                if chk == len(ms) - 1: # then we have a core atom match
                    droplist.append(atoms[m1][at1])
        for i in droplist:
            at=matches[m1].pop(i)
            cores[m1][i]=at
    
    if debug:
        print("AFTER FIRST EMPTY or 1 MATCH CHECKS")
        printListDict(cores)
        printListDict(matches)
        print("")
    
    ## (iiC) Make sure all core lists have the same number of atoms: Correct inconsistencies
    chkcores=[]
    for m1 in range(len(ms)):
        chkcores.append([])
        m1cores=list(cores[m1].keys())
        for m2 in range(len(ms)):
            chkcores[m1].append([])
            if m1 != m2:
                m2cores=list(cores[m2].keys())
                # generate list of m1 core matches from cores[m2] key-values
                m2matches=[]
                for at2 in m2cores:
                    m2matches.extend(cores[m2][at2][m1])
                for at1 in m1cores:
                    if at1 in m2matches:
                        chkcores[m1][m2].append(at1)
                        line=deepcopy(cores[m1][at1])
                        line[m1]=[at1]
                        cores[m2][line[m2][0]]=line
                        cores[m2][line[m2][0]][m2]=[]

    if debug:
        print(chkcores)
        print(cores)
        print(matches)

    ## (iiD) Double check that core atoms do not have the same match. If they do,
    ## then move them out of cores and back into matches
    for m1 in range(len(ms)):
        m1keys=list(cores[m1].keys())
        for m2 in range(len(ms)):
            if m1 != m2:
                droplist=[]
                for k1 in range(len(m1keys)):
                    at=cores[m1][m1keys[k1]][m2]
                    for k2 in range(k1+1,len(m1keys)):
                        if at == cores[m1][m1keys[k2]][m2]:
                            # then we have a match that we need to fix
                            if not (m1keys[k1] in droplist):
                                droplist.append(m1keys[k1])
                            if not (m1keys[k2] in droplist):
                                droplist.append(m1keys[k2])
                if len(droplist) > 0:
                    for k in droplist:
                        matches[m1][k]=cores[m1][k]
                        cores[m1].pop(k)
                    m1keys=list(cores[m1].keys())      
    
    # debug
    if debug:
        print("AFTER THE FIRST CORE CHK")
        printListDict(matches)
        printListDict(cores)
        print("")

    ## (iiE) RMSD analysis between corrent core atoms
    ## (exit with error message if RMSD is > cutoff)
    
    #cutoff = 0.8    # disregard atom pairs with RMSDs above this value ## uncomment to hardcode this
       
    # For each core atom, compute an "average" xyz position, and then calc the RMSD for all
    # real atomic coordinates to this average "reference" position
    rmsd={}
    for m1 in range(1):   # treat first molecule as the "reference molecule"
        for k in cores[m1].keys():
            avg={"X":0.00,"Y":0.00,"Z":0.00}
            rmsd[k]=0.00
            # calc avg position first
            avg['X']+=xyzs[m1][k]["X"]
            avg['Y']+=xyzs[m1][k]["Y"]
            avg['Z']+=xyzs[m1][k]["Z"]
            for m2 in range(1,len(ms)):
                avg['X']+=xyzs[m2][cores[m1][k][m2][0]]["X"]
                avg['Y']+=xyzs[m2][cores[m1][k][m2][0]]["Y"]
                avg['Z']+=xyzs[m2][cores[m1][k][m2][0]]["Z"]
            avg['X']=avg['X']/float(len(ms))
            avg['Y']=avg['Y']/float(len(ms))
            avg['Z']=avg['Z']/float(len(ms))
    
            # then calc rmsd for all atoms pairs to this average position
            tmp=xyzs[m1][k]["X"]-avg['X']
            rmsd[k]+=tmp*tmp
            tmp=xyzs[m1][k]["Y"]-avg['Y']
            rmsd[k]+=tmp*tmp
            tmp=xyzs[m1][k]["Z"]-avg['Z']
            rmsd[k]+=tmp*tmp
            for m2 in range(1,len(ms)):
                tmp=xyzs[m2][cores[m1][k][m2][0]]["X"]-avg['X']
                rmsd[k]+=tmp*tmp
                tmp=xyzs[m2][cores[m1][k][m2][0]]["Y"]-avg['Y']
                rmsd[k]+=tmp*tmp
                tmp=xyzs[m2][cores[m1][k][m2][0]]["Z"]-avg['Z']
                rmsd[k]+=tmp*tmp
            rmsd[k]=np.sqrt(rmsd[k]/float(len(ms)))
    
    chk=0
    for k in rmsd.keys():
        if rmsd[k] > cutoff:
            chk+=1
            print("Large RMSD for CORE atom",k,"- RMSD >",cutoff,"("+str(round(rmsd[k],3))+").")
    if chk == len(cores[0].keys()) and len(cores[0].keys()) != 0:
        print("ERROR: Large RMSDs detected for ALL core atoms.\n     ",
              "This suggests that the molecules are not properly aligned!. Check and resubmit")
        quit()


    ## (iiF) match the atoms based on their spatial distance again
    ## *** !!! This requires the molecules to be ALIGNED to work correctly !!! ***
    ## Move on to using distance to try to deduce if the remaining atoms pairs can be matched and moved into the core
    ## Calc distance between atom pairs and eliminate options above the cutoff value
    matches=[]
    distance=[]
    for m1 in range(len(ms)):
        matches.append({})
        distance.append({})
        for at1 in range(len(atoms[m1])):
            # compare atoms not already in core atoms and heavy atoms only
            if xyzs[m1][atoms[m1][at1]]['Elem'] == 'H':
                continue
            if atoms[m1][at1] in cores[m1].keys():
                continue
            else:
                matches[-1][atoms[m1][at1]]=[]
                distance[-1][atoms[m1][at1]]=[]
                dist = 0.00
                X = xyzs[m1][atoms[m1][at1]]['X']
                Y = xyzs[m1][atoms[m1][at1]]['Y']
                Z = xyzs[m1][atoms[m1][at1]]['Z']
                distance[-1][atoms[m1][at1]].append(dist)
                for m2 in range(len(ms)):
                    matches[-1][atoms[m1][at1]].append([])
                    if m2 == m1:
                        continue
                    if m2 != m1:
                        for at2 in range(len(atoms[m2])):
                            if xyzs[m2][atoms[m2][at2]]['Elem'] != xyzs[m1][atoms[m1][at1]]['Elem']:
                                continue
                            dist = 0.00
                            tmp = X-xyzs[m2][atoms[m2][at2]]['X']
                            dist += tmp*tmp
                            tmp = Y-xyzs[m2][atoms[m2][at2]]['Y']
                            dist += tmp*tmp
                            tmp = Z-xyzs[m2][atoms[m2][at2]]['Z']
                            dist += tmp*tmp
                            dist = np.sqrt(dist)
                            if dist < 0.10:  # require ligands to align very strictfully
                                distance[m1][atoms[m1][at1]].append(dist)
                                matches[m1][atoms[m1][at1]][m2].append(atoms[m2][at2])
                            else:
                                continue
    
    if debug:
        print("SECOND MATCH based on atom_spatial_direction")
        print(matches)

    ## remove entries with more than one empty list       
    for m1 in range(len(ms)):
        droplist=[]
        ibuff=-1
        for at1 in range(len(atoms[m1])):
            empty=0
            if xyzs[m1][atoms[m1][at1]]['Elem'] == 'H':
                continue
            if atoms[m1][at1] in cores[m1].keys():
                continue
            ibuff+=1
            for m2 in range(len(ms)):
                if matches[m1][atoms[m1][at1]][m2] == []:
                    empty+=1
            if empty > 1:
                droplist.append(atoms[m1][at1])
        # pop off items from matches dict
        for i in droplist:
            matches[m1].pop(i)

    ## remove entries with exactly one match and add them to cores
    for m1 in range(len(ms)):
        droplist=[]
        for at1 in range(len(atoms[m1])):
            if atoms[m1][at1] in matches[m1].keys():
                chk=0
                for m2 in range(len(ms)):
                    if m1 != m2:
                        if len(matches[m1][atoms[m1][at1]][m2]) == 1:
                            chk+=1
                if chk == len(ms) - 1: # then we have a core atom match
                    droplist.append(atoms[m1][at1])
        for i in droplist:
            at=matches[m1].pop(i)
            cores[m1][i]=at
    
    if debug:
        for m in range(len(ms)):
            print(cores[m].keys())

    ## remove entries without same number H(s) attached except for the anchor heavy atoms
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):   # check at1 in cores[m1].keys()
            if atoms[m1][at1] in cores[m1].keys():
                hbond1=0
                for h1 in range(len(atoms[m1])):
                    if xyzs[m1][atoms[m1][h1]]['Elem'] == 'H' and bonds[m1][at1][h1] == 1:
                        hbond1+=1
                for m2 in range(len(ms)):
                    if m2 != m1:
                        for at2 in range(len(atoms[m2])):
                            if atoms[m2][at2] in cores[m1][atoms[m1][at1]][m2]:
                                chk=0
                                for h2 in range(len(atoms[m2])):
                                    if xyzs[m2][atoms[m2][h2]]['Elem'] == 'H' and bonds[m2][at2][h2] == 1:
                                        chk+=1
                        if hbond1 != chk:
                            for at12 in range(len(atoms[m1])):  # at12 which attaches to at1(without same H(s) attached)
                                if atoms[m1][at12] in cores[m1].keys() and atoms[m1][at12] != atoms[m1][at1] and bonds[m1][at1][at12] == 1:
                                    bond=0
                                    for at13 in range(len(atoms[m1])):
                                        if atoms[m1][at13] in cores[m1].keys() and atoms[m1][at13] != atoms[m1][at12]:
                                            if bonds[m1][at12][at13] == 1:
                                                bond+=1
                                    if bond == 1:
                                        cores[m1].pop(atoms[m1][at12])
                                        break
                            break
    
    if debug:
        print(cores[0].keys())

    ## remove entries without core-core attachment >= 1
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):   # check at1 in cores[m1].keys()
            if atoms[m1][at1] in cores[m1].keys():
                at1_at2_bondcount=0
                for at2 in range(len(atoms[m1])):
                    if atoms[m1][at2] in cores[m1].keys() and atoms[m1][at2] != atoms[m1][at1]:
                        if bonds[m1][at1][at2] == 1:
                            at1_at2_bondcount+=1
                if at1_at2_bondcount == 0:
                    cores[m1].pop(atoms[m1][at1]) 

    if debug:
        print(cores[0].keys())
    
    ## remove (1)bezene rings with > 1 substituent sites (2)core atoms indicate different conjugated heterocycle
    def GetIndex(atomname,molnum):
        """ given an atomname for a molecule number, return the atom index"""
        for atomindex in range(len(atoms[molnum])):
            if atoms[molnum][atomindex] == atomname:
                break
        return atomindex    
    
    ## (i) attachment situation of conjugated cyclic atoms
    attach=[]
    for m1 in range(len(ms)):
        sub=0
        attach.append({})
        for at1 in range(len(atoms[m1])):
            if atoms[m1][at1] in cores[m1].keys():
                if types[m1][atoms[m1][at1]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                    attach[-1][atoms[m1][at1]]=[]
                    bond=0
                    for at12 in range(len(atoms[m1])):
                        if bonds[m1][at1][at12] == 1:
                            attach[-1][atoms[m1][at1]].append(atoms[m1][at12])
                    
    drop=[] # would be deleted if attach to non-conjugates-cyclic atoms
    for m1 in range(len(ms)):
        drop.append([])
        for at1 in attach[m1].keys():
            chk=0
            for at2 in range(len(attach[m1][at1])):
                if attach[m1][at1][at2] in attach[m1].keys():
                    chk+=1
            if chk == 0:
                drop[-1].append(at1)
    for m1 in range(len(ms)):
        for at1 in drop[m1]:
            attach[m1].pop(at1)
    
    ## (ii) find atoms where would be a substituent site
    site=[]
    attach_type=[]
    type_count=[]
    for m1 in range(len(attach)):
        attach_type.append({})
        for at1 in attach[m1].keys():
            attach_type[-1][at1]=[]
            for at2 in attach[m1][at1]:
                attach_type[-1][at1].append(xyzs[m1][at2]['Elem'])
    for m1 in range(len(ms)):
        type_count.append({})
        for at1 in attach_type[m1].keys():
            type_count[-1][at1]=[]
            countC=0
            countN=0
            countH=0
            countO=0
            countX=0
            for at2 in attach_type[m1][at1]:
                if at2 == 'C':
                    countC+=1
                if at2 == 'N':
                    countN+=1
                if at2 == 'H':
                    countH+=1
                if at2 == 'O':
                    countO+=1
                if at2 == 'Cl' or at2 == 'F' or at2 == 'Br':
                    countX+=1
            type_count[-1][at1].append(countC)
            type_count[-1][at1].append(countN)
            type_count[-1][at1].append(countH)
            type_count[-1][at1].append(countO)
            type_count[-1][at1].append(countX)
    for m1 in range(len(ms)):
        site.append([])
        for at1 in attach_type[m1].keys():
            chk=0
            for m2 in range(len(ms)):
                if m2 != m1:
                    for at2 in attach_type[m2].keys():
                        if at2 in cores[m1][at1][m2]:
                            if type_count[m2][at2] != type_count[m1][at1]:
                                chk+=1
            if chk >= 1:
                site[-1].append(at1)
    
    ## (iii) delete all cyclic atoms related to these sites
    for m in range(len(ms)):
        if len(site[m])>1:
            cycle=[]
            add=[]
            for m1 in range(len(ms)):
                cycle.append([])
                for at1 in range(len(atoms[m1])):
                    if atoms[m1][at1] in site[m1]:
                        if not (attach_type[m1][atoms[m1][at1]] in [['C','C','C'],['C','C','N'],['C','N','N'],['O','C','C'],['N','N','N'],['N','C','O']]):
                            for at2 in range(len(atoms[m1])):
                                if bonds[m1][at1][at2]==1 and types[m1][atoms[m1][at2]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                    if not (atoms[m1][at2] in cycle[-1]):
                                        cycle[-1].append(atoms[m1][at2])
                            if not (atoms[m1][at1] in cycle[-1]):
                                cycle[-1].append(atoms[m1][at1])
                            break
                add.append([])
                for at1 in range(len(atoms[m1])):
                    if atoms[m1][at1] in cycle[m1]:
                        for at2 in range(len(atoms[m1])):
                            if bonds[m1][at1][at2]==1 and types[m1][atoms[m1][at2]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                chk=0
                                for at3 in range(len(atoms[m1])):
                                    if types[m1][atoms[m1][at3]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                        if bonds[m1][at2][at3]==1 and atoms[m1][at3] in cycle[m1]:
                                            chk+=1
                                if chk >= 1:
                                    if not (atoms[m1][at2] in add[m1]):
                                        add[-1].append(atoms[m1][at2])
                for atom in add[m1]:
                    if not (atom in cycle[m1]):
                        cycle[m1].append(atom)
                for at1 in range(len(atoms[m1])):
                    if atoms[m1][at1] in cycle[m1]:
                        chk=0
                        for at2 in range(len(atoms[m1])):
                            if not (atoms[m1][at2] in cycle[m1]):
                                if bonds[m1][at1][at2]==1 and types[m1][atoms[m1][at2]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                    chk+=1
                        if chk > 1:
                            if atoms[m1][at1] in cycle[m1]:
                                cycle[m1].remove(atoms[m1][at1])
                ## check cycle_atoms
                for at1 in range(len(atoms[m1])):
                    if atoms[m1][at1] in cycle[m1]:
                        for at2 in range(len(atoms[m1])):
                            if bonds[m1][at1][at2]==1 and types[m1][atoms[m1][at2]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                chk=0
                                for at3 in range(len(atoms[m1])):
                                    if types[m1][atoms[m1][at3]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                        if bonds[m1][at2][at3]==1 and atoms[m1][at3] in cycle[m1]:
                                            chk+=1
                                if chk >= 1:
                                    if not (atoms[m1][at2] in add[m1]):
                                        add[-1].append(atoms[m1][at2])
                for atom in add[m1]:
                    if not (atom in cycle[m1]):
                        cycle[m1].append(atom)
                ## check cycle_atoms
                for at1 in range(len(atoms[m1])):
                    if atoms[m1][at1] in cycle[m1]:
                        chk=0
                        for at2 in range(len(atoms[m1])):
                            if not (atoms[m1][at2] in cycle[m1]):
                                if bonds[m1][at1][at2]==1 and types[m1][atoms[m1][at2]] in ['CG2R61','NG2R52','CG2R53','CG2R67','NG2R60','CG2RC0']:
                                    chk+=1
                        if chk > 1:
                            if atoms[m1][at1] in cycle[m1]:
                                cycle[m1].remove(atoms[m1][at1])
                break
            for m2 in range(1,len(ms)):
                cycle.append([])
                for atom in cycle[0]:
                    if atom in cores[0].keys():
                        cycle[m2].extend(cores[0][atom][m2])
            break
    for m in range(len(ms)):
        if len(site[m])>1:
            for m in range(len(ms)):
                for atom in cycle[m]:
                    if atom in cores[m].keys():
                        cores[m].pop(atom)
            break

    if debug:
        for m in range(len(ms)):
            print(cores[m].keys())
    
    ## remove entries without core-core attachment >= 1
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):   # check at1 in cores[m1].keys()
            if atoms[m1][at1] in cores[m1].keys():
                at1_at2_bondcount=0
                for at2 in range(len(atoms[m1])):
                    if atoms[m1][at2] in cores[m1].keys() and atoms[m1][at2] != atoms[m1][at1]:
                        if bonds[m1][at1][at2] == 1:
                            at1_at2_bondcount+=1
                if at1_at2_bondcount == 0:
                    cores[m1].pop(atoms[m1][at1]) 

    if debug:
        print("AFTER SECOND EMPTY or 1 MATCH CHECKS")
        printListDict(cores)
        printListDict(matches)
        print("")
    
    ## (iiG) Try a variety of checks to ensure that all atoms in matches are either classified
    ## as a core atom or discarded as a future fragment atom
    def emptyMatches(mollist):
        """ check if the matches structure is empty for all molecules in mollist """
        echk=0
        for mol in range(len(mollist)):
            if len(matches[m]) != 0:
                echk+=1
        return echk
        
    ## (iiH) Make sure all atoms are out of matches - exit with error if not
    chk=emptyMatches(ms)
    if chk > 0:
        print("ERROR: Atoms have not been classified as core or fragment atoms!\n    ",
              "User input is required to continue running. Please specify the matching atom(s)",
              "for each atom in each molecule\n")
        for m1 in range(len(ms)):
            if len(matches[m1]) > 0:
                print("### Molecule",m1,"("+ms[m1]+")")
                print(matches[m1])
            for at1 in range(len(atoms[m1])):
                if atoms[m1][at1] in matches[m1].keys():
                    print("Molecule",m1,"("+ms[m1]+") - ATOM",atoms[m1][at1],"matches with:")
                    for m2 in range(len(ms)):
                        if m1 != m2:
                            if len(matches[m1][atoms[m1][at1]][m2]) > 1:
                                useratom=input("Molecule "+str(m2)+" ("+ms[m2]+"),("+str(matches[m1][atoms[m1][at1]][m2])+"): ")
                                if useratom == 'exit' or useratom == 'EXIT' or useratom == 'Exit':
                                    print("Exiting...")
                                    quit()
                                matches[m1][atoms[m1][at1]][m2]=[useratom]
                    print("")
            print("\n")
        printListDict(matches)
        quit()
    
    ## (iiI) Make sure all core lists have the same number of atoms: Correct inconsistencies
    chkcores=[]
    for m1 in range(len(ms)):
        chkcores.append([])
        m1cores=list(cores[m1].keys())
        for m2 in range(len(ms)):
            chkcores[m1].append([])
            if m1 != m2:
                m2cores=list(cores[m2].keys())
                # generate list of m1 core matches from cores[m2] key-values
                m2matches=[]
                for at2 in m2cores:
                    m2matches.extend(cores[m2][at2][m1])
                for at1 in m1cores:
                    if at1 in m2matches:
                        chkcores[m1][m2].append(at1)
                        line=deepcopy(cores[m1][at1])
                        line[m1]=[at1]
                        cores[m2][line[m2][0]]=line
                        cores[m2][line[m2][0]][m2]=[]
    
    if debug:
        print("AFTER THE SECOND CORE CHK")
        print(matches)
        print("CORE ATOMS")
        printListDict(cores)
        print("")
    
    ## (5) Figure out the number of fragments in each molecule, pair them together between molecules,
    ## and remove redundancies
    
    fragbonds=deepcopy(bonds)         ## no core-core bonds
    for m1 in range(len(ms)):
        for at1 in range(len(atoms[m1])):
            for at2 in range(at1+1,len(atoms[m1])):
                if (atoms[m1][at1] in cores[m1].keys()) and (atoms[m1][at2] in cores[m1].keys()):
                        fragbonds[m1][at1][at2]=0
                        fragbonds[m1][at2][at1]=0
    
    ## figure out the number of core-to-(not core) attachments
    nsites=[]
    Aatoms=[]
    for m in range(len(ms)):
        Aatoms.append({})
        for at1 in range(len(atoms[m])):
            if atoms[m][at1] in cores[m].keys():
                for at2 in range(len(atoms[m])):
                    if xyzs[m][atoms[m][at2]]['Elem'] != 'H':
                        if fragbonds[m][at1][at2] == 1:
                            if atoms[m][at1] in Aatoms[m].keys():
                                Aatoms[m][atoms[m][at1]].append(atoms[m][at2])
                            else:
                                Aatoms[m][atoms[m][at1]]=[atoms[m][at2]]
        chk=len(Aatoms[m].keys())      # multi_attached core_atom counts once only
        nsites.append(chk)

    ## add Anchor-to-H attachments for multisites
    ## a. for Anchor-to-only-H/H(s)
    ## b. for Anchor-to-heavyatoms+H/H(s)
    
    Acount=[]   # number of Anchor atoms
    Aatomsnum=[]    # atom_number of Anchor atoms
    for m in range(len(ms)):
        Acount.append(len(Aatoms[m].keys()))
        Aatomsnum.append([])
        for at in range(len(atoms[m])):
            if atoms[m][at] in Aatoms[m].keys():
                Aatomsnum[m].append(at)
    
    for m in range(len(ms)):
        if len(Aatomsnum[m]) == max(Acount):
            refmol = m
            refatom=[]
            for at in Aatomsnum[m]:
                refatom.append(at)
            break
    
    for m in range(len(nsites)):   # iterate through each molecule
        chk= nsites[m]   
        if nsites[m] < max(Acount):
            for at1 in range(len(atoms[m])):    # iterate through each atom
                    
                for at0 in Aatomsnum[refmol]:
                    if atoms[m][at1] in cores[refmol][atoms[refmol][at0]][m]:

                        for at2 in range(len(atoms[m])):
                            if fragbonds[m][at1][at2] == 1 and xyzs[m][atoms[m][at2]]['Elem'] == 'H':
                                if atoms[m][at1] in Aatoms[m].keys():
                                    continue
                                else:
                                    chk+=1
                                    Aatoms[m][atoms[m][at1]]=[atoms[m][at2]]
                                break
                        
                        break
        nsites[m]=chk
    
    for m in range(len(nsites)):
        if nsites[m] >= max(Acount):
            for at1 in range(len(atoms[m])):
                if atoms[m][at1] in Aatoms[m].keys():
                    for at2 in range(len(atoms[m])):
                        if fragbonds[m][at1][at2] == 1 and xyzs[m][atoms[m][at2]]['Elem'] == 'H':
                            if not (atoms[m][at2] in Aatoms[m][atoms[m][at1]]):
                                Aatoms[m][atoms[m][at1]].append(atoms[m][at2])
        
    if debug:
        print("NSITES and AATOMS")
        print(nsites)
        print(Aatoms)
    
    
    # now take these sites and start expounding on them to list out all fragment atoms attached to them
    # (remember to look for rings (where one fragment branch joins another nsite/Aatom point)
    
    Fatoms=[]
    for m in range(len(ms)):
        Fatoms.append([])
        for k in Aatoms[m].keys():
            Fatoms[m].append([])    ####so that attached multi_atoms wouldnt create multi_lists
            for fragatom in Aatoms[m][k]:
                Fatoms[m][-1].append(fragatom)
                # find the atom index for fragatom
                at1 = GetIndex(fragatom,m)
                # enumerate all bonds to this atom (excluding core atoms)
                chk=1
                atlist=[]
                while chk > 0:
                    i=fragbonds[m][at1,:]
                    for j in range(len(i)):
                        if i[j] == 1:
                            if not (atoms[m][j] in cores[m].keys()):
                                if not (atoms[m][j] in Fatoms[m][-1]):
                                    if atoms[m][j] in Aatoms[m][k]:
                                        continue
                                    else:
                                        atlist.append(j)
                    if atlist==[]:
                        chk=0
                    else:
                        at1=atlist.pop()
                        if not (atoms[m][at1] in Fatoms[m][-1]):
                            if not (atoms[m][at1] in Aatoms[m][k]):
                                Fatoms[m][-1].append(atoms[m][at1])
        # check for redundant/identical fragments & merge if found
        droplist=[]   # list the indices to drop from Fatoms[m][-1]
        tsites=nsites[m]
        for k1 in range(tsites-1):
            for k2 in range(k1+1,tsites):
                chk=0
                tmpfrag=deepcopy(Fatoms[m][k2])
                for at1 in Fatoms[m][k1]:
                    for at2 in range(len(tmpfrag)):
                        if at1 == tmpfrag[at2]:
                            chk+=1
                            tmpfrag.pop(at2)
                            break
                if chk == len(Fatoms[m][k1]):
                    # then we have a match
                    if not (k2 in droplist):
                        nsites[m]=nsites[m]-1
                        droplist.append(k2)
        # merge found matches; modify: (i) Fatoms, (ii) Aatoms (as part of next routine)
        if len(droplist) > 0:
            tmp=[]
            droplist.sort(reverse=True) # to work highest index to lowest
            for i in droplist:
                tmp=Fatoms[m].pop(i)  # pop in reverse direction so we don't remove things we actually want
        # replace Aatom values with frag atom lists
        for frag in Fatoms[m]:
            droplist=[]
            for k in Aatoms[m].keys():
                for f in range(len(Aatoms[m][k])):
                    if Aatoms[m][k][f] in frag:
                        if not (k in droplist):     ####### anchor attached to multi_atoms wouldnt be read in multi_times
                            droplist.append(k)
            if len(droplist) > 1:
                # merge key_names and replace with new frag list
                newname=""
                for at in range(len(droplist)):
                    if at == 0:
                        newname+=droplist[at]
                    else:
                        if droplist[at] != k:  ####### wouldn't make C17+C17 like DUM
                            newname+="+"+droplist[at]
    
                    if len(Aatoms[m][droplist[at]]) == 1:
                        Aatoms[m].pop(droplist[at]) # doesn't work if Aatoms[m][k] has more than one atom in it's nested list
                    else:
                        # try making a new Aatoms[m][at] list without the matching frag atom
                        for f in range(len(Aatoms[m][droplist[at]])):
                            if Aatoms[m][droplist[at]][f] in frag:
                                break
                        newlist=[]
                        for f2 in range(len(Aatoms[m][droplist[at]])):
                            if f != f2:
                                newlist.append(Aatoms[m][droplist[at]][f2])
                        Aatoms[m][droplist[at]]=newlist
                Aatoms[m][newname] = deepcopy(frag)
            elif len(droplist) == 1:
                Aatoms[m][droplist[0]] = deepcopy(frag)
            else: #len(droplist) < 1:
                print("ERROR: frag mismatch between Fatoms list and Aatoms match")
                # quit()
    
    # Do some quick checks:
    for m in range(len(ms)):
        if nsites[m] != len(Aatoms[m]) or nsites[m] != len(Fatoms[m]):
            print("ERROR: mismatch in the number of sites for molecule",m)
            print("MOLECULE",m)
            print("NSITES[M]",nsites[m])
            print("AATOMS[M]",Aatoms[m])
            print("FATOMS[M]",Fatoms[m])
            quit()
    
    if debug:
        print("nsites, Aatoms, and Fatoms after sorting is complete")
        print("NSITES:",nsites,'\n')
        print("AATOMS:")
        printListDict(Aatoms)
        print("FATOMS:")
        printListList(Fatoms)
        print("")
    
    ## (6) Identify the reference ligand (either hardcode it as the first molecule, 
    ## or pick the molecules with the smallest number of fragment atoms
    
    #refnum=0  # uncomment if you want to hardcode this value

    minfragatoms=-1
    for m in range(len(ms)):
        sumchk=0
        for k in Aatoms[m].keys():
            sumchk+=len(Aatoms[m][k])
        if minfragatoms == -1:
            minfragatoms=sumchk
            refnum=m
        else:
            if sumchk < minfragatoms:
                minfragatoms=sumchk
                refnum=m
    
    if debug:
        print("REFNUM =",refnum,"\n\n")
    
    
    ## Make sure all ligands have the same number of sites
    chk=0
    for m in range(len(ms)):
        if nsites[m] != nsites[refnum]:
            chk+=1
    if chk != 0:
        print("ERROR: Not all molecules have the same number of sites defined! Check and resubmit")
        quit()
    
    ## Order the fragments so that everything is consistent when we do redundant checks
    # use the refnum molecule to decide what fragment is first
    Ftemplate=[]
    for k in Aatoms[refnum].keys():
        Ftemplate.append(k)
    
    # generate the expected keys based off what's in cores, then check to make sure everything is correct
    def chkMergedKey(keytocheck):
        """ Check if a key is a merged key (name+name). Split into a list if yes, return value if not """
        tmp=keytocheck.split('+')
        if len(tmp) > 1:
            keyname = tmp
        else:
            keyname=[keytocheck]
        return keyname
    
    Forder=[]
    for m1 in range(len(ms)):
        Forder.append([])
        for temp in Ftemplate:
            tt=chkMergedKey(temp)
            if len(tt) > 1: # then it is a merged key
                tmp=""
                for t3 in range(len(tt)):
                    if m1 == refnum:
                        i=tt[t3]
                    else:
                        i=cores[refnum][tt[t3]][m1][0]
                    if t3 == 0:
                        tmp+=i
                    else:
                        tmp+="+"+i
                Forder[m1].append(tmp)
    
            else:           # it is a single value key
                if m1 == refnum:
                    Forder[m1].append(tt[0])
                else:
                    Forder[m1].append(cores[refnum][tt[0]][m1][0])
        # make sure no two keys are the same
        chk = 0
        for f1 in range(len(Forder[m1])):
            for f2 in range(len(Forder[m1])):
                if f1 != f2:
                    if Forder[m1][f1] == Forder[m1][f2]:
                        chk+=1
        if chk > 0:
            print("ERROR: Non-unique keys found linking core to fragment atoms!")
            quit()
        # make sure the generated key actually matches what's in Aatoms[m1]
        # and account for order differences
        chk = 0
        tmp=list(Aatoms[m1].keys())
        for f in Forder[m1]:
            f1=chkMergedKey(f)
            if len(f1) == 1:
                for f2 in range(len(tmp)):
                    if tmp[f2] == f:
                        chk+=1
                        tmp.pop(f2)
                        break
            else:
                # consider different orders; if found, update Aatoms keyname
                for f2 in range(len(tmp)):
                    f3=chkMergedKey(tmp[f2])
                    if len(f3) == len(f1):
                        chk2=0
                        for f4 in f1:
                            for f5 in range(len(f3)):
                                if f4 == f3[f5]:
                                    chk2+=1
                                    f3.pop(f5)
                                    break
                        if chk2 == len(f1):
                            chk+=1
                            kk=tmp.pop(f2)
                            break
                # update Aatoms key
                if kk != f:
                    Aatoms[m1][f]=Aatoms[m1][kk]
                    Aatoms[m1].pop(kk)
        if chk != len(Forder[m1]):
            print("ERROR: Mismatch between generated and actual Aatom keys!")
            print("Molecule",m1)
            print("AATOMS",Aatoms[m1].keys())
            print("FORDER",Forder[m1])
            quit()
    
    if debug:
        print("FORDER:")
        printListList(Forder)
    
    ## (7) Now look for redundancies in the fragments in each molecule so that we only print out a single fragment per site
    # refnum fragments are added by default
    
    # shortest distance function
    def findShortestDistance(atname1,atname2,molnum):
        """ Explore bonds to find the minimum distance between two atoms """
        at1=GetIndex(atname1,molnum)
        at2=GetIndex(atname2,molnum)
        
        ## at1 = "anchor atom"; start here and tree down until we find the atom we're looking for (at2)
        chk = 0
        atlist=[at1]
        while True:
            chk+=1
            bdlist=[]
            for i in range(len(atlist)):
                ii=bonds[molnum][atlist[i],:]
                for j in range(len(ii)):
                    if ii[j] == 1:
                        if not (j in bdlist):
                            bdlist.append(j)
            if at2 in bdlist:
                mindist = chk
                break
            else:
                atlist=deepcopy(bdlist)
        return mindist
    
    # use a compare function:
    def fragCompare(ufidx,ufkey,qidx,qkey):
        """ Compare two fragments & try to figure out if they are the same or not.
            Return "TRUE" if they match, otherwise, return "FALSE". """
        uflist=Aatoms[ufidx][ufkey]
        qlist=Aatoms[qidx][qkey]
        # reject if list lenghts are different
        if len(uflist) != len(qlist):
            return False
        # reject if lists have atoms with different bonding patterns or unequal #s of atoms with same bonding patterns
        #   i. calc bonding patterns for each atom in uflist and qlist
        ufpatt=[]
        for at in uflist:
            ufpatt.append(getBonded(GetIndex(at,ufidx),ufidx))
        qpatt=[]
        for at in qlist:
            qpatt.append(getBonded(GetIndex(at,qidx),qidx))
        #   ii. find shortest distance to a single core atom (first core atom in key if its a merged key)
        atname1=ufkey.split('+')[0]
        for atname2 in range(len(uflist)):
            mindist=findShortestDistance(atname1,uflist[atname2],ufidx)
            ufpatt[atname2].append(mindist)
        atname1=qkey.split('+')[0]
        for atname2 in range(len(qlist)):
            mindist=findShortestDistance(atname1,qlist[atname2],qidx)
            qpatt[atname2].append(mindist)
        qtemp=deepcopy(qpatt)
        #   iii. now compare types
        amatch=0
        for at1 in range(len(uflist)):
            for at2 in range(len(qlist)):
                if ufpatt[at1][0] == qtemp[at2][0] and ufpatt[at1][2] == qtemp[at2][2]: # atom type matches
                    bdsum=0
                    for bd1 in ufpatt[at1][1]:
                        if bd1 in qtemp[at2][1]:
                            bdsum+=1
                            #remove bd1 atomtype from qtemp[at2][1] to prevent duplicities
                            for i in range(len(qtemp[at2][1])):
                                if qtemp[at2][1][i] == bd1:
                                    qtemp[at2][1].pop(i)
                                    break
                    if bdsum == len(ufpatt[at1][1]): # then they match
                        amatch+=1
                        qtemp[at2][0]='MATCHED' # prevents duplicities
                        break
        #   iiii. now compare distances
        for at1 in uflist:
            XYZ = (xyzs[ufidx][at1]['X'])*(xyzs[ufidx][at1]['X']) 
            + (xyzs[ufidx][at1]['Y'])*(xyzs[ufidx][at1]['Y']) 
            + (xyzs[ufidx][at1]['Z'])*(xyzs[ufidx][at1]['Z'])
            for at2 in qlist:
                XYZ1 = (xyzs[qidx][at2]['X'])*(xyzs[qidx][at2]['X']) 
                + (xyzs[qidx][at2]['Y'])*(xyzs[qidx][at2]['Y']) 
                + (xyzs[qidx][at2]['Z'])*(xyzs[qidx][at2]['Z'])
                dd = np.sqrt(np.abs(XYZ1-XYZ))
                if dd > 0.8:
                    amatch-=1
        if amatch != len(uflist):
            return False
        else:
            return True
    
    UFrag = []  # list of Unique Fragment ms indices
    for site in range(nsites[refnum]):
        UFrag.append([refnum])
    
    for site in range(nsites[refnum]):
        skip=[refnum] # refnum is included in UFrag by default, so we want to skip it in our comparisons
        ifrag=-1
        while len(skip) != len(ms):
            for frag in range(ifrag+1,len(UFrag[site])):
                newfrag=-1
                for m1 in range(len(ms)):
                    if not (m1 in skip):
                        #if COMPARISON says the groups are the same, then add m1 to skip (b/c it's frag is already represented)
                        #otherwise (else), add first non-match to UFrag and do the comparison over again
                        #repeat until all ms indices are in skip (b/c all unique fragments are in UFrag[site]
                        matched=fragCompare(UFrag[site][frag],Forder[UFrag[site][frag]][site],m1,Forder[m1][site])
                        if matched: # matched == True
                            skip.append(m1)
                        else:
                            if newfrag==-1:
                                newfrag = m1
                                UFrag[site].append(newfrag)
                                skip.append(newfrag)
            ifrag=frag
        
    if debug:
        print("UFRAG")
        print(UFrag,"\n")
    
    ## (8) Print A formatted file of information to pass onto CRN (separate script)
    ## This also allows you to modify things by hand 
    fp=open(mcsout,'w')
    fp.write('# Maximum Common Substructure Search for Multisite Lambda Dynamics (JV 2022)\n')
    fp.write('# %d molecules processed\n\n' % (len(ms)))
    # (1) Print nsubs info
    fp.write("NSUBS")
    for site in range(nsites[refnum]):
        fp.write(" %d" % (len(UFrag[site])))
    fp.write("\n\n")
    fp.write("REFLIG %s\n\n" %(ms[refnum]))
    # (2) Print Core Atoms
    fp.write('CORE \n')
    for mol in range(len(ms)):
        fp.write("%s" % (ms[mol]))
        # print heavy atoms in the core
        for k in cores[refnum].keys():
            if mol == refnum:
                at=k
            else:
                at=cores[refnum][k][mol][0]
            fp.write(" %s" % (at))
            # print attached hydrogen atoms too
            at1=GetIndex(at,mol)
            for at2 in range(len(atoms[mol])):
                if bonds[mol][at1][at2] == 1 and xyzs[mol][atoms[mol][at2]]['Elem'] == 'H':
                    repeat=0    # avoid write core_attached_H repeatly
                    for n in range(nsites[mol]):
                        if atoms[mol][at2] in Fatoms[mol][n]:
                            repeat+=1
                        else:
                            continue
                    if repeat == 0:
                        fp.write(" %s" % (atoms[mol][at2]))
        fp.write("\n")
    fp.write("\n")
    # (3) Print Anchor Atoms
    fp.write('ANCHOR ATOMS\n')
    for mol in range(len(ms)):
        fp.write("%s" % (ms[mol]))
        for site in range(nsites[refnum]):
            aatom=Forder[mol][site]
            chk=0
            for c in range(len(aatom)):
                if aatom[c] == '+':
                    # we have a "merged" Aatom
                    chk=1
            if chk == 1:
                fp.write(" %s" % ("DUM"))
            else:
                fp.write(" %s" % (aatom))
        fp.write("\n")
    fp.write("\n")
    # (4) Print Fragments
    for site in range(nsites[refnum]):
        fp.write('SITE '+str(site+1)+' FRAGMENTS\n')
        for uf in UFrag[site]:
            fp.write("%s" % (ms[uf]))
            for at in Aatoms[uf][Forder[uf][site]]:
                fp.write(" %s" % (at))
            fp.write("\n")
        fp.write("\n")
    
    fp.write("END  \n")
    fp.close()
    
    
    # finished
    return ms[refnum]