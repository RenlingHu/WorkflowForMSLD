# Protein and Ligand Preparation

### This folder contains (1) workflow for protein and ligand preparation, (2) Python/Jupter Notebook scripts for aligning ligands according to the MCS

05/04/2023

Author: Renling Hu

E-mails: 22260237@zju.edu.cn

#### Note
 - LIC_align and multi_mol2_split are just part of ligand preparation
 - the overall workflow for protein and ligand preparation would be illustrated as follows

---
#### WorkFlow
1. download protein(so-called reference_complex) structure file(pdb) and ligand structure file(sdf or mol2)

2. optimizing reference_complex structure with **Schrödinger-"prepwizard"**
    - `${path_to_schrodinger}/utilities/prepwizard -c -fillsidechains -fillloops -s -propka_pH 7.4 ${protein}.pdb ${protein}_schrodinger.pdb`
    - fill sidechains and loops, set pKa to 7.4, etc.
    - **make sure** ligand naming is compatible with CHARMM

3. solvate the reference_complex and set the periodic boundary condition with **CHARMM** or **CHARMM-GUI**
    - then get protein.pdb, solv.pdb, ions.pdb and reference_ligand.pdb
    - **keep track** of box_size and pmegrid_size

4. align all ligands with reference_ligand by `./MCS_align.ipynb` (**not strictly**)
    - **make sure** they are in the same multi_sdf file
    - **make sure** to edit the input file name and output file name according to the system
    - **make sure** to assign the number of reference_ligand in the script
    - **make sure** "is it strictly aligned" = False

5. divide ligands into different groups according to MCS and the size of the perturbation R group
    - consider the MCS first is recommended
    - based on the same MCS, consider the size of the perturbation R group

6. within a group, accurately align all other ligands with the one which has the smallest perturbation R group by `./MCS_align.ipynb` (**strictly**)
    - rdkit and a series of relative packages are required
    - **make sure** they are in the same multi_sdf file
    - **make sure** to edit the input file name and output file name according to the system
    - **make sure** to assign the number of reference_ligand in the script
    - **make sure** "is it strictly aligned" = True, so that the MCS of all ligands have the same spatial position. This is **IMPORTANT** in the "Multiple Topology Model Creation" stage!!
    - one could also use **Schrödinger-"Ligand Alignment"**, but a **CAREFUL** check is required
    - it is recommended to check the structure in **Pymol** or **Schrödinger** after alignment

7. generate force field parameter files and topology files for each ligand
    - CgenFF for example here(https://cgenff.umaryland.edu/)
    - mol2 file as the input, and str file as the output

---
#### File information
`./MCS_align.ipynb`
 - script to align ligands together according to their MCS
 - rdkit and a series of relative packages are required


`./split_multi_mol2_file.py`
 - script to convert multi_mol2 file into single mol2 file
 - usage: `python split_multi_mol2_file.py -i ligands.mol`
