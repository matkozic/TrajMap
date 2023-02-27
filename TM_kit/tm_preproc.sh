#!/bin/bash

### Matej Kožić | mkozic@chem.pmf.hr | 2022.11.25 15:00
### Bash script for preprocessing trajectories into a .csv matrix using TM.py
#############################################################################
#INPUTS:

topology="tm_testkit/test.gro"      #your topology file (gro, pdb, prmtop...)
trajectories="tm_testkit/test.xtc"  #your trajectory file (xtc, nc)
stride="1"			     #trajectory reading stride
saved_pdb="savedpdb.pdb"	     #name of the created pdb
saved_csv="savedcsv.csv"            #name of the created .csv matrix
residues="249"                     #real number of residues (excluding missing)

#############################################################################
#############################################################################

python3 TM.py << INPUTS
traj2pdb
$topology
$trajectories
$stride
$saved_pdb
pdb2csv
$saved_pdb
$saved_csv
$residues
q
INPUTS

echo If successful trajectories processed and saved as $saved_csv

