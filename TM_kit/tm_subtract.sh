#!/bin/bash

### Matej Kožić | mkozic@chem.pmf.hr | 2023.11.30
### Bash script for creating difference matrix
#############################################################################
#INPUTS:

### Subtracts matrix B from matrix A: difference = A - B

matrix_A="tm_testkit/test.gro"      #path to matrix A
matrix_B="tm_testkit/test.xtc"  #path to matrix B
saved_csv="savedcsv.csv"            #name of the created .csv matrix

#############################################################################
#############################################################################

python3 TM.py << INPUTS
7
$matrix_A
$matrix_B
$saved_csv
q
INPUTS

echo If successful, subtracted: A-B and saved the difference csv as $saved_csv

