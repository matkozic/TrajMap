#! /bin/bash

python3 ../TM.py << INPUTS
0
test.gro
test.xtc
1
testpdb.pdb
1
testpdb.pdb
testcsv.csv
249
2
testcsv.csv
3
testfig.png
testtitle
25,5,10,1,0,7,249,0,0
4
170,175,0,49
5
testgraph.png
testgraph
residues 170-175
0,10,1,0.1
0,50,5,1
5
fuchsia
h
q
q
INPUTS
#rm testpdb.pdb testcsv.csv 
rm testpdb.pdb testcsv.csv testfig.png testgraph.png
echo If it exited without errors it probably works fine.

