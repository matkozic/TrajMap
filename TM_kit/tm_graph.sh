#!/bin/bash

### Matej Kožić | mkozic@chem.pmf.hr | 2023.02.27 09:25
### Bash script for creating a shift graph from a .csv matrix
#############################################################################
#INPUTS:

saved_csv="csv.csv"  #load the matrix as a previously saved .csv
shift_params="1,40,0,99" # [residue1, residue2, time1, time2]
savefig="savefig" #savename for the graph
title="title"
label="label" #label of the graph legend
y_params="0,5,1,0.1" #y axis ticks
#	[min_value, max_value, major_tick, minor_tick]
x_params="0,100,50,10" #x axis ticks
#	[min_value, max_value, major_tick, minor_tick]
roll_avg="10" #rolling average value; 10/15/20 recomended
color="black" #color of the graph     

########################################################

python3 TM.py << INPUTS

csv2matrix
$saved_csv
matrix2shift
$shift_params
shift2graph
$savefig
$title
$label
$y_params
$x_params
$roll_avg
$color
q
INPUTS

echo If successful trajectories processed and saved as $saved_csv

