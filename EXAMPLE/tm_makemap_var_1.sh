#!/bin/bash

### Matej Kožić | mkozic@chem.pmf.hr | 2022.11.25 15:00
### Bash script for creating trajectory heatmaps from a .csv matrix
#############################################################################
#INPUTS:

saved_csv="var_1.csv"       #load the matrix as a previously saved .csv
savefig="TrajMap_var_1.png"		#savename of the heatmap
title="Variant 1 simulation"			#title of the map
params="25,5,20,5,0,7,308,0,0" #parameters for heatmap. Adjust to fit .
#	params:
#	 [x_major, x_minor, y_major, y_minor, vmin, vmax, residues, cmap, aspect]
#	default for 500 frame 500 aa simulation:
#		params="25,5,20,5,0,7,number of residues,0,0"
#	x/y major/minor : ticks for the axes
#	vmin / vmax : min and max values of the colorbar
#	cmap : colormap, 0 for linear and 1 for divergent (used for diff maps)
#	aspect : 0 for auto 1 for square; aspect of the pixel

#############################################################################
#############################################################################

python3 TM.py << INPUTS
csv2matrix
$saved_csv
matrix2map
$savefig
$title
$params
q
INPUTS

echo If successful trajectory map $title saved to $savefig

