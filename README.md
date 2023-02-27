# TrajMap
MD simulation trajectory visualization tool

######################################
TrajMap.py quick guide for Linux usage with bash scripts
######################################
### Matej Kožić | mkozic@chem.pmf.hr | 2023.02.27 | github.com/matkozic/TrajMap

The program can be used inside of Linux terminal with provided Bash scripts:
	- tm_preproc.sh
	- tm_makemap.sh
	- tm_graph.sh
Parameters of the scripts need to be adjusted, and then executed in order. 


### General Workflow:
	1) Preprocessing: (tm_preproc.sh)
		I) Convert trajectories to an aligned .pdb trajectory (traj2pdb)
		II) Convert .pdb to .csv (pdb2csv)
	2) Make a trajectory map: (tm_makemap.sh)
		III) Load matrix from .csv (csv2matrix)
		IV) Create a map from a loaded matrix (matrix2map)
	3) Make a shift graph: (tm_graph.sh)
		V) Calculate shift data from a loaded matrix (matrix2shift)
		VI) Plot the shift data (shift2graph)
		

##############################################################################

0) Prerequisites
Make sure Python 3 is installed, as well as modules:
	- Numpy
	- Pandas
	- Matplotlib
	- MdTraj
	
The modules can be installed as  following:
	pip install pandas
	pip install matplotlib
	pip install mdtraj
	
Put the directory TM_kit in the desired location. It should contain:
	- TM.py		Main Python script
	- tm_preproc.sh	Preprocessing Bash script
	- tm_makemap.sh	Map-making Bash script
	- tm_graph.sh		Graph-making Bash script
	- tm_quickguide.txt	(optional)
	- tm_testkit		(optional)
		- tm_test.sh	(optional) A script for testing installation
		- test.gro	(optional)
		- test.xtc	(optional)
		
The installation can be tested with a tm_testkit, by running tm_test.sh
	cd tm_testkit
	./tm_test.sh
The script has to be ran inside of the tm_testkit directory.
The script should run for under a few minutes.
Created files will be deleted automatically.
If the code executes and exits successfully, everything should work.


##############################################################################

1) Preprocessing (tm_preproc.sh)
Converts raw trajectories into a .csv matrix.
Edit the inputs of the script accordingly, and execute it.

Variable "stride" reffers to the reading step of trajectories.
 Optimal number of frames for a map is 500-1000.
  Stride should be adjusted accordingly. For an example:
  If the number of frames in a raw trajectory is 10 000, stride should be 10.
 
Variable "residues" reffers to the real number of residues, excluding missing.
 If residues go from 2 to 100, the number of residues is 99.

Firstly the trajectories are saved as an intermediate .pdb file.
After that, the .pdb is converted and saved into a .csv matrix.
The .csv matrix is the main output that is used in later steps.

##############################################################################

2) Making a map (tm_makemap.sh)
Converts a .csv matrix from the previous step to a trajectory map.
Edit the inputs of the script accordingly, and execute it.

The format of "params" variable is as following:
 - [x_major, x_minor, y_major, y_minor, vmin, vmax, residues, cmap, aspect]

	x/y major/minor reffer to respective ticks on respective axes
	vmin/vmax are respective values of the color scale
	residues is real number of residues
	cmap is the used colormap:
		0 for a linear colormap, for a single trajectory map
		1 for a divergent colormap, for difference maps
	aspect is for a pixel aspect, 0 for auto and 1 for square.
		0 (auto) is recommended in most cases
	
It is recomended to optimize parameters such as color scale resolution.
 Several iterations of parameter adjustments might be needed for best results.

When creatign a single trajectory map the colormap should be linear (0).
 A linear map: a smooth transition from dark to light values.
When creating a difference map the colormap should be divergent (1).
 A divergent cmap: blue for negatives, white for zeroes, and red for positives. 
 The values for color scale shold also have the negative part, e.g:
 	vmin/vmax = -/+8

The map can be saved as either .png or .jpg, depending on the extension.


##############################################################################

3) Making a graph (tm_graph.sh)
Converts a .csv matrix from the previous step to a shift graph.
Calculates the shift of a region of residues specified with "shift_params".
	(residue 1 to residue 2, in an interval from time 1 to time 2)

y/x_params reffer to the axes parameters:
	- [min_value, max_value, major_tick, minor_tick]
roll_avg variable is for the rolling average.

The full list of available colors can be found in matplotlib documentation.


##############################################################################

### Considerations and errors:
TrajMap currently supports only proteins with regular 4 atom backbones.
 (5 for residues, which are handled automatically)
 
Multiple chains are supported only if the residue counting is continuous.
 (meaning chain A is residues 1 to N and chain B is N+1 to 2N+2)
To use it with a multimere, one should adjust the numbering so its continuous.
Alternatively, one can create separate trajectory maps for each chain.
 (Requires exporting individual chain trajectories manually) 

The most common error is preprocessing failing or getting stuck in a loop.
 (It's possible to be stuck in a loop repeating step 1)
This could be caused by a wrong residue input, or a faulty .pdb intermediate.
 Intermediate .pdb should contain only four atoms, and "END" or "ENDMDL".
 Remarks, Crystal data, and Connect statements should be automatically removed; 
  but lines other than that could cause it to freeze or loop indefinetily.	  

When adjusting parameters it is keep in mind zero-based indexing.
 For example, when calculating a shift graph, in "shift_params":
 	Nanoseconds 1 to 500 should be refferenced as 0,499
 Additionally, if the error keeps occuring it might be due to preprocessing.
  Last few frames tend to be dropped due to stride rounding.
  To fix it, try adjusting to a shorter time e.g. 0,490

Additional detail, as well as documentation and the original publication,
can be accessed via GitHub repository and by respective links.


### Matej Kožić | mkozic@chem.pmf.hr | 2023.02.27 | github.com/matkozic/TrajMap
