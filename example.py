

""" TrajMap usage example """

import TrajMap_local as tm

###############################################################################
# Single TrajMap from a pdb file ##############################################

file = "data/apo_stride_5.pdb"
save = "data/apo_stride_5.csv"
residues = 289
tm.pdb2csv(file, save, residues)

file = "data/apo_stride_5.csv"
residues = 289
matrix = tm.csv2matrix(file, residues)

savefig = "figs/apo.png"
title = "apo 1/5 frames"
params = [25, 5, 20, 5, 0, 7, 289, 0, 0]
tm.matrix2map(matrix, savefig, title, params)

###############################################################################
# Single TrajMap from an existing csv #########################################

file = "data/wild.csv"
residues = 308
matrix = tm.csv2matrix(file, residues)

savefig = "figs/wild.png"
title = "wild type protein"
params = [25, 5, 20, 5, 0, 7, 308, 0, 0]
tm.matrix2map(matrix, savefig, title, params)

###############################################################################
# Difference TrajMap from csv files ###########################################

file = "data/mut.csv"
mut = tm.csv2matrix(file, residues)

file = "data/wild.csv"
wild = tm.csv2matrix(file, residues)

diff_matrix = mut - wild

savefig = "figs/mut_wild"
title = "mut - wild"
params = [25, 5, 20, 5, -7, 7, 308, 1, 0]
tm.matrix2map(diff_matrix, savefig, title, params)

###############################################################################
# Average shift of a region from a loaded matrix ##############################

params = [140, 150, 0, 500]
avg = tm.matrix2avg(mut, params)

savefig = "figs/shift.png"
title = "Average shift of region 140 to 150"
label = "mut"
y_params = [0, 7, 0.5, 0.1]
x_params = [0, 500, 50, 10]
avg_params = [15, "orchid"]
tm.avg2graph(avg, savefig, title, label,
             y_params, x_params, avg_params)

###############################################################################
# ProtMap from a pdb file #####################################################  

file = "data/split.pdb"
savefig = "figs/protmap_split.png"
title = "split native"
save = ["n", "-"]
load = ["n", "-"]
params = [20, 5, 0, 50, 308]
tm.protmap(file, savefig, title, save, load, params)  

###############################################################################
###############################################################################
