# -*- coding: utf-8 -*-
# Pairwise distance matrix generator v2.1
# mkozic@chem.pmf.hr

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
from matplotlib import colors
import time

###############################################################################
   
file  = title = "dir/prot.pdb"
savefig = str("figs/prot")
params = [100, 5, 0, 40]

binary = False
cutoff = 12

###############################################################################
###############################################################################

start_time = time.time()  
  
xy_major = params[0]
xy_minor = params[1]
vmin = params[2]
vmax = params[3]
########################-Reader-###########################################

data = pd.DataFrame(open(file), dtype = "string")
data = data.replace(" A1"," A 1",regex=True)         #expand for other chains?
data = data.replace(" B1"," B 1",regex=True)         #expand for other chains?
data = data.replace(" C1"," C 1",regex=True)         #expand for other chains?
data = data.replace(" ",",",regex=True)
print(".")
data = data.replace(",,,,,",",",regex=True)
print("..")
data = data.replace(",,,,",",",regex=True)
print("...")
data = data.replace(",,,",",",regex=True)
print("....")
data = data.replace(",,",",",regex=True)
print(".....")
print("Finishing import...")


import_time = np.around((time.time() - start_time), decimals = 3)
print("--- %s seconds ---" %import_time)
###############################################################################


data = data.squeeze()
data = data.str.split(",", expand = True )     
data = data.loc[data[0] == "ATOM" ]
data = data.loc[data[2] == "CA" ]


data_CA = pd.DataFrame(data = None, dtype = "float64", 
                      index = np.arange(0,0,1), 
                      columns = np.arange(0,0,1))

at_typ = 2
chain_index = 4
resid_index = 5
i = 0
chain_counter = 0
residue = data.iloc[0,resid_index]
while i < len(data):
    try:
        line = data.iloc[i,:10]
        if line[at_typ] == "CA":
            residue = line[resid_index]
            chain = line[chain_index]     
            if i > 1 : 
                chain_before =  data.iloc[i-1,chain_index]
                residue_before = data.iloc[i-1, resid_index]
                if chain != chain_before:
                    chain_counter = int(chain_counter) + int(residue_before)
                    print("New chain")
            line.iloc[0] = str(int(residue) + chain_counter)    
            
            data_CA = pd.concat([data_CA, line],
                               axis = 1, join = "outer", ignore_index=True )
            print("Reading line:", i, "of", file, "; Found CA :D")
            i = i + 1
    except TypeError:
        i = i + 1
        pass
print("---Finished reading---")
data_CA = data_CA.transpose()


reindex_time = np.around((time.time() - start_time), decimals = 3)
print("--- %s seconds ---" %reindex_time)
#######################################################################
start_time2 = time.time()

residues = max(np.array(data_CA[0], dtype = "int"))

M_x = pd.DataFrame(data = None, dtype = "float64", 
                      columns = np.arange(0,int(residues) + 1,1), 
                      index = np.arange(0,int(residues) + 1,1))

M_y = pd.DataFrame(data = None, dtype = "float64", 
                      columns = np.arange(0,int(residues) + 1,1), 
                      index = np.arange(0,int(residues) + 1,1))

M_z = pd.DataFrame(data = None, dtype = "float64", 
                      columns = np.arange(0,int(residues) + 1,1), 
                      index = np.arange(0,int(residues) + 1,1))



a = 6
b = 7
c = 8 
i = 0   
while i < len(data_CA):
    res = int(data_CA.iloc[i,0])
    M_x.iloc[res,:] = float(data_CA.iloc[i,a])
    M_y.iloc[res,:] = float(data_CA.iloc[i,b])
    M_z.iloc[res,:] = float(data_CA.iloc[i,c])
    i = i + 1
    
M_x = M_x - M_x.transpose()
M_y = M_y - M_y.transpose()
M_z = M_z - M_z.transpose()

M_x = np.multiply(M_x, M_x)
M_y = np.multiply(M_y, M_y)
M_z = np.multiply(M_z, M_z)
    
M_x = M_x.fillna(0)
M_y = M_y.fillna(0)
M_z = M_z.fillna(0)

PWDM_squared = M_x + M_y + M_z  
PWDM = np.power(PWDM_squared, 0.5)
    
if binary == True:
    vmax = 1
    i = 0
    j = 0
    while i < len(PWDM):
        while j < len(PWDM):
            print(i,j)
            if PWDM.iloc[i,j] == 0:
                 PWDM.iloc[i,j] = 0
            elif PWDM.iloc[i,j] < cutoff:
                PWDM.iloc[i,j] = 1
            else:
                PWDM.iloc[i,j] = 0
            j = j + 1
        i = i + 1
        j = 0

calc_time = np.around((time.time() - start_time2), decimals = 3)
print("calc time is --- %s seconds ---" %calc_time)
###############################################################################
from matplotlib.ticker import (MultipleLocator)
import matplotlib.colors as mcolors
import matplotlib.cm as cm
# Define the modified colormap
viridis = cm.get_cmap('viridis', 256)
cm_colors = viridis(np.linspace(0, 1, 256))
cm_colors[250:] = mcolors.to_rgba('red')
cm_colors[:5] = mcolors.to_rgba('blue')
viridis_capped = mcolors.LinearSegmentedColormap.from_list('viridis_capped', cm_colors)

###############################################################################

label = "Residue number"
map_title = str(title)
if binary == True:
    savefig = str(savefig + "_bnry")
    cmap = colors.ListedColormap(['white', 'black'])
else : cmap = viridis_capped

ax = plt.subplot(111)
plt.imshow(PWDM, cmap = cmap, vmin = vmin, vmax = vmax, aspect = "equal")
ax.invert_yaxis()
plt.xticks(np.arange(0, int(residues), step= xy_major), fontsize = 12)                 
plt.yticks(np.arange(0, int(residues), step=xy_major), fontsize = 12) 
ax.yaxis.set_minor_locator(MultipleLocator(xy_minor)) 
ax.xaxis.set_minor_locator(MultipleLocator(xy_minor))

if binary == False:
    plt.colorbar( label = "Distance / Å", fraction=0.029, pad=0.028)
elif binary == True:
    cbar = plt.colorbar(ticks=[0,1], fraction=0.029, pad=0.04,)
    cbar.ax.set_yticklabels(['No', 'Yes'])
    cbar.set_label( label = str("Distance less than " +str(cutoff) + " Å"), labelpad=-13)
plt.xlabel (label)
plt.ylabel (label)
ax.set_title(map_title)
plt.savefig(savefig, dpi=800, bbox_inches='tight')
print("Map", title, " saved to:", savefig)
