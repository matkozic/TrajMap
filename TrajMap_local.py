#TrajMap
#Author: Matej Kožić
#Email: mkozic@chem.pmf.hr
#Date: July 16, 2022
#Link: https://github.com/matkozic/TrajMap
###############################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)


def pdb2csv (file, savename, residues):

    save_name = savename
   
    print ("Importing" , file)
    data = pd.DataFrame(open(file), dtype = "string")

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

    data = data.squeeze()
    data = data.str.split(",", expand = True )

    graph_data = pd.DataFrame()
    limit = len(data)
    atom = "ATOM"
    end = "END"
    time = 0 
    N = 14
    CA = C = 12
    O = 16
    a = 6
    b = 7 
    c = 8 
    i = int(residues) * 4 + 3
    h = 1
     
    try:
        while i < limit :
            if atom in data.iloc[i,0]:
                
                residue = data.iloc[i,5]
                print("Calculating: Row: ", i ,"; Timestep: " , time,
                      "; Residue: ", residue, "; File:", file)

                if int(residue) == residues :
                    
                    x_i1 = float(data.iloc[i,a]) #N
                    y_i1 = float(data.iloc[i,b])
                    z_i1 = float(data.iloc[i,c])
                    
                    x_i2 = float(data.iloc[i+1,a]) #CA
                    y_i2 = float(data.iloc[i+1,b])
                    z_i2 = float(data.iloc[i+1,c])
                    
                    x_i3 = float(data.iloc[i+2,a]) #C
                    y_i3 = float(data.iloc[i+2,b])
                    z_i3 = float(data.iloc[i+2,c])
                    
                    x_i4 = float(data.iloc[i+3,a]) #O
                    y_i4 = float(data.iloc[i+3,b])
                    z_i4 = float(data.iloc[i+3,c])
                    
                    x_i5 = float(data.iloc[i+4,a]) #O
                    y_i5 = float(data.iloc[i+4,b])
                    z_i5 = float(data.iloc[i+4,c])
                
                    x_h1 = float(data.iloc[h,a])
                    y_h1 = float(data.iloc[h,b])
                    z_h1 = float(data.iloc[h,c])
                    
                    x_h2 = float(data.iloc[h+1,a])
                    y_h2 = float(data.iloc[h+1,b])
                    z_h2 = float(data.iloc[h+1,c])
                    
                    x_h3 = float(data.iloc[h+2,a])
                    y_h3 = float(data.iloc[h+2,b])
                    z_h3 = float(data.iloc[h+2,c])
                    
                    x_h4 = float(data.iloc[h+3,a])
                    y_h4 = float(data.iloc[h+3,b])
                    z_h4 = float(data.iloc[h+3,c])
                    
                    x_h5 = float(data.iloc[h+4,a]) 
                    y_h5 = float(data.iloc[h+4,b])
                    z_h5 = float(data.iloc[h+4,c])
                
                    
                    x_cm_i = (x_i1*N + x_i2*CA + x_i3*C + x_i4*O + x_i5*O) / (O+N+CA+C+O)
                    y_cm_i = (y_i1*N + y_i2*CA + y_i3*C + y_i4*O + y_i5*O) / (O+N+CA+C+O)
                    z_cm_i = (z_i1*N + z_i2*CA + z_i3*C + z_i4*O + z_i5*O) / (O+N+CA+C+O)
                    
                    x_cm_h = (x_h1*N + x_h2*CA + x_h3*C + x_h4*O + x_h5*O) / (O+N+CA+C+O) 
                    y_cm_h = (y_h1*N + y_h2*CA + y_h3*C + y_h4*O + y_h5*O) / (O+N+CA+C+O)
                    z_cm_h = (z_h1*N + z_h2*CA + z_h3*C + z_h4*O + z_h5*O) / (O+N+CA+C+O)
                
                    value =  ( ((x_cm_i - x_cm_h )**2 +
                                (y_cm_i - y_cm_h)**2 + 
                                (z_cm_i - z_cm_h)**2 ) )**(1/2)  
                   
                    i = i + 5
                    h = h + 5
                else:
                    x_i1 = float(data.iloc[i,a]) #N
                    y_i1 = float(data.iloc[i,b])
                    z_i1 = float(data.iloc[i,c])
                    
                    x_i2 = float(data.iloc[i+1,a]) #CA
                    y_i2 = float(data.iloc[i+1,b])
                    z_i2 = float(data.iloc[i+1,c])
                    
                    x_i3 = float(data.iloc[i+2,a]) #C
                    y_i3 = float(data.iloc[i+2,b])
                    z_i3 = float(data.iloc[i+2,c])
                    
                    x_i4 = float(data.iloc[i+3,a]) #O
                    y_i4 = float(data.iloc[i+3,b])
                    z_i4 = float(data.iloc[i+3,c])
                
                    x_h1 = float(data.iloc[h,a])
                    y_h1 = float(data.iloc[h,b])
                    z_h1 = float(data.iloc[h,c])
                    
                    x_h2 = float(data.iloc[h+1,a])
                    y_h2 = float(data.iloc[h+1,b])
                    z_h2 = float(data.iloc[h+1,c])
                    
                    x_h3 = float(data.iloc[h+2,a])
                    y_h3 = float(data.iloc[h+2,b])
                    z_h3 = float(data.iloc[h+2,c])
                    
                    x_h4 = float(data.iloc[h+3,a])
                    y_h4 = float(data.iloc[h+3,b])
                    z_h4 = float(data.iloc[h+3,c])
                    
                    x_cm_i = (x_i1*N + x_i2*CA + x_i3*C + x_i4*O) / (O+N+CA+C)
                    y_cm_i = (y_i1*N + y_i2*CA + y_i3*C + y_i4*O) / (O+N+CA+C)
                    z_cm_i = (z_i1*N + z_i2*CA + z_i3*C + z_i4*O) / (O+N+CA+C)
                    
                    x_cm_h = (x_h1*N + x_h2*CA + x_h3*C + x_h4*O) / (O+N+CA+C) 
                    y_cm_h = (y_h1*N + y_h2*CA + y_h3*C + y_h4*O) / (O+N+CA+C)
                    z_cm_h = (z_h1*N + z_h2*CA + z_h3*C + z_h4*O) / (O+N+CA+C)
                
                    value =  ( ((x_cm_i - x_cm_h )**2 +
                                (y_cm_i - y_cm_h)**2 + 
                                (z_cm_i - z_cm_h)**2 ) )**(1/2)  
                    h = h + 4
                    i = i + 4
                coords= pd.DataFrame([residue, time, value])
                coords = coords.transpose()
                
                graph_data = pd.concat([graph_data,coords],axis=0,join="outer")
            
            elif end in str(data.iloc[i,0]) :
                print("---------Finished step ", time, "----------")
                time = time + 1
                i = i + 1
                h = 1
                
    except IndexError : pass      
    graph_data.to_csv(save_name)
    print("Saved coordinates from pdb to: ",save_name)

###############################################################################

def csv2matrix (file, residues):

    graph_data = open(file)
    graph_data = pd.Series(graph_data)
    graph_data = graph_data.str.split(",", expand = True)
    graph_data = graph_data.iloc[:,1:]
    
    time = graph_data.iloc[len(graph_data)-10,1]
    graph_data = graph_data.set_index(np.arange(0, 
                                     len(graph_data), step = 1))
    matrix = pd.DataFrame(data=None, 
                          index = np.arange(int(time)+1, 0, -1), 
                          columns = np.arange(0,int(residues) + 1,1))
    k = 0                     
    try:
        while k <= len(graph_data) :
            x = int(graph_data.iloc[k,0])
            y = int(graph_data.iloc[k,1])
            z = float(graph_data.iloc[k,2])
            matrix.iloc[y,x] = z
            print("Transcribing row ", k, "of", file)
            k = k + 1
    except IndexError: pass
    
    matrix = np.array(matrix.iloc[:,1:], dtype = "float64")
    matrix = np.transpose(matrix)
    return(matrix)

###############################################################################

# params = [x_major, x_minor, y_major, y_minor, vmin, vmax, residues, cmap, aspect]

def matrix2map (matrix, savefig, title, params):
    
    print("Making the TrajMap...")
    
    x_major = params[0]
    x_minor = params[1]
    y_major = params[2]
    y_minor = params[3]
    vmin = params[4]
    vmax = params[5]
    residues = params[6]
    cmap = params[7]
    aspect = params[8]

    if cmap == 0 : 
        cmap = "magma"
    elif cmap == 1 : 
        cmap = "seismic"
        
    if aspect == 0 :
        aspect = "auto"
    elif aspect == 1 :
        aspect = "equal"

    ylabel = "Residue"
    xlabel = "Frame"
    size = 800
    frames = len(matrix[0])
    map_title = str(str("TrajMap") + " " + str(title) )

    ax = plt.subplot(111)
    plt.imshow(matrix, cmap = cmap, vmin = vmin, vmax = vmax, aspect = aspect)
    ax.set_title(map_title)

    ax.invert_yaxis()
    plt.xticks(np.arange(0, frames + 1, step= x_major), fontsize = 7)                 
    plt.yticks(np.arange(0, residues, step=y_major), fontsize = 7) 
    ax.yaxis.set_minor_locator(MultipleLocator(y_minor)) 
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    plt.ylabel(ylabel)
    plt.xlabel (xlabel)
    plt.colorbar( label = "Shift / A", fraction=0.029, pad=0.028)
    plt.savefig(savefig, dpi = size, bbox_inches='tight')
    print("TrajMap saved to:", savefig)

###############################################################################

# params = [res1, res2, time1, time2]

def matrix2avg (matrix, params):                                                           
    
    res1 = params[0]
    res2 = params[1]
    time1 = params[2]
    time2 = params[3]

    output = pd.DataFrame(data = np.arange(time1,time2,1))    

    y = matrix[res1:res2,time1:time2]
    i = 0
    while i < len(output) : 
        output.iloc[i] = np.average(y[:,i])
        i = i + 1
    return output

###############################################################################

# y_params = [y_min, y_max, y_step, y_tick]
# x_params = [x_min, x_max, x_step, x_tick]
# avg_params = [rolling_average, color]

def avg2graph (avg, savefig, title, label, y_params, x_params, avg_params):
    
    y_min = y_params[0]
    y_max = y_params[1]
    y_step = y_params[2]
    y_tick = y_params[3]
    
    x_min = x_params[0]
    x_max = x_params[1]
    x_step = x_params[2]
    x_tick = x_params[3]
    
    roll_avg = avg_params[0]
    color = avg_params[1]
               
    y = avg    
    size = 800
    ya = y.rolling(roll_avg).mean()
    ax = plt.subplot(111)
    plt.plot(y, linewidth=0.3, color=color, label = label)
    ax.legend(loc='upper left',
              fancybox=True, shadow=True, prop={'size': 7.5})
    plt.plot(ya, linewidth=1, color=color) 

    ax.set_xlim(x_min, x_max) 
    ax.set_ylim(y_min) 
    plt.yticks(np.arange(y_min, y_max, step=y_step)) 
    plt.xticks(np.arange(x_min, x_max + x_step, step=x_step)) 
    plt.title(title)
    plt.xlabel("Frame")
    plt.ylabel ("Shift / A")
    ax.yaxis.set_minor_locator(MultipleLocator(y_tick)) 
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick))

    plt.savefig(savefig, dpi= size)

###############################################################################

# params = [xy_major, xy_minor, vmin, vmax, residues]
# save = ["y"/"n", savename]
# load = ["y"/"n", loadname]

def protmap (file, savefig, title, save, load, params):
    
    xy_major = params[0]
    xy_minor = params[1]
    vmin = params[2]
    vmax = params[3]
    residues = params[4]
    
    if load[0] == "n":
        
        data = pd.DataFrame(open(file), dtype = "string")
        
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
        
        data = data.squeeze()
        data = data.str.split(",", expand = True )
        
        
        matrix = pd.DataFrame(data = None, dtype = "float64", 
                              columns = np.arange(0,int(residues) + 1,1), 
                              index = np.arange(0,int(residues) + 1,1))
        
        data_CA = pd.DataFrame(data = None, dtype = "float64", 
                              index = np.arange(0,int(residues) + 1,1), 
                              columns = np.arange(0,0,1))
        
        at_typ = 2
        i = 0
        
        while i < len(data):
            try:
                line = data.iloc[i,:]
                if line[at_typ] != "CA":
                    print("Reading line:", i, "of", file)
                    i = i + 1
                elif line[at_typ] == "CA":
                    print("Reading line:", i, "of", file, "; Found CA :D")
                    data_CA = pd.concat([data_CA, line],
                                       axis = 1, join = "outer" )
                    i = i + 1
            except TypeError:
                i = i + 1
                pass
        print("---Finished reading---")
        data_CA = data_CA.transpose()
        data_CA = data_CA.reset_index()
        
        
        at_typ = 2
        aa_num = 5
        a = 6
        b = 7
        c = 8
        
        i = 0
        j = 0
        
        while i < len(data_CA):
            line_i = data_CA.iloc[i,:]
            x_i = float(line_i[a])
            y_i = float(line_i[b])
            z_i = float(line_i[c])
            while j < len(data_CA):
                print("Calculating distance of residue", i, "to", j, "for", file)
                line_j = data_CA.iloc[j,:]
                x_j = float(line_j[a])
                y_j = float(line_j[b])
                z_j = float(line_j[c]) 
                value =  ( ((x_i - x_j )**2 +
                            (y_i - y_j)**2 + 
                            (z_i - z_j)**2 ) )**(1/2)
                
                aa_i = int(line_i[aa_num])
                aa_j = int(line_j[aa_num])
                matrix.iloc[aa_i, aa_j] = value
                matrix.iloc[aa_j, aa_i] = value
                
                j = j + 1
            i = i + 1
            j = i
        
        matrix = matrix.iloc[1:,:]
        matrix = matrix.iloc[:,1:]
        
        if save[0] == "y" :
            matrix.to_csv(str(save[1] + ".csv"))
            print("Saved to: ", save[1])
        else:
            print("NOT SAVING")
    else:
        matrix = pd.DataFrame(open(load[1]))
        
    label = "Residue"
    map_title = str(str("ProtMap ") + str(title) )
    
    ax = plt.subplot(111)
    plt.imshow(matrix, cmap = "viridis", vmin = vmin, vmax = vmax, aspect = "equal")
    ax.invert_yaxis()
    plt.xticks(np.arange(0, int(residues), step= xy_major), fontsize = 7)                 
    plt.yticks(np.arange(0, int(residues), step=xy_major), fontsize = 7) 
    ax.yaxis.set_minor_locator(MultipleLocator(xy_minor)) 
    ax.xaxis.set_minor_locator(MultipleLocator(xy_minor))
    plt.colorbar( label = "Distance / A", fraction=0.029, pad=0.028)
    plt.xlabel (label)
    #plt.ylabel (label)
    ax.set_title(map_title)
    plt.savefig(savefig, dpi = 800, bbox_inches='tight')
    print("ProtMap", title, " saved to:", savefig)

###############################################################################
