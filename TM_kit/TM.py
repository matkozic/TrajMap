###############################################################################
# TrajMap_v_f1 ##### Matej Kožić | mkozic@chem.pmf.hr | 2022.11.25 15:00
# MatejKozic, 2023.03.11
###############################################################################

import mdtraj as mdt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
import matplotlib.colors as mcolors
import matplotlib.cm as cm

###############################################################################
# Additional colormaps
###############################################################################

# Viridis capped
viridis = cm.get_cmap('viridis', 256)
colors = viridis(np.linspace(0, 1, 256))
colors[250:] = mcolors.to_rgba('red')
colors[:5] = mcolors.to_rgba('blue')
viridis_capped = mcolors.LinearSegmentedColormap.from_list('viridis_capped', colors)

# Seismic capped
seismic = cm.get_cmap('seismic', 256)
colors = seismic(np.linspace(0, 1, 256))
colors[250:] = mcolors.to_rgba('fuchsia')
colors[:5] = mcolors.to_rgba('cyan')
seismic_capped = mcolors.LinearSegmentedColormap.from_list('seismic_capped', colors)

# ChatGPT AI generated really bad random unusable noice colormap
def random_segmented_cmap(n):
    c = np.random.rand(n+1,3)
    return mcolors.LinearSegmentedColormap.from_list('r', np.repeat(c,2,0))

# Viridis segmented
viridis_segmented = cm.get_cmap('viridis', 8)

###############################################################################
# Functions
###############################################################################

def traj2pdb(topology, trajectories, stride, savename) :

    if ".pdb" not in savename:
        savename = str (savename + ".pdb" )
    
    print("Loading the trajectory file...")
    load_traj = mdt.load(trajectories,top = topology, stride = stride)
    print(load_traj)
    print("Selecting the backbone...")
    select = load_traj.topology.select("backbone")
    load_traj = load_traj.atom_slice(select)
    print("Aligning...")
    print(load_traj)
    load_traj = load_traj.superpose(load_traj,0)
    print("Saving to ", savename)
    load_traj.save_pdb(savename)

###############################################################################

def pdb2csv (file, savename, residues):

    print ("Importing" , file)
    data = pd.DataFrame(open(file), dtype = "string")
    data = data.replace(" A1"," A 1",regex=True)
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
    
    data = data.replace( pd.NA ,"0",regex=True)
    data = data.loc[data[2] != "OXT" ]
    data = data.loc[data[2] != "OT2" ]
    data = data.loc[data[0] != "TER" ]
    data = data.loc[data[0] != "MODEL" ]
    data = data.loc[data[0] != "REMARK" ]
    data = data.loc[data[0] != "CONNECT" ]
    data = data.loc[data[0] != "CONECT" ]

    graph_data = pd.DataFrame()
    graph_data_temp = pd.DataFrame()
    concat_counter = 0
    limit = len(data)
    atom = "ATOM"
    end = "END"
    time = 0 
    N = 14
    CA = C = 12
    O =16
    a = 6
    b = 7 
    c = 8 
    i = int(residues) * 4 + 2
    h = 1
     
    try:
        while i < limit :
            if atom in data.iloc[i,0]:
                
                residue = data.iloc[i,5]
                print("Calculating: Row: ", i ,"; Timestep: " , time,
                      "; Residue: ", residue, "; File:", file)
                    
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
                
                graph_data_temp = pd.concat([graph_data_temp,coords],axis=0,join="outer")

                if concat_counter == 5:
                    graph_data = pd.concat([graph_data,graph_data_temp],axis=0,join="outer")
                    graph_data_temp = pd.DataFrame()
                    concat_counter = 0
            
            elif end in str(data.iloc[i,0]) :
                print("---------Finished step ", time, "----------")
                time = time + 1
                i = i + 1
                h = 1
                concat_counter = concat_counter + 1      
    except IndexError :
        pass      
    
    graph_data = pd.concat([graph_data,graph_data_temp],axis=0,join="outer")  
    
    max_resid = int(max(np.array(graph_data[0], dtype = "int")))
    
    time = graph_data.iloc[len(graph_data)-10,1]
    
    graph_data = graph_data.set_index(np.arange(0, 
                                     len(graph_data), step = 1))
    matrix = pd.DataFrame(data=None, 
                          index = np.arange(0, int(time)+1, 1), 
                          columns = np.arange(0,max_resid + 1,1))
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
    
    matrix = matrix.transpose()
    matrix = matrix.fillna(0)
    matrix.to_csv(savename, index = False)
    
    print("Saved coordinates from pdb to: ",savename)

###############################################################################

def csv2matrix (file):
    matrix = pd.read_csv(file, index_col = None)
    return(matrix)

###############################################################################

# params = [x_major, x_minor, y_major, y_minor, vmin, vmax, residues, cmap, aspect]

def matrix2map (matrix, savefig, title, params):
    
    print("Making the TrajMap...")
    
    matrix = np.array(matrix)
    
    x_major = params[0]
    x_minor = params[1]
    y_major = params[2]
    y_minor = params[3]
    vmin = params[4]
    vmax = params[5]
    residues = params[6]
    cmap = params[7]
    aspect = params[8]

    if cmap == 0 : cmap = "magma"
    elif cmap == 1 :  cmap = "seismic"
    elif cmap == 2 : cmap = seismic_capped
    elif cmap == 3 : cmap = viridis_capped
    elif cmap == 4 : cmap = random_segmented_cmap(np.random.randint(1, 99))
    elif cmap == 5 : cmap = viridis_segmented
    elif cmap == 6 : cmap = "Greys"
    elif cmap == 7 : cmap = "turbo"
     
    if aspect == 0 : aspect = "auto"
    elif aspect == 1 : aspect = "equal"

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
    plt.colorbar( label = "Shift / Å", fraction=0.029, pad=0.028)
    plt.savefig(savefig, dpi = size, bbox_inches='tight')
    print("TrajMap saved to:", savefig)

###############################################################################

# params = [res1, res2, time1, time2]

def matrix2shift (matrix, params):                                                           
    
    matrix_array = np.array(matrix)
    
    res1 = params[0]
    res2 = params[1]
    time1 = params[2]
    time2 = params[3]
    
    output = pd.DataFrame(data = np.arange(time1,time2,1))    
    
    if res1 == res2 :
        output = matrix_array[res1]
    else:
        y = matrix_array[res1:res2,time1:time2]
        
        i = 0
        while i < len(output) : 
            output.iloc[i] = np.average(y[:,i])
            i = i + 1
    return output

###############################################################################

# y_params = [y_min, y_max, y_step, y_tick]
# x_params = [x_min, x_max, x_step, x_tick]

def shift2graph (shift_data, savefig, title, label, y_params, x_params, roll_avg, color):
    
    y_min = y_params[0]
    y_max = y_params[1]
    y_step = y_params[2]
    y_tick = y_params[3]
    
    x_min = x_params[0]
    x_max = x_params[1]
    x_step = x_params[2]
    x_tick = x_params[3]
                    
    size = 800
    shift_data = pd.DataFrame(shift_data, dtype = "float64")
    roll_avg = shift_data.rolling(roll_avg).mean()
    ax = plt.subplot(111)
    plt.plot(shift_data, linewidth=0.3, color=color, label = label)
    ax.legend(loc='upper left',
              fancybox=True, shadow=True, prop={'size': 7.5})
    plt.plot(roll_avg, linewidth=1, color=color) 

    ax.set_xlim(x_min, x_max) 
    ax.set_ylim(y_min) 
    plt.yticks(np.arange(y_min, y_max, step=y_step)) 
    plt.xticks(np.arange(x_min, x_max + x_step, step=x_step)) 
    plt.title(title)
    plt.xlabel("Frame")
    plt.ylabel ("Shift / Å")
    ax.yaxis.set_minor_locator(MultipleLocator(y_tick)) 
    ax.xaxis.set_minor_locator(MultipleLocator(x_tick))

    plt.savefig(savefig, dpi= size)


###############################################################################
# Main program loop
###############################################################################

main = True
entry_counter = 0
matrix_loaded = False
shift_loaded = False

while main == True:
    
    print ("")
    if entry_counter == 0:
        print("                     ___         ") 
        print("         ___        /__/\        ")
        print("        /  /\      |  |::\       ")
        print("       /  /:/      |  |:|:\      ")
        print("      /  /:/     __|__|:|\:\     ")
        print("     /  /::\    /__/::::| \:\    ")
        print("    /__/:/\:\   \  \:\~~\__\/    ")
        print("    \__\/  \:\   \  \:\          ")
        print("         \  \:\   \  \:\         ") 
        print("          \__\/    \  \:\        ") 
        print("                    \__\/        ")
        print("                                 ")
        print("      = Welcome to TrajMap =     ")
    
    print ("_____________________________________________________")
    print("Choose an action by entering its coresponding number.")
    print("List of actions: ")
    print(" 0: traj2pdb | 1: pdb2csv")
    print(" 2: csv2matrix | 3: matrix2map")
    print(" 4: matrix2shift | 5: shift2graph ")
    print(" 6: matrix_average | 7: matrix_substract ")
    print (" h: help | q: quit")
    user_input = input("Enter your input: ")
    
    ########## Main command branch ##############
    if user_input == "0" or user_input == "traj2pdb": 
        print (" = traj2pdb(topology, trajectories, stride, savename) =")
        print(" Ideal number of frames is 500 to 1000, and it's determined by loading stride. ")
        print(" Accepted topologies and trajectories are: pdb + xtc ; gro + xtc")
        print(" If loading from a .pdb trajectory enter it as both topology and trajectory.")
        topology = str(input("Input the topology file: "))
        trajectories = str(input("Input the trajectory file: "))
        stride = int(input("Input loading stride: "))
        savename = str(input("Input the save name: "))
        
        traj2pdb(topology, trajectories, stride, savename)
        print("Great success!")
    
    elif user_input == "1" or user_input == "pdb2csv":
        print (" = pdb2csv(file, savename, residues) = ")
        file = str(input("Input the path to the .pdb file: "))
        savename = str(input("Input the savename for .csv file: "))
        residues = int(input("Input the number of residues: "))
        pdb2csv(file, savename, residues)
        print("Great success!")
        print("From ", file, "made .csv and saved it to ", savename)
        
    elif user_input == "2" or user_input == "csv2matrix":
        print (" = csv2matrix(file) = ")
        file = str(input("Input the path to the .csv file: "))
        loaded_matrix = csv2matrix(file)
        matrix_loaded = True
        print("Great success!")
        print("Loaded the matrix from file ", file)
              
    elif user_input == "3" or user_input == "matrix2map":
        print (" = matrix2map(matrix, savefig, title, params) =")
        
        if matrix_loaded == False:
            input(" !! Error, matrix not loaded. Load the matrix with csv2matrix. Press enter to return.")
        elif matrix_loaded == True:
            print ("Loaded matrix = ", file)
            savefig = str(input("Enter the save name: "))
            title = str(input("Enter the title: "))
            print(" params: [x_major, x_minor, y_major, y_minor, vmin, vmax, residues, cmap, aspect] ")
            print("  cmap- 0: linear | 1: divergent")
            print("  aspect- 0: auto | 1: square")
            params = np.fromstring(input("Enter the parameter list: "),
                                   sep = ",", dtype = "int")
            
            matrix2map(loaded_matrix, savefig, title, params)  
            print("Great success!")
            print("Map from file ", file, "saved to ", savefig)
                  
        
    elif user_input == "4" or user_input == "matrix2shift":
        print (" = matrix2shift (matrix, params) =")
        if matrix_loaded == False: 
            input(" !! Error, matrix not loaded. Load the matrix with csv2matrix. Press enter to return.")
        elif matrix_loaded == True:
            print(" Loads the data for plotting a graph of shifts of a residue/region in time from time 1 to time 2.")
            print (" Region is defined as between and including residue 1 to residue 2")
            print (" For a single residue input it twice.")
            print(" Parameters: [residue 1, residue 2, time 1, time 2]")
            params =np.fromstring(input("Enter the parameter list: "),
                                   sep = ",", dtype = "int")
            shift_data = matrix2shift(loaded_matrix, params)
            
            shift_loaded = True
            
            print("Great success!")
            print("Shifts from residues ", params[0], " to ", params[1], " loaded.")
            
    elif user_input == "5" or user_input == "shift2graph":
        print (" = shift2graph(shift_data, savefig, title, label, y_params, x_params, roll_avg, color) =")
        if shift_loaded == True:
            print (" y_params = [y_min, y_max, y_step, y_tick]")
            print (" x_params = [x_min, x_max, x_step, x_tick]")
            savefig = input("Input the save name: ")
            title = input("Input the title: ")
            label = input("Input a legend lable: ")
            y_params =np.fromstring(input("Enter the y axis parameter list: "),
                                   sep = ",", dtype = "float64")
            x_params =np.fromstring(input("Enter the x axis parameter list: "),
                                   sep = ",", dtype = "float64")
            roll_avg = int(input(("Input the rolling average factor: ")))
            color = str(input("Input the color: "))
            
            shift2graph(shift_data, savefig, title, label, y_params, x_params, roll_avg, color)
            
            print("Great success!")
            print("Shifts graph saved to ", savefig)
            
        elif shift_loaded == False:
            input(" !! Error, shift data not loaded. Load the shfit data with matrix2shift. Press enter to return.")
        
    elif user_input == "6" or user_input == "matrix_average":
        print (" = matrix_average(number_of_matrices, files, savename) =")
        print(" Calculate the average of 2 or 3 matrices.")
        number_of_matrices = int(input("Input the number of matrices to average (2 or 3): "))
        matrix_1 = input("Input the path to the first matrix: ")
        matrix_2 = input("Input the path to the second matrix: ")
        save_name = input("Input the save name: ")
        
        if number_of_matrices == 2:
            matrix_1 = csv2matrix(matrix_1)
            matrix_2 = csv2matrix(matrix_2)
            
            avg_matrix = (matrix_1 + matrix_2) / 2
            
            avg_matrix.to_csv(save_name, index = False)
            
            print("Great success!")
            print("Matrix average saved to ", save_name)
            
            del avg_matrix
            
        elif number_of_matrices == 3:
            matrix_3 = input("Input the path to the third matrix: ")
            
            matrix_1 = csv2matrix(matrix_1)
            matrix_2 = csv2matrix(matrix_2)
            matrix_3 = csv2matrix(matrix_3)

            avg_matrix = (matrix_1 + matrix_2 + matrix_3) / 3
            
            avg_matrix.to_csv(save_name, index = False)
            
            print("Great success!")
            print("Matrix average saved to ", save_name)
            
            del avg_matrix
            
        else:
            input("Command not understood, press enter to return to the main menu.")
        
    elif user_input == "7" or user_input == "matrix_substract":
        print (" = matrix_substract =")
        print("Substract matrix B from matrix A and save matrix(A-B) as .csv")
        matrix_A = str(input("Input .csv for matrix A: "))
        matrix_B = str(input("Input .csv for matrix B: "))
        save_name = str(input("Input the save name for difference matrix: "))
        
        m_A = csv2matrix(matrix_A)
        m_B = csv2matrix(matrix_B)
        diff_matrix = m_A - m_B
        
        diff_matrix.to_csv(save_name, index = False)
        
        print("Saved the difference of ", matrix_A, " - ", matrix_B, " to ", save_name )
        
        del m_A
        del m_B
        del diff_matrix
                    
    elif user_input == "h" or user_input == "help":
        print(" = help = ")
        print (" - Main workflow for creating heatmaps: pdb to csv, csv to matrix, matrix to map.")
        print (" - Main workflow for creating graphs: csv to matrix, matrix to shift data, shift data to graph.")
        print (" - Difference maps: matrix_average to average, and matrix_substract to substract matrices.")
        print (" - For more help visit reffer to the documentation.")
        user_input = input("Press enter to return to the main menu.")
    
    elif user_input == "q" or user_input == "quit":
        main = False
        print("Exiting TrajMap...")
        print("Exited successfully.")
        
    else:
        print("-Command not understood, try again.")
        
    entry_counter = entry_counter + 1
    
