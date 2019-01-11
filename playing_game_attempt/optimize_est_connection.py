'''
@author: Zachary McNulty & Kateka Seth
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import itertools as it
import scipy

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from change_biotapestry import add_biotapestry


"""
Runs differential evolution in an attempt to roughly estimate parameter values using data
Uses "broken" model, so the estimate is based on an incorrect model, but filtering out
likely poor/incorrect species from this estimation can improve its performance. You can
select which species to use in the estimation by altering the "selections" variable
Timepoints is the a vector describing the timepoints present in the data. It is in the form
[start, stop, step_size]. This is needed so the simulated data from the partially complete model
lines up with the experimental data
"""
broken_model = "model_files/biotapestry_broken.csv"
selections = ["mRNA" + str(i) for i in range(1,9)] 
timepoints = [0,200, 40]


'''
Load in experimental data.
'''
# from first RNA seqeunce test
data_table = pd.read_csv('model_files/RNASeq_HiRes.csv') # RNASeq_HiRes has timepoints = [0,200,20]
data_table.set_index('time', inplace=True) 
#data_table[selections].plot(style='.-')


# converts dataframe into numpy 2D array
#data = data_table.as_matrix(columns=selections)-----(depreciated)
data = data_table[selections].values

# Load in our current "best guess" for the model
ant_str = convert_biotapestry_to_antimony(broken_model, 8, [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472, 4.25284957, 0.06687737])
r = te.loada(ant_str)
#r.simulate(timepoints[0],timepoints[1], timepoints[2], selections = ['time'] + selections) 
#r.plot()
#plt.show()




#testing purposes
data_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_network.csv", 8,  [1.0/60, 1, 1.0/60, 1,5.0/60, 5, 1.0/60])
data_model = te.loada(data_str)	
data = data_model.simulate(0,200,40, selections=selections)



"""
the one for gene interaction
gene: source gene you want to look at
data: the experimental data
timepoints: [start,stop, step] for r.simulate
DONT RUN 5 
"""
# TODO: plz optimize
def estimate_connections(gene, data, timepoints, csv_filename, csv_newfile, selections):  
    f_new = open(csv_filename, 'r')
    connections = {i:0 for i in gene} # gene # : in connection count
    for line in f_new:
        line = line.replace("\"", "")
        words = line.split(",")
        if words[0].strip() == "general":
            next_target = words[5].strip()
            next_target = int(next_target[5:])
            if next_target in connections:
                connections[next_target] += 1
    
    f_new.close()
    print ("connection counts: " + str(connections))
   
    #--------------------
    # csv_newfile starts as copy of initial file
    with (open(csv_filename)) as f1:
        with (open(csv_newfile, 'w')) as f2:
            f2.truncate(0)
            for line in f1:
                f2.write(line)
    
    f1.close()
    f2.close()
    temp_file = csv_newfile[:-4] + "_other.csv" 
    #-------------------

    mapping = {}
    permError = [] 
    perms = list(it.permutations(gene)) # added feature so each gene considers option to choose no connections

    ant_str = convert_biotapestry_to_antimony(csv_filename, 8,
        [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
    # simulate
    r = te.loada(ant_str)
    result = r.simulate(timepoints[0],timepoints[1],timepoints[2], selections=selections)
    diff = data - result
    error = np.sum(np.power(diff, 2))
    permError.append(error)
    mapping[str(error)] = [(-1,-1,-1)] # flag for add nothing; current model is already best


    for add_order in perms:
        chosen_connections = []
        for next_gene in add_order:
            for _ in range(2-connections[next_gene]):
                #find best connection to add to next_gene
                add = find_best_connection(8,next_gene, data, timepoints, csv_newfile, selections)
                if not add == None:
                    chosen_connections.extend(add)
                    print ("add: " + str(add))
                    print ( "chosen connections: " + str(chosen_connections) )
                
                    add_biotapestry(add, csv_newfile, temp_file)
                
                    temp = temp_file
                    temp_file = csv_newfile
                    csv_newfile = temp 

         
        ant_str = convert_biotapestry_to_antimony(csv_newfile, 8,
                [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
        # simulate
        r = te.loada(ant_str)
        result = r.simulate(timepoints[0],timepoints[1],timepoints[2], selections=selections)
        diff = data - result
        next_error = np.sum(np.power(diff, 2))
        if str(next_error) not in mapping.keys():
            mapping[str(next_error)] = [chosen_connections]
        else:
            mapping[str(next_error)].append(chosen_connections)
        permError.append(next_error)

        # permError stores 10 lowest erros (which act as keys for the map)
        if len(permError) > 10:
            permError.sort()
            del mapping[str(permError.pop())]
       # assess overall error of new model; add to mapping error_of_connections : chosen_connections 
        
        # csv_newfile starts as copy of initial file
        with (open(csv_filename)) as f1:
            with (open(csv_newfile, 'w')) as f2:
               f2.truncate(0) 
               for line in f1:
                    f2.write(line)

        f1.close()
        f2.close()
        
    permError.sort()
    print("perm error: " + str(permError))
    best_connections = []
    for next_error in permError:   
        best_connections.append(mapping[str(next_error)])
    #return list of best connections to add for each gene
    print("best connections:" + str(best_connections))

    return zip(permError, best_connections)

def find_best_connection(total_gene_count, gene_num, data, timepoints, csv_filename, selections):
    ant_str = convert_biotapestry_to_antimony(csv_filename, 8,
        [1/60, 1, 1/60, 1, 5/60, 5, 1/60])

    r = te.loada(ant_str)
    result = r.simulate(timepoints[0], timepoints[1], timepoints[2], selections = selections)
    diff = data - result
    min_error =  np.sum(np.power(diff,2))
    best_choice = None # add no connection
    for i in range(1, total_gene_count + 1):
        for k in (-1,1):
            next_connection = [(i,gene_num,k)]
            add_biotapestry(next_connection, csv_filename, "temp_biotap.csv")
            # assess error, if better than min, reassign min/best_choice
            ant_str = convert_biotapestry_to_antimony("temp_biotap.csv", 8,
                            [1/60, 1, 1/60, 1, 5/60, 5, 1/60])

            r = te.loada(ant_str)
            result = r.simulate(timepoints[0], timepoints[1], timepoints[2], selections=selections)
            diff = data - result
            error = np.sum(np.power(diff, 2))
            if error < min_error:
                min_error = error
                best_choice = next_connection

    os.remove("temp_biotap.csv")
    return best_choice



'''
Probes for possible connections; we can investigate the feasibility of these connections using further experimental data
'''

#[2,8,7,5]
connection = estimate_connections([2,8,7,5], data, timepoints, "../Biotapestry/8gene_broken.csv", "../Biotapestry/8gene_ie.csv", selections)
output = str(list(connection))
f_new = open("est_connection_output.txt", 'w')
f_new.write(output)
f_new.close()
print("Best connection " + output)
