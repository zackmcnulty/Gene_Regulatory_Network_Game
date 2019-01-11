'''
@author: Zachary McNulty & Kateka Seth
'''


import numpy as np
import tellurium as te
import itertools as it

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from change_biotapestry import add_biotapestry

def estimate_connections(genes, num_genes, data, timepoints, csv_filename, csv_newfile, selections, params):  
    """
    Attempts to find new connections that improve the broken model. Uses a least squared error approach 
    with the given data to evaluate the merits of each possible connection, resimulating the broken model
    each time with added connections. . 

    genes: gene numbers you want to find the in-connections to 
    data: the experimental data
    timepoints: [start,stop, step] for r.simulate; must match timepoints in data
    DONT RUN 5 OR MORE GENES AT ONCE 
    """

    mapping = {}
    permConnections = []
    permError = [1000000]*10
    perms = []
    
    ant_str = convert_biotapestry_to_antimony(csv_filename, 8, params)
    # simulate
    r = te.loada(ant_str)
    start = timepoints[0]
    stop = timepoints[1]
    steps = timepoints[2]
    result = r.simulate(start, stop, steps, selections=selections)
    diff = data - result
    error = np.sum(np.power(diff, 2))
    
    permError[0] = error
    mapping[str(error)] = [(-1,-1,-1)] # flag for add nothing; current model is already best
    
    
    perms = list(it.permutations(genes))
    length = len(perms)
    count = 0
    for ii in range(len(perms)):
        choice = perms[ii]
        numAdded = -1
        connection = [] # stores the chosen gene connection
        for jj in range(len(choice)):
            gene = choice[jj]
            numAdded = -1
            singleConnection = connection[:]
            singleError = 1000000
            for i in range(0, num_genes + 1): # 0 = flag for single connection
                for j in range(1, num_genes + 1):
                    for k in (-1, 1):
                        for m in (-1, 1):   
#                            if m == 0: # don't add connection just add it to permError
#                                permError.sort()
#                                #place = -1
#                                for kk in range(10, -1, -1):
#                                    if singleError <= permError[kk]:
#                                        place = kk
#                                if place != -1:
#                                    permError.insert(kk, singleError)
#                                    mapping[str(singleError)] = connection
                            if (i != j):
                                if (i == 0):
                                    add = [(j, gene, k)]
                                else:
                                    add = [(i, gene, k), (j, gene, m)]
                                # add the new connection
                                singleConnection.extend(add)
                                try:
                                    add_biotapestry(singleConnection, csv_filename, csv_newfile)
                                except ValueError:
                                    break
                                ant_str = convert_biotapestry_to_antimony(csv_newfile, 8, 
                                                  params)
                                os.remove(csv_newfile)
                                # simulate
                                r = te.loada(ant_str)
                                r.Vm2 = 12.5
                                r.Vm7 = 10
                                r.Vm1 = 8
                                result = r.simulate(timepoints[0], timepoints[1], 
                                                    timepoints[2], selections=selections)
                                diff = data - result
                                error = np.sum(np.power(diff, 2))
                                # test error for best error
                                if error < singleError:
                                    # removes old (worst) connection
                                    if numAdded == 1:  
                                        connection.pop()
                                    elif numAdded == 2: 
                                        connection.pop()
                                        connection.pop()
                                    # stores how many new connections are added
                                    if len(add) == 1:
                                        numAdded = 1
                                    else:
                                        numAdded = 2
                                    singleError = error
                                    connection.extend(add)
                                # remove choice for next loop
                                singleConnection.pop()
                                if len(add) == 2:
                                    singleConnection.pop()   
        count += 1
        print("progress: " + str(count) + "/" + str(length))
        permError.sort()
        place = -1
        for kk in range(1,11):
            print(kk)
            print(permError)
            if singleError <= permError[kk]:
                place = kk
                break
        if place != -1:
            print(permError)
            permError.insert(kk, singleError)
            remove=permError[10]
            permError.pop()
            if remove in mapping.keys():
                del mapping[str(remove)]
            mapping[str(singleError)] = connection
            print(connection)
    permError.sort()
    print(permError)
    for i in range(1,11):
        if permError[i] != 1000000:
            permConnections.append(mapping.get(str(permError[i])))
    return permConnections





def objective_func(paramset, r, data, timepoints, selections):
    """
    returns the total squared error between the data and the simulated data generated from
    "best guess" model. Used for parameter estimation with an optimization tool requiring an
    objective function. Only approximates the means for each parameter type (i.e. the mean Vm,
    mean L, mean d_p rather than d_p1, d_p2, etc...)

    paramset: vector of values for each parameter (length of 7)
    r : roadrunner instance storing model
    data: numpy 2D array storing concentrations of each species at each timepoint
    timepoints: [start, stop, steps] gives information on timepoints that our model (r) 
                will simulate species concentrations at; needs to match timepoints provided in data
    selections: which species should be compared; note, you will have to pre-process "data"
                so that it only contains the species you are interested in
    """
    
    
    ids = r.getGlobalParameterIds()
    r.resetToOrigin()
    for i, next_id in enumerate(ids):
        if i % 9 < 6:
            val = paramset[i%9]
        else:
            val = paramset[6]

        exec("r.%s = %f" % (next_id, val)) #runs the given formatted string as if it is a line of code

    start = timepoints[0]
    stop = timepoints[1]
    steps = timepoints[2]
    result = r.simulate(start,stop,steps, selections=selections)
    diff = data - result
    return  np.sum(np.power(diff, 2))


def single_parameter_sweep(paramset, ID,  r, data, timepoints, selections):
    """
    Error function for approximating the individual parameter values rather than the means
    """
    r.resetToOrigin()
    for i, next_param in enumerate(paramset):
        exec("r.%s = %f" % (ID + str(i+1), next_param)) #runs the given formatted string as if it is a line of code

    result = simulate(timepoints[0], timepoints[1], timepoints[2], selections=selections)
    diff = data - result
    return np.sum(np.power(diff,2))
