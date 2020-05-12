

import random
from collections import defaultdict as dd
import numpy as np




# TODO: Conform to PEP8
# • change all variable names to lowercase with underscores instead of camelcase
# • Add docstrings





def acqs_to_imnz(egosNet, acqsNet, numToImnz, timesRefdB4Imnzd):
    
    nodesToImnz = [] 
    numRefsPerNode = dd(int) #A count of the number of times a node has been "marked"
    samplingFrame  = list(egosNet.nodes)
    while (len(nodesToImnz) < numToImnz):
        ego = random.choice(samplingFrame) #randomly choose a node from the available node list
        acqs = [n for n in acqsNet[ego]] #look at the list of neighbors of that node according to the available network data
        if acqs: #if they have neighbors we know of (acqs is not empty)
            acq = random.choice(acqs)
    
            numRefsPerNode[acq] += 1 # "mark" the nodes as important by adding 1 to their number of markings. If they have not been marked before, default dictionary will add them to the dictionary of marked nodes
            if numRefsPerNode[acq] == timesRefdB4Imnzd: #if they have been marked the number of times specified by the threshold
                nodesToImnz.append(acq) #then add them to the list of nodes to immunize
        
    return nodesToImnz






def acqs_to_imnz_updating_sample(egosNet, acqsNet, numToImnz, timesRefdB4Imnzd):
    samplingFrame = list(egosNet.nodes()) 
    allNodes = list(acqsNet.nodes())
    
    nodesToImnz = [] 
    numRefsPerNode = dd(int) #A count of the number of times a node has been "marked"
    while (len(nodesToImnz) < numToImnz):
        ego = random.choice(samplingFrame) #randomly choose a node from the available node list
                 
        acqs =  [n for n in acqsNet[ego]] #We know all of egos neighbors
        
        #In the event that the one of the alters was not already in K, then you just discovered more of the network. Update it
        for alter in acqs:
            if alter not in samplingFrame:
                samplingFrame.append(alter)
                    
        if acqs: #if they have neighbors we know of (acqs is not empty)
            acq = random.choice(acqs)
    
            numRefsPerNode[acq] += 1 # "mark" the nodes as important by adding 1 to their number of markings. If they have not been marked before, default dictionary will add them to the dictionary of marked nodes
            if numRefsPerNode[acq] == timesRefdB4Imnzd: #if they have been marked the number of times specified by the threshold
                nodesToImnz.append(acq) #then add them to the list of nodes to immunize
        
    return nodesToImnz





def subset_network_nodes_missing_at_random(trueNet, numMissing):
    nodesToRemove = np.random.choice(trueNet.nodes(), size=numMissing, replace=False)
    trueNet.remove_nodes_from(nodesToRemove)
        
    return(trueNet, nodesToRemove)



