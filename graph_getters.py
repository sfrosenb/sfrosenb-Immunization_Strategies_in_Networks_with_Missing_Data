
import random
import networkx as nx
import numpy as np


"""
This file contains the python functions used to get or create the networks
 upon which our methods were tested in the main text & supplement of this paper

Note that the Colorado Springs-like networks and the Add Health-like networks-
 were created in R (code and edgelists in data/ERGM_nets). 
The code related to those nets in this file simply reads the edgelists.

The Erdős-Rényi networks were created directly with a built-in networkx function
Thus they don't have any relevant functions in here.

If you would like to add to this file your own network generator or reader,
 note that ALL networks used in the simulation MUST have node labels as ints
 from 0 to N-1. Otherwise the spreading_CR code can behave unreliably.
 If labels are important, you can store a mapping seperately from your labels
 to the 0 through N-1 labeling, 
 or as a seperate node attribute named something other than "label".
 
Additionally, it is possible there are other aspects the network could have
 that would cause issues with the rest of the code.
If you identify one, please file an issue report
  


"""




# TODO: Make all of this conform to PEP8



def read_our_csv(startOfPath, netVers):
    """ Read in the edge list of a graph.
    In our case we created the graphs in R and output the edgelists in a certain manner
     with file names following a certain convention.
    This function simply reads in csvs formatted in that specific way,
     with the file paths formtted in the manner of the edgelists in the data/ERGM_nets folders
    
    If you want to read in your networks formatted in a different way simply modify this function.
    networkx offers a host of read options.
    If none of them work for you can usually find an intermediary option,
     where the network you have is in some format that can be read by some  other software,
     which can then  output it to something readable by networkx
    
    REQUIREMENT: Whatever network is returned by this function MUST have its nodes labeled from 0 to N-1
    Otherwise the SIR simulations can behave unexpectedly
    

    Parameters
    ----------
    startOfPath : str
        DESCRIPTION.
    netVers : TYPE
        DESCRIPTION.

    Returns
    -------
    G
    The graph read in from the csv with node labels from 0 to N-1

    """
    
    pathname = startOfPath + str(netVers)+ '.csv'
        
    GwithRlabels = nx.read_edgelist(pathname, delimiter=',', nodetype=int, comments='V') # the comments ='V' simply ignores the headers in the csv, otherwise they would be counted as edges
    G = nx.convert_node_labels_to_integers(GwithRlabels, first_label=0, ordering='sorted') # doing sorted ordering so that we know for sure that the same labels always correspond to the same nodes, see notes:https://networkx.github.io/documentation/stable/reference/generated/networkx.relabel.convert_node_labels_to_integers.html#networkx.relabel.convert_node_labels_to_integers I dont think we need that but just in case. 

    return G



# The following two functions are modified versions of NetworkX code, which is licensed using the BSD license
# The text of this license can be found at the bottom of this document



def _random_subset(seq, m):
    """ Return m unique elements from seq.

    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.
    
    Note this function was adapated from a networkx function of the same name

    """
    targets = set()
    while len(targets) < m:
        x = random.choice(seq)

        targets.add(x)
    return targets



def make_ba(n, m0=5, m=3):
    """Make a Barabási-Albert network with m_0 starting nodes and m edges added per additional nodes
    
    This function is adapted from https://networkx.github.io/documentation/latest/_modules/networkx/generators/random_graphs.html#barabasi_albert_graph
    to allow for control in function does not allow for control of m_0 (the number of starting nodes)
    This modification was necessary to reproduce the methods of Chen and Lu because they specifices m_0=5 in their paper https://www.nature.com/articles/s41598-017-03379-4#Sec9
    


    """



    # Add m initial nodes (m0 in barabasi-speak)
    G = nx.empty_graph(m0)


    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes = [node for node in range(m0)]
    
    
    # For the first round,  we pick the targets before the loop, the loop always picks the new targets at the end before restarting and adding a new node
    # in the networkx version, they got to just set targets =list(range(m0)) i.e. all the existing nodes because thats all there was to choose from but thats not true in the m_0 > m case
    
    targets = _random_subset(repeated_nodes, m)
    
    
    # Start adding the other n-m nodes. The first node is m0.
    source = m0
    while source < n:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip([source] * m, targets))
        # Add one node to the list for each new edge just created.
        repeated_nodes.extend(targets)
        # And the new node "source" has m edges to add to the list.
        repeated_nodes.extend([source] * m)
        # Now choose m unique nodes from the existing nodes
        # Pick uniformly from repeated_nodes (preferential attachment)
        targets = _random_subset(repeated_nodes, m)

        source += 1
    return G




def random_choice_no_replace_slow(options, numPairs):
    """
    Out of all possible pairs of options, pick numPairs at random without replacement
    TODO: There should be a faster way to do this

    """
    
    data = np.array([np.random.choice(options, 2, replace=False) for pair in range(numPairs)])
    return data


def make_sj(numCommunities=50, communitySize=40, initDegInCom=8, watts_strogatz_p = 0.1, numBtwnCommEdges=2000, numToRewire=0): 
    """
    Make networks in the way outlined in Salathé and Jones (2010),
    "Dynamics and Control of Diseases in Networks with Community Structure".
    
    Generate numCommunities Watts-Strogatz "small-world" networks, of size communitySize
     with initial in-community degree of initDegInCom,
     and in-community re-wiring probability watts_strogatz_p.
    Then place numBtwnCommEdges randomly between nodes in different communities,
     and then rewire numToRewire of those between-community edges to be within-community edges.
        

    Parameters
    ----------
    numCommunities : int
        Number of Watts-Strogatz small-world networks to create as the first step.
        The default is 50 as was used in Salathé and Jones, (2010)
    communitySize : int
        Number of nodes in each of the numCommunities small-world networks.
        The default is 40 as was used in Salathé and Jones (2010).
    initDegInCom : even int (if uneven, the floor will be used) 
        Initial in-community degree of all nodes.
        The default is 8 as was used in Salathé and Jones, (2010).
    watts_strogatz_p : int
        The rewiring probability of the in-community edges.
        For more information see https://networkx.github.io/documentation/stable/reference/generated/networkx.generators.random_graphs.watts_strogatz_graph.html or Duncan J. Watts and Steven H. Strogatz, Collective dynamics of small-world networks, Nature, 393, pp. 440–442, 1998.
        Note that Salathé and Jones do not specify the value used for this parameter in their work, 
        So we just used an arbitrary value for watts_strogatz_p that resulted in networks which were "small-world".
        The purpose of including networks generated in the same manner as Salathé and Jones (2010) in the supplement-
         was not to reproduce the results of that paper,
         but merely to demonstrate that the phenomena we describe in the main text extends to other kinds of networks.

        The default is 0.1 (arbitrary choice but was used in the supplementary material of this paper)
    numBtwnCommEdges : int
        Number of edges to place randomly between communities.
        The default is 2000 as was used in Salathé and Jones, (2010).
    numToRewire : int between 0 and numBtwnCommEdges (inclusive)
        Number of between-community edges to rewire into within-community edges.
        This increases the modularity.
        The default is 0 which is the lowest possible value and one of the two values we used. 

    Returns
    -------
    Gunion : Networkx Graph
        The resulting "network with community structure" as created by the process outlined in Salathé and Jones (2010)

    
    
    TODO: Likely this could all much much faster if we just made the adjacency matrix and then turned that into a network 
     instead of making all the networks individually and unioning them
    
    """
    
    
    
    
    #Make the base network with all the communitites before connecting them
    net_i = range(numCommunities)
    netsToMerge = []
    for i in net_i:
        #net = nx.connected_watts_strogatz_graph(n=communitySize, k=initDegInCom, p=watts_strogatz_p)
        net = nx.watts_strogatz_graph(n=communitySize, k=initDegInCom, p=watts_strogatz_p)
        net= nx.convert_node_labels_to_integers(net, first_label=i*communitySize)
        netsToMerge.append(net)
    Gunion = nx.union_all(netsToMerge)
        

    
    sendCommRecCom = random_choice_no_replace_slow(options=numCommunities, numPairs=numBtwnCommEdges)
    
    #Within the sender and receiver community, pick which nodes are the the sender and receiver
    idsWthnComm = np.random.randint(communitySize, size=(numBtwnCommEdges, 2))
    btwnCommEdges = (sendCommRecCom*communitySize)+idsWthnComm #If this line is confusing, a more explicit version of what it is doing is in the while loop below 
    btwnCommEdgesSet = set([tuple(np.sort(pair)) for pair in btwnCommEdges])
    
    
    # Sometimes we have some duplicate edges 
    # so in order to make sure that there are precisely 2000 between-community edges we fill those last couple in
    # Sampling edges without replacement here across multiple axes with multiple constraints 
    # would be a very complicated process and so while this seems inefficient,
    # it actually couldn't get much better
    
    while len(btwnCommEdgesSet) < numBtwnCommEdges:
        coms = np.random.choice(numCommunities, 2, replace=False)
        idsInCom = np.random.randint(communitySize, size=2)
        btwnCommEdge = (coms[0]*communitySize + idsInCom[0], coms[1]*communitySize + idsInCom[1])
        btwnCommEdgesSet.add(btwnCommEdge)
        
    btwnCommEdges = list(btwnCommEdgesSet)
    
        
    # Rewire Step
    numRewired = 0
    while numRewired < numToRewire:
        
        #The following steps correspond directly to the paper by Salathe and Jones (2010) in the "Generation of network with community structure" subsection of Methods
        #Step (i): randomly choose a between-community edge
        edgeToRewire = btwnCommEdges[np.random.randint(len(btwnCommEdges))]
        
        #Step (ii): randomly choose one of the two communities that the edge connects
        whichNode = np.random.randint(2) #coinflip between 0 or 1
        nodeWhoseCommunityWeChoose = edgeToRewire[whichNode]
        commOfNode = nodeWhoseCommunityWeChoose // communitySize # floor divide the node's id by community size to get the community number
        
        #Step (iii): pick a random node of the chosen community
        randomNodeOfChosenComm = commOfNode*communitySize + np.random.randint(communitySize)
        if randomNodeOfChosenComm == nodeWhoseCommunityWeChoose:
            continue #no self loops
        
        # Step (iv): rewire the edge by detaching it from the node of the community that was not chosen in step (ii), and attaching it to the new node in the community that was chosen in step (iii)               
        if Gunion.has_edge(nodeWhoseCommunityWeChoose, randomNodeOfChosenComm) == False: 
            
            Gunion.add_edge(nodeWhoseCommunityWeChoose, randomNodeOfChosenComm)
            btwnCommEdges.remove(edgeToRewire)
            numRewired+=1
        
        
    

    
    Gunion.add_edges_from(btwnCommEdges)
    
    

    return (Gunion)







"""
License
=======

NetworkX is distributed with the 3-clause BSD license.

::

   Copyright (C) 2004-2020, NetworkX Developers
   Aric Hagberg <hagberg@lanl.gov>
   Dan Schult <dschult@colgate.edu>
   Pieter Swart <swart@lanl.gov>
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

     * Neither the name of the NetworkX Developers nor the names of its
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""