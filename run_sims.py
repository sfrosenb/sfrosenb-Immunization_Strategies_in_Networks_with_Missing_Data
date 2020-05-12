import networkx as nx
import numpy as np
import pandas as pd
import os

from collections import Counter
import time
import sys
import json
import random
from spreading_CR import SpreadingProcess

MAX_CPLUSPLUS_INT = 4294967295


from graph_getters import read_our_csv, make_ba, make_sj
from sim_utils import acqs_to_imnz, acqs_to_imnz_updating_sample, subset_network_nodes_missing_at_random


"""
TODO
----
* Put the actual SirsPerTrial simulations loop into a seperate function in sim_utls.py because its just copy pasted in here twice, once for the no immunization bit and then again for all the immunization strategies
* Make saving the outbreak size for every  individual simulation optional as it would save considerably on storage and on the time it takes to read the files back in later on (and somewhat on run time for this as well)
* Maybe use argparse instead of the flexArgs dict? Heres a good example  of how: https://gitlab.com/compstorylab/hashtag-ecology/-/blob/cnww/parser.py

"""

def run_sims(networkVersions, netType, beta, gamma, fracsMissing, fracsImmunized, SirsPerTrial, flexArgs, numTrialsToEstNoImmOutbreak = 10000, numTrialsToEstR0 = 10000):
    
    obs=[]
    
    
    
        
                
    for netVers in networkVersions: #iterate over the number of times you create a new known network
        print("network version: ", netVers)
        startThisNetworkVers = time.time()
        
        
        
        # The individual network getter/creator functions can take additional arguments if you want to change the params of the networks themselves.
        # This driver simply uses the default params where available to reproduce the analysis from our paper specifically.
        if netType in ['low', 'med', 'high', 'schoolnet']:
            if netType == 'schoolnet':
                startOfPath = "data/ERGM_nets/"+netType+'/edgelist.'+ netType
            elif netType in ['low', 'med', 'high']:
                startOfPath = "data/ERGM_nets/"+netType+'/edgelist.'+ netType + '.clustering'
            G = read_our_csv(startOfPath, netVers)
            
        
        elif netType in ['salatheRewire1800', 'salatheNoRewire']: #in hindsight this was bad naming convention. TODO: organize data to use better naming convention. Maybe salatheRewire{numToRewire}
            G = make_sj(numToRewire=int(flexArgs['numToRewire']))
            
        elif netType == 'ER':
            G = nx.gnm_random_graph(n=int(flexArgs['sizeOfG']), m=int(flexArgs['numberOfEdges']))
            
        elif netType in ['BA1000', 'BA10000']: # in hindsight this was bad naming convention. TODO: organize data to use better naming convention.
            G = make_ba(n= int(flexArgs['sizeOfG']))
        
        sizeOfG = G.number_of_nodes()
        print("sizeOfG: ", sizeOfG)
        listOfGsNodes = list(G.nodes)
        setOfGsNodes = set(listOfGsNodes)
        
        
        listOfGsEdges = list(G.edges())
                
        degDictG = {n: d for n, d in G.degree()}    
            
            #                    sp = SpreadingProcess(listOfGsEdges, transmission_rate, 
#                        recovery_rate, waning_immunity_rate)
                    
                    
        sp = SpreadingProcess(listOfGsEdges, beta, 
                              gamma, 0)

        
        
        #/*******************************************************************************************************\
            # Now do everything we want to do involving contagion on the unimunized network
            # Only need to do this batch once per network version
        
        
        
        
        unimmunizedNodes = listOfGsNodes 
        
        arrayOfPatient0sForEachTrial = np.random.choice(unimmunizedNodes, size=numTrialsToEstNoImmOutbreak  ) 

        patient0s = [[patient0] for patient0 in arrayOfPatient0sForEachTrial] 
        

        Rnode_list = []
        final_size_list = []
        for i in range(numTrialsToEstNoImmOutbreak): # we do more for this one so its stable
            Inode_list = patient0s[i]
            seed = np.random.randint(MAX_CPLUSPLUS_INT+1) #I could make these all at once if they take a long time. I know exactly how many I will need
            sp.initialize(Inode_list, Rnode_list, seed)
            sp.evolve(np.inf)
            
            finalSizeInclTargets =  sp.get_Rnode_number_vector()[-1] 
            final_size_list.append(finalSizeInclTargets) # In this case there are no targets so all recovered nodes were infected, none were immunized
            sp.reset()

        finalSizeDist = Counter(final_size_list)
        finalSizeDist_noImm = dict(finalSizeDist) #We can later use ast.literal_eval to read this back in as a dict when needed
        pctTotalPopInf_noImm = np.mean(final_size_list)/sizeOfG  

        
        
        # Estimate R0: analytically
        
        degs = list(degDictG.values())
        deg_sqrds = [deg**2 for deg in degs]
        
        deg_avg = np.mean(degs)# frequeently written as <k>
        deg_sqr_avg = np.mean(deg_sqrds) # frequently written as <k^2>
        
        R0_formula = (beta/(beta+gamma)) * (deg_sqr_avg/deg_avg -1) #the second term is the "excess degree"
        
        
        # Estimate R0: empirically
        
        seed = np.random.randint(MAX_CPLUSPLUS_INT+1)
        
        
        R0_empirical_mean, R0_empirical_std = sp.estimate_R0(numTrialsToEstR0, seed) 
        sp.reset()
        
        print("R0 formula: ", R0_formula)
        print("R0 emprical: ", R0_empirical_mean)
        print("R0 empirical std: ",R0_empirical_std)
        
        
        
                

        #\_______________________________________________________________________________________________________/                    

        
                    

                    

        
        
        
        
        
        
            
            
            
        for fracImmunized in fracsImmunized:
            print(fracImmunized)

            numNodesRemovedByIntrv = np.round(sizeOfG * fracImmunized).astype(int)
            
            print("numNodesRemovedByIntrv: ", numNodesRemovedByIntrv)

            

            
            
            for fracMissing in fracsMissing:
                startFracMissing = time.time()
                

                
                
                
                
                #create a new known network by removing nodes from the true network (simulating missing data)
                numMissing = np.round(sizeOfG*fracMissing).astype(int)
                K,  nodesOutOfSample = subset_network_nodes_missing_at_random(trueNet=G.copy(), numMissing=numMissing) 
                nodesInSample = K.nodes()

        
        
        #/***************************************************************************************************\

                degDictKK = {n: d for n, d in K.degree()}
                selfReportedDegree = {n: degDictG[n] for n in nodesInSample}
                betDict = nx.betweenness_centrality(K, normalized=False)

        #\___________________________________________________________________________________________________/

        
        
         #/***************************************************************************************************\ 

                scoreDicts = [betDict,  degDictKK, selfReportedDegree]
                stratNames = ["Bet",    "DegKK",        "Srd",         "acqKK2", "acqKG2", "acqGG2", "acqKG2_up", "RandSamp", "No_imm"]

                
                
                
                targetsByStrat = {name:[] for name in stratNames}
                
                acqKK2_can_converge=True
                for strat_idx, strat in enumerate(stratNames): 
                    
                    if strat == "RandSamp":
                        targetsByStrat[strat] = np.random.choice(nodesInSample, numNodesRemovedByIntrv, replace=False) #select the nodes to be immunized via random immunization
                    elif strat == "acqKK2":
                        # if (fracMissing == 0.8 and fracImmunized >= 0.1) or (fracMissing == 0.7 and fracImmunized == 0.15): #TODO: Better way to handle this without checking all networks of that missingness? inflexible but efficient because the alternatives I can think of would involve checking each network to see if it has enough nodes with neighbors to be named an acq and if it doesnt then putting NA or something, or just  immunizing all nodes with degree >0 and randomly immunizing the rest, which is exactly what degree  and other "centrality scores" would do in that instance (because if less than numToImnz even have degree >0 then ANY node with deg>0 will be in the "top" of any measure  of centrality and since  ties are broken randomly then the ties  for last place are  too). Anyway, probably any more flexible way of doing this would either be inefficient because you would do 1999 simulations and one would fail and you couldnt use the whole batch, or it would be incompatible with multiprocessing because you cant check the others when they are  parallel
                        #     continue
                        # else:
                        #     targetsByStrat[strat] = acqs_to_imnz(egosNet=K, acqsNet=K, numToImnz = numNodesRemovedByIntrv, timesRefdB4Imnzd = 2)
                                                # We check whether acqKK2 can converge on this network and if not skip it
                        # If you have the ensemble of networks premade, you could check this beforehand and not even try to run acqKK2 at all which would save you the time spent running those that you wont use, but I wanted this code to be flexible and work without preechecking
                        
                        # I decided to go with flexibility over efficiciency for now
                        
                        degCountKK = Counter(degDictKK.values())
                        numIsolates = degCountKK[0]
                        numNodesWithNonZeroDeg = sizeOfG-numIsolates
                        if numNodesWithNonZeroDeg < numNodesRemovedByIntrv:
                            acqKK2_can_converge=False
                            continue
                        else:
                            targetsByStrat[strat] = acqs_to_imnz(egosNet=K, acqsNet=K, numToImnz =numNodesRemovedByIntrv, timesRefdB4Imnzd =2) 
                    elif strat == "acqKG2":
                        targetsByStrat[strat] = acqs_to_imnz(egosNet=K, acqsNet=G, numToImnz = numNodesRemovedByIntrv,timesRefdB4Imnzd = 2) 
                    elif strat == "acqGG2":
                        targetsByStrat[strat] = acqs_to_imnz(egosNet=G, acqsNet=G, numToImnz =numNodesRemovedByIntrv, timesRefdB4Imnzd =2)
                    elif strat == "No_imm":
                        continue
                    elif strat == "acqKG2_up":
                        targetsByStrat[strat] = acqs_to_imnz_updating_sample(egosNet=K, acqsNet=G, numToImnz=numNodesRemovedByIntrv, timesRefdB4Imnzd=2)
                        
                        
                    else:
                        stratDictAsTuples =  list(scoreDicts[strat_idx].items())  #this is a list of tuples of the form (<node>, <centrality_of_node>)
                        random.shuffle(stratDictAsTuples)
                        stratDictAsTuples.sort(key=lambda x: x[1], reverse=True) #score is the second element in each tuple, node is the first
                        targetsByStrat[strat] = [stratDictAsTuples[target_idx][0] for target_idx in range(numNodesRemovedByIntrv)] #stratDictAsTuples[target_idx] is the target_idx'th element of stratDictAsTuples, which is a tuple where the first element is the target_idx'th scoring node according the node scoring measure (e.g. betweenness), and the second element is its score. We just need the node. 
        #\___________________________________________________________________________________________________/


##############################################################################################################                                    
                    
############################################################################################################## 




           
        #/*******************************************************************************************************\
           
           
               
 

                
                      
                for strat in stratNames:
                    #if strat == "acqKK2" and ((fracMissing == 0.8 and fracImmunized >= 0.1) or (fracMissing == 0.7 and fracImmunized == 0.15)): #TODO: Better way to handle this without checking all networks of that missingness? inflexible but efficient because the alternatives I can think of would involve checking each network to see if it has enough nodes with neighbors to be named an acq and if it doesnt then putting NA or something, or just  immunizing all nodes with degree >0 and randomly immunizing the rest, which is exactly what degree  and other "centrality scores" would do in that instance (because if less than numToImnz even have degree >0 then ANY node with deg>0 will be in the "top" of any measure  of centrality and since  ties are broken randomly then the ties  for last place are  too). Anyway, probably any more flexible way of doing this would either be inefficient because you would do 1999 simulations and one would fail and you couldnt use the whole batch, or it would be incompatible with multiprocessing because you cant check the others when they are  parallel
                    if strat == "acqKK2" and acqKK2_can_converge==False:
                        ob ={'networkVersion:': netVers, 'fracMissing': fracMissing, 'beta':beta, 'gamma':gamma, 'pctImmunized': fracImmunized, 'strategy': strat, 'pctTotalPopInf': np.nan, 'finalSizeDist':np.nan, 'R0_empirical':R0_empirical_mean,  'R0_empirical_std':R0_empirical_std, 'R0_formula': R0_formula, "acqKK2_didntConverge":1}
                        continue
                    if strat == 'No_imm':
                        # For ease of figure making later on we want to have a copy of the no immunization scenario to compare to every other observation, but we only need to run it once, so we just copy it a bunch. It would be more efficient not to do this (save memory, save time reading and writing of files, but it would be complicated to analyze later)
                        ob ={'networkVersion:': netVers, 'fracMissing': fracMissing, 'beta':beta, 'gamma':gamma, 'pctImmunized': fracImmunized, 'strategy': strat, 'pctTotalPopInf': pctTotalPopInf_noImm, 'finalSizeDist':finalSizeDist_noImm, 'R0_empirical':R0_empirical_mean,  'R0_empirical_std':R0_empirical_std, 'R0_formula': R0_formula}
                        obs.append(ob)
                        continue

                    
                    
                    
                    
                    if "imm_effectiveness" in  flexArgs.keys():
                        
                        coin_flips = np.random.rand(len(targetsByStrat[strat]))    
                        
                        immunizedNodes = []
                        for q, target in enumerate(targetsByStrat[strat]):
                            if coin_flips[q] < flexArgs["imm_effectiveness"]:
                                immunizedNodes.append(target)
                    else:
                       immunizedNodes = list(targetsByStrat[strat])     
                    
                    
                    # With the spreading_cr implementation we are simply setting the immunized nodes to the recovered phase. Note that this is mathematically equivalent to removing them from the network as is done in some other network SIR implmentations                  
                    unimmunizedNodes = list(setOfGsNodes - set(immunizedNodes))

                    
                    

                    
                    
                    
                    arrayOfPatient0sForEachTrial = np.random.choice(unimmunizedNodes, size=SirsPerTrial ) #this only takes 2.5 ms per strategy per outer trial(out of totalTrials not per SIR) for 50000 different patient0s, so its not worth changing for a speedup

                    patient0s = [[patient0] for patient0 in arrayOfPatient0sForEachTrial] 
                    
        


                    
                    Rnode_list = targetsByStrat[strat]
                    final_size_list = []
                    for i in range(SirsPerTrial):
                        Inode_list = patient0s[i]
                        seed = np.random.randint(MAX_CPLUSPLUS_INT+1) #I could make these all at once if they take a long time. I know exactly how many I will need
                        sp.initialize(Inode_list, Rnode_list, seed) 
                        sp.evolve(np.inf)
                        
                        finalSizeInclTargets =  sp.get_Rnode_number_vector()[-1] #Still does what its supposed to???
                        adjustedFinalSize = finalSizeInclTargets -  numNodesRemovedByIntrv#                        finalSizeAsFrac =  adjustedFinalSize/sizeOfG
#                        final_size_list.append(finalSizeAsFrac) #dividing one by one only takes about 20 ms per strat per trial for 50000 sir so only a couple seconds a trial or much less. Not worth speeding up regardless of number of sir or trials
                        final_size_list.append(adjustedFinalSize)
                        sp.reset()

                    finalSizeDist = Counter(final_size_list)
                    finalSizeDist = dict(finalSizeDist) #We can later use ast.literal_eval to read this back in as a dict when needed
                    pctTotalPopInf = np.mean(final_size_list)/sizeOfG  
                    

                    
                    ob ={'networkVersion:': netVers, 'fracMissing': fracMissing, 'beta':beta, 'gamma':gamma, 'pctImmunized': fracImmunized, 'strategy': strat, 'pctTotalPopInf': pctTotalPopInf, 'finalSizeDist':finalSizeDist, 'R0_empirical':R0_empirical_mean,  'R0_empirical_std':R0_empirical_std, 'R0_formula': R0_formula}

                    obs.append(ob)
                    
            #\___________________________________________________________________________________________________/
                
                    
                endFracMissing = time.time()
                print('time for loop with ', fracMissing, ' fraction of nodes missing: ', endFracMissing-startFracMissing)
            
            endThisNetworkVers = time.time()
            print('time for single trial, trial ', netVers, ' : ', endThisNetworkVers-startThisNetworkVers)
            







    return obs




def main():
    
    
    # fracsImmunized can be just one value  like "0.10" (i.e. 10% immunized)
    # or it can be multiple values seperated by spaces inside the same string
    # e.g. '0.05 0.1 0.15' -> [0.05, 0.1, 0.15] i.e 5%, 10%, and 15%
    # If you are going to run simulations for multiple values anyway, it saves time to do multiple immunization levels at once, since the %infected with no immunization is the same and so is the estimated R0, together those account for like 5% of total run time
    # However, if you do plan to utilize this, triple check with a stats person that it makes sense for your experimental design because now you are getting caught up with question of independence. 
    # I had some stats people tell me they thought this was fine but I ended up not utilizing for this paper just in case. 
    # I did use it for testing out new stuff though
    
    fracsImmunizedStrs =  sys.argv[1].split()  
    
    
    fracsImmunized  = [float(fracString) for fracString in fracsImmunizedStrs]
    print("\n------------------------------------------------------------------\n")
    print("fractions of N to Immunize:\n", fracsImmunized)
    
    # Though this convention was created to count through the labeled pre-created ergm networks,
    #  it was useful to keep it even for the networks which are created inside the simulation
    # It keeps the results from getting mixed up with each other and lets you keep track of  all the jobs  you have sent out by what results you have received back
    netVersStart = int(sys.argv[2])
    netVersEnd = int(sys.argv[3])
    networkVersions = list(range(netVersStart, netVersEnd))
    
    print("network versions:\n",networkVersions )
    
    Beta  = float(sys.argv[4]) 
    Gamma = float(sys.argv[5])
    
    print("beta: ", Beta, ", gamma: ", Gamma)
    
    
    netType = sys.argv[6]
    
    print("netType: ", netType)
    
    SirsPerTrial = int(sys.argv[7])
    
    print("SirsPerTrial: ", SirsPerTrial)
    
    
    # Ok so clearly I must be re-inventing the wheel here for flexArgs, right
    # But I couldnt something that worked exactly how I wanted it to
    # kwargs was a bust
    # and argparse seems to require I have a whole nother intermediate function which imo would be very confusing for future readers. I had help my first time looking at argparse code and I was still confused. I STILL am
    # TODO: Look into otherways to implement commandline flexible arguments, but I think reinventing the wheel here is more readable than the alternatives by a mile 
    # Sure, argparse lets you see a help message  with -h, and thats all well and good for knowing where to enter stuff, but what if you want to change the code itself. You better strap in and get ready for a lesson in unix systems
    # This way, although I am technically using json, less experienced users can go ahead and treat it like a dictionary all they want as long as they do the quotes right
    if len(sys.argv) >8:
        try:
            flexArgs = json.loads(sys.argv[-1])
        except Exception:
            print("exception type: ", type(Exception))
            print("use single quotes to open and close your json and double quotes for any strings inside. Or the opposite. Just dont mix")
    else:
        flexArgs = dict()
    
    print("flexArgs raw:\n", flexArgs)
    
    
    
    if "fracsMissingMax" in flexArgs.keys():
        # We count down from most missing to least missing because in the debugging  stage there tended to be more issues and things took longer  with more missing data so we would like to know that sooner rather than later
        fracsMissing = np.arange(start=float(flexArgs["fracsMissingMax"]),stop=float(flexArgs["fracsMissingMin"]), step=float(flexArgs["fracsMissingStep"]))
    else:
        fracsMissing = np.arange(start=0.8,stop=0, step=-0.1)
        
        # I switched to using np.arange instead of range because its more clear how to use it, but if you actually look at the elements they can be like 0.00000000002 off from what they are supposed to be. This shouldnt be an issue since we are multiplying against a large integer and rounding but keep this in mind perhaps
        # fracsMissing = [x / 100.0 for x in range(80, -1, -10)]  #Old way
    
    print("fracsMissing:\n", fracsMissing)
    
    
    if "filenameHeader" in flexArgs.keys():
        filenameHeader = flexArgs['filenameHeader']+netType #netType may be redundant but its very important not to forget that in the header
    else:
        filenameHeader = netType
    
    print("filenameHeader:\n", filenameHeader)
    
    
    if "outdir" in flexArgs.keys():
        outdir = flexArgs['outdir']
    else:
        outdir = f'data/results/{filenameHeader}_Beta{str(int(np.round(Beta*100)))}Gamma{str(int(np.round(Gamma*100)))}_SirsPer{str(SirsPerTrial)}_fracsImmunized'+''.join([str(int(np.round(100*fracImmnzd)))+"_" for fracImmnzd in fracsImmunized])
        
    print("output directory: ", outdir)
    
    print("\n\nIf you dont like the way any of those arguments and flex arguments got read in, now would be a good time to stop the processes and try again before  you get your hopes up for a new batch of results and then ahve your dreams dashed...again\n\n")
    
    
    start = time.time()
    print("about to start")
    #obs = run_sims( networkVersions, beta=Beta, gamma=Gamma, fracsMissing=fracsMissing, fracsImmunized=fracsImmunized,  SirsPerTrial=SirsPerTrial, flexArgs=flexArgs)
    obs = run_sims( networkVersions, netType, beta=Beta, gamma=Gamma, fracsMissing=fracsMissing, fracsImmunized=fracsImmunized,  SirsPerTrial=SirsPerTrial, flexArgs=flexArgs)
    
    
    #run_sims(networkVersions, beta, gamma, fracsMissing, pctsImmunized, SirsPerTrial, flexArgs, numTrialsToEstNoImmOutbreak = 10000, numTrialsToEstR0 = 10000 )
    
    obsDF = pd.DataFrame(obs)
    
    
    
    
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #obsDF.to_csv(outdir+"/"+filenameHeader  +'networkVersions'+ str(networkVersions[0])+'thru'+str(networkVersions[-1])+'Beta=' + str(np.round(Beta, decimals=3))+ 'Gamma=' + str(np.round(Gamma, decimals=3))+ 'SirsPerTrial'+str(SirsPerTrial)+'pctsImmunized'+str(fracsImmunized[0]) + '.csv', index = None, header=True) #Don't forget to add '.csv' at the end of the path
    obsDF.to_csv(f"{outdir}/{filenameHeader}networkVersions{str(networkVersions[0])}thru{str(networkVersions[-1])}Beta={str(np.round(Beta, decimals=3))}Gamma={str(np.round(Gamma, decimals=3))}SirsPerTrial{str(SirsPerTrial)}pctsImmunized{str(fracsImmunized[0])}.csv", index = None, header=True) #Don't forget to add '.csv' at the end of the path
    
    end = time.time()
    print("total time: ", end-start)
    


if __name__ == "__main__":
    main()    




