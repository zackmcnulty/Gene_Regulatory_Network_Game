# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:46:56 2018

@author: Yoshi
"""

from RunModel import run_model
from GetModel import get_model
from GetModel import convert_to_biotapestry
#from RunModel_2 import run_model2

model_name='exp5'
numGenes = 8
seednum = 123468

#antStr,biotap_str = get_model(numGenes,model_name=model_name,seed=seednum,export=False)

#antStr= open('C:\Users\Yoshi\Documents\GitHub\DREAM-work\GRN_Game\pathway_antimony.txt','r').read()
noiseLevel = 0.00 # put in a percentage. 0.05 = 5%
tmax=200 # minutes. The complete data will have tmax*5 datapoints
resolution = 1 # The number of minutes between each timepoint you want in your output.
#r,result,resultN = run_model(antStr,noiseLevel,exportData=[0,'P',True,True],inputData=[1,tmax,resolution],showTimePlots=True,bioTap=biotap_str)

antStr = open('pathway_antimony.txt', 'r').read()

#r,res,resN = run_model(antStr,noiseLevel,exportData=[0,'P',True,False,False],inputData=[1,tmax,resolution],showTimePlots=True)

#        'perturb (list: [pertSpecies, pertType, mean, stdev]) :'
genes = [[1,3,5], 'P']
genes = [0,'P']
mean = [50,0]
targ = [4,5]

# no perturb (wild-type)
run_model(antStr, noiseLevel, exportData=[False, False], showTimePlots=True, savePath='test', fileName='no_perturb', genesToExport=genes)

# Knockout
run_model(antStr, noiseLevel, exportData=[False, False], showTimePlots=True, savePath='test', fileName='ko', genesToExport=genes, perturb=[targ,"KO", mean,0])

# Gene upreg
run_model(antStr, noiseLevel, exportData=[False, False], showTimePlots=True, savePath='test', fileName='up', genesToExport=genes, perturb=[targ,"UP", mean,0])

# Gene downreg
run_model(antStr, noiseLevel, exportData=[False, False], showTimePlots=True, savePath='test', fileName='down', genesToExport=genes, perturb=[targ,"DOWN", mean,0])

print('done!')

#%% testing biotap removal

from change_biotapestry import remove_biotapestry

#rem = [(8,1),(1,5),(4,7),(5,3),(6,6)]
#rem = [(5,3),(6,6),(8,1),(1,5),(4,7)]
#rem = [(8,8),(1,5),(5,3)]

#remove_biotapestry(rem,'C:\Users\Yoshi\Documents\GitHub\DREAM-work\zack_kateka\Random_GRNs\exp3\\biotapestry.csv','stestWorking.csv')
print('done')
