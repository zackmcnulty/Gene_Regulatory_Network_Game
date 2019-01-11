# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka & Zachary McNulty
"""

import tellurium as te
import time
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotap2Ant import convert_biotapestry_to_antimony

"""
This code will generate an antimony model for an unbroken network, and output this as a .txt file.
There are a few settings you can mess with to generate a suitable network, and I describe them below.
When you have set them, run this script and a model will be generated in your working directory.
It also plots the timecourse data for the model, so you will have a chance to see if you like the behavior it exhibits.

options:

`seednum` will be the seed for the random processes involved in generating the model.
Some models are not as interesting as others, so you may want to play with the seednum until you get something you want.

`seednum` is saved in `seed.txt` and will be overwritten every time this script is run.
`maxSimTime` is saved in `tmax.txt`. **DO NOT DELETE THIS FILE**. It is used in `run_experiments.py`.

`num_genes` will be the number of genes in the model

`reachability` defines the proportion of  genes in the network you want to be connected to the INPUT, directly or indirectly.
This is really only important if the INPUT is high relative to the other protein concentrations.

`self_feedback_min` is the  min number of self feedback loops that the model will have.
If you want self feedback, increasing this number will do the trick.

Do NOT change the model name from the default of `pathway`; the other code requires the name to be `pathway`.
"""

#seednum = np.random.randint(1,999999999) # You can use this line when trying to generate new networks
seednum = 425300122
num_genes = 8
reachability = 0.9
self_feedback_min = 0
maxSimTime = 1200

ant_str, biotap_str = get_model(num_genes, seed = seednum, model_name ="pathway" , reachability=reachability, self_feedback_min=self_feedback_min, export=True)

r = te.loada(ant_str)
result = r.simulate(0,maxSimTime,maxSimTime*5)
#r.plot()

### Plotting

df = pd.DataFrame.from_dict(result)
selections = r.getFloatingSpeciesIds()
selections.insert(0,'time')
df.columns = selections

genes = map(str, range(1,num_genes+1))
protein = ['P' + s for s in genes]
protein.insert(0,'time')
df1 = df[df.columns.intersection(protein)]
df1.set_index('time', inplace=True)
df1.plot()
plt.ylabel("Protein count")
plt.xlabel("time")
plt.show()

mRNA = ['mRNA' + s for s in genes]
mRNA.insert(0,'time')
df1 = df[df.columns.intersection(mRNA)]
df1.set_index('time', inplace=True)
df1.plot()
plt.ylabel("mRNA count")
plt.xlabel("time")
plt.show()

# Saving text
f = open("seed.txt", "w")
f.write('Seed: ' + str(seednum))
f.close()

f = open("tmax.txt", "w")
f.write(str(maxSimTime))
f.close()

print("Done! Model has been generated")