from run_experiments import export_experiments
import os
import numpy as np


"""
This section will create a team score file if it doesn't exist.
***Be sure to change team names in the `team_names` array.***
***Specify starting amount of credits below.***
"""

maxcredits = 14000

if os.path.isfile("team_scores.csv") == False:
    # change team names here. Make sure it matches the google form team names
    team_names = np.array(["Quintus", "Marcus", "Gaius", "Claudius", "Lucius", "Scipio"])
    money = np.array([str(maxcredits)]*len(team_names))
    result = np.vstack((team_names, money))
    np.savetxt("team_scores.csv", result, fmt="%s",delimiter=",")
    print("Team score file made!")

"""
!!Make sure sendEmail and updateMoney is set to False if debugging!!

This file will send out the next batch of experimental data.
Simply ensure that num_genes is the number of genes in the network, and run the script.

!!Make sure sendEmail and updateMoney is set to False if debugging!!
"""
num_genes = 8
filename = "BIOEN 498 Experiment Request Form.csv" # the file where all the experimental data orders are stored

export_experiments(num_genes, csv_file=filename, sendEmail=True, updateMoney=True)
