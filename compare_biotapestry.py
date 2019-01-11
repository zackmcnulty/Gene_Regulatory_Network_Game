# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 13:48:20 2018

@author: Kateka Seth
"""

from change_biotapestry import add_biotapestry
import os

"""
given:
* student's connections
* broken csv file
* name of new student's file
* true csv file
"""
def compare_biotapestry(add, csv_broken, csv_network, csv_student):
    add_biotapestry(add, csv_broken, csv_student)
    network = set()

    #strip useless stuff, add to set
    f = open(csv_network)
    for line in f:
        line = line.replace("\"", "")
        line = line.replace(" ","")
        line = line.strip()
        words = line.split(",")
        for i in range(len(words)):
            words[i] = words[i].strip()
        if words[0] == "general" or words[0] == "nodeOnly":
            network.add(line)
    f.close()
    f = open(csv_student)
    for line in f:
        line = line.replace("\"", "")
        line = line.replace(" ","")
        line = line.strip()
        if line not in network and len(line) > 0:
            words = line.split(",")
            if words[0] == "general":
                print("wrong connection: " + words[3] + " to " + words[5])
            elif words[0] == "nodeOnly":
                print("wrong connection: " + words[3] + " not a node only.")
    f.close()

#debug code, ignore below
#add=[(6,8,1),(3,5,-1),(6,7,-1),(7,7,-1)]
#csv_broken="Biotapestry/8gene_broken.csv"
#csv_network="Biotapestry/8gene_network.csv"
#csv_student="Biotapestry/8gene_student.csv"
#compare_biotapestry(add, csv_broken, csv_network, csv_student)