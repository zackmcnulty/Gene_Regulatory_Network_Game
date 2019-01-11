# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 09:38:56 2018

@author: Kateka Seth
"""


"""
Given a biotapestry csv file, removes connections specified by remove. Prints
a new csv file to given location.

remove: A list containing int tuples of connections to be remove. The tuple
        should be formatted as (source,target). Ex: (1,3) to remove connection
        starting at Gene 1 going to Gene 3. Use "INPUT" for "INPUT" box.
csv_filename: Location with name of the original csv file.
csv_newfile: Location with name of new file. Ex: `file_name.csv`.
"""
def remove_biotapestry(remove, csv_filename, csv_newfile):

    escapeChars = ['t','r','n','b','f','0']
    if csv_newfile[0] in escapeChars:
        raise ValueError('The first character of your model name cannot start with an escape character [\'t\',\'r\',\'n\',\'b\',\'f\',\'0\'].')

    f = open(csv_filename)
    f_new = open(csv_newfile,'w')

    og_nodes = set()
    new_sources = set()
    new_targets = set()
    canPrint = True

    for line in f:
        line = line.replace("\"", "")
        words = line.split(",")
        if "#" not in words[0].strip() and words[0].strip() != "general":
            f_new.write(line)
        elif words[0].strip() == "general":
            gene_source = words[3].strip()
            gene_target = words[5].strip()
            og_nodes.add(gene_source)
            og_nodes.add(gene_target)
            for i in range(len(remove)):
                given_source = "Gene " + str(remove[i][0])
                if remove[i][0] == "INPUT":
                    given_source ="INPUT"
                given_target = "Gene " + str(remove[i][1])
                if (given_source == gene_source) and (given_target == gene_target):
                    canPrint = False
            if canPrint:
                f_new.write(line)
                new_sources.add(gene_source)
                new_targets.add(gene_target)
            canPrint = True
    # print the node only
    node_only=og_nodes.difference(new_sources).difference(new_targets)
    for node in node_only:
        f_new.write("nodeOnly, root, gene, " + node + "\n")
    f.close()
    f_new.close()


"""
Given a biotapestry csv file, adds new connections specified by add. Prints
a new csv file to given location.

:param add: A list containing tuples of connections to be added. The tuple
             should be formatted as (source,target,type), where type is 1 for
             a positive connection and -1 for a negative connection.
             Use "INPUT" for "INPUT" box.
:param csv_filename: Location with name of the original csv file
:param csv_newfile: Location with new file name
"""
def add_biotapestry(add, csv_filename, csv_newfile):
    f = open(csv_filename)
    f_new = open(csv_newfile,'w')
    new_nodes = set()

    # If a nodeOnly is inside new_nodes, erase the nodeOnly line
    for i in add:
        if i[0] == "INPUT" or i[1] == "INPUT":
            new_nodes.add("INPUT");
        else:
            new_nodes.add("Gene " + str(i[0]))
            new_nodes.add("Gene " + str(i[1]))

    # erase nodeOnly
    for line in f:
        line = line.replace("\"", "")
        words = line.split(",")
        if "Instance" not in words[0].strip() and words[0].strip() != "nodeOnly":
            write_fencepost(f_new, words)
        elif words[0].strip() == "nodeOnly":
            if words[3].strip() not in new_nodes:
                write_fencepost(f_new, words)

    # write new nodes
    for i in range(0, len(add)):
        f_new.write("general, root, ")
        for j in range(0,2):
            if add[i][j] == "INPUT":
                f_new.write("box, ")
            else:
                f_new.write("gene, ")
            f_new.write("Gene " + str(add[i][j]) + ", ")
        if add[i][2] == 1:
            f_new.write("positive\n")
        else:
            f_new.write("negative\n")
    f.close()
    f_new.close()

    # check if too many connections have been formed; if so, remove csv_newfile and throw exception
    f_new = open(csv_newfile, 'r')
    connections = {}
    for line in f_new:
        line = line.replace("\"", "")
        words = line.split(",")
        if words[0].strip() == "general":
            next_target = words[5].strip()
            if next_target not in connections:
                connections[next_target] = 1
            else:
                connections[next_target] += 1

    #print (connections)
    if any(x > 2 for x in list(connections.values())):
            #os.remove(csv_newfile)
            raise ValueError("You have added too many connections to one of the genes.\n A gene can have at most 2 connections.")
    else:
        f_new.close()


"""
Writes the given words array to the given file with a comma between each array element.
"""
def write_fencepost(f_new, words):
    f_new.write(words[0].strip())
    for i in range(1, len(words)):
        f_new.write(", " + words[i].strip())
    f_new.write("\n")
