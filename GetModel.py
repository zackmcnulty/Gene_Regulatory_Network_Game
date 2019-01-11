'''
@author: Zachary McNulty
'''


import numpy as np
import time
import math
import os


def get_model(num_genes, reg_probs = [0.2, 0.2, 0.2, 0.2, 0.2], model_name="pathway",init_params=[1.0/60,1,1.0/60,1,5.0/60,5,1.0/60],
              param_std = 0.25, seed = 0, reachability=0.9, self_feedback_min = 0, max_builds = 1000, export = False):
    """
    Generates and returns an antimony string for a random biological pathway involving num_genes genes.
    The pathway is fully connected, and contains no orphans.
    Also generates a CSV format compatible with the Biotapestry program.
    Both are returned in a tuple `(antimony, biotapestry)`.


    num_genes = number of genes included in the network (an INPUT is also included, but not counted as part of num_genes)
    reg_probs = probability a gene is  each regulation type -> [prob(SA), prob(SR), prob(DA), prob(SA+SR), prob(DR)]
        * SA/SR = single activator/repressor, DA/DR = doubel activator/repressor, SA+SR = one repressor one activator
    model_name = what the antimony model will be named;file with antimony string will be created with name "model_name_antimony.txt" in working directory (if export == True)
    init_params = mean values used to randomly generate parameter values for our reactions ['d_protein', 'd_mRNA' , 'L' , 'Vm' , 'a_protein' , 'H', 'K']
        * d_protein/d_mRNA = degradation rate of protein/mRNA, L = leak rate (of transcription), Vm = maximum rate of gene expression (transcription rate),
        * a_protein = translation rate, H = Hill coefficient, K = rate constant
    param_std = controls std used to randomly generate parameter values from normal dist;
                proportion of mean that standard deviation should be
                parameters will be generated using distribution ~N(mean in init_params, param_std * mean)
    seed = controls seed for randomly generated portions of this method
    reachability = lower bound for proportion of network that should be reachable by INPUT (allows filtering out of less active networks)
    self_feedback_min = lower bound for the number of self feedback loops that are desired in the network
    max_builds = maximum number of attempts to build model that matches desired specifications
    export = if True, files model_name_antimony.txt and model_name_biotapestry.csv will be made in working directory
    """

    ### Invalid parameter handling

    #catch for illegal file names
    forbiddenChar = ['<','>',':','"','/','\'','|','?','*']
    try:
        assert list(set(model_name).intersection(forbiddenChar))==[]
    except AssertionError:
        raise AssertionError('You used an illegal character in your model name!')

    escapeChars = ['t','r','n','b','f','0']
    if model_name[0] in escapeChars:
        raise ValueError('The first character of your model name cannot start with an escape character [\'t\',\'r\',\'n\',\'b\',\'f\',\'0\'].')

    if not type(num_genes) == int or num_genes < 2:
        raise ValueError("num_genes is invalid: it must be a single integer greater than 1")

    if not type(reg_probs) == list or any([not (type(x) == int or type(x) == float) for x in reg_probs]) or not sum(reg_probs) == 1 or not len(reg_probs) == 5 or any(x < 0 for x in reg_probs):
        raise ValueError("reg_probs is invalid: reg_probs must be a list of 5 positive decimals that sum to 1")

    if not type(init_params) == list or any([not (type(x) == int or type(x) == float) for x in init_params]) or not len(init_params) == 7 or any(x < 0 for x in init_params):
        raise ValueError("init_params is invalid: init_params must be a list of 7 positive numbers")

    if not type(model_name) == str or model_name == None or len(model_name) == 0:
        raise ValueError("model name is invalid: must be a valid string with no illegal characters")

    if not (type(reachability) == int or type(reachability) == float) or reachability < 0 or reachability > 1:
        raise ValueError("reachability must be a number between 0 and 1")

    if not (type(param_std) == float or type(param_std) == int) or param_std < 0:
        raise ValueError("param_std must be a positive number")

    if not type(self_feedback_min == int) or self_feedback_min < 0 or self_feedback_min > num_genes:
        raise ValueError("self_feedback_min must be a non-negative integer")

    if not type(max_builds == int) or max_builds < 0:
        raise ValueError("max_builds must be a positive integer")



    if seed == 0:
       np.random.seed()
    else:
        print('Using seed: ' + str(seed))
        np.random.seed(seed)
    if not type(export) == bool:
        raise ValueError("export option must be True/False")

    '''
     Algorithm: look through each gene. Based on its type, assign other proteins/genes
     to act as its activators/repressors. Before assigning however, check if this process
     will create an orphan (within each disjoint set, keep track of how many inputs are
     left between all genes within the set; if the previous action connects two genes within
     a set AND the # inputs remaining in whole set is only 1, do not connect the two). The INPUT messes with
     this algorithm because it acts as a gene with no in-connections. Thus, if it initially connects
     to a single-regulation type (SA or SR) it may create an orphan, as there is no guarantee
     this SA/SR will connect to other genes, and since it used up its only in connection with
     INPUT, it might not be connected at all. This is likely only an issue if the proportion of
     the double input genes (DA, DR, SA+SR) is significantly low. To avoid orphans, we
     regenerate any models that contain orphans or that do not satisfy reachability/feedback conditions
    '''

    current_build = 1
    while current_build <= max_builds:
        gene_types = ["SA", "SR", "DA", "SA+SR", "DR"]

        # chooses a random collection of gene types of the required size with the given probabilities
        random_reg_types = np.random.choice(gene_types, size=num_genes, p=np.asarray(reg_probs), replace=True)

        # assigns genes their names and regulation types
        all_genes = []
        for i, reg_type in enumerate(random_reg_types):
            protein_name = "P" + str(i+1)
            all_genes.append(Gene(protein_name, reg_type))

        # DisjointSets helps keeps track of which genes are connected, directly or indirectly.
        # Used to prevent the formation of "orphans"
        gene_sets = DisjointSets()
        for gene in all_genes:
            gene_sets.make_set(gene.protein_name, -1*(gene.remaining_connections + 1))

        # number of self feedback loops created
        feedback_count = assign_connections(all_genes, gene_sets)

        # Handles the case where INPUT causes algorithm to fail; this is only likely when the proportion
        # of double input genes (DA, DR, SA+SR) is low. Second condition handles case where graph
        # is poorly reachable from INPUT (leading to a low activity network)
        quality = check_input_quality(all_genes)
        if not gene_sets.get_set_count() == 1 or quality < reachability or feedback_count < self_feedback_min:
            print("Model " + str(current_build) +  " does not meet desired specifications. Rebuilding model...")
            current_build += 1
        else:
            ant_str = convert_to_antimony(all_genes, model_name, init_params, param_std)
            biotap_str = convert_to_biotapestry(all_genes)

            if export:
                f = open(model_name + "_antimony.txt", 'w')
                f.write(ant_str)
                f.close()

                f2 = open(model_name + "_biotapestry.csv", 'w')
                f2.write(biotap_str)
                f2.close()

            print("\n\nThe reachability of this network from INPUT is " + str(round(quality*100)) +
                    "% and there are " + str(feedback_count) + " self feedback loops")
            return (ant_str, biotap_str)

    raise TimeoutError("A model could not be generated with the desired specifications in the allotted number of builds." +
            "Try increasing max_builds, or lowering the constraints on the model (such as the reachability or self_feedback_min)")

def assign_connections(all_genes, gene_sets):
    '''
    Assigns all in connections for each gene (i.e. assigns proteins which act as regulators for each gene)
    params
      all_genes = a list of all the genes in the network
      gene_sets = a collection of disjoint sets, keeping track of which genes are already connected

    Assumptions:
      * a protein cannot be connected to same gene more than once
      * connections should not be formed if they result in the formation of an orphan
    '''

    # randomly choose a gene to have INPUT as one of its regulators
    input_gene = np.random.choice(all_genes)
    input_gene.add_in_connection(Gene("INPUT"))
    gene_sets.decrement_value(input_gene.protein_name)

    #keep track of the total # of unassigned connections within entire network
    total_connections_left = sum([gene.remaining_connections for gene in all_genes])

    # keeps track of the number of self feedback loops created
    feedback_count = 0

   # go gene to gene and fill in any empty regulatory sites
    for gene in all_genes:

        # while the gene still has unoccupied/unassigned regulatory sites
        while gene.remaining_connections > 0:

            # to avoid connecting same protein to both regulatory sites on a gene, do not consider
            # genes whose proteins are already acting as a regulator for this current gene
            available_genes = [x for x in all_genes if x not in gene.in_connections]

            gene_to_add = np.random.choice(available_genes)

            # this section checks if new connection is valid (does not result in orphan)
            name1 = gene.protein_name
            name2 = gene_to_add.protein_name

            set1 = gene_sets.find_set(name1)
            set2 = gene_sets.find_set(name2)

            # If the two genes are not already connected, connect them
            if not set1 == set2:
                gene_sets.union(name1, name2)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1

            # else genes are part of same set (already connected), but there may
            # still be enough connections remaining to form another without creating an orphan.

            # Second condition handles edge case where final connection is trying to be made.
            # In this case, the situation is similar to when an orphan is being formed. You are
            # connecting two genes in same network, and only one connection remains. However,
            # the "orphan" is the entire, complete network
            elif gene_sets.get_total_connections(name1) > 1 or total_connections_left == 1:
                gene_sets.decrement_value(name1)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1
                if name1 == name2:
                    feedback_count += 1

    return feedback_count


def convert_to_antimony(all_genes, model_name, init_params, std_dev_perc):
    """
    Converts the generated network to an antimony string. Uses init_params and std_dev_perc to randomly generate parameter
    values for each reaction with a mean given in init_params and a standard deviation of mean * std_dev_perc. Gives the model
    the name model_name. To load, use tellurium.loada(ant_str).

    """


    ant_str = ""
    ant_str += "model *" + model_name + "()\n\n"
    ant_str += "\t// Compartments and Species:\n"

    species = "".join([gene.protein_name + ", " + "mRNA" + str(i+1) + ", " for i, gene in enumerate(all_genes)])
    ant_str += "\tspecies INPUT, " + species[:-2] + ";\n\n"

    rules = {}
    ant_str += "\t// Assignment Rules (production rates used in reactions):\n"
    for i, gene in enumerate(all_genes):
        reg_type = gene.reg_type
        inputs = gene.in_connections
        if (len(inputs) > 0):
            P1 = inputs[0].protein_name
        if len(inputs) == 2:
            P2 = inputs[1].protein_name

        # Ki_j corresponds to rate constant i for protein j
        K1 = "K1_" + str(i+1)
        K2 = "K2_" + str(i + 1)
        K3 = "K3_" + str(i + 1)
        H1 = "H" + str(i + 1)

        num = ""
        denom = ""
        if reg_type == "SA":
            num = '(' + K1 + '*' + P1 + '^' + H1 + ')'
            denom = '(1 +' + K1 + '*' + P1 + '^' + H1 + ')'

        elif reg_type == "SR":
            num = ' 1 '
            denom = '(1 +' + K1 + '*' + P1 + '^' + H1 + ')'

        elif reg_type == "DA":
            num = '(' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + '+ K1 + '*' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'
            denom = '(1 + ' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + ' + K1 + '*' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'

        elif reg_type == "DR":
            num = ' 1 '
            denom = '(1 + ' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + ' + K1 + '*' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'

        elif reg_type == "SA+SR":
            num = '(' + K1 + '*' + P1 + '^' + H1 + ')'
            denom = '(1 + ' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + '+ K1 + '*' + K2 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'

        # else:
        #    raise ValueError("gene type does not match any of the 5 standard varieties")

        expression = "Vm" + str(i+1) + '*(' + num  + '/' + denom + ')'

        if reg_type == "N/A":
            rules["v"+str(i+1)] = "0"
        else:
            rules["v" + str(i+1)] = expression
        ant_str += "\t// transcription" + str(i+1) + "("+str(gene.reg_type)+" : in connections = " + str(gene.in_connections)+ ")" + " uses production rate := " + expression +  ";\n"

    var_names = ["d_protein", "d_mRNA", "L", "Vm", "a_protein", "H", "K1_", "K2_", "K3_"]
    ant_str += "\n\tconst "
    for i in range(len(all_genes)):
        for var in var_names:
            ant_str += var + str(i+1) + ", "
    ant_str = ant_str[:-2] + ";\n"


    ant_str += "\n\t// Reactions:\n"
    for i in range(len(all_genes)):
        ant_str += "\t//transcription" + str(i+1) + "\n"
        ant_str += "\tJ" + str(i+1) + ": => mRNA" + str(i+1) + " ; L" + str(i+1) + " + "+ rules["v" + str(i+1)] + " - d_mRNA" + str(i+1) + " * mRNA" + str(i+1) + ";\n"
        ant_str += "\t//translation" + str(i+1) + "\n"
        ant_str += "\tF" + str(i+1) + ": => P" + str(i+1) + " ; " + "a_protein" + str(i+1) + " * mRNA" + str(i+1) + " - d_protein" + str(i+1) + " * P" + str(i+1) + ";\n"

    ant_str += "\n\t// Species initializations:\n"
    ant_str += "\tINPUT = 1;\n"
    for i in range(len(all_genes)):
        ant_str += "\tmRNA" + str(i+1) + " = 0;\n"
        ant_str += "\tP" + str(i + 1) + " = 0;\n"

    var_names = ["d_protein", "d_mRNA", "L", "Vm", "a_protein", "H", "K1_", "K2_", "K3_"]
    for i in range(len(all_genes)):
        for k, var in enumerate(var_names):
            mean = 0
            if k < 5:
                mean = init_params[k]
            elif k == 5:
                mean = math.ceil(init_params[k])
            else:
                mean = init_params[6]
            std = mean * std_dev_perc
            value = np.random.normal(loc=mean, scale=std)
            ant_str += "\t" + var + str(i+1) + " = " + str(value) + ";\n"

    ant_str += "\n\nend"
    return ant_str


def check_input_quality(all_genes):
    """
    Checks how many genes can be activated directly or indirectly by INPUT
    As input is the only source of stimulus besides the leak rate, this determines how much of the network
    can be activated. Returns proportion of genes that can be reached. This is only important if the INPUT
    is relatively high compared to the leak rates.
    """

    # out_connections -> {protein_name : out connects [P1, P2, ...]}
    out_connections = {"INPUT": []}
    for gene in all_genes:
        out_connections[gene.protein_name] = []

    for gene in all_genes:
        for in_connect in gene.in_connections:
            out_connections[in_connect.protein_name].append(gene.protein_name)

    # Standard Depth first search starting at INPUT
    reachable_genes = 0
    to_vist = ["INPUT"]
    genes_counted = []
    while not len(to_vist) == 0:
        curr_gene = to_vist.pop()
        for adjacent_gene in out_connections[curr_gene]:
            if adjacent_gene not in genes_counted:
                reachable_genes += 1
                to_vist.append(adjacent_gene)
                genes_counted.append(adjacent_gene)

    return 1.0*reachable_genes/len(all_genes)


def convert_to_biotapestry(all_genes):
    """
    Converts the given network to a CSV format that is readable by the program Biotapestry
    and returns this as a string.
    To load into Biotapestry, from the program choose:  File > Import > Import Full Model Hierarchy from CSV
    """

    biotap = ""
    biotap += "model,root,,,,\n"
    reg_patterns = {"SR": ["Repressor"], "SA": ["Activator"], "SA+SR": ["Activator","Repressor"],
                    "DA": ["Activator", "Activator"], "DR": ["Repressor", "Repressor"]}
    for gene in all_genes:
        for i, in_connect in enumerate(gene.in_connections):
            biotap += "general,root,"
            if in_connect.gene_name == "INPUT":
                biotap += "box,"
            else:
                biotap += "gene,"
            biotap += in_connect.gene_name + ","
            biotap += "gene," + gene.gene_name + ","

            reg_type = reg_patterns.get(gene.reg_type)[i]
            if reg_type == "Activator":
                biotap += "positive\n"
            elif reg_type == "Repressor":
                biotap += "negative\n"
            else:
                raise ValueError("reg_type is not recognized")

    return biotap


class Gene:
    """
    Keeps track of the gene regulatory type, the protein name associated to this gene,
    the other genes this gene has already been connected to, and the remaining connections to be made
    """

    def __init__(self, protein_name, reg_type=None):
        # name of protein created by this gene
        self.protein_name = protein_name

        # name of gene
        self.gene_name = ""
        if self.protein_name == "INPUT":
            self.gene_name = "INPUT"
        else:
            self.gene_name = "Gene " + self.protein_name[1:]

        # regulation type this gene undergoes
        self.reg_type = reg_type

        # list of genes connected to this one (genes whose proteins are activators/repressors)
        self.in_connections = []
        if self.protein_name == "INPUT":
            self.remaining_connections = 0
        elif self.reg_type == "SA" or self.reg_type == "SR":
            self.remaining_connections = 1
        else:
            self.remaining_connections = 2

    # connect another gene to this one (add another regulator)
    def add_in_connection(self, other):
        if self.remaining_connections <= 0:
            raise ValueError("Too many connections have been added to a gene of this type")
        self.in_connections += [other]
        self.remaining_connections -= 1

    def equals(self, other):
        return self.protein_name == other.protein_name

    # toString method mostly used for debugging
    def __repr__(self):
        #return str(self.protein_name)
        return str(self.reg_type) + " (" + str(self.protein_name) + "): "
        # str(self.remaining_connections) + " connection(s) remaining")


class DisjointSets:
    """
    DisjointSets keeps a collection of disjoint sets. In this code, these sets represent collections of
    genes that are connected, directly or indirectly, to each other (when the network is considered as an
    undirected graph). Helps keep track of conditions important for avoiding the creation of orphans.

    The following link gives some more information on disjoint sets (pg 32-57)
    https://courses.cs.washington.edu/courses/cse373/18su/files/slides/lecture-19.pdf
    """

    def __init__(self):
        # stores values that point towards index of head of disjoint set, or stores a sentinel value representing
        # total # of remaining connections within that set if that element is the head of a disjoint set
        self.pointers = []

        # allows us to find index in array that stores value associated to given gene (keys are protein names)
        self.converter = {}

        # number of total elements in all disjoint sets
        self.size = 0

        # count of the number of disjoint sets in the collection
        self.set_count = 0

    # makes a new set of a single element, and stores the given value as a sentinel value
    def make_set(self, item, value):
        if value >= 0:
            raise ValueError("The sentinel value must be less than zero")
        if self.contains(item):
            raise ValueError("This item is already present in the DisjointSet")

        self.pointers.append(value)
        self.converter[item] = self.size
        self.size += 1
        self.set_count += 1

    # returns total # of connections remaining in set that contains the given item
    def get_total_connections(self, item):
        head_index = self.find_set(item)
        return -1 * (self.pointers[head_index] + 1)

    # returns the # of sets in the collection
    def get_set_count(self):
        return self.set_count

    # returns the index (in self.pointers) of the representative element of the disjoint
    # set that the given protein_name is in
    def find_set(self, protein_name):
        index = self.converter.get(protein_name)
        return self.find_helper(index)

    def find_helper(self, index):
        if self.pointers[index] < 0:
            return index
        else:
            representative_index = self.find_helper(self.pointers[index])
            self.pointers[index] = representative_index
            return representative_index

    # returns true if any disjoint set in the collection contains the given item
    def contains(self, item):
        return item in self.converter.keys()

    # combines the two sets of item1 and item2; item1 becomes the new head (arbitrary in this case)
    def union(self, item1, item2):
        if not self.contains(item1) or not self.contains(item2):
            raise ValueError("One of those items given is not part of this disjoint set collection")
        head1 = self.find_set(item1)
        head2 = self.find_set(item2)
        if head1 == head2:
            raise ValueError("These two items are already part of the same set")

        # +2 is due to 1 connection lost as new connection forms and an extra +1 is added to account for the artificial
        #  extra -1 we added at the beginning to keep all our sentinel values < 0
        self.pointers[head1] += self.pointers[head2] + 2
        self.pointers[head2] = head1
        self.set_count -= 1

    # used for the case where we want to form a connection between two genes that are already indirectly
    # connected. This will not create an orphan as long as that set has some remaining unassigned connections
    def decrement_value(self, item):
        # reduces total # of connections remaining in this set (uses +1 vs -1 as sentinal values are negative)
        head_index = self.find_set(item)
        self.pointers[head_index] += 1

    # toString mostly used for debugging purposes
    def __repr__(self):
        return str(self.pointers) + str(self.converter)
