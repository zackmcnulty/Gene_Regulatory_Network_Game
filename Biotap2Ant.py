'''
@author: Zachary McNulty
'''
import tellurium as te
from GetModel import convert_to_antimony
from GetModel import Gene
# From Biotapestry: Export > Export to SBML
def convert_biotapestry_to_antimony(csv_filename, num_genes, init_params, model_name = "pathway"):
    """
    converts a biotapestry CSV format to an antimony string.
    Uses GetModel as a dependency.

    init_params = parameter values (i.e. rate constants) to be used in the system.
    csv_filename = must be a gene network in the csv format specified by Biotapestry

    """

    f = open(csv_filename)

    # name : [(P1, positive)(P2, negative)]
    all_in_connects = {}
    for i in range(num_genes):
        all_in_connects["P"+ str(i+1)] = []

    all_genes = []
    converter = {"INPUT":Gene("INPUT")} # protein_name : Gene()
    for line in f:
        line = line.replace("\"", "")
        words = line.split(",")
        if words[0].strip() == "general":
            gene_source = words[3].strip()
            source_name = ""
            if gene_source == "INPUT":
                source_name = "INPUT"
            else:
                source_name = "P" + gene_source[5:].strip()

            gene_target = words[5].strip()
            target_name = "P" + gene_target[5:].strip()

            reg_type = words[6].strip()

            all_in_connects[target_name].append((source_name, reg_type))

    #print (all_in_connects)

    for protein_name in all_in_connects.keys():
        in_connects = all_in_connects[protein_name]
        if len(in_connects) == 0:
            next_gene = Gene(protein_name, "N/A")
        elif len(in_connects) == 1:
            if  in_connects[0][1] == "positive":
                next_gene = Gene(protein_name, "SA")
            else:
                next_gene = Gene(protein_name, "SR")
        else:
            if in_connects[0][1] == "positive" and in_connects[1][1] == "positive":
                next_gene = Gene(protein_name, "DA")
            elif in_connects[0][1] == "negative" and in_connects[1][1] == "negative":
                next_gene = Gene(protein_name, "DR")
            else:
                next_gene = Gene(protein_name, "SA+SR")

        all_genes.append(next_gene)
        converter[protein_name] = next_gene

    for gene in all_genes:
        in_connects = all_in_connects[gene.protein_name]
        for other in in_connects:
            gene.add_in_connection(converter[other[0]])

    return convert_to_antimony(all_genes, model_name, init_params, 0)



