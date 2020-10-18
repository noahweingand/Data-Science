from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import re

# print("seq %s is %i bases long" % (test, len(test)))
# print("reverse complement is %s" % test.reverse_complement())
# print("protein translation is %s" % test.translate())

lengths = []
weights = []

training_data = pd.read_csv("./8319945.txt")

for protein in training_data.itertuples():
    print(protein[1])
    protein_length = len(str(protein[1])) # length of protein sequence
    lengths.append(protein_length)
    analyzed_protein = ProteinAnalysis(str(protein[1]))
    ambigious_match = re.findall("X+|Z+", protein[1])
    if ambigious_match:
        molecular_weight = "?"
    else:
        molecular_weight = analyzed_protein.molecular_weight()
    weights.append(molecular_weight)
    #amino_acid_count = analyzed_protein.count_amino_acids()

    # print("Index: %s" % protein[0])
    # print("Protein Length: " + str(protein_length))
    # print("Molecular Weight: %s" % molecular_weight)
    # print(" ")

training_data.columns = ["PROTEIN SEQUENCE", "CLASS"]
training_data["LENGTH"] = lengths
training_data["MOLECULAR WEIGHT"] = weights

training_data.to_csv('training_data.csv')