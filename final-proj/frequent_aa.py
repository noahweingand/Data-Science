
import pandas
import numpy

def main():
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    sequences = pandas.read_csv('protein_data.csv', header=None)
    # remove bad amino acids from sequences
    for i in range(len(sequences)):
        sequences[0][i] = sequences[0][i].replace('B', '')
        sequences[0][i] = sequences[0][i].replace('U', '')
        sequences[0][i] = sequences[0][i].replace('X', '')
        sequences[0][i] = sequences[0][i].replace('Z', '')
    pandas.DataFrame(sequences).to_csv('updated_protein_data.csv', index_label=None, header=None)

    # use amino acid composition results from pfeature to generate most common amino acid
    data = pandas.read_csv('final_amino_acid_result.csv', header=None, skiprows=1)
    data = data.drop([0], axis=1)
    data = numpy.asarray(data)
    most_frequent = []
    for i in range(len(data)):
        max = 0
        col = 0
        for j in range(len(data[i])):
            if((data[i][j]) > max):
                max = data[i][j]
                col = j
        most_frequent.append(aa[col])
    pandas.DataFrame(most_frequent).to_csv('frequent_aa.csv', index_label=None)
    
    

if __name__ == "__main__":
    main()