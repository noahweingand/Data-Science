
import pandas
import numpy
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def main():
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # sequences = pandas.read_csv('protein_data.csv', header=None)
    # # remove bad amino acids from sequences
    # for i in range(len(sequences)):
    #     sequences[0][i] = sequences[0][i].replace('B', '')
    #     sequences[0][i] = sequences[0][i].replace('U', '')
    #     sequences[0][i] = sequences[0][i].replace('X', '')
    #     sequences[0][i] = sequences[0][i].replace('Z', '')
    # pandas.DataFrame(sequences).to_csv('updated_protein_data.csv', index_label=None, header=None)

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
    
    # more features
    amino_acid = {}
    first_aa = []
    last_aa = []
    arom = []
    ii = []
    ip = []
    mec_rc = []
    mec_db = []
    ssf_helix = []
    A = []
    C = []
    D = []
    E = []
    F = []
    G = []
    H = []
    I = []
    K = []
    L= []
    M = [] 
    N = []
    P = []
    Q = []
    R = []
    S = []
    T = []
    V = []
    W = []
    Y = []
    data = pandas.read_csv('updated_protein_data.csv', header=None)
    for protein in data.itertuples():
        analyzed_protein = ProteinAnalysis(str(protein[1]))
        amino_acid = (analyzed_protein.count_amino_acids())
        A.append(amino_acid.get('A'))
        C.append(amino_acid.get('C'))
        D.append(amino_acid.get('D'))
        E.append(amino_acid.get('E'))
        F.append(amino_acid.get('F'))
        G.append(amino_acid.get('G'))
        H.append(amino_acid.get('H'))
        I.append(amino_acid.get('I'))
        K.append(amino_acid.get('K'))
        L.append(amino_acid.get('L'))
        M.append(amino_acid.get('M'))
        N.append(amino_acid.get('N'))
        P.append(amino_acid.get('P'))
        Q.append(amino_acid.get('Q'))
        R.append(amino_acid.get('R'))
        S.append(amino_acid.get('S'))
        T.append(amino_acid.get('T'))
        V.append(amino_acid.get('V'))
        W.append(amino_acid.get('W'))
        Y.append(amino_acid.get('Y'))

        first_aa.append(str(protein[1])[0])
        last_aa.append(str(protein[1])[-1])
        arom.append(analyzed_protein.aromaticity())
        ii.append(analyzed_protein.instability_index())
        ip.append(analyzed_protein.isoelectric_point())
        mec_rc.append(analyzed_protein.molar_extinction_coefficient()[0])
        mec_db.append(analyzed_protein.molar_extinction_coefficient()[1])
        ssf_helix.append(analyzed_protein.secondary_structure_fraction()[0])
    

    data.columns = ["PROTEIN SEQUENCE", "CLASS"]
    data["most frequent aa"] = most_frequent
    data["first amino acids"] = first_aa
    data["last amino acid"] = last_aa
    data["aromaticity"] = arom
    data["instability index"] = ii
    data["isolectric point"] = ip
    data["molecular extinction coefficient - reduced cysteines"] = mec_rc
    data["molecular extinction coefficient - disulfid bridges"] = mec_db
    data["secondary structure fraction helix"] = ssf_helix
    data['A'] = A
    data['C'] = C
    data['D'] = D
    data['E'] = E
    data['F'] = F
    data['G'] = G
    data['H'] = H
    data['I'] = I
    data['K'] = K
    data['L'] = L
    data['M'] = M
    data['N'] = N
    data['P'] = P
    data['Q'] = Q
    data['R'] = R
    data['S'] = S
    data['T'] = T
    data['V'] = V
    data['W'] = W
    data['Y'] = Y
    data.to_csv('features.csv')

    

if __name__ == "__main__":
    main()