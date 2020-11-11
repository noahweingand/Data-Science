
import pandas
import numpy
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
import re



def main():
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    dipeptide = ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL', 'AM', 'AN', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AV', 'AW', 'AY', 
    'CA', 'CC', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT', 'CV', 'CW', 'CY', 
    'DA', 'DC', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN', 'DP', 'DQ', 'DR', 'DS', 'DT', 'DV', 'DW', 'DY', 
    'EA', 'EC', 'ED', 'EE', 'EF', 'EG', 'EH', 'EI', 'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER', 'ES', 'ET', 'EV', 'EW', 'EY', 
    'FA', 'FC', 'FD', 'FE', 'FF', 'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT', 'FV', 'FW', 'FY', 
    'GA', 'GC', 'GD', 'GE', 'GF', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS', 'GT', 'GV', 'GW', 'GY', 
    'HA', 'HC', 'HD', 'HE', 'HF', 'HG', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP', 'HQ', 'HR', 'HS', 'HT', 'HV', 'HW', 'HY', 
    'IA', 'IC', 'ID', 'IE', 'IF', 'IG', 'IH', 'II', 'IK', 'IL', 'IM', 'IN', 'IP', 'IQ', 'IR', 'IS', 'IT', 'IV', 'IW', 'IY', 
    'KA', 'KC', 'KD', 'KE', 'KF', 'KG', 'KH', 'KI', 'KK', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR', 'KS', 'KT', 'KV', 'KW', 'KY', 
    'LA', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LK', 'LL', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT', 'LV', 'LW', 'LY', 
    'MA', 'MC', 'MD', 'ME', 'MF', 'MG', 'MH', 'MI', 'MK', 'ML', 'MM', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY', 
    'NA', 'NC', 'ND', 'NE', 'NF', 'NG', 'NH', 'NI', 'NK', 'NL', 'NM', 'NN', 'NP', 'NQ', 'NR', 'NS', 'NT', 'NV', 'NW', 'NY', 
    'PA', 'PC', 'PD', 'PE', 'PF', 'PG', 'PH', 'PI', 'PK', 'PL', 'PM', 'PN', 'PP', 'PQ', 'PR', 'PS', 'PT', 'PV', 'PW', 'PY', 
    'QA', 'QC', 'QD', 'QE', 'QF', 'QG', 'QH', 'QI', 'QK', 'QL', 'QM', 'QN', 'QP', 'QQ', 'QR', 'QS', 'QT', 'QV', 'QW', 'QY', 
    'RA', 'RC', 'RD', 'RE', 'RF', 'RG', 'RH', 'RI', 'RK', 'RL', 'RM', 'RN', 'RP', 'RQ', 'RR', 'RS', 'RT', 'RV', 'RW', 'RY', 
    'SA', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SK', 'SL', 'SM', 'SN', 'SP', 'SQ', 'SR', 'SS', 'ST', 'SV', 'SW', 'SY', 
    'TA', 'TC', 'TD', 'TE', 'TF', 'TG', 'TH', 'TI', 'TK', 'TL', 'TM', 'TN', 'TP', 'TQ', 'TR', 'TS', 'TT', 'TV', 'TW', 'TY', 
    'VA', 'VC', 'VD', 'VE', 'VF', 'VG', 'VH', 'VI', 'VK', 'VL', 'VM', 'VN', 'VP', 'VQ', 'VR', 'VS', 'VT', 'VV', 'VW', 'VY', 
    'WA', 'WC', 'WD', 'WE', 'WF', 'WG', 'WH', 'WI', 'WK', 'WL', 'WM', 'WN', 'WP', 'WQ', 'WR', 'WS', 'WT', 'WV', 'WW', 'WY', 
    'YA', 'YC', 'YD', 'YE', 'YF', 'YG', 'YH', 'YI', 'YK', 'YL', 'YM', 'YN', 'YP', 'YQ', 'YR', 'YS', 'YT', 'YV', 'YW', 'YY']

    sequences = pandas.read_csv('protein_data.csv', header=None)

    lengths = []
    weights = []
    for protein in sequences.itertuples():
        protein_length = len(str(protein[1])) # length of protein sequence
        lengths.append(protein_length)
        analyzed_protein = ProteinAnalysis(str(protein[1]))
        ambigious_match = re.findall("X+|Z+", protein[1])
        if ambigious_match:
            molecular_weight = "?"
        else:
            molecular_weight = analyzed_protein.molecular_weight()
        weights.append(molecular_weight)
    # remove bad amino acids from sequences
    for i in range(len(sequences)):
        sequences[0][i] = sequences[0][i].replace('B', '')
        sequences[0][i] = sequences[0][i].replace('U', '')
        sequences[0][i] = sequences[0][i].replace('X', '')
        sequences[0][i] = sequences[0][i].replace('Z', '')
    pandas.DataFrame(sequences).to_csv('updated_protein_data.csv', index_label=None, header=None, index=None)

    # use amino acid composition results from pfeature to generate most common amino acid and dipeptide
    data = pandas.read_csv('updated_protein_data.csv', header=None)
    data = numpy.asarray(data)
    most_frequent_di = []
    most_frequent = []
    for i in range(len(data)):
        max = 0
        col = 0
        for j in range(len(dipeptide)):
            c = data[i][0].count(dipeptide[j])
            if(c > max):
                max = c
                col = j
        most_frequent_di.append(dipeptide[col])
        for j in range(len(aa)):
            c = data[i][0].count(aa[j])
            if(c > max):
                max = c
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
    ssf_turn = []
    ssf_sheet = []
    gravy = []
    ph_0 = []
    ph_7 = []
    ph_14 = []
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
    classes = []
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
        ssf_turn.append(analyzed_protein.secondary_structure_fraction()[1])
        ssf_sheet.append(analyzed_protein.secondary_structure_fraction()[2])
        gravy.append(analyzed_protein.gravy())
        ph_0.append(analyzed_protein.charge_at_pH(0.0))
        ph_7.append(analyzed_protein.charge_at_pH(7.0))
        ph_14.append(analyzed_protein.charge_at_pH(14.0))
        classes.append(protein[2])
    
    features = pandas.DataFrame()
    features["LENGTH"] = lengths
    features["MOLECULAR WEIGHT"] = weights
    features["most frequent aa"] = most_frequent
    features["first amino acids"] = first_aa
    features["last amino acid"] = last_aa
    features["most frequence dipeptide"] = most_frequent_di
    features["aromaticity"] = arom
    features["instability index"] = ii
    features["isolectric point"] = ip
    features["molecular extinction coefficient - reduced cysteines"] = mec_rc
    features["molecular extinction coefficient - disulfid bridges"] = mec_db
    features["secondary structure fraction helix"] = ssf_helix
    features["secondary structure fraction turn"] = ssf_turn
    features["secondary structure fraction sheet"] = ssf_sheet
    features["gravy"] = gravy
    features["charge at ph 0"] = ph_0
    features["charge at ph 7"] = ph_7
    features["charge at ph 14"] = ph_14
    features['A'] = A
    features['C'] = C
    features['D'] = D
    features['E'] = E
    features['F'] = F
    features['G'] = G
    features['H'] = H
    features['I'] = I
    features['K'] = K
    features['L'] = L
    features['M'] = M
    features['N'] = N
    features['P'] = P
    features['Q'] = Q
    features['R'] = R
    features['S'] = S
    features['T'] = T
    features['V'] = V
    features['W'] = W
    features['Y'] = Y
    features["CLASS"] = classes
    features.to_csv('features.csv', index=None)

    

if __name__ == "__main__":
    main()