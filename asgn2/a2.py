import csv
import copy
import numpy as np


# ------------------- I/O logic ------------------- 

# read dataset files and load into numpy
dataset_complete_file = open("dataset_complete.csv", "r")
dataset_complete = np.genfromtxt("dataset_complete.csv", delimiter=',',skip_header=1,usecols=[0,1,2,3,4,5,6,7,8,9,10],dtype=None,encoding=None)
dataset_missing01_file = open("dataset_missing01.csv", "r")
dataset_missing01 = np.genfromtxt("dataset_missing01.csv", delimiter=',',skip_header=1,usecols=[0,1,2,3,4,5,6,7,8,9,10],dtype=None,encoding=None)
dataset_missing20_file = open("dataset_missing20.csv", "r")
dataset_missing20 = np.genfromtxt("dataset_missing20.csv", delimiter=',',skip_header=1,usecols=[0,1,2,3,4,5,6,7,8,9,10],dtype=None,encoding=None)

imputedHd01_file = open("V00763681_missing01_imputed_hd.csv", "w")
imputedHd01_cond_file = open("V00763681_missing01_imputed_hd_conditional.csv", "w")
imputedHd20_file = open("V00763681_missing20_imputed_hd.csv", "w")
imputedHd20_cond_file = open("V00763681_missing20_imputed_hd_conditional.csv", "w")

# ------------------- functions ------------------- 
def calc_mean(dataset, i):
    total = 0
    count = 0
    for row in dataset:
        if(row[i] != '?'):
            total += float(row[i])
            count += 1
    mean = total/count
    return round(mean, 5)

def calc_mean_cond(dataset, i):
    yesTotal = 0
    noTotal = 0
    yesCount = 0
    noCount = 0
    for row in dataset:
        if(row[i] != '?' and row[10] == 'Yes'):
            yesTotal += float(row[i])
            yesCount += 1
        if(row[i] != '?' and row[10] == 'No'):
            noTotal += float(row[i])
            noCount += 1
    yesMean = yesTotal/yesCount
    noMean = noTotal/noCount
    return round(yesMean, 5), round(noMean, 5)

def findMAE(completeData, imputedData):
    completeErrorValues = []
    imputedErrorValues = []
    maeValues = []
    total = 0
    for i, row in enumerate(completeData):
        for j, num in enumerate(row):
            if (num == "Yes" or num == "No"):
                continue
            if (num != float(imputedData[i][j])):
                completeErrorValues.append(completeData[i][j])
                imputedErrorValues.append(float(imputedData[i][j]))
    for x, t in zip(completeErrorValues, imputedErrorValues):
        maeValues.append(abs(x - t))
    for i in range(0, len(maeValues)):
        total += maeValues[i]
    mae = round((total/len(maeValues)), 4)
    return mae
    
# ------------------- mean logic ------------------- 

# imputed mean logic
    # missing values are ? in numpy
    # binary label has no missing values

# create lists and hold each datset's column's mean in them
missing01_means = []
missing20_means = []
i = 0
while i < 10: # ignore binary column since no missing values
    missing01_means.append(calc_mean(dataset_missing01, i))
    missing20_means.append(calc_mean(dataset_missing20, i))
    i += 1
i = 0

# copy the datasets and replace each column's missing value with the respective mean
imputedDataset01 = copy.deepcopy(dataset_missing01)
imputedDataset20 = copy.deepcopy(dataset_missing20)

for i, row in enumerate(imputedDataset01):
    for j, num in enumerate(row):
        if (num == '?'):
            imputedDataset01[i][j] = missing01_means[j]

for i, row in enumerate(imputedDataset20):
    for j, num in enumerate(row):
        if (num == '?'):
            imputedDataset20[i][j] = missing20_means[j]

# ------------------- mean conditional logic ------------------- 
missing01_cond_yes_means = []
missing20_cond_yes_means = []
missing01_cond_no_means = []
missing20_cond_no_means = []
i = 0
while i < 10: # ignore binary column since no missing values
    yes_01, no_01 = calc_mean_cond(dataset_missing01, i)
    yes_20, no_20 = calc_mean_cond(dataset_missing20, i)
    missing01_cond_yes_means.append(yes_01)
    missing01_cond_no_means.append(no_01)
    missing20_cond_yes_means.append(yes_20)
    missing20_cond_no_means.append(no_20)
    i += 1

imputedCondDataset01 = copy.deepcopy(dataset_missing01)
imputedCondDataset20 = copy.deepcopy(dataset_missing20)

for i, row in enumerate(imputedCondDataset01):
    for j, num in enumerate(row):
        if (num == '?' and row[10] == 'Yes'):
            imputedCondDataset01[i][j] = missing01_cond_yes_means[j]
        if (num == '?' and row[10] == 'No'):
            imputedCondDataset01[i][j] = missing01_cond_no_means[j]
for i, row in enumerate(imputedCondDataset20):
    for j, num in enumerate(row):
        if (num == '?' and row[10] == 'Yes'):
            imputedCondDataset20[i][j] = missing20_cond_yes_means[j]
        if (num == '?' and row[10] == 'No'):
            imputedCondDataset20[i][j] = missing20_cond_no_means[j]

# ------------------- MAE calls ------------------- 

np.savetxt('V00763681_missing01_imputed_mean.csv', imputedDataset01, fmt='%s')
np.savetxt('V00763681_missing20_imputed_mean.csv', imputedDataset20, fmt='%s')
np.savetxt('V00763681_missing01_imputed_mean_conditional.csv', imputedCondDataset01, fmt='%s')
np.savetxt('V00763681_missing20_imputed_mean_conditional.csv', imputedCondDataset20, fmt='%s')

dataset01_mean_mae = findMAE(dataset_complete, imputedDataset01)
dataset20_mean_mae = findMAE(dataset_complete, imputedDataset20)
dataset01_mean_cond_mae = findMAE(dataset_complete, imputedCondDataset01)
dataset20_mean_cond_mae = findMAE(dataset_complete, imputedCondDataset20)

# 4 decimal places
print("MAE_01_mean =", dataset01_mean_mae) 
print("MAE_01_mean_conditional =",dataset01_mean_cond_mae) 
print("MAE_20_mean =", dataset20_mean_mae) 
print("MAE_20_mean_conditional =", dataset20_mean_cond_mae) 