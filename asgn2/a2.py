import csv
import numpy as np

# read dataset files
dataset_complete_file = open("dataset_complete.csv", "r")
dataset_complete = np.genfromtxt("dataset_complete.csv", delimiter=',',skip_header=1,usecols=[0,1,2,3,4,5,6,7,8,9,10],dtype=None,encoding=None)
dataset_missing01_file = open("dataset_missing01.csv", "r")
dataset_missing01 = np.genfromtxt("dataset_missing01.csv", delimiter=',',skip_header=1,usecols=[0,1,2,3,4,5,6,7,8,9,10],dtype=None,encoding=None)
dataset_missing20_file = open("dataset_missing20.csv", "r")
dataset_missing20 = np.genfromtxt("dataset_missing20.csv", delimiter=',',skip_header=1,usecols=[0,1,2,3,4,5,6,7,8,9,10],dtype=None,encoding=None)


# create result files
imputedMean01 = open("V00763681_missing01_imputed_mean.csv", "w")
imputedMean01_cond = open("V00763681_missing01_imputed_mean_conditional.csv", "w")
imputedMean20 = open("V00763681_missing20_imputed_mean.csv", "w")
imputedMean20_cond = open("V00763681_missing20_imputed_mean_conditional.csv", "w")
imputedHd01 = open("V00763681_missing01_imputed_hd.csv", "w")
imputedHd01_cond = open("V00763681_missing01_imputed_hd_conditional.csv", "w")
imputedHd20 = open("V00763681_missing20_imputed_hd.csv", "w")
imputedHd20_cond = open("V00763681_missing20_imputed_hd_conditional.csv", "w")

# print MAE values


# ------------------- main logic ------------------- 

# imputed mean logic



# missing values are ? in numpy
for row in dataset_missing02:
    print(row)
