import numpy as np 
import csv
import pickle as pkl
import sys


filename = "temp.csv"
filename = "collated_128_p128_v3.csv"
#filename = "collated_32_p32_v3.csv"
#filename = sys.argv[1]

import pandas as pd
print("start reading")
csvreader = pd.read_csv(filename, engine="pyarrow",header=None)
print("data read in")
dataframe = csvreader.iloc[:, 3:]
names = csvreader.iloc[:, :3] 
#names["names"] = names[0].astype(str) + names[1].astype(str) + names[2].astype(str) 
print("data split")

dataframe = dataframe.to_numpy()
print("to numpy")
dataframe = dataframe.astype(float)
print("asfloat")
row_sums = dataframe.sum(axis=1)
new_matrix = dataframe / row_sums[:, np.newaxis]
print("normalised")
dataframe = np.log(dataframe)
print("logged all the values in the array")
dataframe = dataframe.transpose()
print("transposed")

# save the column index as a dictionary and as a list
filename_dict = filename.split(".csv")[0]+"_dict.pickle"
with open(filename_dict, 'wb') as handle:
    pkl.dump(names, handle, protocol=pkl.HIGHEST_PROTOCOL)
print("names saved")

num_rows, num_cols = dataframe.shape
for i in range(num_rows):
    # SAVE THE NUMPY ARRAY
    filename_np = filename.split(".csv")[0]+"_"+str(i)+".npy"
    #np.savez_compressed(filename_np, myarray, allow_pickle=True, fix_imports=True)
    np.save(filename_np, dataframe[i,], allow_pickle=True, fix_imports=True)
    print(str(i)+" done")



