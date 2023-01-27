
import pickle
import numpy
import pandas

output_filename = "collated_128_p128_v3_dict.pickle"
input_filename = "collated_128_p128_v3_list.pickle"
data = pickle.load(open(input_filename,'rb'))


mydict = {}
for index, row in data.iterrows():
	if row[0] not in mydict:
		mydict[row[0]] = {}
	if row[1] not in mydict[row[0]]:
		mydict[row[0]][row[1]] = {}
	mydict[row[0]][row[1]][row[2]]=index

pickle.dump(mydict,open(output_filename,'wb'))
