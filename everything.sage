# precompute
import sys
import os

args = sys.argv
i=int(args[1])
j=int(args[2])
n=int(args[3])
dimension=n
l=int(args[4])
p=args[5]


collated_output_file = "precomputed/collated_"+str(i)+"_"+str(j)+"_"+str(dimension)+"_"+p+".csv"
print("START")
print(collated_output_file)

import pickle
import numpy as np
from numpy import linalg as LA

if not os.path.isfile(collated_output_file):
  outfile = "precomputed/pre_mat"+str(n)+"_u"+str(i)+"_d"+str(j)+".pickle"
  if not os.path.isfile(outfile):
    m=load("matrix_"+str(n)+".sobj")
    m2 = m.subs(u=i/100,d=j/100)
    m2 = m2.apply_map(RR)
    #m2 = m2.np(dtype='longdouble')
    m2 = m2.numpy(dtype='double')
    with open(outfile, 'wb') as m_data:
      pickle.dump(m2, m_data)
    #save(m2,outfile)


  #for i2 in range(101):
  #  for j2 in range(101):
  #    if m2[i2,j2] > 1:
  #      m2[i2,j2] = 1
  #    if m2[i2,j2] < 0:
  #      m2[i2,j2] = 0
  #    assert m2[i2,j2] >= 0 
  #    assert m2[i2,j2] <= 1

  # all_paths
  # need to also edit this for constrained paths that can't go to zero.

  l=sys.argv[4]
  path_combinations = sys.argv[5]
  file = "precomputed/pre_mat"+str(n)+"_u"+sys.argv[1]+"_d"+sys.argv[2]+".pickle"

  with open(file,'rb') as m_data:
    m = pickle.load(m_data)
  #m = load(file)
  with open("GD.pickle",'rb') as GD_data:
    G = pickle.load(GD_data)
  #G = load("GD.pickle")
  
  dimension = int(n) # int(file.split("pre_mat")[1].split("_")[0])
  m = m[:dimension,:dimension]
  G = G[:dimension,:dimension]

  # need to normalise all the rows still!!!!!

  for i in range(dimension):
    total = sum(sum(m[i,:]))
    for j in range(dimension):
      m[i,j] = m[i,j]/total

  all_paths = load("all_path_combinations_"+path_combinations+".sobj")
  #with open("all_path_combinations_"+path_combinations+".pickle",'rb') as path_combinations_data:
  #  all_paths = pickle.load(path_combinations_data)
  #
  single_paths = [x for x in all_paths if "G" not in x]

  powers_filename = file.split(".pickle")[0] + ".powers.pickle"
  if os.path.isfile(powers_filename):
    with open(powers_filename,'rb') as infile:
      powers = pickle.load(infile)
  else:
    powers = {}
  for path in single_paths:
    if int(path) not in powers:
      res = LA.matrix_power(m,int(path))
      powers[int(path)] = res
  with open(powers_filename,'wb') as infile:
    pickle.dump(powers,infile)
  
  #output = file.split(".sobj")[0] + ".precomputed.sobj"
  output = file.split(".pickle")[0] + ".precomputed.pickle"
  count = 0
  if os.path.isfile(output):
    try:
      with open(output,'rb') as precomputed_data:
        myd = pickle.load(precomputed_data)
      #myd = load(output)
    except:
      # delete the file
      os.remove(output)
      myd = {}
  else:
    myd = {}

  for path in all_paths:
    if path in myd:
      continue
    splits = path.split("G")
    G1 = 0
    G2 = 0
    if len(splits) == 1:
      pre = int(path)
      mid = 0
      post = 0
    if len(splits) == 2:
      pre = int(splits[0])
      mid = 0
      post = int(splits[1])
      G1 = 1
    if len(splits) == 3:
      pre = int(splits[0])
      mid = int(splits[1])
      post = int(splits[2])
      G1 = 1
      G2 = 1

    # compute the probabilities for this path:
    res = powers[pre] #LA.matrix_power(m,pre)
    if G1 > 0:
      res = np.matmul(res,G)
    res = np.matmul(res,powers[mid]) #LA.matrix_power(m,mid))
    if G2 > 0:
      res = np.matmul(res,G)
    res = np.matmul(res,powers[post]) #LA.matrix_power(m,post))

    #res = m^pre * G^(G1>0) * m^mid * G^(G2>0) * m^post

    myd[path] = res
    if count % 200 == 0:
      with open(output,'wb') as precomputed_data:
        pickle.dump(myd,precomputed_data)
    count = count + 1

  with open(output,'wb') as precomputed_data:
    pickle.dump(myd,precomputed_data)

  # collate
  i = int(sys.argv[1])
  j = int(sys.argv[2])
  dimension = int(sys.argv[3])
  length = int(sys.argv[4])
  p = sys.argv[5]
  # iterate over all combinations of u and d and all paths and collate them into one data frame
  input_file = file.split(".pickle")[0] + ".precomputed.pickle"
  with open(input_file,'rb') as myd_data:
    myd = pickle.load(myd_data)
  out_text = ""
  for path in myd:
    m = myd[path] 
    line = str(i)+","+str(j)
    line = line+","+path+","
    line = line+",".join([str(x) for x in m[1,:]])+"\n"
    out_text += line
  
  with open(collated_output_file,'w') as output_file:
    output_file.write(out_text)
