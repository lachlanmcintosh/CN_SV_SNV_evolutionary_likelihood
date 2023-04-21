# precompute
import sys
import os

args = sys.argv
p_up = args[1]
p_down = args[2]
max_CN = args[3]
path_length = args[4]
path_description = args[5]


collated_output_file = "MATRICES/collated_u"+p_up+"_d"+p_down+"_"+path_description+".csv"
print("START")
print(collated_output_file)

import pickle
import numpy as np
from numpy import linalg as LA

if not os.path.isfile(collated_output_file):
  base_filename = "MATRICES/subbed_mat_u"+p_up+"_d"+p_down+"_"+path_description+".pickle"
  if not os.path.isfile(base_filename):
    m = load("MATRICES/matrix_"+path_description+".sobj")
    m2 = m.subs(u=int(p_up)/100,d=int(p_down)/100)
    m2 = m2.apply_map(RR)
    m2 = m2.numpy(dtype='double')
    with open(base_filename, 'wb') as m_output:
      pickle.dump(m2, m_output)

  with open(base_filename,'rb') as m_data:
    m = pickle.load(m_data)


  if not os.path.isfile("GD.pickle"):
    # Create an empty matrix of zeros
    G = np.zeros((int(max_CN)+1, int(max_CN)+1))

    # Fill in the appropriate entries with ones
    for i in range(round(int(max_CN)/2)):
      G[i, 2*i] = 1

  else:
    with open("GD.pickle",'rb') as GD_data:
      G = pickle.load(GD_data)

  m = m[:(int(max_CN)+2),:(int(max_CN)+2)]
  G = G[:(int(max_CN)+2),:(int(max_CN)+2)]

  for row in range(int(max_CN)+1):
    total = sum(sum(m[row,:]))
    for col in range(int(max_CN)+1):
      if total != 0:
        m[row,col] = m[row,col]/total

  all_paths = load("all_path_combinations_"+path_description+".sobj")
  
  single_paths = [x for x in all_paths if "G" not in x]

  powers_filename = base_filename.split(".pickle")[0] + ".powers.pickle"
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

  precomputed_paths_filename = base_filename.split(".pickle")[0] + ".precomputed_paths.pickle"
  count = 0
  if os.path.isfile(precomputed_paths_filename):
    try:
      with open(precomputed_paths_filename,'rb') as precomputed_data:
        myd = pickle.load(precomputed_data)
    except:
      os.remove(precomputed_paths_filename)
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


    myd[path] = res
    if count % 200 == 0:
      with open(precomputed_paths_filename,'wb') as precomputed_data:
        pickle.dump(myd,precomputed_data)
    count = count + 1

  with open(precomputed_paths_filename,'wb') as precomputed_data:
    pickle.dump(myd,precomputed_data)

  # collate
  # iterate over all combinations of u and d and all paths and collate them into one data frame

  with open(precomputed_paths_filename, 'rb') as myd_data:
    myd = pickle.load(myd_data)

  out_text = ""

  for path in myd:
    m = myd[path]
    line = str(p_up)+","+str(p_down)
    line = line+","+path+","
    line = line+",".join([str(x) for x in m[1,:]])+"\n"
    out_text += line

  with open(collated_output_file,'w') as output_file:
    output_file.write(out_text)
