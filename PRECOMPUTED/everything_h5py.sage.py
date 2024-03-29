

# This file was *autogenerated* from the file everything_h5py.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3); _sage_const_4 = Integer(4); _sage_const_5 = Integer(5); _sage_const_100 = Integer(100); _sage_const_0 = Integer(0); _sage_const_200 = Integer(200)# precompute
import sys
import os
import h5py
import numpy as np
from numpy import linalg as LA

args = sys.argv
p_up = args[_sage_const_1 ]
p_down = args[_sage_const_2 ]
max_CN = args[_sage_const_3 ]
path_length = args[_sage_const_4 ]
path_description = args[_sage_const_5 ]

collated_output_file = "MATRICES/collated_u"+p_up+"_d"+p_down+"_"+path_description+".csv"
print("START")
print(collated_output_file)

if not os.path.isfile(collated_output_file):
    base_filename = "MATRICES/subbed_mat_u"+p_up+"_d"+p_down+"_"+path_description+".hdf5"
    if not os.path.isfile(base_filename):
        m = load("MATRICES/matrix_"+path_description+".sobj")
        m2 = m.subs(u=int(p_up)/_sage_const_100 , d=int(p_down)/_sage_const_100 )
        m2 = m2.apply_map(RR)
        m2 = m2.numpy(dtype='double')
        with h5py.File(base_filename, 'w') as m_output:
            m_output.create_dataset('m2', data=m2)

    with h5py.File(base_filename, 'r') as m_data:
        m = np.array(m_data['m2'])

    if not os.path.isfile("GD.hdf5"):
        # Create an empty matrix of zeros
        G = np.zeros((int(max_CN)+_sage_const_1 , int(max_CN)+_sage_const_1 ))

        # Fill in the appropriate entries with ones
        for i in range(round(int(max_CN)/_sage_const_2 )):
            G[i, _sage_const_2 *i] = _sage_const_1 
    else:
        with h5py.File("GD.hdf5", 'r') as GD_data:
            G = np.array(GD_data['G'])

    m = m[:(int(max_CN)+_sage_const_2 ),:(int(max_CN)+_sage_const_2 )]
    G = G[:(int(max_CN)+_sage_const_2 ),:(int(max_CN)+_sage_const_2 )]

    for row in range(int(max_CN)+_sage_const_1 ):
        total = sum(sum(m[row, :]))
        for col in range(int(max_CN)+_sage_const_1 ):
            if total != _sage_const_0 :
                m[row, col] = m[row, col] / total

    all_paths = load("all_path_combinations_"+path_description+".sobj")

    single_paths = [x for x in all_paths if "G" not in x]

    powers_filename = base_filename.split(".hdf5")[_sage_const_0 ] + ".powers.hdf5"
    if os.path.isfile(powers_filename):
        with h5py.File(powers_filename, 'r') as infile:
            powers = {int(k): np.array(v) for k, v in infile.items()}
    else:
        powers = {}
    for path in single_paths:
        if int(path) not in powers:
            res = LA.matrix_power(m, int(path))
            powers[int(path)] = res
        with h5py.File(powers_filename, 'w') as infile:
            for k, v in powers.items():
                infile.create_dataset(str(k), data=v)

    precomputed_paths_filename = base_filename.split(".hdf5")[_sage_const_0 ] + ".precomputed_paths.hdf5"
    count = _sage_const_0 
    if os.path.isfile(precomputed_paths_filename):
        try:
            with h5py.File(precomputed_paths_filename, 'r') as precomputed_data:
                myd = {k: np.array(v) for k, v in precomputed_data.items()}
        except:
            os.remove(precomputed_paths_filename)
            myd = {}
    else:
        myd = {}


    for path in all_paths:
        if path in myd:
            continue
        splits = path.split("G")
        G1 = _sage_const_0 
        G2 = _sage_const_0 
        if len(splits) == _sage_const_1 :
            pre = int(path)
            mid = _sage_const_0 
            post = _sage_const_0 
        if len(splits) == _sage_const_2 :
            pre = int(splits[_sage_const_0 ])
            mid = _sage_const_0 
            post = int(splits[_sage_const_1 ])
            G1 = _sage_const_1 
        if len(splits) == _sage_const_3 :
            pre = int(splits[_sage_const_0 ])
            mid = int(splits[_sage_const_1 ])
            post = int(splits[_sage_const_2 ])
            G1 = _sage_const_1 
            G2 = _sage_const_1 

        # compute the probabilities for this path:
        res = powers[pre]

        if G1 > _sage_const_0 :
            res = np.matmul(res, G)

        res = np.matmul(res, powers[mid])

        if G2 > _sage_const_0 :
            res = np.matmul(res, G)

        res = np.matmul(res, powers[post])

        myd[path] = res
        if count % _sage_const_200  == _sage_const_0 :
            with h5py.File(precomputed_paths_filename, 'w') as precomputed_data:
                for k, v in myd.items():
                    precomputed_data.create_dataset(k, data=v)
        count = count + _sage_const_1 

    with h5py.File(precomputed_paths_filename, 'w') as precomputed_data:
        for k, v in myd.items():
            precomputed_data.create_dataset(k, data=v)

    # collate
    # iterate over all combinations of u and d and all paths and collate them into one data frame

    with h5py.File(precomputed_paths_filename, 'r') as myd_data:
        myd = {k: np.array(v) for k, v in myd_data.items()}

    out_text = ""

    for path in myd:
        m = myd[path]
        line = str(p_up) + "," + str(p_down)
        line = line + "," + path + ","
        line = line + ",".join([str(x) for x in m[_sage_const_1 , :]]) + "\n"
        out_text += line

    with open(collated_output_file, 'w') as output_file:
        output_file.write(out_text)


