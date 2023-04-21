# precompute
import sys
import os
import h5py
import numpy as np
from numpy import linalg as LA

args = sys.argv
p_up = args[1]
p_down = args[2]
max_CN = args[3]
path_length = args[4]
path_description = args[5]

collated_output_file = "MATRICES/collated_u"+p_up+"_d"+p_down+"_"+path_description+".csv"
print("START")
print(collated_output_file)

if not os.path.isfile(collated_output_file):
    base_filename = "MATRICES/subbed_mat_u"+p_up+"_d"+p_down+"_"+path_description+".hdf5"
    if not os.path.isfile(base_filename):
        m = load("MATRICES/matrix_"+path_description+".sobj")
        m2 = m.subs(u=int(p_up)/100, d=int(p_down)/100)
        m2 = m2.apply_map(RR)
        m2 = m2.numpy(dtype='double')
        with h5py.File(base_filename, 'w') as m_output:
            m_output.create_dataset('m2', data=m2)

    with h5py.File(base_filename, 'r') as m_data:
        m = np.array(m_data['m2'])

    if not os.path.isfile("GD.hdf5"):
        # Create an empty matrix of zeros
        G = np.zeros((int(max_CN)+1, int(max_CN)+1))

        # Fill in the appropriate entries with ones
        for i in range(round(int(max_CN)/2)):
            G[i, 2*i] = 1
    else:
        with h5py.File("GD.hdf5", 'r') as GD_data:
            G = np.array(GD_data['G'])

    m = m[:(int(max_CN)+2),:(int(max_CN)+2)]
    G = G[:(int(max_CN)+2),:(int(max_CN)+2)]

    for row in range(int(max_CN)+1):
        total = sum(sum(m[row, :]))
        for col in range(int(max_CN)+1):
            if total != 0:
                m[row, col] = m[row, col] / total

    all_paths = load("all_path_combinations_"+path_description+".sobj")

    single_paths = [x for x in all_paths if "G" not in x]

    powers_filename = base_filename.split(".hdf5")[0] + ".powers.hdf5"
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

    precomputed_paths_filename = base_filename.split(".hdf5")[0] + ".precomputed_paths.hdf5"
    count = 0
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
        res = powers[pre]

        if G1 > 0:
            res = np.matmul(res, G)

        res = np.matmul(res, powers[mid])

        if G2 > 0:
            res = np.matmul(res, G)

        res = np.matmul(res, powers[post])

        myd[path] = res
        if count % 200 == 0:
            with h5py.File(precomputed_paths_filename, 'w') as precomputed_data:
                for k, v in myd.items():
                    precomputed_data.create_dataset(k, data=v)
        count = count + 1

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
        line = line + ",".join([str(x) for x in m[1, :]]) + "\n"
        out_text += line

    with open(collated_output_file, 'w') as output_file:
        output_file.write(out_text)

