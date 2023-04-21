base="/export/share/prkfs2/shared/bioinf-data/Papenfuss_lab/projects/lachlan/PHD/GD_matrices3/precomputed"

from glob import glob
import sys
i = int(sys.argv[1])
j = int(sys.argv[2])
dimension = int(sys.argv[3])
length = int(sys.argv[4])
p = sys.argv[5]
# iterate over all combinations of u and d and all paths and collate them into one data frame
output = open("/export/share/prkfs2/shared/bioinf-data/Papenfuss_lab/projects/lachlan/PHD/GD_matrices/collated_"+str(i)+"_"+str(j)+"_"+str(dimension)+"_"+p+".csv",'w')
print(i)
print(j)
print(dimension)
print(p)
file_search = base+"/pre_mat"+str(dimension)+"_"+p+"_u"+str(i)+"_d"+str(j)+"/*"
print(file_search)

all_paths = 

for file in glob(file_search,recursive=False):
    print(file)
    file = str(file)
    history = file.split(".precomputed.sobj")[0].split(base)[1].split("/")[2]
    print(file)
    print(history)
    m = load(file)
    line = str(i)+","+str(j)
    line = line+","+history+","
    line = line+",".join([str(RR(x)) for x in m[1,:][0] ])+"\n"
    output.write(line)
