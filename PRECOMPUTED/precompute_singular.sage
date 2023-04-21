

import sys

args = sys.argv
i=int(args[1])
j=int(args[2])
n=int(args[3])
m=load("matrix_"+str(n)+".sobj")
m2 = m.subs(u=i/100,d=j/100)
outfile = "precomputed/pre_mat"+str(n)+"_u"+str(i)+"_d"+str(j)
save(m2,outfile)
