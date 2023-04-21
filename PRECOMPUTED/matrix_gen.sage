import sys

var('u,d')
R.<s> = ZZ[]
def f(s): return d+(1-u-d)*s+u*s^2

# m[i,j] tells you the probability of going from copy number i to copy number j in 1 generation
# m[i,j]^k tells you the probability of going from copy number i to copy number j in k generations

poly = f(s)
path_length = int(sys.argv[1])
dimension = 1+2**(path_length) 
m = matrix(QQ['u,d'],dimension,dimension)

# base case
m[0,0] = 1

# other cases
# WE ONLY NEED TO GO TO ONE POWER LESS THAT THE MAX COPY NUMBER STATE BECAUSE
for i in range(1,2**(path_length-1)+1):
  print(i)
  my_coefs = poly.coefficients(s,sparse=False)
  for j in range(len(my_coefs)):
    m[i,j] = my_coefs[j]
  poly = poly * f(s) 

outfile = "MATRICES/matrix_p"+str(path_length)+"_v4"
save(m,outfile)


