var('u,d')
R.<s> = ZZ[]
def f(s): return d+(1-u-d)*s+u*s^2

poly = f(s)
dimension = 400
m = matrix(QQ['u,d'],dimension,2*dimension)
m[0,0] = 1
for i in range(1,dimension):
  print(i)
  my_coefs = poly.coefficients(s,sparse=False)
  for j in range(len(my_coefs)):
    m[i,j] = my_coefs[j]
  poly = poly * f(s)
  outfile = "matrix_"+str(i)
  save(m,outfile)


# find eigenvalues and eigenvectors
#square_m = m[:,0:dimension]

#print(square_m)
#eigenvectors = square_m.eigenvectors_right()
#print(eigenvectors)

#import numpy as np
#outfile = "matrix.npy"
#np.save(outfile, m)

# strategies
# we could precompute all of the matrices to a certain precision
# and store the precomputed values of the matrices
# thats actually a great idea!
# say that we did a grid of 0.01 precision for the u and d values, then there would be 0.5/0.01^2=5000 matrices to store
# for each of these suppose we used a timeline of 100 epochs with up to two rounds of genome doubling, then there would be 100 + 100 + 100*99 = 10100 timelines for each of these matrices
# then there would be a search space of 50.5 million matrices to iterate over. 
# it might be fine to iterate over all these matrices when building them but it would be a lot better to search through them faster
# is there a way to do that?
# can we do assume a unimodal distribution on any of the data? it is definitely unimodal in the u/d grid space. Seen that over and over.
# we can do gradient descent in the u-d space if we can write down the gradient, which we can very easily
# I think we definately precompute the u-d things and then just make a grid like we already did. that is easy. 

# plan is to precompute all of the 5000 matrices AND find their JNF/ EIGENVALUE FORM. Then compute all the timelines for each u/d combo
# 20 epochs probably is enough surely, 20+20+20*19=420 
# is 0.01 square grid the right grid to use?
# a probability of 0.3 and 0.31 will not result in very different values BUT a probability change in 0.01 to 0.02 will, surely?
# apparently ggplot can print a million points quite quickly... in a matter of seconds. 

# ok method right now is to get a large initial matrice
# find the 5000 matrices at the evaluation points
# find the eigenvalue forms of each of these matrices. 
# optimise the likelihood and timelines for specific examples....
# plot them for specific examples. 

# the right grid really does matter... remember that as the number of epochs go up the values for u and d get smaller and smaller... kind of as there is a continuous version of the epochs. 


# do we want to remove the sage stuff "QQ" from this matrice? Turn it into a purely numpy thing?
# need to find the best way to save and store these matrices.

# the most important thing to do here is to come up with all possible products, try and simplify them and save them all
# NOT TRUE - it will be faster to multiply them etc when 
# the integration of these polynomials gives something very similar to a binomial structure. It is very cool, but not really relevant. 


#m_shape = (dimension,dimension*2) 
#diff1 = matrix(QQ,m_shape[0],m_shape[1])
#for i in range(m_shape[0]):
#  print("integrated " +str(i))
#  for j in range(m_shape[1]):
#    diff1[i,j] = diff(diff(m[i,j],u),d)
#
#print(diff1)
#
#from sage.symbolic.integration.integral import definite_integral
#int1 = matrix(QQ,m_shape[0],m_shape[1])
#
#for i in range(m_shape[0]):
#  print("integrated " +str(i))
#  for j in range(m_shape[1]):
#    int1[i,j] = definite_integral(definite_integral(m[i,j],u,0,1-d),d,0,1)
#
#print(int1)
#
## now integrate the product
#int2 = matrix(QQ,m_shape[0],m_shape[1])
#m2 = m[:,0:dimension] * m[:,0:dimension]
#for i in range(m_shape[0]):
#  print("integrated " +str(i))
#  for j in range(m_shape[0]):
#    int2[i,j] = definite_integral(definite_integral(m2[i,j],u,0,1-d),d,0,1)
#
#print(int2)
#
