from random_methods import *

a = 2000 * np.ones(10000,dtype=np.int)
d = poisson_mutation(a, 30, 33, 1e-6, tag="UID")

print(d.head(10))