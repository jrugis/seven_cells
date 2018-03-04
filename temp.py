import numpy as np

l = []

A = np.array(range(20),dtype=int)
B = np.array(range(2,22),dtype=int)
L = np.empty((3,0), dtype=int)

T = np.concatenate([[A],np.full([1,A.shape[0]],2,dtype=int),[B]])
print T

L = np.concatenate((L,T), axis=1)
print L
L = np.concatenate((L,T), axis=1)
print L

l.append(L)
print l

print l[0][0]
