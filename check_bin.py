# -*- coding: utf-8 -*-
#
# check_bin.py

import numpy as np
import subprocess
import struct
import sys
from evtk.hl import pointsToVTK
from evtk.hl import unstructuredGridToVTK
from evtk.vtk import VtkTriangle
from evtk.vtk import VtkTetra
import read_write as rw  # read_write.py
import utils as ut       # utils.py

###########################################################################
# main program
###########################################################################

fname = '4sim_out_N4_p3-p2-p4-'
cnum = 1
c2cnum = 2

verts, v2a, tris, tets, t2a, t2b, atris, btris, ctris \
  = rw.read_bin(fname + str(cnum) + 'tet')

print('vertices', end=' ')
print(verts.shape)
#print(np.amin(verts, axis=0))
#print(np.amax(verts, axis=0))

print('v2a', end=' ')
print(v2a.shape)
#print(min(v2a))
#print(max(v2a))

print('triangles', end=' ')
print(tris.shape)
#print(np.amin(tris))
#print(np.amax(tris))

print('tetrahedrons', end=' ')
print(tets.shape)
#print(np.amin(tets))
#print(np.amax(tets))

print('t2a', end=' ')
print(t2a.shape)
#print(min(t2a))
#print(max(t2a))

print('t2b', end=' ')
print(t2b.shape)
#print(min(t2b))
#print(max(t2b))

print('apical triangles', end=' ')
print(atris.shape)
#print(np.amin(atris))
#print(np.amax(atris))

print('basal triangles', end=' ')
print(btris.shape)
#print(np.amin(btris))
#print(np.amax(btris))

print('common triangles', end=' ')
print(ctris.shape)
#print(np.amin(ctris, axis=0))
#print(np.amax(ctris, axis=0))

t1 = ctris[:,1]
t2 = np.nonzero(t1 == c2cnum)
t3 = ctris[t2]
t4 = t3[:,2]
prev = t4[0]
for i in range(1,len(t4)):
  current = t4[i]
  if (current == prev + 1): print('.', end='')
  else: print('x', end='')
  prev = current
print()

print(t3[range(20)])