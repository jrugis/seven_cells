#! /usr/bin/env python

from evtk.hl import linesToVTK
import numpy as np

###########################################################################
# get verts and lines from data file
def getVertsLines():
  f = open("lumen_lines.txt", 'r')
  n = int(f.next().split()[0]) # vertex count
  v = np.empty((n, 3), dtype=np.float)

  for i in range(n): # verticies
    v[i] = map(float, f.next().split())

  n = int(f.next().split()[0]) # line count
  l = np.empty((n, 2), dtype=np.int)
  for i in range(n): # lines
    l[i] = map(int, f.next().split())
  f.close()
  return v, l  # return (exit verticies, lines)

###########################################################################
# main program
###########################################################################

verts, lines = getVertsLines()
npoints = 2 * lines.shape[0]

x = np.zeros(npoints)
y = np.zeros(npoints)
z = np.zeros(npoints)

for i, l in enumerate(lines):
  for j in range(2):
    x[2*i+j] = verts[l[j]-1,0]
    y[2*i+j] = verts[l[j]-1,1]
    z[2*i+j] = verts[l[j]-1,2]

#pressure = np.random.rand(npoints)
#temp = np.random.rand(npoints)
#vel = np.zeros(2)
#vel[0] = 1.0
#vel[1] = 5.0
#

print x.shape
linesToVTK("lumen_lines", x, y, z, cellData = None, pointData = None)


