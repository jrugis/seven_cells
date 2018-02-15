###########################################################################
# Create seven cells lumen tree data
# J. Rugis
# 13.11.17
#
# input file "tree.dat" from blender "model_20d.blend"
# output file "tree.txt" contains vertices, lines, cell counts per line
# output file "tree.vtu" for checking with paraview
#
###########################################################################

import math
import numpy as np
from evtk.hl import pointsToVTK

###########################################################################
# functions
###########################################################################
def getDistQ(A, B, P): # get distance from point P to line segment AB
  AB = B-A
  AP = P-A
  dAP = np.linalg.norm(AP)
  if dAP > 5.0: return 100.0 # cull (too far away)
  AQ = (np.dot(AB, AP)) * AB / (AB**2).sum() 
  Q = AQ + A                 # QP is perpendicular to line AB
  d2AQ = (AQ**2).sum()
  d2BQ = ((Q-B)**2).sum()
  d2AB = (AB**2).sum()
  if((d2AQ > d2AB) or (d2BQ > d2AB)): # Q not on segment AB?
    d = min(dAP, np.linalg.norm(P-B)) # use distance to A or B
  else:
    d = np.linalg.norm(P-Q)           # use perpendicular distance
  return d

###########################################################################
def getCloseSeg(tv, tl, P):
  cd = 10.0 # closest line distance
  cs = 0    # closest segment index
  for i, l in enumerate(tl):
    d = getDistQ(tv[l[0]-1], tv[l[1]-1], P)
    if d < cd:
      cs = i
      cd = d
  return cs

###########################################################################
def getCellCounts(tv, tl, sn): # tverts, tlines, snodes
  cc = np.zeros([len(tl), 8]) # counts for each of seven cells per line
  nl = np.zeros(len(sn[0]))
  for i, p in enumerate(sn[0]):
    j = getCloseSeg(tv, tl, p[0:3])
    cc[j, int(p[3])] += 1.0 # add cell count to line
    nl[i] = j + 1
  return cc, nl  

###########################################################################
def getSurfNodes(): # get idealised lumen surface nodes
  f = open("tubes_20c.msh", 'r') # open the mesh file

  for line in f: # get the node coordinates
    if line.startswith("$Nodes"): break
  n = int(f.next())
  xyz = np.zeros([n, 5])
  for t in range(n):
    v = f.next().split()
    xyz[t,0:3] = map(float, v[1:4])

  for line in f:   # get the surface nodes
    if line.startswith("$Elements"): break
  n = int(f.next())
  for t in range(n): # identify the surface nodes
    v = f.next().split()
    if(int(v[1]) == 2): # surface triangle?
      for i in range(5,8):
        xyz[int(v[i])-1, 4] = 1.0 # mark the node
    else: break
  mk = np.where(xyz[:,4]==1.0) # surface vertex indices

  for line in f: # get the nodes cell label
    if line.startswith('"nearest cell number"'): break
  for t in range(5): f.next()
  n = int(f.next())
  for t in range(n): # get the cell number
    v = f.next().split()
    xyz[t,3] = float(v[1])
  f.close()
  return xyz[mk,0:4] # surface vertices with cell number

###########################################################################
# calculate rotation matrix from euler angles
def euler2rotation(t):
  x = np.array([[1, 0, 0],
                [0, math.cos(t[0]), -math.sin(t[0])],
                [0, math.sin(t[0]), math.cos(t[0])]])
  y = np.array([[math.cos(t[1]), 0, math.sin(t[1])],
                [0, 1, 0],
                [-math.sin(t[1]), 0, math.cos(t[1])]])
  z = np.array([[math.cos(t[2]), -math.sin(t[2]), 0],
                [math.sin(t[2]), math.cos(t[2]), 0],
                [0, 0, 1]])
  R = np.dot(z, np.dot(y, x))
  return R.transpose()

###########################################################################
# get index of closest vert to point
def getiVert(p, vs): # point, vertices
  d = 100.0 # lowest distance dummy value
  for i, v in enumerate(vs):
    nd = np.linalg.norm(p - v)
    if(nd < d): # lower distance?
      id = i 
      d = nd 
  return id+1 # index of closest vertex 

###########################################################################
# get verts and calculate lines from data file
def getVertsLines():
  f = open("tree.dat", 'r')
  n = int(f.next()) # vertex count
  v = np.empty((n, 3), dtype=np.float)
  xv = 0
  for i in range(n): # verticies
    vals = map(float, f.next().split()[0:4])
    v[i] = vals[0:3]
    if(vals[3] >= 0.05): # the lumen exit vertex, HARD CODED!!!
      xv = i+1
  n = int(f.next()) # line count
  l = np.empty((n, 2), dtype=np.int)
  for i in range(n): # lines
    vals = map(float, f.next().split()[0:7])
    c = (vals[0], vals[1], vals[2])   # center of line segment
    rm = euler2rotation(vals[3:6])    # line segment rotation
    v1 = c + np.dot((0, 0, vals[6]/2.0), rm)  # line segment endpoints
    v2 = c + np.dot((0, 0, -vals[6]/2.0), rm)
    iv = (getiVert(v1, v), getiVert(v2, v)) # verticies closest to endpoints
    l[i] = (iv)
  f.close()
  return xv, v, l # return (exit vertex index, verticies, lines)

###########################################################################
def saveVtk(v, nln):
  v = np.transpose(v)
  x = np.array(v[0,:]) # deep copy required
  y = np.array(v[1,:])
  z = np.array(v[2,:])
  ncn = np.array(v[3,:])
  d = {}
  d["ncn"] = ncn  # nearest cell number
  d["nln"] = nln  # nearest line number
  pointsToVTK("./tree", x, y, z, data = d) # write out vtk file
  return

###########################################################################
def saveTree(tv, tl, cc):
  f = open("tree.txt", 'w')
  f.write("%d %d (vertices, exit vertex)\n" % (len(tverts), txvert))
  for v in tv:
    f.write("%5.3f %5.3f %5.3f\n" % (v[0], v[1], v[2]))
  f.write("%d (lines)\n" % len(tlines))
  for v in tl:
    f.write("%d %d\n" % (v[0], v[1]))
  f.write("%d (cell counts per line)\n" % len(tlines))
  for v in cc:
    f.write("%d %d %d %d %d %d %d %d\n" %\
      (v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7]))
  f.close()
  return

###########################################################################
# main program
###########################################################################

print "get tree data"
txvert, tverts, tlines = getVertsLines()

print "get lumen surface nodes"
snodes = getSurfNodes()

print "get tree cell counts and node lines"
cellcnts, nlines = getCellCounts(tverts, tlines, snodes)

print "save tree"
saveTree(tverts, tlines, cellcnts)

print "save vtk file"
saveVtk(snodes[0], nlines)

###########################################################################

