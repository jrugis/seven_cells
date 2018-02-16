import numpy as np
import subprocess
import struct

###########################################################################
# functions
###########################################################################
def getDistQ(A, B, P): # get distance from point P to line segment AB
  AB = B-A
  AP = P-A
  dAP = np.linalg.norm(AP)
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
def read_lumen():
  f1 = open('lumen_lines.txt', 'r') # open the lumen file

  nlverts = int(f1.next().split()[0]) 
  lverts = np.empty([nlverts, 3])
  for i in range(nlverts): # get the lumen verts
    lverts[i] = map(float, f1.next().split()[0:3])

  nlines = int(f1.next().split()[0])
  lines = np.empty([nlines, 2], dtype=int)
  for i in range(nlines): # get the lumen line segments
    lines[i] = map(int, f1.next().split()[0:2])

  lsegs = np.empty([nlines, 6])
  for i in range(nlines): # get the lumen line segments
    v1 = map(float, lverts[lines[i,0]-1])
    v2 = map(float, lverts[lines[i,1]-1])
    lsegs[i] = np.array([v1, v2]).flatten()

  f1.close()
  return lsegs

###########################################################################
def get_dnl(tris, verts, lsegs):
  ntris = tris.shape[0]
  dnl = np.empty([ntris])
  for i in range(ntris):
    C = np.average(verts[tris[i]-1], axis=0) # center of tri
    d = 100.0 # large dummy initial distance
    for seg in lsegs:
      ds = getDistQ(seg[0:3], seg[3:6], C)
      if ds < d: 
        d = ds
    dnl[i] = d # save the minimum distance
  return dnl

###########################################################################
def read_mesh(fname):
  f1 = open(fname + '.msh', 'r') # open the mesh file

  # get the mesh coordinates
  for line in f1: 
    if line.startswith("$Nodes"): break
  nverts = int(f1.next())
  verts = np.empty([nverts, 3])
  for i in range(nverts):
    verts[i] = map(float, f1.next().split()[1:4])

  # get the mesh surface elements
  for line in f1: 
    if line.startswith("$Elements"): break
  nelements = int(f1.next())
  tris = np.empty([nelements, 3], dtype=int) # overkill for now
  for i, line in enumerate(f1):
    v = map(int, line.split())
    if v[1] == 4:
      tris = tris[0:i,:]  # trim the tris
      break
    tris[i] = v[5:8]
  f1.close # close the mesh file
  return verts, tris

###########################################################################
def write_bin(fname, verts, tris, dnl):
  f1 = open(fname + '.bin', 'wb') # create the binary file
  nverts = verts.shape[0]
  f1.write(struct.pack('i', nverts))
  for v in verts:
    f1.write(struct.pack('fff', v[0], v[1], v[2]))
  ntris = tris.shape[0]
  f1.write(struct.pack('i', ntris))
  for i,t in enumerate(tris):
    f1.write(struct.pack('iiif', t[0], t[1], t[2], dnl[i]))

  f1.close # close the binary file 
  return

###########################################################################
# main program
###########################################################################

lsegs = read_lumen()
mesh_names = subprocess.check_output("ls *tet.msh", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname
  verts, tris = read_mesh(fname)
  dnl = get_dnl(tris, verts, lsegs)
  write_bin(fname, verts, tris, dnl)

