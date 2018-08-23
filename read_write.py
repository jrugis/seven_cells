import numpy as np
import struct
from evtk.hl import pointsToVTK
from evtk.hl import unstructuredGridToVTK
from evtk.vtk import VtkTriangle
from evtk.vtk import VtkTetra

###########################################################################
# functions
###########################################################################

def read_mesh(fname):
  f1 = open(fname + '.msh', 'r') # open the mesh file
  for line in f1: # get the mesh coordinates
    if line.startswith("$Nodes"): break
  nverts = int(f1.next())
  verts = np.empty([nverts, 3])
  for i in range(nverts):
    verts[i] = map(float, f1.next().split()[1:4])
  for line in f1: # get the mesh surface elements
    if line.startswith("$Elements"): break
  nelements = int(f1.next())
  tris = np.empty([nelements, 3], dtype=int) # overkill for now
  tets = np.empty([nelements, 4], dtype=int) # overkill for now
  ntris = 0
  ntets = 0
  for i in range(nelements):
    v = map(int, f1.next().split())
    if v[1] == 2:
      tris[ntris] = v[5:8]
      ntris += 1
    if v[1] == 4:
      tets[ntets] = v[5:9]
      ntets += 1
  tris = tris[0:ntris,:]  # trim the tris
  tets = tets[0:ntets,:]  # trim the tets
  f1.close # close the mesh file
  return verts, tris, tets

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
def read_bin(fname):   ### NEEDS UPDATE!!!!
  f1 = open(fname + '.bin', 'rb') # open the binary file
  # get the vertices
  nverts = struct.unpack('i', f1.read(4))[0]
  verts = np.empty([nverts, 3])
  for i in range(nverts):
    verts[i] = struct.unpack('fff', f1.read(12))
  # get the tris
  ntris = struct.unpack('i', f1.read(4))[0]
  tris = np.empty([ntris, 3], dtype=int)
  dnl = np.empty([ntris])
  for i in range(ntris):
    tris[i] = struct.unpack('iii', f1.read(12))
    dnl[i] = struct.unpack('f', f1.read(4))[0]
  # get the tets
  ntets = struct.unpack('i', f1.read(4))[0]
  tets = np.empty([ntets, 4], dtype=int)
  for i in range(ntets):
    tets[i] = struct.unpack('iiii', f1.read(16))

  f1.close # close the binary file 
  return verts, tris, dnl, tets

###########################################################################
def write_bin(fname, verts, dfa, tris, tets, dfb, apical, basal, common):
  f1 = open(fname + '.bin', 'wb') # create the binary file

  f1.write(struct.pack('i', verts.shape[0]))
  for i,x in enumerate(verts):
    f1.write(struct.pack('ffff', x[0], x[1], x[2], dfa[i]))

  f1.write(struct.pack('i', tris.shape[0]))
  for i,x in enumerate(tris):
    f1.write(struct.pack('iii', x[0], x[1], x[2]))

  f1.write(struct.pack('i', tets.shape[0]))
  for i,x in enumerate(tets):
    f1.write(struct.pack('iiiif', x[0], x[1], x[2], x[3], dfb[i]))

  f1.write(struct.pack('i', apical.shape[0]))
  for i,x in enumerate(apical):
    f1.write(struct.pack('i', x))

  f1.write(struct.pack('i', basal.shape[0]))
  for i,x in enumerate(basal):
    f1.write(struct.pack('i', x))

  f1.write(struct.pack('i', common.shape[1]))
  for i,x in enumerate(np.transpose(common)):
    f1.write(struct.pack('iii', x[0], x[1], x[2]))

  f1.close # close the binary file 
  return

###########################################################################
def write_points(fname, verts, pdata=None):
  nverts = verts.shape[0]
  xyz = np.empty([3, nverts]) # needs re-ordering
  for i in range(nverts):
    for j in range(3):
      xyz[j,i] = verts[i,j]  
  pointsToVTK(fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    data=pdata)  # write out vtu file
  return

###########################################################################
def write_tris(fname, verts, tris, data=None):
  nverts = verts.shape[0]
  xyz = np.empty([3, nverts]) # because it needs re-ordering
  for i in range(nverts):
    for j in range(3):
      xyz[j,i] = verts[i,j]  
  ntris = tris.shape[0]
  conn = np.empty(3*ntris)
  for i in range(ntris):
    conn[3*i] = tris[i,0]-1 # index from zero
    conn[3*i+1] = tris[i,1]-1
    conn[3*i+2] = tris[i,2]-1
  offset = np.zeros(ntris, dtype=int)
  for i in range(ntris):
    offset[i] = 3*(i+1) 
  ctype = np.zeros(ntris)
  for i in range(ntris):
    ctype[i] = VtkTriangle.tid 
  unstructuredGridToVTK(fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    connectivity=conn, offsets=offset, cell_types=ctype, \
    cellData=data, pointData=None)  # write out vtu file
  return

###########################################################################
def write_tets(fname, verts, tets, data=None):
  nverts = verts.shape[0]
  xyz = np.empty([3, nverts]) # because it needs re-ordering
  for i in range(nverts):
    for j in range(3):
      xyz[j,i] = verts[i,j]  
  ntets = tets.shape[0]
  conn = np.empty(4*ntets)
  for i in range(ntets):
    conn[4*i] = tets[i,0]-1 # index from zero
    conn[4*i+1] = tets[i,1]-1
    conn[4*i+2] = tets[i,2]-1
    conn[4*i+3] = tets[i,3]-1
  offset = np.zeros(ntets, dtype=int)
  for i in range(ntets):
    offset[i] = 4*(i+1) 
  ctype = np.zeros(ntets)
  for i in range(ntets):
    ctype[i] = VtkTetra.tid 
  unstructuredGridToVTK(fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    connectivity=conn, offsets=offset, cell_types=ctype, \
    cellData=data, pointData=None)  # write out vtu file
  return

###########################################################################


