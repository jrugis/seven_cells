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
  nverts = int(next(f1))
  verts = np.empty([nverts, 3])
  for i in range(nverts):
    verts[i] = [float(x) for x in next(f1).split()[1:4]]
  for line in f1: # get the mesh surface elements
    if line.startswith("$Elements"): break
  nelements = int(next(f1))
  tris = np.empty([nelements, 3], dtype=int) # overkill for now
  tets = np.empty([nelements, 4], dtype=int) # overkill for now
  ntris = 0
  ntets = 0
  for i in range(nelements):
    v = [int(x) for x in next(f1).split()]
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
  f1 = open('lumen_lines.txt', 'r+') # open the lumen file

  nlverts = int(next(f1).split()[0]) 
  lverts = np.empty([nlverts, 3])
  for i in range(nlverts): # get the lumen verts
    lverts[i] = [float(x) for x in next(f1).split()[0:3]]

  nlines = int(next(f1).split()[0])
  lines = np.empty([nlines, 2], dtype=int)
  for i in range(nlines): # get the lumen line segments
    lines[i] = [int(x) for x in next(f1).split()[0:2]]

  lsegs = np.empty([nlines, 6])
  for i in range(nlines): # get the lumen line segments
    v1 = [float(x) for x in lverts[lines[i,0]-1]]
    v2 = [float(x) for x in lverts[lines[i,1]-1]]
    lsegs[i] = np.array([v1, v2]).flatten()

  f1.close()
  return lsegs

###########################################################################
def read_bin(fname):   ### NOTE: one based indexing !!!
  f1 = open(fname + '.bin', 'rb') # open the binary file
  # get the vertices
  nverts = struct.unpack('i', f1.read(4))[0]
  verts = np.empty([nverts, 3])
  v2a = np.empty([nverts, 1])
  for i in range(nverts):
    verts[i] = struct.unpack('fff', f1.read(12))
    v2a[i] = struct.unpack('f', f1.read(4))
  # get the tris
  ntris = struct.unpack('i', f1.read(4))[0]
  tris = np.empty([ntris, 3], dtype=int)
  for i in range(ntris):
    tris[i] = struct.unpack('iii', f1.read(12))
  # get the tets
  ntets = struct.unpack('i', f1.read(4))[0]
  tets = np.empty([ntets, 4], dtype=int)
  t2a = np.empty([ntets, 1])
  t2b = np.empty([ntets, 1])
  for i in range(ntets):
    tets[i] = struct.unpack('iiii', f1.read(16))
    t2a[i] = struct.unpack('f', f1.read(4))
    t2b[i] = struct.unpack('f', f1.read(4))
  # get the apical tris
  ntris = struct.unpack('i', f1.read(4))[0]
  atris = np.empty([ntris, 1], dtype=int)
  for i in range(ntris):
    atris[i] = struct.unpack('i', f1.read(4))
  # get the bascal tris
  ntris = struct.unpack('i', f1.read(4))[0]
  btris = np.empty([ntris, 1], dtype=int)
  for i in range(ntris):
    btris[i] = struct.unpack('i', f1.read(4))
  # get the common tris
  ntris = struct.unpack('i', f1.read(4))[0]
  ctris = np.empty([ntris, 3], dtype=int)
  for i in range(ntris):
    ctris[i] = struct.unpack('iii', f1.read(12))

  f1.close # close the binary file 
  return verts, v2a, tris, tets, t2a, t2b, atris, btris, ctris

###########################################################################
def read_basic_bin(fname):   ### NOTE: one based indexing !!!
  f1 = open(fname + '.bin', 'rb') # open the binary file
  # get the vertices
  nverts = struct.unpack('i', f1.read(4))[0]
  verts = np.empty([nverts, 3])
  for i in range(nverts):
    verts[i] = struct.unpack('ddd', f1.read(24))
  # get the tris
  ntris = struct.unpack('i', f1.read(4))[0]
  tris = np.empty([ntris, 3], dtype=int)
  for i in range(ntris):
    tris[i] = struct.unpack('iii', f1.read(12))
  # get the tets
  ntets = struct.unpack('i', f1.read(4))[0]
  tets = np.empty([ntets, 4], dtype=int)
  for i in range(ntets):
    tets[i] = struct.unpack('iiii', f1.read(16))

  f1.close # close the binary file 
  return verts, tris, tets

###########################################################################
def write_bin(fname, verts, dfa_vert, tris, tets, dfa_tet, dfb, apical, basal, common):
  f1 = open(fname + '.bin', 'wb') # create the binary file

  f1.write(struct.pack('i', verts.shape[0]))
  for i,x in enumerate(verts):
    f1.write(struct.pack('ffff', x[0], x[1], x[2], dfa_vert[i]))

  f1.write(struct.pack('i', tris.shape[0]))
  for i,x in enumerate(tris):
    f1.write(struct.pack('iii', x[0], x[1], x[2]))

  f1.write(struct.pack('i', tets.shape[0]))
  for i,x in enumerate(tets):
    f1.write(struct.pack('iiiiff', x[0], x[1], x[2], x[3], dfa_tet[i], dfb[i]))

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
def write_basic_bin(fname, verts, tris, tets):
  f1 = open(fname + '.bin', 'wb') # create the binary file

  f1.write(struct.pack('i', verts.shape[0]))
  for i,x in enumerate(verts):
    f1.write(struct.pack('ddd', x[0], x[1], x[2],))

  f1.write(struct.pack('i', tris.shape[0]))
  for i,x in enumerate(tris):
    f1.write(struct.pack('iii', x[0], x[1], x[2]))

  f1.write(struct.pack('i', tets.shape[0]))
  for i,x in enumerate(tets):
    f1.write(struct.pack('iiii', x[0], x[1], x[2], x[3]))

  f1.close # close the binary file 
  return

###########################################################################
def write_indicies(fname, idx):
  f1 = open(fname + '.bin', 'wb') # create the binary file
  f1.write(struct.pack('i', idx.shape[0]))

  for n in idx:
    f1.write(struct.pack('i', n))

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


