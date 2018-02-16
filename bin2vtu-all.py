import numpy as np
import subprocess
import struct
from evtk.hl import pointsToVTK
from evtk.hl import unstructuredGridToVTK
from evtk.vtk import VtkTriangle

##################################################################
# functions
##################################################################

def read_bin(fname):
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

  #for i in range(nverts):
  #  verts[i] = map(float, f1.next().split()[1:4])
  #  f2.write(struct.pack('fff', verts[i,0], verts[i,1], verts[i,2]))

  # get "distance to nearest lumen" data
  #for line in f1: 
    #if line.startswith('"distance to nearest lumen"'): break
  #  if line.startswith('"nearest cell number"'): break
  #for t in range(6): # skip 6 lines
  #  f1.next()
  #dnl = np.empty(pcount)
  #for t in range(pcount):
  #  v = f1.next().split()
  #  dnl[t] = float(v[1])

  f1.close # close the binary file 
  return verts, tris, dnl

def write_points(fname, verts):
  nverts = verts.shape[0]
  xyz = np.empty([3, nverts]) # needs re-ordering
  for i in range(nverts):
    for j in range(3):
      xyz[j,i] = verts[i,j]  
  pointsToVTK("nodes_"+fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    data=None)  # write out vtu file
  return

def write_tris(fname, verts, tris, dnl):
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

  unstructuredGridToVTK("surface_"+fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    connectivity=conn, offsets=offset, cell_types=ctype, \
    cellData= {"dnl" : dnl}, pointData=None)  # write out vtu file

  return

##################################################################
# main program
##################################################################

mesh_names = subprocess.check_output("ls *tet.bin", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname
  verts, tris, dnl = read_bin(fname)
  write_points(fname, verts)
  write_tris(fname, verts, tris, dnl)
