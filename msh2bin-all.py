import numpy as np
import subprocess
import struct

##################################################################
# functions
##################################################################

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

def write_bin(fname, verts, tris):
  f1 = open(fname + '.bin', 'wb') # create the binary file
  nverts = verts.shape[0]
  f1.write(struct.pack('i', nverts))
  for v in verts:
    f1.write(struct.pack('fff', v[0], v[1], v[2]))
  ntris = tris.shape[0]
  f1.write(struct.pack('i', ntris))
  for t in tris:
    f1.write(struct.pack('iii', t[0], t[1], t[2]))

  f1.close # close the binary file 
  return

##################################################################
# main program
##################################################################

mesh_names = subprocess.check_output("ls *tet.msh", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname
  verts, tris = read_mesh(fname)
  write_bin(fname, verts, tris)

