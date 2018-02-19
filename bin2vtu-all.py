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

  f1.close # close the binary file 
  return verts, tris, dnl

##################################################################
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

##################################################################
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
def write_apical(fname, verts, tris, dnl, radius):
  conn = np.empty(3*tris.shape[0]) # more than needed...
  c = 0
  for i, t in enumerate(tris):
    if dnl[i] < radius:
      conn[c] = t[0]-1; c+=1 # index from zero
      conn[c] = t[1]-1; c+=1
      conn[c] = t[2]-1; c+=1
  conn = conn[0:c]  # trim the array to actually used

  ntris = c/3
  offset = np.zeros(ntris, dtype=int)
  for i in range(ntris):
    offset[i] = 3*(i+1) 

  ctype = np.zeros(ntris)
  for i in range(ntris):
    ctype[i] = VtkTriangle.tid 

  nverts = verts.shape[0]
  xyz = np.empty([3, nverts])
  for i in range(nverts): # all verts (not just the apical)
    for j in range(3):
      xyz[j,i] = verts[i,j]  

  unstructuredGridToVTK("apical_"+fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    connectivity=conn, offsets=offset, cell_types=ctype, \
    cellData=None, pointData=None)  # write out vtu file

  return

##################################################################
def common_tris(lvertsA, trisA, lvertsB):
  cvi = []
  for v in lvertsB:
    ans = np.where(lvertsA==v)
    if len(ans[0]) == 3:
      cvi.append(ans[0][0]+1)

  cvi = np.array(cvi, dtype=int)
  cti = []
  for t in trisA:
    ans = np.where(cvi == t[0])
    if len(ans[0]) == 0: continue
    ans = np.where(cvi == t[1])
    if len(ans[0]) == 0: continue
    ans = np.where(cvi == t[2])
    if len(ans[0]) == 0: continue
    cti.append(t)

  return np.array(cti, dtype=int)

##################################################################
def write_common(fname, verts, tris):
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

  unstructuredGridToVTK("common_"+fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    connectivity=conn, offsets=offset, cell_types=ctype, \
    cellData=None, pointData=None)  # write out vtu file

  return

##################################################################
# main program
##################################################################

lverts = []  # vertices for each cell
ltris = []       # tris for each cell

mesh_names = subprocess.check_output("ls *tet.bin", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname
  verts, tris, dnl = read_bin(fname)
  lverts.append(verts)
  ltris.append(tris)
  #write_points(fname, verts)
  #write_tris(fname, verts, tris, dnl)
  write_apical(fname, verts, tris, dnl, 0.8)

ncells = len(mesh_names)
for c1 in range(ncells):
  for c2 in range(c1+1, ncells):
    ctris = common_tris(lverts[c1], ltris[c1], lverts[c2])
    if ctris.shape[0] != 0:
      fname = str(c1+1) + str(c2+1)
      print fname
      write_common(fname, lverts[c1], ctris)
print "DONE."


