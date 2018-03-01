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
  pointsToVTK(fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    data=None)  # write out vtu file
  return

##################################################################
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

##################################################################
def get_apical(tris, dnl, small_radius, large_radius):
  atrisS = []
  atrisL = []
  for i, d in enumerate(dnl):
    if d < small_radius: atrisS.append(i+1)
    if d < large_radius: atrisL.append(i+1)
  return(np.array(atrisS, dtype=int), np.array(atrisL, dtype=int))

##################################################################
def get_basal(tris, ctrisi, atrisi):
  btrisi = np.array([range(tris.shape[0])], dtype=int) + 1 # all indices
  btrisi = np.setdiff1d(btrisi, ctrisi) # remove common
  btrisi = np.setdiff1d(btrisi, atrisi) # remove apical Large
  return(btrisi)

##################################################################
def find_tris(vi, tris): # vertex indices, tris
  ti = []
  for i, t in enumerate(tris): # find tris associated with verticies in list
    ans = np.where(vi == t[0]) # first vertex in list?
    if len(ans[0]) == 0: continue
    ans = np.where(vi == t[1]) # second vertex in list?
    if len(ans[0]) == 0: continue
    ans = np.where(vi == t[2]) # third vertex in list?
    if len(ans[0]) == 0: continue
    ti.append(i+1)
  return np.array(ti, dtype=int) # return tri indicies

##################################################################
def get_common(vertsA, trisA, vertsB):
  cviA = []
  for v in vertsB: # first find common vertices 
    ans = np.where(np.linalg.norm(vertsA-v,axis=1) < 0.001)
    if len(ans[0]) != 0:
      if len(ans[0]) > 1:
        print("ERROR: duplicate vertices")
      cviA.append(ans[0][0]+1)
  cviA = np.array(cviA, dtype=int)
  return(find_tris(cviA,trisA)) # return tri indicies

##################################################################
# main program
##################################################################

# lists of data for each cell 
#  NOTE: one-based indexing
lverts   = []  # vertices
ltris    = []  # tris
latrisiS = []  # apical tri indices Small
latrisiL = []  # apical tri indices Large
lctrisi  = []  # common tri indices
lbtrisi  = []  # basal tri indices

print("individual cell calcs")
mesh_names = subprocess.check_output("ls *tet.bin", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname[-4]
  verts, tris, dnl = read_bin(fname)
  atrisiS, atrisiL = get_apical(tris, dnl, 0.8, 1.4) # apical distances
  write_points("nodes_"+fname, verts)
  write_tris("surface_"+fname, verts, tris, data={"dnl" : dnl})
  write_tris("apical_"+fname, verts, tris[atrisiS-1]) # has all the verts
  lverts.append(verts)
  ltris.append(tris)
  latrisiS.append(atrisiS)
  latrisiL.append(atrisiL)

print("cell-to-cell calcs")
ncells = len(mesh_names)
for c1 in range(ncells):
  ctrisi = np.array([], dtype=int) # common tri indices per cell
  for c2 in range(ncells):
    if(c1==c2): continue
    print c1+1, c2+1
    cti = get_common(lverts[c1], ltris[c1], lverts[c2]) # pair-wise common
    if cti.shape[0] != 0:
      ctrisi = np.concatenate([ctrisi,cti])
      if(c2 > c1):
        fname = "common_" + str(c1+1) + "c" + str(c2+1) + "c"
        write_tris(fname, lverts[c1], ltris[c1][cti-1]) # has all the verts
  lctrisi.append(ctrisi)

print("additional cell calcs")
for i, mesh in enumerate(mesh_names):
  print str(i+1)
  fname = mesh.split('.')[0]
  btrisi = get_basal(ltris[i], lctrisi[i], latrisiL[i])
  write_tris("basal_"+fname, lverts[i], ltris[i][btrisi-1]) # has all the verts
  lbtrisi.append(btrisi)

print "DONE."

##################################################################
##################################################################

