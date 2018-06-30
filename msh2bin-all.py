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

# lists of data for each cell 
#  NOTE: one-based indexing
lverts   = []  # vertices
ltris    = []  # tris
ltets    = []  # tets
ldfb     = []  # dfb
latrisiS = []  # apical tri indices Small
latrisiL = []  # apical tri indices Large
lctrisi  = []  # common tri indices

print("read mesh, calc dnl/apical, write nodes/surface/apical")
lsegs = rw.read_lumen()
mesh_names = subprocess.check_output("ls *.msh", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname[-4],; sys.stdout.flush()
  verts, tris, tets = rw.read_mesh(fname)
  dnl = ut.get_dnl(tris, verts, lsegs)
  atrisiS, atrisiL = ut.get_apical(tris, dnl, 0.8, 1.4) # apical distances
#  rw.write_points("nodes_"+fname, verts)
#  rw.write_tris("surface_"+fname, verts, tris, data={"dnl" : dnl})
#  rw.write_tris("apical_"+fname, verts, tris[atrisiS-1]) # has all the verts
  lverts.append(verts)
  ltris.append(tris)
  ltets.append(tets)
  latrisiS.append(atrisiS)
  latrisiL.append(atrisiL)

print("\ncalc/write common")
ncells = len(mesh_names)
for c1 in range(ncells):
  ctrisi = np.empty((3,0), dtype=int) # common tri indices per cell
  for c2 in range(ncells):
    if(c1==c2): continue
    print str(c1+1)+'/'+str(c2+1),; sys.stdout.flush()
    ctiA,ctiB = ut.get_common(lverts[c1], ltris[c1], lverts[c2], ltris[c2]) # pair-wise common
    if ctiA.shape[0] != 0:
      T = np.concatenate([[ctiA],np.full([1,ctiA.shape[0]],c2+1,dtype=int),[ctiB]])
      ctrisi = np.concatenate((ctrisi,T), axis=1)
      if(c2 > c1):
        fname = "common_" + str(c1+1) + "c" + str(c2+1) + "c"
#        rw.write_tris(fname, lverts[c1], ltris[c1][ctiA-1]) # has all the verts
  lctrisi.append(ctrisi)

print("\ncalc basal, calc dfb, write basal/elements/binary")
for i, mesh in enumerate(mesh_names):
  print str(i+1),; sys.stdout.flush()
  fname = mesh.split('.')[0]
  btrisi = ut.get_basal(ltris[i], lctrisi[i][0], latrisiL[i])
  dfb = ut.get_dfb(btrisi, lverts[i], ltets[i])
#  ldfb.append(dfb)
#  rw.write_tris("basal_"+fname, lverts[i], ltris[i][btrisi-1]) # has all the verts
  rw.write_tets("elements_"+fname, lverts[i], ltets[i])
#  rw.write_bin("4sim_"+fname, lverts[i], ltris[i], 
#    ltets[i], ldfb[i], latrisiS[i], btrisi, lctrisi[i]) 

print "\nDONE."

###########################################################################
###########################################################################

