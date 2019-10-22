# -*- coding: utf-8 -*-
#
# msh2bin-all.py

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

#  NOTE: one-based indexing
print("read mesh, write basic bin, write vtu's")
mesh_names = subprocess.check_output("ls *.msh", shell=True).split()
mesh_names = [s.decode("utf-8") for s in mesh_names] # bytes -> string
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print(fname[-4], end=''); sys.stdout.flush()
  verts, tris, tets = rw.read_mesh(fname)
  #rw.write_points("nodes_"+fname, verts)
  #rw.write_tris("surface_"+fname, verts, tris)
  #rw.write_tets("elements_"+fname, verts, tets)
  rw.write_basic_bin("4bsim_"+fname, verts, tris, tets) 

###########################################################################
###########################################################################

