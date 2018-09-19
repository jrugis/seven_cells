import numpy as np
import subprocess
import struct
import sys
import read_write as rw  # read_write.py

###########################################################################
# main program
###########################################################################

# lists of data for each cell 
#  NOTE: one-based indexing
lverts   = []  # vertices
ltris    = []  # tris
ltets    = []  # tets

mesh_names = subprocess.check_output("ls *.msh", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname[-4]

  verts, tris, tets = rw.read_mesh(fname)

  rvertsi = np.array([range(verts.shape[0])], dtype=int) + 1 # all vert indices
  rvertsi = np.setdiff1d(rvertsi, tris) # remove surface tri indices
  rverts = verts[rvertsi-1] # all the non-surface verts

  nrverts = rverts - np.min(rverts, axis=0) # normalise all verts to non-negative 
  max = (np.max(nrverts, axis=0)) # get the range of vertex values

  # create a 4D grid (as an array) for extracting a uniform spatial vertex subset
  # - stores distance to nearest vertex (dnv) and the associated vertex index at each grid point
  # - the integer parts of every vertex coordinate are used to index the grid 
  vgrid = np.zeros((np.concatenate((np.floor(max+1),[2])).astype(int)))
  toohigh = 1000000 # high dummy values for dnv and index
  vgrid.fill(toohigh) 

  ifverts = np.modf(nrverts) # get the integer and fractional part of all nrverts

  # iterate through all verts and store the vert that is closest to each (integer) grid point
  for i in range(nrverts.shape[0]):
    dist = np.linalg.norm(ifverts[0][i])
    vgridi = (ifverts[1][i]).astype(int) # grid index is simply the vertex location integer part
    if dist > 1: # don't bother with grid points that have no close vertex
      continue

    # is this vertex closest to the grid point?
    noise = (2.0 * np.random.ranf()) - 1.0 # some spatial dithering to break up an aligned visual
    if (dist + noise) < vgrid[vgridi[0]][vgridi[1]][vgridi[2]][0]: 
      vgrid[vgridi[0]][vgridi[1]][vgridi[2]][0] = dist # store the new closer distance
      vgrid[vgridi[0]][vgridi[1]][vgridi[2]][1] = i # update the associated vertex index

  # extract the close vertex indices
  cvi = vgrid[:,:,:,1]
  cvi = np.extract(cvi < toohigh, cvi)
  print "Vertex count reduction:", verts.shape[0], "->", cvi.shape[0]

  rw.write_points("reduced-nodes_"+fname, rverts[cvi.astype(int)])
  rw.write_indicies("reduced-indices_"+fname, cvi)

###########################################################################
###########################################################################

