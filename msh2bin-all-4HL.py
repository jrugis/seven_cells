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
  rvertsi = np.array([range(1, verts.shape[0]+1)], dtype=int) # all vert indices
  rvertsi = np.setdiff1d(rvertsi, tris) # remove surface tri indices
  rvertsi -= 1; # change to zero indexed
  nverts = 0.6 * verts # node reduction factor

  nverts = nverts - np.min(nverts, axis=0) # normalise all verts to non-negative 
  max = (np.max(nverts, axis=0)) # get the range of vertex values

  # create a 4D grid (as an array) for extracting a uniform spatial vertex subset
  # - stores distance to nearest vertex (dnv) and the associated vertex index at each grid point
  # - the integer parts of every vertex coordinate are used to index the grid 
  vgrid = np.zeros((np.concatenate((np.floor(max+1),[2])).astype(int)))
  toohigh = 1000000 # high dummy values for dnv and index
  vgrid.fill(toohigh) 

  ifverts = np.modf(nverts) # get the integer and fractional part of all nverts

  # iterate through rvertsi and store the vert that is closest to each (integer) grid point
  for i in rvertsi:
    dist = np.linalg.norm(ifverts[0][i])
    vgridi = (ifverts[1][i]).astype(int) # grid index is simply the vertex location integer part
    if dist > 0.5: # don't bother with grid points that have no close vertex
      continue
    noise = 0.3 * np.random.ranf() # some spatial dithering to break up an aligned visual
    dist += noise

    # is this vertex closest to the grid point?
    if dist < vgrid[vgridi[0]][vgridi[1]][vgridi[2]][0]: 
      vgrid[vgridi[0]][vgridi[1]][vgridi[2]][0] = dist # store the new closer distance
      vgrid[vgridi[0]][vgridi[1]][vgridi[2]][1] = i # update the associated vertex index

  # extract the close vertex indices
  cvi = vgrid[:,:,:,1]
  cvi = np.extract(cvi < toohigh, cvi)
  print "Vertex count reduction:", verts.shape[0], "->", cvi.shape[0]

  rw.write_points("reduced-nodes_"+fname, verts[cvi.astype(int)])
  rw.write_indicies("reduced-indices_"+fname, cvi)

###########################################################################
###########################################################################

