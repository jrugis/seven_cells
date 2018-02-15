# -*- coding: utf-8 -*-
'''
Created on Mon Oct  6 14:17:23 2014
@author: jrugis
'''

import sys
import math
import numpy as np

for p in range(1,8):
  cfname = "out_N4_p3-p2-p4-" + str(p) + "tet.msh"
  tfname = "tubes_20c.msh"
  print cfname, tfname

  # read in tubes coordinates
  f = open(tfname)
  while f.readline().strip() != '$Nodes': None # skip to the $Nodes section
  n = int(f.readline().strip()) # get the number of nodes
  print 'tubes node count:', n
  tubes = np.zeros((n,3))
  for i in range(n):  # get the node coordinate values
    tubes[i] = f.readline().split()[1:4]
  f.close()

  # read in cell coordinates
  f = open(cfname)
  while f.readline().strip() != '$Nodes': None # skip to the $Nodes section
  n = int(f.readline().strip()) # get the number of nodes
  print 'cell node count:', n
  cell = np.zeros((n,3))
  for i in range(n):  # get the node coordinate values
    cell[i] = f.readline().split()[1:4]
  f.close()

  # calculate a 'stretched' cell bounding box
  cmin = np.amin(cell,0)
  cmax = np.amax(cell,0)
  d = (cmax - cmin) / 10  # stretch factor
  cmin -= d
  cmax += d

  # create a trimmed tubes array using the cell bounding box
  tb = []
  for j in range(tubes.shape[0]):
    if (tubes[j,0] >= cmin[0]) and (tubes[j,0] <= cmax[0]):
      if (tubes[j,1] >= cmin[1]) and (tubes[j,1] <= cmax[1]):
        if (tubes[j,2] >= cmin[2]) and (tubes[j,2] <= cmax[2]):
          tb.append([tubes[j,0], tubes[j,1], tubes[j,2]])
  tb = np.array(tb)
  print 'trimmed tubes node count:', tb.shape[0]

  # write out the distance from each cell node 
  #   to the nearest tube node within the cell's bounding box
  f = open(cfname, 'a')
  f.write('$NodeData\n')
  f.write('1\n')                             # one string tag
  f.write('"distance to nearest lumen"\n')   #   name
  f.write('1\n')                             # one real tag
  f.write('0.0\n')                           #   time value
  f.write('3\n')                             # three integer tags
  f.write('0\n')                             #   time step
  f.write('1\n')                             #   one component (scalar)
  f.write('%s\n' % str(cell.shape[0]))       #   node count
  for i in range(cell.shape[0]):
    dmin = sys.float_info.max
    for j in range(tb.shape[0]):
      dxyz = cell[i] - tb[j]
      d = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2]
      if d < dmin: dmin = d
    f.write('%d %f\n' % (i+1, math.sqrt(dmin)))  # node, value
    if (i % 100) == 0: 
      print i,
      sys.stdout.flush()
  f.write('$EndNodeData\n')
  f.close()
  print


