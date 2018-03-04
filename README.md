Cell: binary data file format (4sim_*.bin)
  NOTE: all indicies start from 1

number of verticies (int32)
  vertices x,y,z (3x float32)
  .
  .
  .
number of surface triangles (int32)
  triangle vertex indicies v1,v2,v3 (3x int32), 
    distance to nearest lumen (1x float32)
  .
  .
  .
number of element tetrahedrons (int32)
  tetrahedron vertex indicies t1,t2,t3,t4 (4x int32)
  .
  .
  .
number of apical surface triangles (int32)
  apical surface triangle index (1x int32)
  .
  .
  .
number of basal surface triangles (int32)
  basal surface triangle index (1x int32)
  .
  .
  .
Number of cell-to-cell common surface triangles (int32)
  this cell common surface triangle index (int32), 
    other cell index (int32), 
    other cell common surface triangle index (int32)
  .
  .
  .

