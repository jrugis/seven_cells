import numpy as np

###########################################################################
# functions
###########################################################################

def getDistQ(A, B, P): # get distance from point P to line segment AB
  AB = B-A
  AP = P-A
  dAP = np.linalg.norm(AP)
  AQ = (np.dot(AB, AP)) * AB / (AB**2).sum() 
  Q = AQ + A                 # QP is perpendicular to line AB
  d2AQ = (AQ**2).sum()
  d2BQ = ((Q-B)**2).sum()
  d2AB = (AB**2).sum()
  if((d2AQ > d2AB) or (d2BQ > d2AB)): # Q not on segment AB?
    d = min(dAP, np.linalg.norm(P-B)) # use distance to A or B
  else:
    d = np.linalg.norm(P-Q)           # use perpendicular distance
  return d

###########################################################################
def get_dnl(tris, verts, lsegs):
  ntris = tris.shape[0]
  dnl = np.empty([ntris])
  for i in range(ntris):
    C = np.average(verts[tris[i]-1], axis=0) # center of tri
    d = 100.0 # large dummy initial distance
    for seg in lsegs:
      ds = getDistQ(seg[0:3], seg[3:6], C)
      if ds < d: 
        d = ds
    dnl[i] = d # save the minimum distance
  return dnl

###########################################################################
def get_apical(tris, dnl, small_radius, large_radius):
  atrisS = []
  atrisL = []
  for i, d in enumerate(dnl):
    if d < small_radius: atrisS.append(i+1)
    if d < large_radius: atrisL.append(i+1)
  return(np.array(atrisS, dtype=int), np.array(atrisL, dtype=int))

###########################################################################
def get_basal(tris, ctrisi, atrisi):
  btrisi = np.array([range(tris.shape[0])], dtype=int) + 1 # all indices
  btrisi = np.setdiff1d(btrisi, ctrisi) # remove common
  btrisi = np.setdiff1d(btrisi, atrisi) # remove apical Large
  return(btrisi)

###########################################################################
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

###########################################################################
def get_common(vertsA, trisA, vertsB, trisB):
  cviA = []
  cviB = []
  for i,v in enumerate(vertsB): # first find common vertices
    ans = np.where(np.linalg.norm(vertsA-v,axis=1) < 0.001)
    if len(ans[0]) != 0:
      if len(ans[0]) > 1:
        print("ERROR: duplicate vertices")
      cviA.append(ans[0][0]+1)
      cviB.append(i+1)
  cviA = np.array(cviA, dtype=int)
  cviB = np.array(cviB, dtype=int)

  tiA = find_tris(cviA,trisA)
  tiB = np.empty((0,1))
  return(tiA, tiB) # return tri indicies
###########################################################################


