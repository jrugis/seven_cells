import numpy as np
import sys

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
def get_dfb(btris, tris, verts, tets):
  nbtris = btris.shape[0]
  ntets = tets.shape[0]
  ctris = np.empty([nbtris,3]) # center of basal tris
  dfb = np.empty([ntets]) # distance from basal to center of tet
  for i in range(nbtris):
    ctris[i] = np.average(verts[tris[btris[i]-1]-1], axis=0) # center of basal tri
  for i in range(ntets):
    if i%1000 == 0: print(".", end='.'); sys.stdout.flush() # indication of progress 
    T = np.average(verts[tets[i]-1], axis=0) # center of tet
    d = 100.0 # initial large dummy distance
    for B in ctris:
      ds = np.linalg.norm(T-B)
      if ds < d: 
        d = ds
    dfb[i] = d # save the minimum distance
  return dfb

###########################################################################
#def get_dfa(atris, tris, verts):
#  nverts = verts.shape[0]
#  natris = atris.shape[0]
#  catris = np.empty([natris,3]) # center of apical tris
#  dfa = np.empty([nverts]) # distance from apical (node-wise)
#  for i in range(natris):
#    catris[i] = np.average(verts[tris[atris[i]-1]-1], axis=0) # center of apical tri
#  for i in range(nverts):
#   if i%100 == 0: print".",; sys.stdout.flush() # indication of progress 
#   d = 100.0 # initial large dummy distance
#    for A in catris:
#      ds = np.linalg.norm(verts[i]-A)
#      if ds < d: 
#        d = ds
#    dfa[i] = d # save the minimum distance
#  return dfa

def get_dfa(atris, tris, verts, tets):
  nverts = verts.shape[0]
  natris = atris.shape[0]
  ntets = tets.shape[0]
  catris = np.empty([natris,3]) # center of apical tris
  dfa_vert = np.empty([nverts]) # distance from apical to vertices
  dfa_tet = np.empty([ntets]) # distance from apical to center of tet
  for i in range(natris):
    catris[i] = np.average(verts[tris[atris[i]-1]-1], axis=0) # center of apical tri
  for i in range(nverts):
   if i%1000 == 0: print(".", end=''); sys.stdout.flush() # indication of progress 
   d = 100.0 # initial large dummy distance
   for A in catris:
     ds = np.linalg.norm(verts[i]-A)
     if ds < d: 
       d = ds
   dfa_vert[i] = d # save the minimum distance
  for i in range(ntets):
    dfa_tet[i] = np.average(dfa_vert[tets[i]-1]) # average dfa of verts
  return dfa_vert, dfa_tet

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
def find_tri(vi, tris): # vertex indices, tris
  for i,t in enumerate(tris): # find tri containing all three vertex indices
    if len(np.where(t == vi[0])[0]) == 0: continue # first vertex in tri?
    if len(np.where(t == vi[1])[0]) == 0: continue # second vertex in tri?
    if len(np.where(t == vi[2])[0]) == 0: continue # third vertex in tri?
    return i+1 # found
  return 0     # not found

###########################################################################
def find_tris(viA, trisA, viB, trisB): # vertex indices, tris
  tiA = [] # matching tri pairs in cell A & B
  tiB = [] #
  v = np.array([0,0,0], dtype=int)
  for i, t in enumerate(trisA): # find tris associated with verticies in list
    ans = np.where(viA == t[0]) # first vertex in list?
    if len(ans[0]) == 0: continue
    v[0] = ans[0][0]
    ans = np.where(viA == t[1]) # second vertex in list?
    if len(ans[0]) == 0: continue
    v[1] = ans[0][0]
    ans = np.where(viA == t[2]) # third vertex in list?
    if len(ans[0]) == 0: continue
    v[2] = ans[0][0]
    t = find_tri(viB[v], trisB) # find matching tri
    if t != 0: # found?
      tiA.append(i+1)
      tiB.append(t) # t is already one indexed
  return(np.array(tiA, dtype=int),np.array(tiB, dtype=int)) # return tri indicies

###########################################################################
def get_common(vertsA, trisA, vertsB, trisB):
  cviA = [] # matching vertex pairs in cell A & B
  cviB = [] #
  for i,v in enumerate(vertsB): # first find common vertices
    ans = np.where(np.linalg.norm(vertsA-v,axis=1) < 0.001)
    if len(ans[0]) != 0:
      if len(ans[0]) > 1:
        print("ERROR: duplicate vertices")
      cviA.append(ans[0][0]+1)
      cviB.append(i+1)
  cviA = np.array(cviA, dtype=int)
  cviB = np.array(cviB, dtype=int)
  tiA, tiB = find_tris(cviA,trisA,cviB,trisB) # find associated tris
#  for n in range(len(tiA)): # check that triangles actually match
#    t1 = np.average(vertsA[trisA[tiA[n]-1]-1],axis=0)
#    t2 = np.average(vertsB[trisB[tiB[n]-1]-1],axis=0)
#    d = np.linalg.norm(t1-t2)
#    if(d > 0.001):
#      print(tiB[n])
#      print(vertsA[trisA[tiA[n]-1]-1])
#      print(t1)
#      print(vertsB[trisB[tiB[n]-1]-1])
#      print(t2)
#      print(d)
  return(tiA, tiB) # return tri indicies
###########################################################################


