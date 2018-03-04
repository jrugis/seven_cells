import numpy as np

cti0 = lctrisi[0][:,110]
triAi = cti0[0]-1
celli = cti0[1]-1
triBi = cti0[2]-1
triA = ltris[0][triAi]
triB = ltris[celli][triBi]
vertsA = lverts[0][triA-1]
vertsB = lverts[celli][triB-1]
print vertsA
print vertsB


tris = np.array([[5,7,6],[2,4,3]], dtype=int)
tri = np.array([5,7,6], dtype=int)
ans = np.in1d(tris, [7,6,5])
print ans
#print np.all(ans)


ctrisi = np.empty((3,0), dtype=int) # common tri indices per cell
ctiA,ctiB = ut.get_common(lverts[0], ltris[0], lverts[1], ltris[1]) # pair-wise common
T = np.concatenate([[ctiA],np.full([1,ctiA.shape[0]],c2+1,dtype=int),[ctiB]])
ctrisi = np.concatenate((ctrisi,T), axis=1)
print ctrisi

cviA = []
cviB = []
for i,v in enumerate(lverts[1]): # first find common vertices
  ans = np.where(np.linalg.norm(lverts[0]-v,axis=1) < 0.001)
  if len(ans[0]) != 0:
    if len(ans[0]) > 1:
      print("ERROR: duplicate vertices")
    cviA.append(ans[0][0]+1)
    cviB.append(i+1)
cviA = np.array(cviA, dtype=int)
cviB = np.array(cviB, dtype=int)

tiA = []
tiB = []
v = np.array([0,0,0], dtype=int)
for i, t in enumerate(ltris[0]): # find tris associated with verticies in list
  ans = np.where(cviA == t[0]) # first vertex in list?
  if len(ans[0]) == 0: continue
  v[0] = ans[0][0]
  ans = np.where(cviA == t[1]) # second vertex in list?
  if len(ans[0]) == 0: continue
  v[1] = ans[0][0]
  ans = np.where(cviA == t[2]) # third vertex in list?
  if len(ans[0]) == 0: continue
  v[2] = ans[0][0]
  tiA.append(i+1)
  break
  tiB.append(ut.find_tri(cviB[v], ltris[1]))

print lverts[0][cviA[0]-1]
print lverts[1][cviB[0]-1]
print lverts[0][cviA[100]-1]
print lverts[1][cviB[100]-1]

print tiA[0]
print ltris[0][tiA[0]-1]
print ltris[0][tiA[0]-1][0]-1
print lverts[0][ltris[0][tiA[0]-1][0]-1]

print lverts[0][ltris[0][tiA[0]][0]-1]
print lverts[1][ltris[1][tiB[0]][0]-1]

























