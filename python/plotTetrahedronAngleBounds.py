from discretized4dSphere import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#paper
mpl.rc('font',size=30) 
mpl.rc('lines',linewidth=3.)
figSize = (14, 5.5)
figSize = (14, 12)

def ToRad(deg):
  return deg*np.pi/180.
def ToDeg(rad):
  return rad/np.pi*180.

s3 = S3Grid(0)
tetras = s3.GetTetrahedra(0)

lvls = 20
dotMaxPred = [-1*np.ones(1)]*20
dotMaxPred[0] = tetras[0].GeMinMaxVertexDotProduct()[1]
for lvl in range(1,lvls):
  dotMaxPred[lvl] = dotMaxPred[lvl-1]*2./(1.+dotMaxPred[lvl-1])

dotMaxPredSqrt = [-1*np.ones(1)]*20
dotMaxPredSqrt[0] = tetras[0].GeMinMaxVertexDotProduct()[1]
for lvl in range(1,lvls):
  dotMaxPredSqrt[lvl] = np.sqrt((dotMaxPredSqrt[lvl-1]+1.)/2.)

dotMaxPred2 = [-1*np.ones(1)]*20
dotMaxPred2[0] = tetras[0].GeMinMaxVertexDotProduct()[1]
for lvl in range(1,lvls):
  dotMaxPred2[lvl] = (3.*dotMaxPred2[lvl-1]+1.)/(2.*(1. +dotMaxPred2[lvl-1]))

dotMax = [-1.]*20

for i in range(600):
  tetra = tetras[np.random.randint(0,len(tetras),1)]
  dotMax[0] = max(dotMax[0], tetra.GeMinMaxVertexDotProduct()[1])
  for lvl in range(1,lvls):
    tetra = tetra.Subdivide()[np.random.randint(0,6,1)]
    dotMax[lvl] = max(dotMax[lvl], tetra.GeMinMaxVertexDotProduct()[1])

dotMaxPred = [ToDeg(np.arccos(dotMaxPred_i)) for dotMaxPred_i in dotMaxPred]
dotMaxPredSqrt = [ToDeg(np.arccos(dotMaxPredSqrt_i)) for dotMaxPredSqrt_i in dotMaxPredSqrt]
dotMaxPred2 = [ToDeg(np.arccos(dotMaxPred2_i)) for dotMaxPred2_i in dotMaxPred2]
dotMax = [ToDeg(np.arccos(dotMax_i)) for dotMax_i in dotMax]

print dotMaxPred
print dotMax

fig = plt.figure(figsize = figSize, dpi = 80, facecolor="w",
    edgecolor="k")
plt.plot(dotMaxPred, label="upper bound")
plt.plot(dotMax, label="actual")
plt.xlabel("subdivision level of the tetrahedron")
plt.ylabel("min. angle between any two vertices [deg]")
plt.legend()
#plt.tight_layout()
plt.savefig("../subdivisionVsMinAngle_ActualAndBound.png", figure=fig)

fig = plt.figure(figsize = figSize, dpi = 80, facecolor="w",
    edgecolor="k")
plt.plot(dotMaxPredSqrt,'-', label=r"$\sqrt{\frac{1+\gamma}{2}}$")
plt.plot(dotMaxPred2, '--', label=r"$\frac{1+3\gamma}{2(1+\gamma)}$")
plt.plot(dotMaxPred,'-', label=r"$\frac{2\gamma}{1+\gamma}$")
plt.plot(dotMax,"-.", label="actual")
plt.xlabel("subdivision level of the tetrahedron")
plt.ylabel("min. angle between any two vertices [deg]")
plt.legend()
#plt.tight_layout()
plt.savefig("../subdivisionVsMinAngle_ActualAndBound_allBounds.png", figure=fig)
plt.show()
