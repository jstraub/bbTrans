import numpy as np
from js.geometry.rotations import *
from discretized4dSphere import S3Grid

class vMF(object):
  def __init__(self, mu, tau):
    self.mu = np.copy(mu)
    self.tau = tau
    self.Z = self.ComputePartitionFunction()
  def GetTau(self):
    return self.tau
  def GetMu(self):
    return self.mu
  def GetZ(self):
    return self.Z
  def ComputePartitionFunction(self):
    return self.tau / (2.*np.pi*(np.exp(self.tau) - np.exp(-self.tau)))

class vMFMM(object):
  def __init__(self, pis, vMFs):
    self.pis = pis / pis.sum()
    self.vMFs = vMFs
  def GetvMF(self, k):
    return  self.vMFs[k]
  def GetPi(self, k):
    return  self.pis[k]
  def GetK(self):
    return self.pis.size

def ComputevMFtovMFcost(R, vMFMM_A, vMFMM_B, j, k):
  C = 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
      vMFMM_A.GetvMF(j).GetZ() * vMFMM_A.GetvMF(k).GetZ() 
  z_jk = np.sqrt(((vMFMM_A.GetvMF(j).GetTau() *
    vMFMM_A.GetvMF(j).GetMu() + vMFMM_A.GetvMF(k).GetTau() *
    R.dot(vMFMM_A.GetvMF(k).GetMu()))**2).sum())
  C *= (np.exp(z_jk) - np.exp(-z_jk)) / z_jk
  return C

def LowerBound(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  Compute a lowerbound on the objective by evaluating it at the center
  point of the tetrahedron in 4D
  '''
  center = 0.25*vertices[tetra[0],:] + 0.25*vertices[tetra[1],:] + \
    0.25*vertices[tetra[2],:] +  0.25*vertices[tetra[3],:]
  center /= np.sqrt((center**2).sum())
  q = Quaternion(vec=center)
  lb = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      lb += ComputevMFtovMFcost(q.toRot().R, vMFMM_A,
          vMFMM_B, j, k)
  return lb

def UpperBoundConvexity(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  TODO
  '''
  ub = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      ub += 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
          vMFMM_A.GetvMF(j).GetZ() * vMFMM_A.GetvMF(k).GetZ() 
  return ub


def ClosestMu(muA, muB, qs):
  mus = np.array([q.rotate(muB) for q in qs])
  a = np.linalg.solve(mus, muA)
  if np.all(a > 0.):
    return muA
  else:
    #TODO
   

def UpperBound(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  '''
  qs = [Quaternion(vec=vertices[tetra[i],:]) for i in range(4)]
  ub = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      mu_star = ClosestMu(vMFMM_A.GetvMF(j).GetMu(), vMFMM_B.GetvMF(k).GetMu(), qs)
  return ub

if __name__ == "__main__":
  s3 = S3Grid(0)
  print s3.tetra_levels


  vMFs_A = [vMF(np.array([1.,0.,0.]), 1.)]
  vMFs_B = [vMF(np.array([1.,0.,0.]), 1.)]
  vMFMM_A = vMFMM(np.array([1.]), vMFs_A)
  vMFMM_B = vMFMM(np.array([1.]), vMFs_B)

  tetras = s3.GetTetras(0)
  print tetras.shape
  lb = np.zeros((tetras.shape[0], 1))
  for i in range(s3.tetra_levels[-1]):
    lb[i] = LowerBound(vMFMM_A, vMFMM_B, s3.vertices, tetras[i,:])
  print lb.T

  print np.argmax(lb)
  print np.max(lb)
  
  vertices = s3.vertices
  tetraMax = tetras[np.argmax(lb),:] 
  center = 0.25*vertices[tetraMax[0],:] + 0.25*vertices[tetraMax[1],:] + \
    0.25*vertices[tetraMax[2],:] +  0.25*vertices[tetraMax[3],:]
  center /= np.sqrt((center**2).sum())
  q = Quaternion(vec=center)
  print center
  print q.toRot().R
