import numpy as np
from js.geometry.rotations import *
from js.geometry.sphere import Sphere
from discretized4dSphere import S3Grid
import mayavi.mlab as mlab

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

def ComputevMFtovMFcost(vMFMM_A, vMFMM_B, j, k, nu):
  C = 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
      vMFMM_A.GetvMF(j).GetZ() * vMFMM_A.GetvMF(k).GetZ() 
  z_jk = np.sqrt(((vMFMM_A.GetvMF(j).GetTau() *
    vMFMM_A.GetvMF(j).GetMu() + vMFMM_B.GetvMF(k).GetTau() *
    nu)**2).sum())
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
      lb += ComputevMFtovMFcost(vMFMM_A,
          vMFMM_B, j, k, q.toRot().R.dot(vMFMM_B.GetvMF(k).GetMu()))
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

def near(a, b):
  return np.abs(a-b) < 1e-6

def ComputeClosestLocationOnGeodesic(mu1, mu2, nu):
#  print mu1, mu2, nu
  theta12 = np.arccos(mu1.dot(mu2))
  theta1 = np.arccos(mu1.dot(nu))
  theta2 = np.arccos(mu2.dot(nu))
  if near(theta1, np.pi*0.5) and near(theta2, np.pi*0.5):
    # any t is good.
    t = 0.5
  t = (np.arctan2(np.cos(theta2) - np.cos(theta12)*np.cos(theta1),np.cos(theta1)*np.sin(theta12))) / theta12
  t = min(1.0, max(0.0, t))
  mu_star = (mu1*np.sin((1.-t)*theta12) + mu2*np.sin(t*theta12))/np.sin(theta12)
#  print t, mu_star, nu.dot(mu_star)
  return mu_star

def ClosestMu(muA, muB, qs, figm = None):
  mus = np.array([q.toRot().R.dot(muB) for q in qs]).T
#  print mus.shape, muA.shape
#  print mus.T.dot(mus), mus.T.dot(muA)
#  a = np.linalg.solve(mus.T.dot(mus), mus.T.dot(muA))
#  print a
#  if np.all(a > 0.):
#    return muA
#  else:
  mus = mus.T
  closestLocations = []
#  print mus
  for i, mu1 in enumerate(mus):
    for mu2 in mus[i+1::]:
      closestLocations.append(ComputeClosestLocationOnGeodesic(mu1,
        mu2, muA)) 
  dists = np.array([mu.dot(muA) for mu in closestLocations])
  mu_star = closestLocations[np.argmax(dists)]
#  print mu_star, dists
  if not figm is None:
    for mu in mus:
      mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.1, opacity
          = 0.5)
    mlab.points3d([mu_star[0]],[mu_star[1]],[mu_star[2]],
        scale_factor=0.1, color=(0,1,0), opacity = 0.5, mode="2dcross")
    mlab.points3d([muA[0]],[muA[1]],[muA[2]],
        scale_factor=0.1, color=(1,0,0), opacity = 0.5, mode="2dcircle")
    mlab.show(stop=True)
  return mu_star
   
def UpperBound(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  '''
  qs = [Quaternion(vec=vertices[tetra[i],:]) for i in range(4)]
  ub = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
#      figm = mlab.figure(bgcolor=(1,1,1))
#      s = Sphere(2)
#      s.plotFanzy(figm, 1)
      mu_star = ClosestMu(vMFMM_A.GetvMF(j).GetMu(),
          vMFMM_B.GetvMF(k).GetMu(), qs, None)
      ub += ComputevMFtovMFcost(vMFMM_A, vMFMM_B, j, k, mu_star)
  return ub

if __name__ == "__main__":
  s3 = S3Grid(0)
  print s3.tetra_levels


  vMFs_A = [vMF(np.array([1.,0.,0.]), 1.), vMF(np.array([0.,1.,0.]), 10.)]
  vMFs_B = [vMF(np.array([1.,0.,0.]), 1.), vMF(np.array([0.,0.,1.]), 10.)]
  vMFMM_A = vMFMM(np.array([0.5, 0.5]), vMFs_A)
  vMFMM_B = vMFMM(np.array([0.5, 0.5]), vMFs_B)

  tetras = s3.GetTetras(0)
  print tetras.shape
  lb = np.zeros((tetras.shape[0], 1))
  ub = np.zeros((tetras.shape[0], 1))
  for i in range(s3.tetra_levels[-1]):
    lb[i] = LowerBound(vMFMM_A, vMFMM_B, s3.vertices, tetras[i,:])
    ub[i] = UpperBound(vMFMM_A, vMFMM_B, s3.vertices, tetras[i,:])
  print lb.T
  print ub.T
  print np.all(ub > lb)
  print np.sum(ub > lb)

  print np.argmax(lb), np.argmax(ub)
  print np.max(lb), np.max(ub)
  
  vertices = s3.vertices
  tetraMax = tetras[np.argmax(lb),:] 
  center = 0.25*vertices[tetraMax[0],:] + 0.25*vertices[tetraMax[1],:] + \
    0.25*vertices[tetraMax[2],:] +  0.25*vertices[tetraMax[3],:]
  center /= np.sqrt((center**2).sum())
  q = Quaternion(vec=center)
  print "UB", center
  print q.toRot().R

  tetraMax = tetras[np.argmax(ub),:] 
  center = 0.25*vertices[tetraMax[0],:] + 0.25*vertices[tetraMax[1],:] + \
    0.25*vertices[tetraMax[2],:] +  0.25*vertices[tetraMax[3],:]
  center /= np.sqrt((center**2).sum())
  q = Quaternion(vec=center)
  print "LB", center
  print q.toRot().R

  import matplotlib.pyplot as plt
  plt.figure()
  plt.plot(lb)
  plt.plot(ub)
  plt.show()
