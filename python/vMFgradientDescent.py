import numpy as np
from js.geometry.rotations import *
from js.geometry.sphere import Sphere
from vMFMM import *
import mayavi.mlab as mlab
import matplotlib.pyplot as plt

def ToDeg(theta):
  return theta*180./np.pi
def ComputeF(z):
  if np.abs(z) < 1e-6:
    return 2.
  else:
    return (np.exp(z) - np.exp(-z)) / z
def ComputedF(z):
  return ((np.exp(z)+np.exp(-z))/z) - ((np.exp(z)-np.exp(-z))/z**2)

def skew(a):
  return np.array([
    [0., -a[2], a[1]],
    [a[2], 0., -a[0]],
    [-a[1],a[0], 0.]])

def ComputeGradient(vMFMM_A, vMFMM_B, R):
  J = np.zeros(3)
  C = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      tau_A = vMFMM_A.GetvMF(j).GetTau()
      tau_B = vMFMM_B.GetvMF(k).GetTau()
      mu_A = vMFMM_A.GetvMF(j).GetMu()
      mu_B = vMFMM_B.GetvMF(k).GetMu()
      z = np.sqrt(((tau_A * mu_A + tau_B * R.dot(mu_B))**2).sum())
      D = 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
        vMFMM_A.GetvMF(j).GetZ() * vMFMM_B.GetvMF(k).GetZ() 
      mu_B_x = skew(R.dot(mu_B))
      J += -mu_A.T.dot(mu_B_x) * D*(-(tau_A * tau_B)/(np.sqrt(z))) * ComputedF(z)
      C += D * ComputeF(z)
  return J, C

def ComputeCostFunction(vMFMM_A, vMFMM_B, R):
  C = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      tau_A = vMFMM_A.GetvMF(j).GetTau()
      tau_B = vMFMM_B.GetvMF(k).GetTau()
      mu_A = vMFMM_A.GetvMF(j).GetMu()
      mu_B = vMFMM_B.GetvMF(k).GetMu()
      z = np.sqrt(((tau_A * mu_A + tau_B * R.dot(mu_B))**2).sum())
      D = 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
        vMFMM_A.GetvMF(j).GetZ() * vMFMM_B.GetvMF(k).GetZ() 
      C += D * ComputeF(z)
  return C

def ComputeArmijoStepSize(vMFMM_A, vMFMM_B, R0, J):
  beta = 0.5
  alpha = 100000.0 / beta
  sigma = 0.00001
  fx = ComputeCostFunction(vMFMM_A, vMFMM_B, R0)
  fxPlus = fx
  dR = Rot3(np.eye(3))
  normJ = np.sqrt((J**2).sum())
  while fxPlus - fx <= sigma*alpha*normJ and alpha>1e-12:
#  while fx - fxPlus <= sigma*alpha*np.sqrt((J**2).sum()) and alpha>1e-12:
    alpha *= beta 
    dR.expMap(J[:, np.newaxis]*alpha/normJ)
#    dR.expMap(-J[:, np.newaxis]*alpha)
    R = dR.R.dot(R0)
    fxPlus = ComputeCostFunction(vMFMM_A, vMFMM_B, R)
#    print "    ", -fx+fxPlus, sigma*alpha*normJ
#    print R0
#    print R
  return alpha

class GradientDescent:
  def __init__(self, vMFMM_A, vMFMM_B):
    self.vMFMM_A = vMFMM_A
    self.vMFMM_B = vMFMM_B
  def Compute(self, R0, maxIt, R_gt, figm=None ):
    R = np.copy(R0)
    counter = 0
    eps = np.zeros((maxIt, self.vMFMM_A.GetK()+2))
    J = np.ones(3)
    while counter < maxIt: # and J.sum() > 1e-6:
      J, C = ComputeGradient(self.vMFMM_A, self.vMFMM_B, R)
      J *= -1. # Gradient ascent!
      alpha = ComputeArmijoStepSize(self.vMFMM_A, self.vMFMM_B, R, J)
      dR = Rot3(np.eye(3))
#      alpha = 1.
      dR.expMap(J[:, np.newaxis]*alpha/np.sqrt((J**2).sum()))
      #########
      dAng = np.zeros(self.vMFMM_A.GetK())
      for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
        dAngs = np.zeros(self.vMFMM_B.GetK())
        for k, vMF_B in enumerate(self.vMFMM_B.vMFs):
          dAngs[k] = ToDeg(np.arccos(vMF_A.GetMu().dot(R.dot(vMF_B.GetMu()))))
        dAng[j] = np.min(dAngs)
      eps[counter,:self.vMFMM_A.GetK()] = dAng 
      eps[counter,self.vMFMM_A.GetK()] = C
#      eps[counter,3] = ToDeg(Rot3(R.T.dot(R_gt)).toQuat().toAxisAngle()[0])
#      eps[counter,4] = ToDeg(Rot3(R).toQuat().angleTo(Rot3(R_gt).toQuat()))
      eps[counter,self.vMFMM_A.GetK()+1] = ToDeg(np.sqrt(((Rot3(R.dot(R_gt)).logMap())**2).sum()))
      print counter, C, J, alpha
      for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
        print "A", j, vMF_A.GetMu()
      for k, vMF_B in enumerate(self.vMFMM_B.vMFs):
        print "B", k, R.dot(vMF_B.GetMu())
      counter += 1
      if not figm is None and counter == maxIt-1:
        for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
          mu = vMF_A.GetMu()
          mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.05, opacity
              = 0.5, color=(0,1,0))
        for j, vMF_B in enumerate(self.vMFMM_B.vMFs):
          mu = vMF_B.GetMu()
          mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.05, opacity
              = 0.5)
        for j, vMF_B in enumerate(self.vMFMM_B.vMFs):
          mu = R.dot(vMF_B.GetMu())
          mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.05, opacity
              = 0.5, color=(1,0,0))
        mlab.show(stop=True)
      R = dR.R.dot(R)
    return R, eps


if __name__ == "__main__":
  q = Quaternion()
  q.setToRandom()
  R_gt = q.toRot().R

  print "q True: ", q.q, np.sqrt((q.q**2).sum())
  print R_gt

  vMFs_A = [vMF(np.array([1.,0.,0.]), 1.), 
      vMF(np.array([0.,1.,0.]), 10.),
      vMF(np.array([0.,0.,1.]), 100.)]
  vMFs_B = [vMF(R_gt.dot(np.array([1.,0.,0.])), 1.),
      vMF(R_gt.dot(np.array([0.,1.,0.])), 10.), 
      vMF(R_gt.dot(np.array([0.,0.,1.])), 100.)]
  vMFMM_A = vMFMM(np.array([0.3, 0.3, 0.4]), vMFs_A)
  vMFMM_B = vMFMM(np.array([0.3, 0.3, 0.4]), vMFs_B)

  gd = GradientDescent(vMFMM_A, vMFMM_B)

  fig = plt.figure()
  plt.ylabel("Sum over angular deviation between closest vMF means.")
  for i in range(1):
    q.setToRandom()
    R0 = q.toRot().R
    R0 = np.copy(R_gt)
#    R0 = np.eye(3)
    figm = mlab.figure(bgcolor=(1,1,1))
    s = Sphere(3)
    s.plotFanzy(figm, 1.0)
    Rt, eps = gd.Compute(R0, 1000, R_gt, figm)
    plt.subplot(3,1,1)
    plt.plot(eps[:,:3], label="UpperBound")
    plt.subplot(3,1,2)
    plt.plot(eps[:,3], label="UpperBound")
    plt.subplot(3,1,3)
    plt.plot(eps[:,4], label="UpperBound")
  plt.legend()

  print R_gt
  print Rt
#  plt.savefig("./vMFMM_GradientComparison_residuals.png")
#  fig = plt.figure()
#  plt.ylabel("Angle between GT and inferred rotation.")
#  plt.plot(eps[:,2], label="UpperBound")
#  plt.legend()
#  plt.savefig("./vMFMM_GradientComparison_angleToGT.png")
  plt.show()

