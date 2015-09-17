import numpy as np

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
  if np.abs(z_jk) < 1e-6:
    C *= 2.
  else:
    C *= (np.exp(z_jk) - np.exp(-z_jk)) / z_jk
#  C *= (np.exp(z_jk) - np.exp(-z_jk)) / z_jk
  return C
