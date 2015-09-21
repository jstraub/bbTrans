import numpy as np
from scipy.linalg import det, eig, inv, solve
import scipy
from itertools import combinations
from js.geometry.rotations import *
from js.geometry.sphere import Sphere
from discretized4dSphere import S3Grid
from vMFMM import *
import mayavi.mlab as mlab
import matplotlib.pyplot as plt

def near(a, b):
  return np.abs(a-b) < 1e-6

def ToDeg(theta):
  return theta*180./np.pi

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

def ComputeF(z):
  if np.abs(z) < 1e-6:
    return 2.
  else:
    return 2.*np.sinh(z) / z
#    return (np.exp(z) - np.exp(-z)) / z

def FindMaximumQAQ(A, vertices, tetra):
  lambdas = []
  Q = np.zeros((4,4))
  for i in range(4):
    Q[:,i] = vertices[tetra[i]]
  # Only one q: 
  for i in range(4):
    lambdas.append((Q[:,i]).T.dot(A).dot(Q[:,i]))
  # Full problem:
  A_ = Q.T.dot(A).dot(Q) 
  B_ = Q.T.dot(Q)  
  try: 
    e, V = eig(A_, B_)
    alpha = np.real(V[:,np.argmax(e)])
    if np.all(alpha >= 0.) or np.all(alpha <= 0.):
      lambdas.append(np.max(np.real(e)))
  except ValueError:
    return np.max(np.array(lambdas))
  # Only three qs: 
  for comb in combinations(range(4), 3):
    A__ = np.array([[A_[i,j] for j in comb] for i in comb])
    B__ = np.array([[B_[i,j] for j in comb] for i in comb])
    try: 
      e, V = eig(A__, B__)
      alpha = np.real(V[:,np.argmax(e)])
      if np.all(alpha >= 0.) or np.all(alpha <= 0.):
        lambdas.append(np.max(np.real(e)))
    except ValueError:
      pass
  # Only two qs: 
  for comb in combinations(range(4), 2):
    A__ = np.array([[A_[i,j] for j in comb] for i in comb])
    B__ = np.array([[B_[i,j] for j in comb] for i in comb])
    try: 
      e, V = eig(A__, B__)
      alpha = np.real(V[:,np.argmax(e)])
      if np.all(alpha >= 0.) or np.all(alpha <= 0.):
        lambdas.append(np.max(np.real(e)))
    except ValueError:
      pass
  return np.max(np.array(lambdas))

def BuildM(u,v):
  ui, uj, uk = u[0], u[1], u[2]
  vi, vj, vk = v[0], v[1], v[2]
  M = np.array([
    [u.dot(v),    uk*vj-uj*vk,       ui*vk-uk*vi,       uj*vi-ui*vj],
    [uk*vj-uj*vk, ui*vi-uj*vj-uk*vk, uj*vi+ui*vj,       ui*vk+uk*vi],
    [ui*vk-uk*vi, uj*vi+ui*vj,       uj*vj-ui*vi-uk*vk, uj*vk+uk*vj],
    [uj*vi-ui*vj, ui*vk+uk*vi,       uj*vk+uk*vj,       uk*vk-ui*vi-uj*vj]])
  return M

def LowerBound(vMFMM_A, vMFMM_B, vertices, tetra, returnBestRotation = False):
  ''' 
  Compute a lowerbound on the objective by evaluating it at the center
  point of the tetrahedron in 4D
  '''
  center = 0.25*vertices[tetra[0]] + 0.25*vertices[tetra[1]] + \
    0.25*vertices[tetra[2]] +  0.25*vertices[tetra[3]]
  center /= np.sqrt((center**2).sum())
  qs = [Quaternion(vec=center),
      Quaternion(vec=vertices[tetra[0]]),
      Quaternion(vec=vertices[tetra[1]]),
      Quaternion(vec=vertices[tetra[2]]),
      Quaternion(vec=vertices[tetra[3]])]
  lb = np.zeros(5)
  for i in range(5):
    for j in range(vMFMM_A.GetK()):
      for k in range(vMFMM_B.GetK()):
        lb[i] += ComputevMFtovMFcost(vMFMM_A,
            vMFMM_B, j, k, qs[i].toRot().R.dot(vMFMM_B.GetvMF(k).GetMu()))
  if not returnBestRotation:
    return np.max(lb)
  else:
    return np.max(lb), qs[np.argmax(lb)]

def UpperBoundConvexity(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  '''
  qs = [Quaternion(vec=vertices[tetra[i]]) for i in range(4)]
  A = np.zeros((4,4))
  B = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      tau_A = vMFMM_A.GetvMF(j).GetTau()
      tau_B = vMFMM_B.GetvMF(k).GetTau()
      mu_U = ClosestMu(vMFMM_A.GetvMF(j).GetMu(),
          vMFMM_B.GetvMF(k).GetMu(), qs, None) #figm)
      mu_L = FurtestMu(vMFMM_A.GetvMF(j).GetMu(),
          vMFMM_B.GetvMF(k).GetMu(), qs, None) # figm)
      U = np.sqrt(((tau_A * vMFMM_A.GetvMF(j).GetMu() + tau_B *
        mu_U)**2).sum())
      L = np.sqrt(((vMFMM_A.GetvMF(j).GetTau() *
        vMFMM_A.GetvMF(j).GetMu() + tau_B * mu_L)**2).sum())
      fUfLoU2L2 = 0.
      L2fUU2fLoU2L2 = 0.
      if np.abs(U-L) < 1e-6:
        # Assymptotics for U-L -> 0
        fUfLoU2L2 = (1. + U - np.exp(2.*U) + U * np.exp(2.*U))/(2.*U**3*np.exp(U))
        L2fUU2fLoU2L2 = -(3+U-3*np.exp(2.*U) + U*np.exp(2.*U))/(2.*U*np.exp(U))
      else:
        f_U = ComputeF(U)
        f_L = ComputeF(L)
        fUfLoU2L2 = ((f_U-f_L)/(U**2 - L**2))
        L2fUU2fLoU2L2 = ((U**2*f_L - L**2*f_U)/(U**2-L**2))
      M = BuildM(vMFMM_A.GetvMF(j).GetMu(), vMFMM_B.GetvMF(k).GetMu())
      D = 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
        vMFMM_A.GetvMF(j).GetZ() * vMFMM_B.GetvMF(k).GetZ() 
      A += 2.*tau_A*tau_B*D*fUfLoU2L2 * M
      B += D*(tau_A**2*fUfLoU2L2 + tau_B**2*fUfLoU2L2 + L2fUU2fLoU2L2)
  lambda_max = FindMaximumQAQ(A, vertices, tetra)
  return B + lambda_max #, B, lambda_max

def ComputeExtremaLocationOnGeodesic(mu1, mu2, nu):
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

def FurtestMu(muA, muB, qs, figm = None):
  mus = [q.toRot().R.dot(muB) for q in qs]
  A = np.zeros((3,3))
  for tri in combinations(range(4),3):
    A[:,0] = mus[tri[0]]
    A[:,1] = mus[tri[1]]
    A[:,2] = mus[tri[2]]
    try:
      # If the triangle of mus is not degenrate (i.e. they ly on a
      # line)
      a = np.linalg.solve(A, -muA)
      if np.all(a > 0.):
        # Closest mu is interior point.
        mu_star = -muA
        return mu_star
    except scipy.linalg.LinAlgError:
      pass
  furthestLocations = [mu for mu in mus]
  for i, mu1 in enumerate(mus):
    for mu2 in mus[i+1::]:
      furthestLocations.append(ComputeExtremaLocationOnGeodesic(mu1,
        mu2, muA)) 
  dists = np.array([mu.dot(muA) for mu in furthestLocations])
  mu_star = furthestLocations[np.argmin(dists)]
  return mu_star

def ClosestMu(muA, muB, qs, figm = None):
  mus = [q.toRot().R.dot(muB) for q in qs]
  A = np.zeros((3,3))
  for tri in combinations(range(4),3):
    A[:,0] = mus[tri[0]]
    A[:,1] = mus[tri[1]]
    A[:,2] = mus[tri[2]]
    try:
      # If the triangle of mus is not degenrate (i.e. they ly on a
      # line)
      a = np.linalg.solve(A, muA)
      if np.all(a > 0.):
        # Closest mu is interior point.
        mu_star = muA
        return muA
    except scipy.linalg.LinAlgError:
      pass
  closestLocations = [mu for mu in mus]
  for i, mu1 in enumerate(mus):
    for mu2 in mus[i+1::]:
      closestLocations.append(ComputeExtremaLocationOnGeodesic(mu1,
        mu2, muA)) 
  dists = np.array([mu.dot(muA) for mu in closestLocations])
  mu_star = closestLocations[np.argmax(dists)]
  return mu_star
   
def UpperBound(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  '''
  qs = [Quaternion(vec=vertices[tetra[i]]) for i in range(4)]
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

class Node(object):
  def __init__(self, tetrahedron):
    self.tetrahedron = tetrahedron

class BB:
  def __init__(self, vMFMM_A, vMFMM_B, lowerBound, upperBound):
    self.vMFMM_A = vMFMM_A
    self.vMFMM_B = vMFMM_B
    self.upperBound_ = upperBound
    self.lowerBound_ = lowerBound

  def Branch(self, node):
    tetrahedra = node.tetrahedron.Subdivide()
    return [Node(tetrahedron) for tetrahedron in tetrahedra]

  def UpperBound(self, node):
    return self.upperBound_(self.vMFMM_A, self.vMFMM_B,
        node.tetrahedron.vertices, node.tetrahedron.tetra)

  def LowerBound(self, node, returnBestRotation=False):
    return self.lowerBound_(self.vMFMM_A, self.vMFMM_B,
        node.tetrahedron.vertices, node.tetrahedron.tetra, 
        returnBestRotation)

  def Compute(self, nodes, maxIt, q_gt):
    lb = -1e6
    ub = 1e6
    node_star = None
    counter = 0
    lbs = [self.LowerBound(node) for node in nodes]
    ubs = [self.UpperBound(node) for node in nodes]
    eps = np.zeros((maxIt, 4))
    while counter < maxIt and ub-lb > 1e-6:
      lbs = [lbs[i] for i, ubn in enumerate(ubs) if ubn > lb]
      nodes = [nodes[i] for i, ubn in enumerate(ubs) if ubn > lb]
      ubs = [ubn for ubn in ubs if ubn > lb]

      i = len(nodes) - 1
      i = np.argmax(np.array(ubs))
      #i = np.argmax(np.array(lbs))
      #i = 0
      node = nodes.pop(i)
      lbn = lbs.pop(i)
      ubn = ubs.pop(i)
      if lbn > lb:
        lb = lbn
        ub = ubn
        node_star = node

        nodes.append(node_star)
        lbs.append(lb)
        ubs.append(ub)
      else:
        new_nodes = self.Branch(node)
        for n in new_nodes:
          ubn = self.UpperBound(node)
          if ubn > lb:
            # Upper bound of the node is greater than the global lower
            # bound. So keep the node arround.
            nodes.append(n)
            lbs.append(self.LowerBound(node))
            ubs.append(ubn)

      q_star = self.LowerBound(node_star, True)[1]
      dAng = np.zeros(self.vMFMM_A.GetK())
      for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
        dAngs = np.zeros(self.vMFMM_B.GetK())
        for k, vMF_B in enumerate(self.vMFMM_B.vMFs):
          dAngs[k] = ToDeg(np.arccos(vMF_A.GetMu().dot(q_star.rotate(vMF_B.GetMu()))))
        dAng[j] = np.min(dAngs)
      print lb, ub, counter, len(nodes), node_star.tetrahedron.lvl, \
        dAng, ToDeg(q_gt.angleTo(q_star))
#      for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
#        print "A", j, vMF_A.GetMu()
#      for k, vMF_B in enumerate(self.vMFMM_B.vMFs):
#        print "B", k, q_star.rotate(vMF_B.GetMu())
      eps[counter,:2] = dAng 
      eps[counter,2] = ComputeCostFunction(self.vMFMM_A, self.vMFMM_A, q_star.toRot().R)
      eps[counter,3] = ToDeg(q_gt.angleTo(q_star))
      counter += 1
    return eps

if __name__ == "__main__":
  s3 = S3Grid(0)
  print s3.tetra_levels
  
  q = Quaternion()
  q.setToRandom()
  R = q.toRot().R
  print "q True: ", q.q, np.sqrt((q.q**2).sum())

  vMFs_A = [vMF(np.array([1.,0.,0.]), 1.), vMF(np.array([0.,1.,0.]), 10.)]
  vMFs_B = [vMF(R.dot(np.array([1.,0.,0.])), 1.),
      vMF(R.dot(np.array([0.,1.,0.])), 10.)]
  vMFMM_A = vMFMM(np.array([0.5, 0.5]), vMFs_A)
  vMFMM_B = vMFMM(np.array([0.5, 0.5]), vMFs_B)

  tetras = s3.GetTetras(0)
  tetrahedra = s3.GetTetrahedra(0)

  nodes = [Node(tetrahedron) for tetrahedron in tetrahedra]
  bb = BB(vMFMM_A, vMFMM_B, LowerBound, UpperBound)
  eps = bb.Compute(nodes, 2000, q)

  nodes = [Node(tetrahedron) for tetrahedron in tetrahedra]
  bb = BB(vMFMM_A, vMFMM_B, LowerBound, UpperBoundConvexity)
  epsC = bb.Compute(nodes, 2000, q)

  fig = plt.figure()
  plt.ylabel("Sum over angular deviation between closest vMF means.")
  plt.plot(eps[:,:2].sum(axis=1), label="UpperBound")
  plt.plot(epsC[:,:2].sum(axis=1), label="UpperBoundConvexity")
  plt.legend()
  plt.savefig("./vMFMM_BoundComparison_residuals.png")
  fig = plt.figure()
  plt.ylabel("Angle between GT and inferred rotation.")
  plt.plot(eps[:,2], label="UpperBound")
  plt.plot(epsC[:,2], label="UpperBoundConvexity")
  plt.legend()
  plt.savefig("./vMFMM_BoundComparison_angleToGT.png")
  plt.show()

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
