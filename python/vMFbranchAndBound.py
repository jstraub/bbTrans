import numpy as np
from js.geometry.rotations import *
from js.geometry.sphere import Sphere
from discretized4dSphere import S3Grid
from vMFMM import *
import mayavi.mlab as mlab

def near(a, b):
  return np.abs(a-b) < 1e-6

def ToDeg(theta):
  return theta*180./np.pi

def LowerBound(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  Compute a lowerbound on the objective by evaluating it at the center
  point of the tetrahedron in 4D
  '''
  center = 0.25*vertices[tetra[0]] + 0.25*vertices[tetra[1]] + \
    0.25*vertices[tetra[2]] +  0.25*vertices[tetra[3]]
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
  def __init__(self, vMFMM_A, vMFMM_B):
    self.vMFMM_A = vMFMM_A
    self.vMFMM_B = vMFMM_B

  def Branch(self, node):
    tetrahedra = node.tetrahedron.Subdivide()
    return [Node(tetrahedron) for tetrahedron in tetrahedra]

  def UpperBound(self, node):
    return UpperBound(self.vMFMM_A, self.vMFMM_B, node.tetrahedron.vertices, 
          node.tetrahedron.tetra)

  def LowerBound(self, node):
    return LowerBound(self.vMFMM_A, self.vMFMM_B, node.tetrahedron.vertices, 
          node.tetrahedron.tetra)

  def Compute(self, nodes, q_star):
    lb = -1e6
    ub = -1e6
    node_star = None
    counter = 0
    lbs = [self.LowerBound(node) for node in nodes]
    ubs = [self.UpperBound(node) for node in nodes]
    while counter < 50000:
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
#        if node.tetrahedron.lvl > 20:
#          continue
        new_nodes = self.Branch(node)
        for n in new_nodes:
          ubn = self.UpperBound(node)
          if ubn > lb:
            # Upper bound of the node is greater than the global lower
            # bound. So keep the node arround.
            nodes.append(n)
            lbs.append(self.LowerBound(node))
            ubs.append(ubn)

      q_star = Quaternion(vec=node_star.tetrahedron.Center())
      dAng = np.zeros(self.vMFMM_A.GetK())
      for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
        dAngs = np.zeros(self.vMFMM_B.GetK())
        for k, vMF_B in enumerate(self.vMFMM_B.vMFs):
          dAngs[k] = ToDeg(np.arccos(vMF_A.GetMu().dot(q_star.rotate(vMF_B.GetMu()))))
        dAng[j] = np.min(dAngs)
      print lb, ub, counter, len(nodes), node_star.tetrahedron.lvl, \
        dAng
      for j, vMF_A in enumerate(self.vMFMM_A.vMFs):
        print "A", j, vMF_A.GetMu()
      for k, vMF_B in enumerate(self.vMFMM_B.vMFs):
        print "B", k, q_star.rotate(vMF_B.GetMu())
      counter += 1

if __name__ == "__main__":
  s3 = S3Grid(0)
  print s3.tetra_levels
  
  q = Quaternion()
  q.setToRandom()
  R = q.toRot().R
  print "q True: ", q.q, np.sqrt((q.q**2).sum())

  vMFs_A = [vMF(np.array([1.,0.,0.]), 10.), vMF(np.array([0.,1.,0.]), 10.)]
  vMFs_B = [vMF(R.dot(np.array([1.,0.,0.])), 10.),
      vMF(R.dot(np.array([0.,1.,0.])), 10.)]
  vMFMM_A = vMFMM(np.array([0.5, 0.5]), vMFs_A)
  vMFMM_B = vMFMM(np.array([0.5, 0.5]), vMFs_B)

  tetras = s3.GetTetras(0)
  tetrahedra = s3.GetTetrahedra(0)

  bb = BB(vMFMM_A, vMFMM_B)
  nodes = [Node(tetrahedron) for tetrahedron in tetrahedra]
  bb.Compute(nodes, q)

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
