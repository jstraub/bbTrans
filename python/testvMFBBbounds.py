import numpy as np
from scipy.linalg import det, eig, inv, solve
import scipy
from itertools import combinations
from js.geometry.rotations import *
from js.geometry.sphere import Sphere
from discretized4dSphere import S3Grid
from vMFMM import *
import mayavi.mlab as mlab

def near(a, b):
  return np.abs(a-b) < 1e-6

def ComputeF(z):
  if np.abs(z) < 1e-6:
    return 2.
  else:
    return (np.exp(z) - np.exp(-z)) / z

def LowerBound(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  Compute a lowerbound on the objective by evaluating it at the center
  point of the tetrahedron in 4D
  '''
  center = 0.25*vertices[tetra[0],:] + 0.25*vertices[tetra[1],:] + \
    0.25*vertices[tetra[2],:] +  0.25*vertices[tetra[3],:]
  center /= np.sqrt((center**2).sum())
  qs = [Quaternion(vec=center),
      Quaternion(vec=vertices[tetra[0],:]),
      Quaternion(vec=vertices[tetra[1],:]),
      Quaternion(vec=vertices[tetra[2],:]),
      Quaternion(vec=vertices[tetra[3],:])]
  lb = np.zeros(5)
  for i in range(5):
    for j in range(vMFMM_A.GetK()):
      for k in range(vMFMM_B.GetK()):
        lb[i] += ComputevMFtovMFcost(vMFMM_A,
            vMFMM_B, j, k, qs[i].toRot().R.dot(vMFMM_B.GetvMF(k).GetMu()))
        if np.isnan(lb[i]):
          import ipdb
          ipdb.set_trace()
           
  return np.max(lb)

def FindMaximumQAQ(A, vertices, tetra):
  lambdas = []
  Q = np.zeros((4,4))
  for i in range(4):
    Q[:,i] = vertices[tetra[i],:]
  # Full problem:
  A_ = Q.T.dot(A).dot(Q) 
  B_ = Q.T.dot(Q)  
  e, V = eig(A_, B_)
  alpha = np.real(V[:,np.argmax(e)])
  if np.all(alpha >= 0.) or np.all(alpha <= 0.):
    lambdas.append(np.max(np.real(e)))
#  print V
#  print alpha
  # Only three qs: 
  for comb in combinations(range(4), 3):
    A__ = np.array([[A_[i,j] for j in comb] for i in comb])
    B__ = np.array([[B_[i,j] for j in comb] for i in comb])
    e, V = eig(A__, B__)
    alpha = np.real(V[:,np.argmax(e)])
    if np.all(alpha >= 0.) or np.all(alpha <= 0.):
      lambdas.append(np.max(np.real(e)))
  # Only two qs: 
  for comb in combinations(range(4), 2):
    A__ = np.array([[A_[i,j] for j in comb] for i in comb])
    B__ = np.array([[B_[i,j] for j in comb] for i in comb])
    e, V = eig(A__, B__)
    alpha = np.real(V[:,np.argmax(e)])
    if np.all(alpha >= 0.) or np.all(alpha <= 0.):
      lambdas.append(np.max(np.real(e)))
  # Only one q: 
  for i in range(4):
    lambdas.append((Q[:,i]).T.dot(A).dot(Q[:,i]))
#  print lambdas
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

def UpperBoundConvexity(vMFMM_A, vMFMM_B, vertices, tetra):
  ''' 
  TODO
  '''
  qs = [Quaternion(vec=vertices[tetra[i],:]) for i in range(4)]
  A = np.zeros((4,4))
  B = 0.
  for j in range(vMFMM_A.GetK()):
    for k in range(vMFMM_B.GetK()):
      tau_A = vMFMM_A.GetvMF(j).GetTau()
      tau_B = vMFMM_B.GetvMF(k).GetTau()
#      s = Sphere(3)
#      figm = mlab.figure(bgcolor=(1,1,1))
#      s.plotFanzy(figm, 1.)
      mu_U = ClosestMu(vMFMM_A.GetvMF(j).GetMu(),
          vMFMM_B.GetvMF(k).GetMu(), qs, None) #figm)
      mu_L = FurtestMu(vMFMM_A.GetvMF(j).GetMu(),
          vMFMM_B.GetvMF(k).GetMu(), qs, None) # figm)
#      mlab.show(stop=True)
      U = np.sqrt(((vMFMM_A.GetvMF(j).GetTau() *
        vMFMM_A.GetvMF(j).GetMu() + vMFMM_B.GetvMF(k).GetTau() *
        mu_U)**2).sum())
      L = np.sqrt(((vMFMM_A.GetvMF(j).GetTau() *
        vMFMM_A.GetvMF(j).GetMu() + vMFMM_B.GetvMF(k).GetTau() *
        mu_L)**2).sum())
#      print U, L
      fUfLoU2L2 = 0.
      L2fUU2fLoU2L2 = 0.
      if np.abs(U-L) < 1e-6:
        # Assymptotics for U-L -> 0
        fUfLoU2L2 = (1. + U - np.exp(2.*U) + U * np.exp(2.*U))/(2.*U**3*np.exp(U))
        L2fUU2fLoU2L2 = -(3+U-3*np.exp(2.*U) + U*np.exp(2.*U))/(2.*U*np.exp(U))
      else:
        f_U = ComputeF(U)
        f_L = ComputeF(L)
#        print f_U, f_L
        fUfLoU2L2 = ((f_U-f_L)/(U**2 - L**2))
        L2fUU2fLoU2L2 = ((U**2*f_L - L**2*f_U)/(U**2-L**2))
#      M = ToRightQuaternionProductMatrix(vMFMM_A.GetvMF(j).GetMu()).T.dot(\
#          ToLeftQuaternionProductMatrix(vMFMM_B.GetvMF(k).GetMu()))
      M = BuildM(vMFMM_A.GetvMF(j).GetMu(), vMFMM_B.GetvMF(k).GetMu())
      D = 2. * np.pi * vMFMM_A.GetPi(j) * vMFMM_B.GetPi(k) * \
        vMFMM_A.GetvMF(j).GetZ() * vMFMM_A.GetvMF(k).GetZ() 
      A += 2.*tau_A*tau_B*D*fUfLoU2L2 * M
      B += D*(tau_A**2*fUfLoU2L2 + tau_B**2*fUfLoU2L2 + L2fUU2fLoU2L2)
  #print A, B
  lambda_max = FindMaximumQAQ(A, vertices, tetra)
#  print "--", B, lambda_max
  return B + lambda_max, B, lambda_max

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
#    if np.abs(det(A)) > 1e-6:
      # If the triangle of mus is not degenrate (i.e. they ly on a
      # line)
      a = np.linalg.solve(A, -muA)
      if np.all(a > 0.):
        # Closest mu is interior point.
        mu_star = -muA
        if not figm is None:
          for mu in mus:
            mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.1, opacity
                = 0.5)
          mlab.points3d([mu_star[0]],[mu_star[1]],[mu_star[2]],
              scale_factor=0.1, color=(0,1,1), opacity = 1.0,
              mode="2ddiamond")
          mlab.points3d([muA[0]],[muA[1]],[muA[2]],
              scale_factor=0.1, color=(1,0,0), opacity = 0.5, mode="2dcircle")
#          mlab.show(stop=True)
        return mu_star
    except scipy.linalg.LinAlgError:
      print "singular matrix in FurthesstMu"
#      pass
#      if np.sum(a == 1.) == 1 and np.sum(a == 0.) == 2:
#        return A.dot(a)
  furthestLocations = [mu for mu in mus]
  for i, mu1 in enumerate(mus):
    for mu2 in mus[i+1::]:
      furthestLocations.append(ComputeExtremaLocationOnGeodesic(mu1,
        mu2, muA)) 
  dists = np.array([mu.dot(muA) for mu in furthestLocations])
  mu_star = furthestLocations[np.argmin(dists)]
#  print mu_star, dists
  if not figm is None:
    for mu in mus:
      mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.1, opacity
          = 0.5)
    mlab.points3d([mu_star[0]],[mu_star[1]],[mu_star[2]],
        scale_factor=0.1, color=(0,1,1), opacity = 1.0, mode="2ddiamond")
    mlab.points3d([muA[0]],[muA[1]],[muA[2]],
        scale_factor=0.1, color=(1,0,0), opacity = 0.5, mode="2dcircle")
#    mlab.show(stop=True)
  return mu_star

def ClosestMu(muA, muB, qs, figm = None):
  mus = [q.toRot().R.dot(muB) for q in qs]
  A = np.zeros((3,3))
  for tri in combinations(range(4),3):
    A[:,0] = mus[tri[0]]
    A[:,1] = mus[tri[1]]
    A[:,2] = mus[tri[2]]
    try:
#    if np.abs(det(A)) > 1e-6:
      # If the triangle of mus is not degenrate (i.e. they ly on a
      # line)
      a = np.linalg.solve(A, muA)
      if np.all(a > 0.):
        # Closest mu is interior point.
        mu_star = muA
        if not figm is None:
          for mu in mus:
            mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.1, opacity
                = 0.5)
          mlab.points3d([mu_star[0]],[mu_star[1]],[mu_star[2]],
              scale_factor=0.1, color=(0,1,0), opacity = 1., mode="2dcross")
          mlab.points3d([muA[0]],[muA[1]],[muA[2]],
              scale_factor=0.1, color=(1,0,0), opacity = 0.5, mode="2dcircle")
#          mlab.show(stop=True)
#        print a, muA, A, 
#        import ipdb
#        ipdb.set_trace()
        return muA
    except scipy.linalg.LinAlgError:
      print "singular matrix in ClosestMu"
#      pass
#      if np.sum(a == 1.) == 1 and np.sum(a == 0.) == 2:
#        return A.dot(a)
  closestLocations = [mu for mu in mus]
  for i, mu1 in enumerate(mus):
    for mu2 in mus[i+1::]:
      closestLocations.append(ComputeExtremaLocationOnGeodesic(mu1,
        mu2, muA)) 
  dists = np.array([mu.dot(muA) for mu in closestLocations])
  mu_star = closestLocations[np.argmax(dists)]
#  print mu_star, dists
  if not figm is None:
    print "plotting closest"
    for mu in mus:
      mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.1, opacity
          = 0.5)
    mlab.points3d([mu_star[0]],[mu_star[1]],[mu_star[2]],
        scale_factor=0.1, color=(0,1,0), opacity = 1., mode="2dcross")
    mlab.points3d([muA[0]],[muA[1]],[muA[2]],
        scale_factor=0.1, color=(1,0,0), opacity = 0.5, mode="2dcircle")
#    mlab.show(stop=True)
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
  if False:
    a, b = 0.1, 2.
    z = np.linspace(a, b,100)
    fig = plt.figure()
    f = lambda (z): (np.exp(z) - np.exp(-z)) / z
    plt.plot(z, f(z), label="f")
    g = lambda z,a,b: z**2*((ComputeF(b) - ComputeF(a))/(b**2-a**2)) +\
      ((b**2*ComputeF(a)-a**2*ComputeF(b))/(b**2-a**2))
    plt.plot(z, g(z,a,b), label="g")
    plt.legend()
    plt.show()


  if False:
    a, b = 0.1, 0.11
    z = np.linspace(a, b,100)
    fig = plt.figure()
    f = lambda (z): (np.exp(z) - np.exp(-z)) / z
    plt.plot(z, f(z), label="f")
    g = lambda z,a,b: z**2*((ComputeF(b) - ComputeF(a))/(b**2-a**2)) +\
      ((b**2*ComputeF(a)-a**2*ComputeF(b))/(b**2-a**2))
    plt.plot(z, g(z,a,b), label="g")

    h = lambda z,a,b: z**2*((1. + a - np.exp(2.*a) + a \
      * np.exp(2.*a))/(2.*a**3*np.exp(a))) -\
      ((3+a-3*np.exp(2.*a) + a*np.exp(2.*a))/(2.*a*np.exp(a)))
    plt.plot(z, h(z,a,b), label="h")
    plt.legend()
    plt.show()
#        fUfLoU2L2 = (1. + U - np.exp(2.*U) + U * np.exp(2.*U))/(2.*U**3*np.exp(U))
#        L2fUU2fLoU2L2 = (3+U-3*np.exp(2.*U) + U*np.exp(2.*U))/(2.*U*np.exp(U))

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
  ubC = np.zeros((tetras.shape[0], 1))
  ubCB = np.zeros((tetras.shape[0], 1))
  ubCL = np.zeros((tetras.shape[0], 1))
  for i in range(s3.tetra_levels[-1]):
    lb[i] = LowerBound(vMFMM_A, vMFMM_B, s3.vertices, tetras[i,:])
    ub[i] = UpperBound(vMFMM_A, vMFMM_B, s3.vertices, tetras[i,:])
    ubC[i], ubCB[i], ubCL[i] = UpperBoundConvexity(vMFMM_A, vMFMM_B, s3.vertices, tetras[i,:])
    if lb[i] > ubC[i] + 1e-6:
      print i, lb[i], ub[i], ubC[i]
      tetra = tetras[i,:]
      vertices = s3.vertices
      qs = [Quaternion(vec=vertices[tetra[0],:]),
          Quaternion(vec=vertices[tetra[1],:]),
          Quaternion(vec=vertices[tetra[2],:]),
          Quaternion(vec=vertices[tetra[3],:])]
      for j in range(vMFMM_A.GetK()):
        for k in range(vMFMM_B.GetK()):
          figm = mlab.figure(bgcolor=(1,1,1))
          s = Sphere(3)
          s.plotFanzy(figm,1.)
          muA = vMFMM_A.GetvMF(j).GetMu()
          muB = vMFMM_B.GetvMF(k).GetMu()
          print ClosestMu(muA, muB, qs, figm)
          print FurtestMu(muA, muB, qs, figm)
          center = 0.25*vertices[tetra[0],:] + 0.25*vertices[tetra[1],:] + \
            0.25*vertices[tetra[2],:] +  0.25*vertices[tetra[3],:]
          center /= np.sqrt((center**2).sum())
          mus = [q.toRot().R.dot(muB) for q in qs]
          for mu in mus:
            mlab.points3d([mu[0]],[mu[1]],[mu[2]], scale_factor=0.1, opacity
                = 0.5)
#          mlab.points3d([mu_star[0]],[mu_star[1]],[mu_star[2]],
#              scale_factor=0.1, color=(0,1,0), opacity = 0.5, mode="2dcross")
          mlab.points3d([muA[0]],[muA[1]],[muA[2]],
              scale_factor=0.1, color=(1,0,0), opacity = 0.5, mode="2dcircle")
          mlab.show(stop=True)

  print lb.T
  print ub.T
  print ubC.T
  print np.all(ub > lb)
  print np.sum(ub > lb)
  print np.sum(ubC > lb)

  print np.argmax(lb), np.argmax(ub), np.argmax(ubC)
  print np.max(lb), np.max(ub), np.max(ubC)
  
#  vertices = s3.vertices
#  tetraMax = tetras[np.argmax(lb),:] 
#  center = 0.25*vertices[tetraMax[0],:] + 0.25*vertices[tetraMax[1],:] + \
#    0.25*vertices[tetraMax[2],:] +  0.25*vertices[tetraMax[3],:]
#  center /= np.sqrt((center**2).sum())
#  q = Quaternion(vec=center)
#  print "UB", center
#  print q.toRot().R
#
#  tetraMax = tetras[np.argmax(ub),:] 
#  center = 0.25*vertices[tetraMax[0],:] + 0.25*vertices[tetraMax[1],:] + \
#    0.25*vertices[tetraMax[2],:] +  0.25*vertices[tetraMax[3],:]
#  center /= np.sqrt((center**2).sum())
#  q = Quaternion(vec=center)
#  print "LB", center
#  print q.toRot().R

  import matplotlib.pyplot as plt
  plt.figure()
  plt.plot(lb, label = "lb")
  plt.plot(ub, label = "ub")
  plt.plot(ubC, label="ub convex")
#  plt.plot(ubCB, label="B")
#  plt.plot(ubCL, label="lambda_max")
  plt.legend()

  plt.figure()
  plt.hist((ubC-lb)/(ub-lb), 100)
  plt.show()

