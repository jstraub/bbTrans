from vMFgradientDescent import *
from vMFbranchAndBound import *

def LoadvMFMM(path):
  param = np.loadtxt(path)
  vMFs = []
  pi = param[4,:] / param[4,:].sum()
  print pi
  for k in range(param.shape[1]):
    vMFs.append(vMF(param[:3,k], param[3,k]))
    print param[:3,k]
    print param[3,k]
  return vMFMM(pi, vMFs)

if __name__ == "__main__":
  s3 = S3Grid(0)
  q = Quaternion()
  q.setToRandom()
  R_gt = q.toRot().R
  print "q True: ", q.q, np.sqrt((q.q**2).sum())
  
  path = ["../data/middle_cRmf.csv", "../data/left_cRmf.csv"]
  if path is None:
    vMFs_A = [vMF(np.array([1.,0.,0.]), 1.), vMF(np.array([0.,1.,0.]), 10.)]
    vMFs_B = [vMF(R_gt.dot(np.array([1.,0.,0.])), 1.),
        vMF(R_gt.dot(np.array([0.,1.,0.])), 10.)]
    vMFMM_A = vMFMM(np.array([0.5, 0.5]), vMFs_A)
    vMFMM_B = vMFMM(np.array([0.5, 0.5]), vMFs_B)
  else:
    vMFMM_A = LoadvMFMM(path[0])
    vMFMM_B = LoadvMFMM(path[1])

  tetras = s3.GetTetras(0)
  tetrahedra = s3.GetTetrahedra(0)

  maxIter = 1000
  fig = plt.figure()

  print "UpperBound"
  nodes = [Node(tetrahedron) for tetrahedron in tetrahedra]
  bb = BB(vMFMM_A, vMFMM_B, LowerBoundLog, UpperBoundLog)
  eps = bb.Compute(nodes, maxIter, q)

  print "UpperBoundConvexity"
  nodes = [Node(tetrahedron) for tetrahedron in tetrahedra]
  bb = BB(vMFMM_A, vMFMM_B, LowerBoundLog, UpperBoundConvexityLog)
  epsC = bb.Compute(nodes, maxIter, q)

  plt.subplot(2,1,1)
  plt.ylabel("Sum over angular deviation between closest vMF means.")
  plt.plot(eps[:,:2].sum(axis=1), 'r-', label="UpperBound")
  plt.plot(epsC[:,:2].sum(axis=1), 'g-', label="UpperBoundConvexity")
  plt.legend()
  plt.subplot(2,1,2)
  plt.ylabel("Angle between GT and inferred rotation.")
  plt.plot(eps[:,3], 'r-', label="UpperBound")
  plt.plot(epsC[:,3], 'g-', label="UpperBoundConvexity")
  plt.legend()
#  plt.subplot(3,1,3)
#  plt.ylabel("Cost Function.")
#  plt.plot(eps[:,2], 'r-', label="UpperBound")
#  plt.plot(epsC[:,2], 'g-', label="UpperBoundConvexity")
#  plt.legend()

  for i in range(1):
    gd = GradientDescent(vMFMM_A, vMFMM_B)
    q.setToRandom()
    R0 = q.toRot().R
    Rt, epsGC = gd.Compute(R0, maxIter, R_gt)

    plt.subplot(2,1,1)
    plt.plot(epsGC[:,:2].sum(axis=1), 'b-' , label="Gradient")
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(epsGC[:,3], 'b-', label="Gradient")
    plt.legend()
#    plt.subplot(3,1,3)
#    plt.plot(epsGC[:,2], 'b-', label="Gradient")
#    plt.legend()
  plt.savefig("./vMFMM_GradientAndBoundComparison_angleToGT.png")
  plt.savefig("./vMFMM_GradientAndBoundComparison_residuals.png")
  plt.show()

