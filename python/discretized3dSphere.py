import numpy as np
import mayavi.mlab as mlab
import time

class SphereGrid(object):
    def __init__(self, depth):
        self.depth = depth
        self.vertices = np.zeros((12, 3))
        self.tri = np.zeros((20, 3), dtype = np.int)
        self.InitGrid()
        for lvl in range(depth):
          self.SubdivideOnce()
        self.figm = None

    def InitGrid(self):
        ''' 
        Initialize a Icosahedron
        '''
        a = (1. + np.sqrt(5.0)) * 0.5

        self.vertices[0, :] = [-1, a, 0]
        self.vertices[1, :] = [1, a, 0]
        self.vertices[2, :] = [-1, -a, 0]
        self.vertices[3, :] = [1, -a, 0]

        self.vertices[4, :] = [0, -1, a]
        self.vertices[5, :] = [0, 1, a]
        self.vertices[6, :] = [0, -1, -a]
        self.vertices[7, :] = [0, 1, -a]

        self.vertices[8, :] = [a, 0, -1]
        self.vertices[9, :] = [a, 0, 1]
        self.vertices[10, :] = [-a, 0, -1]
        self.vertices[11, :] = [-a, 0, 1]

        self.tri[0, :] = [0, 11, 5]
        self.tri[1, :] = [0, 5, 1]
        self.tri[2, :] = [0, 1, 7]
        self.tri[3, :] = [0, 7, 10]
        self.tri[4, :] = [0, 10, 11]

        self.tri[5, :] = [1, 5, 9]
        self.tri[6, :] = [5, 11, 4]
        self.tri[7, :] = [11, 10, 2]
        self.tri[8, :] = [10, 7, 6]
        self.tri[9, :] = [7, 1, 8]

        self.tri[10, :] = [3, 9, 4]
        self.tri[11, :] = [3, 4, 2]
        self.tri[12, :] = [3, 2, 6]
        self.tri[13, :] = [3, 6, 8]
        self.tri[14, :] = [3, 8, 9]

        self.tri[15, :] = [4, 9, 5]
        self.tri[16, :] = [2, 4, 11]
        self.tri[17, :] = [6, 2, 10]
        self.tri[18, :] = [8, 6, 7]
        self.tri[19, :] = [9, 8, 1]
        self.tri_levels = [0, 20]
        self.vertices /= np.sqrt((self.vertices**2).sum(axis=1))[:, np.newaxis]

    def SubdivideOnce(self):
        '''
        Subdivide each of the existing triangles of the current
        SphereGrid into 4 triangles and pop out the inner corners of
        the new triangles to the sphere.
        '''
        n_vertices = self.vertices.shape[0]
        n_tri = self.tri.shape[0]
        print n_tri, n_vertices
        self.vertices.resize((n_vertices + n_tri * 3, 3), refcheck=False)
        self.tri.resize((n_tri * 5, 3))
        self.tri_levels.append(n_tri * 5)
        for i in range(n_tri):
            i0 = self.tri[i, 0]
            i1 = self.tri[i, 1]
            i2 = self.tri[i, 2]
            i01 = n_vertices + i*3 
            i12 = n_vertices + i*3 + 1
            i20 = n_vertices + i*3 + 2
            self.vertices[i01, :] = 0.5 * (self.vertices[i0, :] + self.vertices[i1, :]) 
            self.vertices[i12, :] = 0.5 * (self.vertices[i1, :] + self.vertices[i2, :]) 
            self.vertices[i20, :] = 0.5 * (self.vertices[i2, :] + self.vertices[i0, :]) 
            self.tri[n_tri + i*4, :] = [i0, i01, i20]
            self.tri[n_tri + i*4 + 1, :] = [i01, i1, i12]
            self.tri[n_tri + i*4 + 2, :] = [i12, i2, i20]
            self.tri[n_tri + i*4 + 3, :] = [i01, i12, i20]
        self.vertices[n_vertices::, :] /= np.sqrt((self.vertices[n_vertices::,:]**2).sum(axis=1))[:, np.newaxis]
        return 

    def GetTrianglID(self, direction):
      j = -1 
      id = []
      for i in range(20):
        i0 = self.tri[i, 0]
        i1 = self.tri[i, 1]
        i2 = self.tri[i, 2]
        if self.Intersects(self.vertices[i0, :], self.vertices[i1, :],
            self.vertices[i2, :], direction):
          j = i 
          print j
          #mlab.show(stop=True)
          break
      if j > 0:
        id.append(j)
        for lvl, n_lvl in enumerate(self.tri_levels[1:-1]):
          found_triangle = False
          for i in range(4):
            i_tri = n_lvl + j*4 + i
            i0 = self.tri[i_tri, 0]
            i1 = self.tri[i_tri, 1]
            i2 = self.tri[i_tri, 2]
            print '---', n_lvl, j*4, i
            if self.Intersects(self.vertices[i0, :], self.vertices[i1, :],
                self.vertices[i2, :], direction):
              id.append(i)
              j = n_lvl + (j*4 + i) 
              found_triangle = True
              #mlab.show(stop=True)
              break
          if not found_triangle:
            print "not good ... Did not find a intersection at level {}".format(lvl)
            break
          else:
            print "intersecting triangle {} at lvl {}".format(i, lvl)
      else:
        print "No intersection with the base triangles found."
      return id, j - self.tri_levels[-1]

    def Intersects(self, p0, p1, p2, direction):
      '''
      http://geomalgorithms.com/a06-_intersect-2.html
      '''
      if False:
        if self.figm is None:
          self.figm = mlab.figure()
        mlab.triangular_mesh(np.array([p0[0], p1[0], p2[0]]),
            np.array([p0[1], p1[1], p2[1]]), np.array([p0[2], p1[2],
              p2[2]]), np.array([0, 1, 2])[np.newaxis,:])
        mlab.plot3d([0,direction[0]], [0,direction[1]], [0,direction[2]])
        mlab.points3d([0], [0], [0], scale_factor=0.1)
        #yymlab.show(stop=True)
        #mlab.close(figm)

      n = np.cross(p1 - p0, p2 - p0)
      denom = n.dot(direction)
      if np.abs(denom) < 1e-6:
        # the direction is parallel to the plane of the triangle
        print "Direction and triangle are parallel. {}".format(denom)
        return False
      r = n.dot(p0) / denom
      if r < 0.:
        # We are intersecting on the other side
        print "Intersection is in the opposite direction. {}".format(r)
        return False
      intersection = r * direction
      v = p2 - p0
      u = p1 - p0
      w = intersection - p0
      uv = u.dot(v)
      uu = u.dot(u)
      vv = v.dot(v)
      wv = w.dot(v)
      wu = w.dot(u)
      denom = uv**2 - uu * vv
      s = ((uv * wv) - (vv * wu)) / denom
      t = ((uv * wu) - (uu * wv)) / denom
      if s >= 0 and t >= 0 and s+t <= 1.:
        return True
      else:
        print "Intersection is outside the triangle {}, {}".format(s,t)
        return False

    def GetNumTrianglesAtLevel(self, level):
      if level < len(self.tri_levels) - 1:
        return self.tri_levels[level+1]
      else:
        print "Do not have that many levels ({}).".format(level)
        return -1

    def GetNumLevels(self):
      return len(self.tri_levels) - 1

    def Plot(self, level, figm):
        if level + 1 >= len(self.tri_levels):
            print "Cannot plot this level in the SphereGrid"
            return 
        mlab.triangular_mesh(self.vertices[:, 0], self.vertices[:, 1],
            self.vertices[:, 2],
            self.tri[self.tri_levels[level]:self.tri_levels[level+1], :])
        return

class SphereHistogram:
  def __init__(self, sphereGrid, level = None):
    self.sphereGrid = sphereGrid
    if level is None:
      self.level = sphereGrid.GetNumLevels()
    else:
      self.level = level
    self.hist = np.zeros(self.sphereGrid.GetNumTrianglesAtLevel(level))

  def Compute(self, pts):
    for pt in pts:
      self.hist[self.sphereGrid.GetTrianglID(pt)[1]] += 1.
    print self.hist

  def Plot(self, level, figm):
    self.sphereGrid.Plot(self.level, figm)
    

run_tests = False
sphere = SphereGrid(2)
######################################### Intersection tests
if run_tests:
  p0 = np.array([0,0,1])
  p1 = np.array([1,0,1])
  p2 = np.array([0,1,1])
  t0 = time.clock()
  q = np.array([0.5, 0.5, -0.5])
  if sphere.Intersects(p0, p1, p2, q):
    print "Wrong intersection detected {}".format(q)
  q = np.array([0.5, 0.5, 1.])
  if not sphere.Intersects(p0, p1, p2, q):
    print "intersection not detected {}".format(q)
  q = np.array([0.5, 0., 1.])
  if not sphere.Intersects(p0, p1, p2, q):
    print "intersection not detected {}".format(q)
  q = np.array([0., 0.5, 1.])
  if not sphere.Intersects(p0, p1, p2, q):
    print "intersection not detected {}".format(q)
  q = np.array([1., 0.0, 0.])
  p0 = np.array([ 0.96193836,  0., -0.27326653])
  p1 = np.array([ 1.,  0.,  0.])
  p2 = np.array([ 0.95105652,  -0.26286556, -0.16245985])
  if not sphere.Intersects(p0, p1, p2, q):
    print "intersection not detected {}".format(q)
  p1 = np.array([ 0.96193836,  0., -0.27326653])
  p0 = np.array([ 1.,  0.,  0.])
  p2 = np.array([ 0.95105652,  -0.26286556, -0.16245985])
  if not sphere.Intersects(p0, p1, p2, q):
    print "intersection not detected {}".format(q)
  t1 = time.clock()
  print "tests tock {} ms".format((t1-t0)*1e3)

figm = mlab.figure(bgcolor=(1,1,1))
for lvl in range(3):
  sphere.Plot(lvl, figm)
  mlab.show(stop=True)

#q = np.array([0.5, 0.5, -0.5])
#q /= np.sqrt((q**2).sum())
#sphere.GetTrianglID(q)

q = np.array([[0.5, 0.5, 0.5],
  [0.5, 0.5, -0.5],
  [0.5, -0.5, -0.5],
  [0.5, -0.5, -0.5],
  [0.5, -0.5, -0.5],
  ])
sphereHist = SphereHistogram(sphere, 2)
sphereHist.Compute(q)

figm = mlab.figure()
sphereHist.Plot(2, figm)
mlab.show(stop=True)
