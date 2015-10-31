import numpy as np
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from js.utils.plot.colors import colorScheme

mpl.rc('font',size=25) 
mpl.rc('lines',linewidth=3.)
figSize = (14, 5.5)
figSize = (14, 10)
figSize = (14, 12)

c1 = colorScheme("labelMap")["turquoise"]
c2 = colorScheme("labelMap")["orange"]
c3 = colorScheme("labelMap")["green"]

#bs = np.loadtxt('./testBound.csv').T
for path in ['./bb_bounds_S3_t0.csv', './bb_bounds_R3_t0.csv']:
  bs = np.loadtxt(path).T
  bs[:3,:] = np.log(bs[:3,:])/np.log(10)
  ids = np.argsort(bs[0,:])

  fig = plt.figure(figsize = figSize, dpi = 80, facecolor="w",
      edgecolor="k")
  ax = plt.subplot(1,1,1)
#  ax.set_yscale("log")
  plt.plot(bs[1,ids],label="independent component upper bound",
      color=c2)
  plt.plot(bs[2,ids],label="joint convex upper bound",color=c1)
  plt.plot(bs[0,ids],label="lower bound",color=c3)
  plt.legend(loc="best")
  plt.ylabel("log$_{10}$(bounds)")
  plt.ylim([bs.min(), bs.max()])
  plt.xlim([0, bs.shape[1]])
  plt.xlabel("nodes in first level of B&B tree")
  #plt.subplot(2,1,2)
  #plt.plot(bs[3,:],label="sorting criterium")
  plt.savefig(re.sub("csv","png",path), figure = fig)
  plt.show()
