import numpy as np
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from js.utils.plot.colors import colorScheme

mpl.rc('font',size=45) 
mpl.rc('lines',linewidth=4.)
figSize = (14, 5.5)
figSize = (14, 12)
figSize = (9, 6)
figSize = (7, 5)
figSize = (14, 10)
figSize = (14, 16)

c1 = colorScheme("labelMap")["turquoise"]
c2 = colorScheme("labelMap")["orange"]
c3 = colorScheme("labelMap")["green"]

#bs = np.loadtxt('./testBound.csv').T
fig = plt.figure(figsize = figSize, dpi = 80, facecolor="w",
    edgecolor="k")
for i,path in enumerate(['./bb_bounds_S3_t0.csv','./bb_bounds_TpS3_t0.csv',
  './bb_bounds_R3_t0.csv']):
  bs = np.loadtxt(path).T
  bs[:3,:] = np.log(bs[:3,:])/np.log(10)
  ids = np.argsort(bs[0,:])

  ax = plt.subplot(3,1,i+1)
#  ax.set_yscale("log")
  plt.plot(bs[1,ids],label="indep. UB",
      color=c3)
  plt.plot(bs[2,ids],label="convex UB",color=c2)
  plt.plot(bs[0,ids],label="LB",color=c1)
  plt.ylabel("log$_{10}$(bound)")
  plt.ylim([bs.min(), bs.max()])
  plt.xlim([0, bs.shape[1]])
  #plt.subplot(2,1,2)
  #plt.plot(bs[3,:],label="sorting criterium")

  if i==2:
    plt.legend(loc="best")
    plt.xlabel("nodes in first level of BB tree")
  plt.tight_layout(0.4)
plt.savefig("./bb_bounds_S3_R3_t0.png", figure = fig)
#  plt.show()
