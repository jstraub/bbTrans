import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#paper
mpl.rc('font',size=30) 
mpl.rc('lines',linewidth=3.)
figSize = (14, 5.5)
figSize = (14, 10)
figSize = (14, 12)

def ToRad(deg):
  return deg*np.pi/180.
def ToDeg(rad):
  return rad/np.pi*180.

eta0 = ToRad(72.)
eta = ToRad(np.exp(np.linspace(np.log(0.001), np.log(180.),1000)))

a = np.cos(eta0)
b = np.cos(eta*0.5)

N = np.ceil(np.log((1./a - 1.)/(1./b - 1.))/np.log(2.))
N[N<0.] = 0.

fig = plt.figure(figsize = figSize, dpi = 80, facecolor="w",
    edgecolor="k")
ax = plt.subplot(111)
plt.plot(ToDeg(eta), N, color="orange")
#plt.plot(N, ToDeg(eta))
#ax.set_yscale("log")
ax.set_xscale("log")
plt.xlabel("desired angular precision $\eta$ [deg]")
plt.ylabel("number of required subdivisions")
plt.savefig("../angularPrecisionVsSubdivisionLvl.png", figure=fig)
plt.show()
