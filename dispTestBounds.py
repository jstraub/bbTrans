import numpy as np

bs = np.loadtxt('./testBound.csv').T
import matplotlib.pyplot as plt

plt.figure()
ax = plt.subplot(2,1,1)
plt.plot(bs[0,:],label="lb")
plt.plot(bs[1,:],label="ub")
plt.plot(bs[2,:],label="ubC")
ax.set_yscale("log")
plt.legend()
plt.subplot(2,1,2)
plt.plot(bs[3,:],label="sorting criterium")
plt.show()
