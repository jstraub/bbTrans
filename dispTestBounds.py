import numpy as np

bs = np.loadtxt('./testBound.csv').T
import matplotlib.pyplot as plt

plt.figure()
plt.plot(bs[0,:],label="lb")
plt.plot(bs[1,:],label="ub")
plt.plot(bs[2,:],label="ubC")
plt.legend()
plt.show()
