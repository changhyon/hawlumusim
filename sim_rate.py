import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("rate_profile.txt")

s     = d[:,0]
mwe   = d[:,1]
ratio = d[:,3]
rate  = d[:,4]

plt.figure(figsize=(12,8))

plt.subplot(2,1,1)
plt.plot(s, mwe, linewidth=2)
plt.ylabel("Overburden [m.w.e.]")
plt.grid()

plt.subplot(2,1,2)
plt.yscale('log')
plt.plot(s, rate, linewidth=2)
plt.xlabel("Tunnel position [m]")
plt.ylabel("Predicted coincidence rate [Hz]")
plt.grid()

plt.tight_layout()
plt.show()
