from pypes.plot import plot
import matplotlib.pyplot as plt
import numpy as np
col = (0,1)
fig = plt.figure()
ax = fig.add_subplot(111)
data = np.loadtxt('v2b_ab_pes.dat', unpack=True)
for i in col:
    if i == 0:
        continue
    plt.scatter(data[0], data[i])
plt.plot(range(-500,200), range(-500,200), 'r--')
#plt.plot(data[0],data[0])
# ax.set_title(title)
ax.set_xlabel('MP2-aVTZ (cm$^{-1}$)')
ax.set_ylabel('PES$_{\mathrm{2b}}$ (cm$^{-1}$)')

plot().save(fig,'v2b_ab_pes.eps')