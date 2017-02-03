import matplotlib.pyplot as plt
import numpy as np
xtitle='angle'
ytitle='angle2'
ztitle='energy'
title='hex'
(x,y,z) = np.loadtxt('log_cart_compare', unpack=True)

xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

fig = plt.figure()
ax = fig.add_subplot(111)
#fig, ax = plt.subplots( sharey=True, figsize=(7, 4))
#fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

hb = ax.hexbin(x, y, z,gridsize=100, cmap='inferno')

ax.axis([xmin, xmax, ymin, ymax])
ax.set_title(title)
ax.set_xlabel(xtitle)
ax.set_ylabel(ytitle)
cb = fig.colorbar(hb, ax=ax)
cb.set_label(ztitle)

plt.show()
