#!/usr/bin/env python
class plot():
    def heatmap(self,file='log_cart_compare',gridsize=100,xtitle='xtitle',ytitle='ytitle',ztitle='ztitle',title=None):
        """This is to construct 3d heat map. """

        import matplotlib.pyplot as plt
        import numpy as np

        (x, y, z) = np.loadtxt(file, unpack=True)

        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # fig, ax = plt.subplots( sharey=True, figsize=(7, 4))
        # fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

        hb = ax.hexbin(x, y, z, gridsize, cmap='inferno')

        ax.axis([xmin, xmax, ymin, ymax])
        ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label(ztitle)

        plt.show()
        return None