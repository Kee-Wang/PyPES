#!/usr/bin/env python
class plot():
    def version(self):
        print('Plot Version: 0.0.1')
        return None

    def heatmap(self,file,bin=100,xtitle='xtitle',ytitle='ytitle',ztitle='ztitle',title=None):
        """This is to construct 3d heat map. """

        import matplotlib.pyplot as plt
        import numpy as np
        self.version()

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
    def scatter(self,file,xtitle='xtitle',ytitle='ytitle',delimiter=' ',cols=(0,1)):
        """This is to plot simple scatter plot for 2 columns."""
        import matplotlib.pyplot as plt
        self.version()

        plt.plotfile(file, delimiter=delimiter, cols=cols,names=(xtitle, ytitle))
        plt.show()
