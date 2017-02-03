#!/usr/bin/env python
class plot():
    def version(self):
        print('Plot Version: 0.0.3')
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

        hb = ax.hexbin(x, y, z, bin, cmap='inferno')

        ax.axis([xmin, xmax, ymin, ymax])
        ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label(ztitle)

        plt.show()
        return None
    def scatter(self,file,col=[0,1],xtitle='xtitle',ytitle='ytitle',title=' '):
        """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
        import matplotlib.pyplot as plt
        import numpy as np
        self.version()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        data = np.loadtxt(file, unpack=True)
        for i in col:
            if i == 0:
                continue
            plt.scatter(data[0],data[i])
        #ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)

        plt.show()
        return None
    def line(self,file,col=[0,1],xtitle='xtitle',ytitle='ytitle',title=' '):
        """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
        import matplotlib.pyplot as plt
        import numpy as np

        self.version()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        data = np.loadtxt(file, unpack=True)
        for i in col:
            if i == 0:
                continue
            plt.plot(data[0],data[i])
        #ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)

        plt.show()
        return None