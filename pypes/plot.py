#!/usr/bin/env python
class plot():

    def __init__(self):
        """Global setting so that it is naturally paper level plot"""
        import matplotlib.pyplot as plt


        SMALL_SIZE = 12
        MEDIUM_SIZE = 14
        BIGGER_SIZE = 16

        plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    def version(self):
        print('Plot Version: 0.0.5 --By Kee')
        return None
    def heatmap(self,file=None,bin=100,xtitle='xtitle',ytitle='ytitle',ztitle='ztitle',title=None,save=None):
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

        hb = ax.hexbin(x, y, z, bin, cmap='inferno')

        ax.axis([xmin, xmax, ymin, ymax])
        ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label(ztitle)

        self.save(fig,save)

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
    def line_back(self,file,col=[0,1],xtitle='xtitle',ytitle='ytitle',title=' ',save=None,linewidth=2):
        """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
        import matplotlib.pyplot as plt
        import numpy as np

        self.version()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        data = np.loadtxt(file, unpack=True) #Read columns
        for i in col:
            if i == 0:
                continue
            plt.plot(data[0],data[i],linewidth=linewidth)
        #ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        self.save(fig, save)

        return None
    def line(self,file=None,data=None,col=[0,1],xtitle='xtitle',ytitle='ytitle',title=' ',save=None,linewidth=2):
        """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
        import matplotlib.pyplot as plt
        import numpy as np

        self.version()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if file is not None:
            data = np.loadtxt(file, unpack=True) #Read columns
        for i in col:
            if i == 0:
                continue
            plt.plot(data[0],data[i],linewidth=linewidth)
        #ax.set_title(title)
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        self.save(fig, save)

        return None
    def save(self,fig,save=None):
        """The module is to show and save figure"""

        import matplotlib.pylab as plt

        plt.tight_layout()
        fig = plt.gcf()
        plt.show()
        if save is None:
            decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
            if decision is 'y':
                filename = input('Please specify .eps (1200 dpi) filename: ').strip()

                fig.savefig(filename, format='eps', dpi=1200)
                print('Plot saved to {}.'.format(save))
        else:
            fig.savefig(save, format='eps', dpi=1200)
            print('Plot saved to {}.'.format(save))


        return None