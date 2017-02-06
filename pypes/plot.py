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
        print('Plot Version: 0.0.6 --By Kee')
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
    def line_backup(self,file,col=[0,1],xtitle='xtitle',ytitle='ytitle',title=' ',save=None,linewidth=2):
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
    def line(self,file=None,file2=None,data=None,xmin=None,xmax=None,ymin=None,ymax=None,col=[0,1],col2=[0,1],xtitle='xtitle',ytitle='ytitle',title=' ',save=None,linewidth=2,label=None,label2=None,legendloc=None,sizex=6,sizey=4):
        """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
        import matplotlib.pyplot as plt
        import numpy as np

        self.version()
        fig = plt.figure(figsize=(sizex, sizey))
        ax = fig.add_subplot(111)
        if xmin is not None:
            axes = plt.gca()
            axes.set_xlim([xmin, xmax])
            axes.set_ylim([ymin, ymax])

        if file is not None:
            data = np.loadtxt(file, unpack=True) #Read columns
        for i in col2:
            if i == 0:
                continue
            plt.plot(data[0], data[i], linewidth=linewidth,label = label)


        if file2 is not None:
            data2 = np.loadtxt(file2, unpack=True)
            for i in col2:
                if i == 0:
                    continue
                plt.plot(data2[0],data2[i],'--',linewidth=linewidth,label = label2)


        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        ax.legend(loc=legendloc)
        """Loc:
        best -- 0
        upper right -- 1
        upper left -- 2
        lower left -- 3
        lower right -- 4
        right -- 5
        center left -- 6
        center right -- 7
        lower center -- 8
        upper center -- 9
        center -- 10
        """
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