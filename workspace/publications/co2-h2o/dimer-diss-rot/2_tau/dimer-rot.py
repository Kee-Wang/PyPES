#from pypes.plot import plot


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
        print('Plot Version: 0.0.7 --By Kee')
        return None


    def line(self,file=None,file2=None,data=None,xmin=None,xmax=None,ymin=None,ymax=None,col=(0,1),col2=(0,1),xtitle='xtitle',ytitle='ytitle',title=' ',save=None,linewidth=2,label=None,label2=None,legendloc=None,sizex=6,sizey=4):
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
            data = np.loadtxt(file,usecols=col, unpack=True) #Read columns
            count = 0
            np.delete(data[1],0)
            np.delete(data[0],0)
            for coord in data[0]:
                data[0][count] = coord - 90
                if coord <= 90:
                    data[0][count] += 360
                count = count + 1



              #print(coord)
        for i in col2:
            if i == 0:
                continue
            plt.plot(data[0], data[i], linewidth=linewidth,label = label)


        if file2 is not None:
            data2 = np.loadtxt(file2,usecols=col2, unpack=True)
            count = 0

            data2[1][0] = data2[1][-1]
            data2[1][1] = data2[1][-2]
            #print(len(data[0]),len(data2))
            for coord in data2[0]:
                print(data2[1][count])
                data2[0][count] = coord - 90
                if coord <= 90 and coord >=0:
                    data2[0][count] += 360

                #if coord >=89 and coord <=91:
                #    print(data2[1][count])
                #if data2[0][count]>=250 and data2[0][count] <=300:
                #    data2[1][count] = 0
                count = count + 1

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
        #fig.set_tight_layout(True)
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


a = plot()
a.line(file='rigid',
       file2='rlx',
       xtitle=r'Internal rotation angle $\tau$ (degree)',
       ytitle='CO$_2$-H$_2$O potential (cm$^{-1}$)',
       label='rigid',
       label2='relaxed',
       legendloc=9,
       xmin = 0,
       xmax=360,
       ymin=0,
       ymax=400,
       save='dimer-rot.eps')
