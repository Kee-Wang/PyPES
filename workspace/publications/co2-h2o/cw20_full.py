import matplotlib.pyplot as plt
from pypes.plot import plot




def line(file=None, file2=None, data=None, xmin=None, xmax=None, ymin=None, ymax=None, col=(0, 1), col2=(0, 1),
         xtitle='xtitle', ytitle='ytitle', title=' ', save=None, linewidth=2, label=None, label2=None, legendloc=None,
         sizex=6, sizey=4):
    """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
    import matplotlib.pyplot as plt
    import numpy as np
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16

    plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    co2_bend1 = 654
    co2_bend2 = 655
    co2_sym = 1352
    co2_asym = 2418
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)
    if xmin is not None:
        axes = plt.gca()
        axes.set_xlim([xmin, xmax])
        axes.set_ylim([ymin, ymax])

    if file is not None:
        data = np.loadtxt(file, usecols=col, unpack=True)  # Read columns
    for i in col2:
        if i == 0:
            continue
        plt.plot(data[0], data[i], linewidth=1, label=label)

    if file2 is not None:
        data2 = np.loadtxt(file2, usecols=col2, unpack=True)
        for i in col2:
            if i == 0:
                continue
            plt.plot(data2[0], data2[i], '--', linewidth=linewidth, label=label2)

    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.legend(loc=legendloc)
    ax.annotate(r'${\nu}_{\mathrm{b1}}$ (CO$_2$), ${\nu}_{\mathrm{b2}}$ (CO$_2$)', xy=(co2_bend1, 4.7), xytext=(500-135,6.5),
                arrowprops=dict(arrowstyle="->"))
    ax.annotate(r'$\nu_{\mathrm{sym}}$(CO$_2$)', xy=(co2_sym, 1.1), xytext=(co2_sym-300,2.7),
                arrowprops=dict(arrowstyle="->"))
    ax.annotate(r'$\nu_{\mathrm{asym}}$(CO$_2$)', xy=(co2_asym, 1.1), xytext=(co2_asym-300,2.7),
                arrowprops=dict(arrowstyle="->"))
    #ax.annotate(' ', xy=(co2_bend2, 4.7), xytext=(500,7),
    #            arrowprops=dict(arrowstyle="->"))
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

line('cw20_full.dat',xtitle='Frequency (cm$^{-1}$)',ytitle='Density of states (arb. unit)',save='cw20-full-nma.eps')