
"""This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
import matplotlib.pyplot as plt
import numpy as np
import math

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
sizex=6
sizey=4
fig = plt.figure()
#figsize=(sizex, sizey)
ax = fig.add_subplot(111)

axes = plt.gca()
#axes.set_xlim([xmin, xmax])
#axes.set_ylim([ymin, ymax])

file = 'nbas.prn'
data = np.loadtxt(file,  unpack=True)  # Read columns
label=['NMAX=2','NMAX=3', 'NMAX=4', 'NMAX=5']
for i in range(4):
    plt.semilogy(data[0], data[i+1], linewidth=2, label=label[i])
    plt.scatter(data[0], data[i+1],s=30,linewidth=1,facecolors='none')

ax.legend(loc=4)
#loc=legendloc

ax.set_xlabel('MAXBAS')
ax.set_ylabel('Size of VCI matrix')
plt.tight_layout()
fig = plt.gcf()
fig.savefig('vci_matrix_size.eps', format='eps', dpi=1200)
plt.show()
