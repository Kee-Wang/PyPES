import matplotlib.pyplot as plt
from pypes.plot import plot
import numpy as np
from scipy.interpolate import spline

q1_4mode = '4-mode-q1.dat'
q2_4mode = '4-mode-q2.dat'
q3_4mode = '4-mode-q3.dat'
q4_4mode = '4-mode-q4.dat'
q1_9mode = '9-mode-q6.dat'
q2_9mode = '9-mode-q7.dat'
q3_9mode = '9-mode-q8.dat'
q4_9mode = '9-mode-q9.dat'
q1_co2 = 'co2-q1.dat'
q2_co2 = 'co2-q2.dat'
q3_co2 = 'co2-q3.dat'
q4_co2 = 'co2-q4.dat'

q1 = list()
q2 = list()
q3 = list()
q4 = list()
q1_4 = np.loadtxt(q1_4mode, usecols=(0,2),unpack=True)
q1.append(q1_4)
q2_4 = np.loadtxt(q2_4mode, usecols=(0,2),unpack=True)
q2.append(q2_4)
q3_4 = np.loadtxt(q3_4mode, usecols=(0,2),unpack=True)
q3.append(q3_4)
q4_4 = np.loadtxt(q4_4mode, usecols=(0,2),unpack=True)
q4.append(q4_4)
q1_9 = np.loadtxt(q1_9mode, usecols=(0,2),unpack=True)
q1.append(q1_9)
q2_9 = np.loadtxt(q2_9mode, usecols=(0,2),unpack=True)
q2.append(q2_9)
q3_9 = np.loadtxt(q3_9mode, usecols=(0,2),unpack=True)
q3.append(q3_9)
q4_9 = np.loadtxt(q4_9mode, usecols=(0,2),unpack=True)
q4.append(q4_9)
q1_c = np.loadtxt(q1_co2, usecols=(0,2),unpack=True)
q1.append(q1_c)
q2_c = np.loadtxt(q2_co2, usecols=(0,2),unpack=True)
q2.append(q2_c)
q3_c = np.loadtxt(q3_co2, usecols=(0,2),unpack=True)
q3.append(q3_c)
q4_c = np.loadtxt(q4_co2, usecols=(0,2),unpack=True)
q4.append(q4_c)

q = [q1,q2,q3,q4]
#For mode 1

xtitle='Q(A0)'
ytitle=r'Potential (cm$^{-1}$)'
legendloc=1
#fig, axs = plt.subplots(2,2, figsize=(15, 6), facecolor='w', edgecolor='k')
#axs = axs.ravel()
#for j in range(0,5):
    #for i in range(1,3):

#        axs[j].set_xlabel(xtitle)
 #       axs[j].set_ylabel(ytitle)
  #      axs[j].legend(loc=legendloc)
   #     axs[j].plot(q[j][i][0],q[j][i][1])
#plot().save(fig,'q1.eps')


fig = plt.figure(figsize=(15,10))
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
legendloc = 0
label = [r'$V_{\mathrm{CO}_2}^{(1)}+V_{{\mathrm{CO}_2}-{\mathrm{H}_{\mathrm{2}}\mathrm{O}}}^{(2)}}$',r'$V_{\mathrm{CO}_2}^{(1)}$']
title = [r'${\mathrm{CO}_2} \mathrm{{\ }bend{\ }1} $',r'${\mathrm{CO}_2} \mathrm{\ bend{\ }2}$',r'${\mathrm{CO}_2} \mathrm{\ sym{\ }stretch}$',r'${\mathrm{CO}_2} \mathrm{\ asym{\ }stretch}$']
for i in range(1,3):
    ax = fig.add_subplot(221)
    plt.plot(q1[i][0],q1[i][1],label=label[i-1])
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.set_title(title[0])
    ax.legend(loc=legendloc)

for i in range(1,3):
    ax = fig.add_subplot(222)
    plt.plot(q2[i][0],q2[i][1],label=label[i-1])
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.set_title(title[1])
    ax.legend(loc=legendloc)

for i in range(1,3):
    ax = fig.add_subplot(223)
    plt.plot(q3[i][0],q3[i][1],label=label[i-1])
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.set_title(title[2])
    ax.legend(loc=legendloc)

for i in range(1,3):
    ax = fig.add_subplot(224)
    plt.plot(q4[i][0],q4[i][1],label=label[i-1])
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.set_title(title[3])
    ax.legend(loc=legendloc)
plot().save(fig,'q_all.eps')