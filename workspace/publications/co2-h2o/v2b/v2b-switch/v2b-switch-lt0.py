import matplotlib.pyplot as plt
from pypes.configs import configs


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

#This script is used for making 3-panel plot for CO2-H2O PES.

a=configs('diss_a7_long.dat')
a1=a.list()
b=configs('diss_bas_whole.dat')
b1=b.list()
c=configs('diss_100A.abE')
c1 = c.list()

d = configs('b_diss_long.dat')
d1 = d.list()
e = configs('b_diss_bas.dat')
e1 = e.list()
f = configs('b_diss_ab.dat')
f1 = f.list()

i = configs('c_diss_ab.dat')
i1=i.list()
h = configs('c_diss_bas.dat')
h1=h.list()
g = configs('c_diss_long.dat')
g1 = g.list()


configss = [a1,b1,c1,d1,e1,f1,g1,h1,i1]

xmin =4;xmax = 15 ;ymin = -200;ymax = 10 ;xmin2 =4;xmax2 =15;ymin2 = -10;ymax2 = 7
atomA=3
atomB=6
s=100
title=''
label='default'
#fontsize=None

aucm = 219474.63

n = len(configss)
dis = list()
e = list()

for i in range(0, n):
    configss[i] = a.configs_check(configss[i])
    configss[i] = a.sort(configss[i], key='distance', subkey1=3, subkey2=6)

fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(311)
axes = plt.gca()
axes.set_xlim([4, 15])
plt.ylim(ymin, ymax)

ax2 = fig.add_subplot(312)  # ,sharex=True)
axes = plt.gca()
axes.set_xlim([4, 15])
plt.ylim(-200, 10)

ax3 = fig.add_subplot(313)  # ,sharex=True)
axes = plt.gca()
axes.set_xlim([4, 15])
plt.ylim(-200, 10)

for n in range(0, len(configss)):
    if (n + 1) % 3 == 0:
        i = 0
        j = 0
        e_temp = list()
        e_temp_ref = list()
        dis_temp = list()
        for config in configss[n]:
            i = i + 1
            if i % 8 == 0:
                j = j + 1
                print(i, j)
                dis_temp.append(a.distance(config, atomA, atomB))
                e_temp.append(config[1][0][0] * aucm)
        dis.append(dis_temp)
        e.append(e_temp)
    else:
        e_temp = list()
        e_temp_ref = list()
        dis_temp = list()
        for config in configss[n]:
            dis_temp.append(a.distance(config, atomA, atomB))
            e_temp.append(config[1][0][0] * aucm)

        dis.append(dis_temp)
        e.append(e_temp)
label1='fit-LR'
label2='fit-SR'
label3='$ab\ initio$'

ax.plot(dis[0], e[0], label=label1, linewidth=2.0)
ax.plot(dis[1], e[1], '--', label=label2, linewidth=2.0)
ax.scatter(dis[2], e[2], s=200, color='k', marker='.', label=label3,
           facecolors='none',
           linewidth='2')
ax.set_title(title)
ax.legend(loc=4)
ax.set_ylabel('$V_{\mathrm{2b}}$ (cm$^{-1}$)')
ax.text(4.2, -10, '(a)')

ax2.plot(dis[3], e[3], label=label1, linewidth=2.0)
ax2.plot(dis[4], e[4], '--', label=label2, linewidth=2.0)
ax2.scatter(dis[5], e[5], s=200, color='k', marker='.', label=label3,
            facecolors='none',
            linewidth='2')
ax2.set_title(title)
ax2.legend(loc=4)
ax2.set_ylabel('$V_{\mathrm{2b}}$ (cm$^{-1}$)')
ax2.text(4.2, -10, '(b)')

ax3.plot(dis[6], e[6], label=label1, linewidth=2.0)
ax3.plot(dis[7], e[7], '--', label=label2, linewidth=2.0)
ax3.scatter(dis[8], e[8], s=200, color='k', marker='.', label=label3,
            facecolors='none',
            linewidth='2')
ax3.set_title(title)
ax3.legend(loc=4)
ax3.set_ylabel('$V_{\mathrm{2b}}$ (cm$^{-1}$)')
ax3.text(4.2, -10, '(c)')
ax3.set_xlabel('$R_{\mathrm{CO}}$ ($\mathrm{\AA}$)')

#for tick in ax.xaxis.get_major_ticks():
    #tick.label.set_fontsize(fontsize)

#for tick in ax2.xaxis.get_major_ticks():
   # tick.label.set_fontsize(fontsize)

#for tick in ax3.xaxis.get_major_ticks():
    #tick.label.set_fontsize(fontsize)

plt.tight_layout()
plt.show()
# decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
decision = 'y'
if decision is 'y':
    # filename = input('Please specify .eps (1200 dpi) filename: ').strip()
    filename = 'v2b-switch-lt0.eps'
    fig.savefig(filename, format='eps', dpi=1200)
    print('Plot saved to {}.'.format(filename))

