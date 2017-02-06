import matplotlib.pyplot as plt
from pypes.configs import configs
"""This script was used for v2b >= 0 dissociation cut plot for CO2-H2O hybrid PES"""

a = configs('v2bgt0_pes.dat')
aa = configs('v2bgt0.abE')
bb = aa.list()
b = a.list()
E = list()
EE = list()
E0 = list()
dis = list()
dis2 = list()
c = list()
i=0
for config in b:
    #print(config)
    E.append(config[1][0][0]*219474.63)
    if aa.distance(bb[i], atom_A=3, atom_B=6) <=3.3 and i%2 ==0:
        EE.append(bb[i][1][0][0]*219474.63)
        dis2.append(aa.distance(bb[i], atom_A=3, atom_B=6))
    elif aa.distance(bb[i], atom_A=3, atom_B=6) <= 4 and i%4==0:
        EE.append(bb[i][1][0][0] * 219474.63)
        dis2.append(aa.distance(bb[i], atom_A=3, atom_B=6))
    elif i%8 == 0:
        EE.append(bb[i][1][0][0]*219474.63)
        dis2.append(aa.distance(bb[i], atom_A=3, atom_B=6))
    E0.append(config[1][0][0]*0)
    dis.append(a.distance(config,atom_A=3,atom_B=6))

    i = i+1
    if a.distance(config,atom_A=3,atom_B=6)>=10 and config[1][0][0]>0:
        c.append(config)
#a.molden(c)

fig = plt.figure(figsize=(6,4))
#plt.gcf().subplots_adjust(bottom=0.15)
ax = fig.add_subplot(111)
ax.set_ylabel('V$_{\mathrm{2b}}$ (cm$^{-1}$)')
ax.set_xlabel('r$_{\mathrm{C--O}}$ ($\mathrm{\AA}$)')
axes = plt.gca()

ax.plot(dis, E,  color='k',linewidth='2',label='PES')
ax.scatter(dis2, EE, s=50, color='m', marker='o',linewidth='2 ',label='$ab\ initio$',facecolors='None')
ax.legend()
#ax.plot(dis, E0,  color='r', label='y=0',
                   # linewidth='2')
plt.gcf().subplots_adjust(bottom=0.15)
axes = plt.gca()
#axes.set_xlim([1, 15])
#axes.set_ylim([-10, 1200])
plt.tight_layout()
plt.show()
fig.savefig('v2bgt0.eps', format='eps', dpi=1200)
