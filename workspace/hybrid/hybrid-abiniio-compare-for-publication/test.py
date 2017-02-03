import matplotlib.pyplot as plt
from pypes.configs import configs
a = configs('gt_6.xyz')
b = a.list()
E = list()
dis = list()
c = list()
for config in b:
    #print(config)
    E.append(config[1][0][0]*219474.63)
    dis.append(a.distance(config,atom_A=3,atom_B=6))
   # if a.distance(config,atom_A=3,atom_B=6)>=10 and config[1][0][0]>0:
   #     c.append(config)
#a.molden(c)

fig = plt.figure(figsize=(12,12))
plt.gcf().subplots_adjust(bottom=0.15)
ax = fig.add_subplot(111)
ax.set_ylabel('V$_{\mathrm{2b}}$ (cm$^{-1}$)')
ax.set_xlabel('r$_{\mathrm{Si--O}}$ ($\AA$)')
axes = plt.gca()

ax.scatter(dis, E, s=2, color='k', marker='.', label='$\mathrm{V}^{\mathrm{2b}}_{ab\ initio}$',
                    linewidth='2')
plt.gcf().subplots_adjust(bottom=0.15)
ax2 = fig.add_subplot(312)
ax2.set_ylabel('V$_{\mathrm{2b}}$ (cm$^{-1}$)')
ax2.set_xlabel('r$_{\mathrm{Si--H}}$ ($\AA$)')
axes = plt.gca()
axes.set_xlim([6, 15])
axes.set_ylim([-200, 100])
ax2.scatter(dis, E, s=2, color='k', marker='.', label='$\mathrm{V}^{\mathrm{2b}}_{ab\ initio}$',
                    linewidth='2')
ax2 = fig.add_subplot(313)
ax2.set_ylabel('V$_{\mathrm{2b}}$ (cm$^{-1}$)')
ax2.set_xlabel('r$_{\mathrm{Si--H}}$ ($\AA$)')
axes = plt.gca()
axes.set_xlim([9, 25])
axes.set_ylim([-6, 5])
ax2.scatter(dis, E, s=2, color='k', marker='.', label='$\mathrm{V}^{\mathrm{2b}}_{ab\ initio}$',
                    linewidth='2')
plt.show()
fig.savefig('v2bdis.eps', format='eps', dpi=1200)
