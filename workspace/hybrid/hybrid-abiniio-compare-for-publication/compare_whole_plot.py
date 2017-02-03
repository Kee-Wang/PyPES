import matplotlib.pyplot as plt
from pypes.configs import configs
#This script is used for making 3-panel plot
c=configs('diss_100A.abE')
c1 = c.list()
a=configs('diss_a7_long.dat')
a1=a.list()
b=configs('diss_bas_whole.dat')
b1=b.list()
i = configs('c_diss_ab.dat')
i1=i.list()
h = configs('c_diss_bas.dat')
h1=h.list()
g = configs('c_diss_long.dat')
g1 = g.list()
d = configs('b_diss_long.dat')
d1 = d.list()
e = configs('b_diss_bas.dat')
e1 = e.list()
f = configs('b_diss_ab.dat')
f1 = f.list()


# "Compare long with whole range"
color = ['m','c','r','k']
marker = ['.', 'd', 's', 'o']
label = [r'${\mathrm{PES}}^{\mathrm{lr}}_{\mathrm{sr}}$','$\mathrm{PES}^{\mathrm{2b}}_\mathrm{sr}$', '$\mathrm{V}^{\mathrm{2b}}_{ab\ initio}$']

compare = [a1,b1,c1,d1,e1,f1,g1,h1,i1]
xmin =4;xmax = 15 ;ymin = -200;ymax = 10 ;xmin2 =4;xmax2 =15;ymin2 = -10;ymax2 = 7
a.co2h2o(compare,atomA=3,atomB=6,s=200,color=color,marker=marker,label=label,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,xmin2=xmin2,xmax2=xmax2,ymin2=ymin2,ymax2=ymax2)
