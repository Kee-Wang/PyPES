from pypes.plot import plot

manual_locations = [(20, 20),(20,40),(50,20),(40,60),(60,20),(80,60),(80,100),(40,110),(30,130),(30,150),(30,160),(55,140),(45,150),(110,80),(160,110),(140,120),(140,150),(160,150),(160,160),(120,50),(140,60),(130,20),(150,40),(145,20),(160,20)]
plot().contour('body-center-rotation.dat',manual_locations=manual_locations,sizex=6,sizey=4,xtitle=r'${\theta}$',ytitle='$\phi$',ztitle='CW20 potential (cm$^{-1}$)',title='',save='cw20-body-center-rot.eps',gridsize=150,linewidths=0.1)
