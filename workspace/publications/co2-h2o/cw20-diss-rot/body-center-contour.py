from pypes.plot import plot
plot().contour('body-center-rotation.dat',sizex=6,sizey=4,xtitle=r'${\theta}$',ytitle='$\phi$',ztitle='CW20 potential (cm$^{-1}$)',title='',save='cw20-body-center-rot.eps',gridsize=150,linewidths=0.1)
