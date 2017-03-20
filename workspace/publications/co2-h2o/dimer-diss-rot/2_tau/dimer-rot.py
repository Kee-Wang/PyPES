from pypes.plot import plot
a = plot()
a.line(file='rigid',
       file2='rlx',
       xtitle=r'Internal rotation angle $\tau$ (degree)',
       ytitle='CO$_2$-H$_2$O potential (cm$^{-1}$)',
       label='rigid',
       label2='relaxed',
       legendloc=8,
       xmin = 0,
       xmax=360,
       ymin=0,
       ymax=400,
       save='dimer-rot.eps')