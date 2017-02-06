from pypes.plot import plot
a = plot()
a.line(file='rigid',
       file2='rlx',
       xtitle='Internal rotation angle $\phi$',
       ytitle='energy (cm$^{-1}$)',
       label='rigid',
       label2='relaxed',
       legendloc=8,
       xmin = 0,
       xmax=360,
       ymin=0,
       ymax=400,
       save='dimer-rot.eps')