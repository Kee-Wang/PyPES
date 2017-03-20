from pypes.plot import plot
a = plot()
a.line(file='r5_rigid_cut.dat',
       file2='r5_rlx_cut.dat',
       xmin = 2.2,
       xmax = 15,
       ymin = 0,
       ymax = 1200,
       xtitle=r'$R_{\mathrm{CO}}\ (\mathrm{\AA})$',
       ytitle='CO$_2$-H$_2$O potential (cm$^{-1}$)',
       label='rigid',
       label2='relaxed',
       legendloc=4,
       save='dimer-diss.eps')