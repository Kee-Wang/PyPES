from pypes.plot import plot
a = plot()
a.line(file='rigid',
       file2='rlx',
       xtitle=r'$R_{\mathrm{CO}}\ (\mathrm{\AA})$',
       ytitle='CO$_2$-H$_2$O potential (cm$^{-1}$)',
       label='rigid',
       label2='relaxed',
       legendloc=4,
       save='dimer-diss.eps')