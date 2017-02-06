from pypes.plot import plot
a = plot()
a.line(file='rigid',
       file2='rlx',
       xtitle=r'$\mathrm{r}_{\mathrm{C--O}}\ (\mathrm{\AA})$',
       ytitle='energy (cm$^{-1}$)',
       label='rigid',
       label2='relaxed',
       legendloc=4,
       save='dimer-diss.eps')