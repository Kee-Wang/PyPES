from pypes.plot import plot

#plot().line('hwww_h2o_channel.cut',xmin=0,xmax=100,ymin=-9000,ymax=-100,ytitle=r'V (cm$^{-1}$)',xtitle = r'r$_\mathrm{OO}$ ($\mathrm{\AA}$)',save = 'hwww_h2o_channel.eps')

#plot().line('hwww_rlx_channel.cut',xmin=5,xmax=30,ymin=-5114,ymax=-5090,ytitle=r'V (cm$^{-1}$)',xtitle = r'r$_\mathrm{OO}$ ($\mathrm{\AA}$)',save = 'hwww_rlx_channel_detail.eps')

#plot().line('hwww_rlx_channel.cut',ytitle=r'V (cm$^{-1}$)',xtitle = r'r$_\mathrm{OO}$ ($\mathrm{\AA}$)',save = 'hwww_rlx_channel.eps')



plot().line('hwww_h2o_channel.cut','hwww_rlx_channel.cut',xmin=0,xmax=80,ymin=-9000,ymax=-2000,label=r'HWWW(GM)$\rightarrow$ HWW + W',label2=r'HWWW$\rightarrow$ HWW(GM) + W(GM)',legendloc=8,ytitle=r'V (cm$^{-1}$)',xtitle = r'r$_\mathrm{OO}$ ($\mathrm{\AA}$)',save = 'hwww_compare.eps')

#plot().line('hwww_h2o_channel.cut','hwww_rlx_channel.cut', ytitle=r'V (cm$^{-1}$)',xtitle = r'r$_\mathrm{OO}$ ($\mathrm{\AA}$)',save = 'hwww_compare.eps')


