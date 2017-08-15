from pypes.plot import plot

plot().hist('hcl.zpe_data',col=[0], xtitle='m/s', ytitle='Count',title='HCl in 2 HCl Channels',bin=100, save='hcl_h_speed_zpe.eps')
plot().hist('hcl.zpe_data',col=[1], xtitle='J value', ytitle='Count',title='HCl in 2 HCl Channels',bin=100,save='hcl_h_j_zpe.eps')
plot().hist('hcl.zpe_j4',col=[0], xtitle='m/s', ytitle='Count',title='HCl in 2 HCl Channels, J=4', save='hcl_h_speed_j4_zpe.eps',bin = 50)
plot().hist('hcl.zpe_j6',col=[0], xtitle='m/s', ytitle='Count',title='HCl in 2 HCl Channels, J=6', save='hcl_h_speed_j6_zpe.eps',bin = 50)

plot().hist('hcl.data',col=[0], xtitle='m/s', ytitle='Count',title='HCl in 2 HCl Channels',bin=100, save='hcl_h_speed.eps')
plot().hist('hcl.data',col=[1], xtitle='J value', ytitle='Count',title='HCl in 2 HCl Channels',bin=100,save='hcl_h_j.eps')
plot().hist('hcl.j4',col=[0], xtitle='m/s', ytitle='Count',title='HCl in 2 HCl Channels, J=4', save='hcl_h_speed_j4.eps',bin = 50)
plot().hist('hcl.j6',col=[0], xtitle='m/s', ytitle='Count',title='HCl in 2 HCl Channels, J=6', save='hcl_h_speed_j6.eps',bin = 50)

