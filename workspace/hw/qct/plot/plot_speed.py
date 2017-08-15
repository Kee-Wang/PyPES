from pypes.plot import plot

plot().hist('water.zpe_data',col=[0], xtitle='m/s', ytitle='Count',title='water in 2 water Channels, correct ZPE',bin=100, save='water_h_speed_zpe.eps')
plot().hist('water.zpe_data',col=[1], xtitle='Rotational Energy (cm$^{-1}$)', ytitle='Count',title='water in 2 water Channels, correct ZPE',bin=100,save='water_h_j_zpe.eps')
plot().hist('water.zpe_j221',col=[0], xtitle='m/s', ytitle='Count',title='water in 2 water Channels, JKaKc=221, correct ZPE', save='water_h_speed_j221_zpe.eps',bin = 50)
plot().hist('water.zpe_j321',col=[0], xtitle='m/s', ytitle='Count',title='water in 2 water Channels, JKaKc=321, correct ZPE', save='water_h_speed_j321_zpe.eps',bin = 50)

plot().hist('water.data',col=[0], xtitle='m/s', ytitle='Count',title='water in 2 water Channels',bin=100, save='water_h_speed.eps')
plot().hist('water.data',col=[1], xtitle='Rotational Energy (cm$^{-1}$)', ytitle='Count',title='water in 2 water Channels',bin=100,save='water_h_j.eps')
plot().hist('water.j221',col=[0], xtitle='m/s', ytitle='Count',title='water in 2 water Channels, JKaKc=221', save='water_h_speed_j221.eps',bin = 50)
plot().hist('water.j321',col=[0], xtitle='m/s', ytitle='Count',title='water in 2 water Channels, JKaKc=321', save='water_h_speed_j321.eps',bin = 50)
