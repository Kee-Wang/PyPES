from pypes.configs import configs
a = configs('minima',xlim=[-1,-4000/219474.63])#,first_n_configs=10000)
a.info()
a.plot(binwidth=0.5,ref = True,ylim=[0,100])