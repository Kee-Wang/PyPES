from pypes.configs import configs

a = configs('all.xyz')
#a.plot(ref=True,binwidth=1)
b=a.sort()
a.write('sorted.xyz',b)