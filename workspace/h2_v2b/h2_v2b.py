from pypes.configs import configs
a = configs('dimer.abE')
a1 = a.list()
b = configs('monomerA.abE')
b1 = b.list()
c = configs('monomerB.abE')
c1 = c.list()
d = a.v2b(a1,b1,c1)
a.write('h2_v2b.abE',d)