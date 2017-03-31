from pypes.configs import configs
a = configs('v2b.abE')
a1 = a.list()
a2 = a.sort_unique(a1)
a.write('unique_v2b.xyz',a2)