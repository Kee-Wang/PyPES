from pypes.configs import configs
import numpy as np
import copy
a= configs('v2b_47558_abinitio.dat')
#b = configs('v2b-pes-lr.dat')
b = configs('v2b-pes-lr-almost-unique.dat')
a1 = a.list()
b1 = b.list()

#a1new = list()
#for config in a1:
#    dis = a.distance(config,3,6)
#    if  dis>= 6:
#        a1new.append(config)

#count = 0
#b1same = list()
#count1 = 0
#count2 = 0
#bindex = list()
#for config in b1:

#    for a1list in a1new:
#        count2 = count2 + 1
#        if (abs(config[1][0][0] - a1list[1][0][0])<0.000001 and abs(a.distance(config,3,6)-a.distance(a1list,3,6))<0.001):
#            bindex.append(count1)

#    count1 = count1 + 1

#count1 = 0
#b1new= list()
#for num in bindex[::-1]:
    #b1new.append(b1[num])
#    del b1[num]
#a.write('v2b-pes-lr-almost-unique.dat',b1)
count2 = 0
d = a1+b1
for config in d:
    if config[1][0][0]*219474.68 >=3100:
        count2 = count2 + 1
print((count2/len(d)))
a.plotv2b(d)
