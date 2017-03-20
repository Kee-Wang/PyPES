from pypes.configs import configs
a= configs('v2b-pes-lr.dat')#,first_n_configs=1000)
a1 = a.list()
count = 0
for config in a1:
    dis = a.distance(config,3,6)
    if  dis>= 21:
        count = count + 1

        print(count,dis)