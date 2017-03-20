from pypes.configs import configs
a= configs('v2b_47558_abinitio.dat')#,first_n_configs=1000)
a1 = a.list()#first_n_configs=1000)
b = configs('v2b-pes-lr.dat')
b1 = b.list()
b.plotevsr_for_publication(a1,b1)
