import numpy as np
file = '1023_local.eng'
f = open(file)
energy = list()

count = 0
for e in f:
    count = count + 1
    if count > 100000:
        break
    try:

        s = float(e)
        print(count)
        if s > 21240:
            continue

        energy.append(float(e))

        #print(s)
    except:
        pass
energy = np.array(energy)
#print(energy,s)
print(energy.mean(), np.std(energy))


#data = np.loadtxt(, unpack=True)
#print(data)