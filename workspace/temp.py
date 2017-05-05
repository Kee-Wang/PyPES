f = open('mode101_full.dat')
a = list()
b = 0
for line in f:
    for a in line.split():
        b = b + float(a)**2
print(b)
f.close()


f = open('mode101.dat')
a = list()
b = 0
for line in f:
    for a in line.split():
        b = b + float(a)**2
print(b)
f.close()

f = open('mode102.dat')
a = list()
b = 0
for line in f:
    for a in line.split():
        b = b + float(a)**2
print(b)
f.close()

f = open('mode103.txt')
a = list()
b = 0
for line in f:
    for a in line.split():
        b = b + float(a)**2
print(b)
f.close()