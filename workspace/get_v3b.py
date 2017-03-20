#!/usr/bin/env python
file = input('Input filename: ')
#file = 'mp.out'
outfile = file.split('.')[0]
outfile = outfile + '.xyz'
f = open(file)
out = open(outfile,'w')
natm = 9
count = 0
line1 = 30
line2 = 39
config = list()
for line in f:
    count = count + 1
    newline = line.split()
    if len(newline)<1:
        continue
    for i in range(line1,line2):
        if count == i:
            config.append(line)
    if newline[0] == 'V3B':
        e = float(newline[2])
out.write( """{:<2d}
{:<14.8f}
"""
.format(natm,e) )
for line in config:
    out.write(str(line))

out.close()