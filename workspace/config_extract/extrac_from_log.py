f = open('result.log')
g = open('result.txt','w')
line_num = list()
text = list()
count = 0
for line in f:
    text.append(line)
    count = count + 1
    line = line.split()
    if len(line) is 0:
        continue

    if line[0] == 'Current':
        print(line)
        line_num.append(count)
print(line_num)
for num in range(line_num[-1]+1,line_num[-1]+15):
    print(text[num])
    g.write(text[num])

f.close()
g.close()
