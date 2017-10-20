f = open('/home/kee/github/PyPES/workspace/que/nodetemp')
cores = list()
nodes = list()
for iter in range(26):
    nodes.append([-1])

'node[0],node[1]...node[25], total length 26, index is the node'

for line in f:
    if len(line) < 1:
        continue
    line = line.split()
    for word in line:
        word=word.split('+')
        if (word[0].find('node')) is 0:
            for core in word:
                core = core.split('node')[1].split('/')
                core[0] = int(core[0])
                core[1] = int(core[1])
                if int(core[0]) < 30: #Those are all xeon nodes
                    #print(core)
                    nodes[core[0]].append(core[1])

print('Xeon16 Node availability: (Node18 is broken)\n')
node_count= list()
node_count=[0,0,0]
core_count = 0
for iter in range(26):
    if iter is 0 or iter is 1 or iter is 18:
        continue
    ava = 17-len(nodes[iter])
    if ava is 16:
        print('Node{:>2d}: {:<2d} (Full node <<<)'.format(iter, ava))
        node_count[0] = node_count[0] + 1
    elif ava is 0:
        print('Node{:>2d}: {:<2d} (NA)'.format(iter, ava))
        node_count[1] = node_count[1] + 1

    else:
        print('Node{:>2d}: {:<2d}'.format(iter, ava))
        node_count[2] = node_count[2] + 1
        core_count = core_count + ava


print('Availability Summary: \n')
print('Full nodes: {:>2d}  (={:3d} cores)'.format(node_count[0],node_count[0]*16 ))
print('Partial nodes: {:>2d} (={:3d} cores)'.format(node_count[2], core_count))
print('Total available cores: {:>3d}'.format(node_count[0]*16+core_count))

