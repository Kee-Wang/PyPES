from pypes.configsqct import configs_hwww_qct
import numpy as np
from pypes.configs import configs
from itertools import chain
import matplotlib.pyplot as plt
import sys


b = str(sys.argv)
path = b.split("'")[3]
print('path:', path)
bond_threshold = 2#  #Angstrom, The max covalent bond length for a molecule
rmax = 6.5 #Angstrom The max distance between two molecules to be interacted
#stop

def inlist(i,dis):
    """i is the number of the atom, dis is list contains all recent connect pairs
        return -1 means i is not in list, return other means which list is it in"""

    #print(dis)
    count = -1
    for sublist in dis:
        count = count + 1
        for subsublist in sublist:
            #print(i,subsublist)
            if i is subsublist:
                #result = True
                return count
    count = -1
    return count


#Use this to first sort out duplicates
#b = configs('./0918/0918.xyz', col=9,duplicate=True)
#b = configs('./1008/1008_all.xyz', col=9,duplicate=True)
#blist = b.list()

#b.write('1003_all_no_dup',blist)

# Then use this modified module to do calculation.
print(path,'test')
a = configs_hwww_qct(path, col=9)

print('\n Sorting Begin: \n')

#raise Exception('stopped')


temp = a.list()
al = temp[0]
step_size=2.5
au_s = 2.418884326509e-5 #unit ps

steps = temp[1] #*  step_size * au_s#All the steps associate with each trajectory
h_steps = list()
w_steps = list()


HWWW=list()
H_WWW=list()
HWW_W=list()
HW_WW=list()
H_W_WW=list()
HW_W_W=list()
H_W_W_W=list()
broken = list()



count = 0
count_step = -1
for config in al:
    count_step = count_step + 1
    broke = False

    bonds = a.bond_length(config)#,show=True,full=True)
    #for b in bonds:
    for atom in bonds:
        new_atom = list()
        for r in atom:
            if r < bond_threshold:
                new_atom.append(r)
            #print(r)
        if len(new_atom) > 1:
            continue
        else:

            print('Molecule broken! on molecule: ', count)
            #print(dis)
            broken.append(config)
            broke = True

            break#Break current for look
    if broke is True:
        continue

    count = count + 1
    dis = list()
    for i in range(7,12):
        if i is 10:
            continue
        for j in range(i+1,12):
            if j is 10:
                continue
            #print(i,j)
            r = a.distance(config,atom_A=i,atom_B=j)
            #print(i,j,r)
            if len(dis) < 1:

                if r > rmax:
                    dis.append([i])
                    dis.append([j])
                else:
                    dis.append([i, j])

            else: #Notic atom i will always in the list, either in the group or by itself, true forever after first step
                loci = inlist(i,dis)
                locj = inlist(j,dis) #Location of j in dis. If -1, not in list, otherwise location

                if r > rmax: #Two molecules are too far

                    if locj < 0: #atom i is already in the list. if j is not in the list, put it in new list
                        dis.append([j])

                else:#i and j are close and are grouped together
                    if locj < 0: #if j haven't showed up before, put j into i group
                        dis[loci].append(j)

                    else: #If j already existed, put j group and i group together
                        if loci is locj:  # If in the same group already, do nothing
                            continue
                        else:
                            for ele in dis[locj]:
                                dis[loci].append(ele)#Merge and then delete group j

                            dis.pop(locj) #Notice i<j
                        #print('t')

    lst = list(chain.from_iterable(dis))
    lst.sort()


    if lst != [7,8,9,11]:

        print('Error for list', count)

    else:

        if len(dis) is 1:
            print('Not dissociated! on molecule: ', count)
            HWWW.append(config)
        elif len(dis) is 2:
            for sublist in dis:
                if len(sublist) is 1:#the `1` could be H or W
                    if sublist[0] is 11:
                        #H_WWW.append(config)
                        t = steps[count_step] * step_size * au_s
                        h_steps.append(t)
                        config[1][0][0] = t

                        H_WWW.append(config)


                    else: #Else is `W`
                    #"This part is important because we have to put the fist water as the dissociated wter"
                        if sublist[0] is 8:#`W` is not the first `W` Bubble up that water
                            config = a.switch(config,neworder=[2,3,0,1,4,5,7,6,8,9,10])[0]

                        elif sublist[0] is 9:#`W` is not the first `W`
                            #print('before')
                            #a.prt(config)
                            config = a.switch(config,neworder=[4,5,0,1,2,3,8,6,7,9,10])[0]
                            #print('after')
                            #a.prt(config)

                        #print(sublist[0])

                        t = steps[count_step] * step_size * au_s
                        w_steps.append(t)
                        config[1][0][0] = t
                        HWW_W.append(config)
                        #print(count_step, steps[count_step])

                elif len(sublist) is 2 and (sublist[0] is 11 or sublist[1] is 11):
                    #print(dis)
                    #print(count,'HW + WW')
                    HW_WW.append(config)
        elif len(dis) is 3:
            for sublist in dis:
                if len(sublist) is 2:
                    if sum(sublist) >= (11 + 7): #Meaning 11 is inside 2-element list

                        #print(count,'HW + W + W')
                        HW_W_W.append(config)
                    else:
                        if sum(sublist) == (7+8):# water 9 is out
                            config = a.switch(config, neworder=[4, 5, 0, 1, 2, 3, 8, 6, 7, 9, 10])[0]
                        elif sum(sublist) == (7+9): #water 8 is out
                            config = a.switch(config, neworder=[2, 3, 0, 1, 4, 5, 7, 6, 8, 9, 10])[0]
                        #print(count,'WW + W + H')
                        H_W_WW.append(config)

        else: #This is all separated
            #print(dis)
            #print(count, 'H + W + W + W')
            H_W_W_W.append(config)






print('\n ---Soring finished! \n')
print('\n ---Result for dissociation:')
print('')
All = [HWWW, H_WWW, HWW_W, HW_WW, H_W_WW, HW_W_W, H_W_W_W]#, broken]
assign = ['HWWW','H + WWW','HWW + W','HW + WW','H + W + WW','HW + W + W','H + W + W + W']#, 'Blow-up!']
count = 0
num = list()
all_configs= list()

for lst in All:
    num.append(len(lst))
    for config in lst:
        all_configs.append(config)
    print('{:14s}: {:d}'.format(assign[count],num[count]))
    count = count + 1
#a.write('result_all.xyz',all_configs)
#print(a.broke, len(broken))
num.append((len(broken)+a.broke)) #Add configs that cannot be read in the initial file
print('{:14s}: {:d}'.format('Blow-up!',num[-1]))
print('')
print('{:14s}: {:d}'.format('Total', sum(num)))


#epoch1 = 1497045960
#epoch2 = 1497966025
#time = epoch2 - epoch1
#0817time = 1502126137

#total_step=sum(a.traj)
#evals = total_step*67 #+ len(a.traj)*11*11*4 #Per core
#print('Number of total evals: {:d}'.format(evals))
#print('Time per eval: {:f} ms'.format(time/evals*1000*((64-13))))


h_steps = np.array(h_steps)# * 2.5*au_s
w_steps = np.array(w_steps) #* 2.5*au_s
#print('Average step for H_WWW: {:f} pm {:f}'.format(h_steps.mean(), h_steps.std()))
print('Average time for H_WWW: {:f}  pm {:f} (ps)'.format(h_steps.mean() , h_steps.std()))
#print('Average step for W_WWW: {:f}'.format(w_steps.mean()))
print('Average time for W_WWW: {:f}  pm {:f} (ps)'.format(w_steps.mean() , w_steps.std() ))
#print(len(h_steps),len(w_steps))#traj = np.array(a.traj)


#plt.hist(traj,bins='auto')
#plt.xlabel('Time for dissociation (ps) \n (At least one of atom pairs with distance > 50 a.u)')
#plt.ylabel('Count')
#plt.tight_layout()
#fig = plt.gcf()
#plt.show()
#save = 'time_to_diss.eps'
#fig.savefig(save, format='eps', dpi=1200)
#print('Plot saved to {}.'.format(save))

assign1 = ['HWWW','H_WWW','HWW_W','HW_WW','H_W_WW','HW_W_W','H_W_W_W']

a.write('../2_final_condition/water/result_HWW_W.xyz', All[2])
a.write('../2_final_condition/hcl/result_H_WWW.xyz', All[1])


"""Problem: Can't recognize HW + WW correctly?"""


