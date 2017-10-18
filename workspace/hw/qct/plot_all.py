import matplotlib.pyplot as plt
import numpy as np
from pypes.plot import plot

def j_plot_back(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note):
    import pylab as p
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.set_title(ax_title)
    data1 = np.loadtxt(file1, unpack=True)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    plt.plot(data1[0],data1[2],'-',linewidth=linewidth,label='IR ON')
    ax.legend()
    ax2 = ax.twinx()
    data2 = np.loadtxt(file2, unpack=True)
    ax.text(note[1], note[2], note[0], fontsize=18)

    #y, binEdges=plt.hist(data2, bins=bin,fill=False,linewidth=3)
    y, binEdges = np.histogram(data2, bins=bin)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    #p.plot(bincenters, y, '-')
    plt.plot(bincenters,y,'r-')
    ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)

    a = plot()
    a.save(fig=fig, save=save)


def j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label,index1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.text(note[1], note[2], note[0], fontsize=18)
    #ax.set_title(ax_title)
    data1 = np.loadtxt(file1, unpack=True)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    plt.plot(data1[0],data1[index1],'k-',linewidth=linewidth,label='Exp')
    ax.legend()
    ax2 = ax.twinx()
    #data2 = list()
    ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)
    for i in range(len(file2)):
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2, bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=3,label = label[i])

    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)

"Overall speed distribution for HCl"
file1 = 'j=4.prn'
file2 = ['./result_hcl/result_H_WWW.s1_j4','./result_hcl/result_H_WWW.s2_j4','./result_hcl/result_H_WWW.s3_j4','./result_hcl/result_H_WWW.s4_j4']

save = 'j=4.eps'
ymax1 = 150
ymax2 = 90
xtitle='Velocity (m/s)'
ytitle1='Experiment signal intensity'
ytitle2='Trajectory count'
linewidth=1
bin = [20,15,15,10]
label = ['No Constraint','HCl only ZPE','Soft ZPE','Hard ZPE']
note = ['J=4',900,140]
#j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label, index1=2)

file1 = 'j=6.prn'
file2 = ['./result_hcl/result_H_WWW.s1_j6','./result_hcl/result_H_WWW.s2_j6','./result_hcl/result_H_WWW.s3_j6','./result_hcl/result_H_WWW.s4_j6']

save = 'j=6.eps'
ymax1 = 150
ymax2 = 60
xtitle='Velocity (m/s)'
ytitle1='Experiment signal intensity'
ytitle2='Trajectory count'
linewidth=1
bin = [20,10,10,10]
label = ['No Constraint','HCl only ZPE','Soft ZPE','Hard ZPE']
note = ['J=6',900,140]
#j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label,index1=1)


def j_distribution(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,label):
    import numpy as np
    """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    x = [3, 4, 5, 6, 7, 8]
    y = [1, 0.6, 0.4, 0.09, 0.034, 0.02]
    yerr = [0.07, 0.072, 0.1, 0.08, 0.02, 0.012]
    ax = ax.twinx()
    #ax.errorbar(x, y, yerr=yerr, fmt='ko', capsize=5, elinewidth=3, markeredgewidth=2,label='Exp')
    #plt.plot(x,y,'k-',linewidth=linewidth)


    ax2 = ax.twinx()
    ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)
    multi=[1,5,10,100]
    unique = [0,0,0,0]
    counts = [0,0,0,0]
    loc=list()

    max = list()
    for i in range(len(file2)):
        data2 = np.loadtxt(file2[i], unpack=True)
        x1 = np.array(data2[0])
        unique[i], counts[i] = np.unique(x1, return_counts=True)

        count = 0
        for j in unique[i]:

            if j == 3:
                loc.append(count)
            count = count + 1
        print(loc)



        #print(counts[0][loc[0]],counts[i][loc[i]])
        #counts[i] = counts[i] * counts[0][loc[0]] / counts[i][loc[i]]
        #print(counts[i])

        plt.plot(unique[i],counts[i]/counts[i][loc[i]],'-',linewidth=2,label = label[i])
        #plt.scatter(unique[i], counts[i]/counts[i][loc[i]])

    ax2.legend(loc=7)
    ax.errorbar(x, y, yerr=yerr, fmt='ko', capsize=5, elinewidth=3, markeredgewidth=2, label='Exp')
    ax.legend()

    a = plot()
    a.save(fig=fig, save=save)

    return None

file2 = ['./result_hcl/result_H_WWW.s1','./result_hcl/result_H_WWW.s2','./result_hcl/result_H_WWW.s3','./result_hcl/result_H_WWW.s4']

save = 'j_distribution.eps'
ymax1 = 1.2
ymax2 = 1.2
xtitle='J(HCl)'
ytitle1='Experiment relative population'
ytitle2='Theory (aligned at J=3)'
linewidth=1
bin = 10
label = ['No Constraint','HCl only ZPE','Soft ZPE','Hard ZPE']
#j_distribution(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,label)

def j_plot2(file1,file2,save,xlim,xtitle,ytitle1,ytitle2,label,index1,index2):
    fig = plt.figure(figsize=(6, 9))


    #ax.set_title(ax_title)
    #data1 = np.loadtxt(file1, unpack=True)
    #
    #ax.set_ylabel(ytitle1)
    #ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    #plt.plot(data1[0],data1[index1],'k-',linewidth=linewidth,label='Exp')
    #ax.legend()
    #ax2 = ax.twinx()
    #data2 = list()
    #ax2.set_ylim([0, ymax2])
    num = 410
    for i in range(len(file2)):
        num = num+1
        ax = fig.add_subplot(num)
        ax.set_ylabel(ytitle2)
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[index1], bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=3)
        ax.set_xlim(xlim)
        ymax = ax.get_ylim()
        ax.text(index2, ymax[1]*0.8, label[i], fontsize=15)
        ax.legend()




    #ax.legend(loc=7)
    ax.set_xlabel(xtitle)
    a = plot()
    a.save(fig=fig, save=save)

"Overall speed distribution for HCl"
file1 = 'j=4.prn'
file2 = ['./result_hcl/result_H_WWW.s1','./result_hcl/result_H_WWW.s2','./result_hcl/result_H_WWW.s3','./result_hcl/result_H_WWW.s4']

save = 'HCl_internal.eps'
xlim=[0,7000]
xtitle='Velocity (m/s)'
ytitle1='Experiment signal intensity'
ytitle2='Trajectory count'
linewidth=1
bin = [100, 100, 100, 20]
xtitle='HCl Internal Energy (cm$^{-1}$)'
label = ['No Constraint','Hard ZPE on HCl','Soft ZPE on both','Hard ZPE on both']
#j_plot2(file1,file2,save,xlim,xtitle,ytitle1,ytitle2,label, index1=2, index2 = 4000)

xlim=[10000,20000]
bin = [100, 100, 100, 50]
save = 'WWW_internal.eps'
xtitle='(H$_2$O)$_3$ Internal Energy (cm$^{-1}$)'
#j_plot2(file1,file2,save,xlim,xtitle,ytitle1,ytitle2,label, index1=3, index2=10150)



a =  np.loadtxt('./result_water/result_HWW_W.s1', unpack=True)
def j_w_plot(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label):
    plt.rc('xtick', labelsize=20)  # fontsize of the axes title
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    #ax.text(note[1], note[2], note[0], fontsize=18)
    #ax.set_title(ax_title)
    #data1 = np.loadtxt(file1, unpack=True)
    ax2.set_xlabel(xtitle)
    #ax.set_ylabel(ytitle1)
    #ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    #plt.plot(data1[0],data1[index1],'k-',linewidth=linewidth,label='Exp')
    #ax.legend()
    #ax2 = ax.twinx()
    #data2 = list()
    #ax2.set_ylim([0, ymax2])
    bin = [50,50,50,20]
    ax2.set_ylabel(ytitle2)

    for i in range(len(file2)):
        print(file2[i])
        f = open('water_j_distribution_s{:d}.txt'.format(i+1), 'w')
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[1], bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=3,label = label[i])
        count = 0
        f.write('#Bin center    #Number of count \n')
        for num in y:
            f.write('{:f},{:d} \n'.format(binEdges[count],num))
            count = count + 1


        f.close()

    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)
file2 = ['./result_water/result_HWW_W.s1','./result_water/result_HWW_W.s2','./result_water/result_HWW_W.s3','./result_water/result_HWW_W.s4']
save = 'water_j_distribution.eps'
label = ['No Constraint','Hard ZPE on monomer','Soft ZPE on both','Hard ZPE on both']
xtitle = 'J(water monomer)'
#j_w_plot(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label)


def w_erot_plot(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label):
    plt.rc('xtick', labelsize=20)  # fontsize of the axes title
    fig = plt.figure()
    ax2 = fig.add_subplot(111)

    ax2.set_xlabel(xtitle)

    bin = [50,50,50,20]
    ax2.set_ylabel(ytitle2)
    ax2.set_xlim([0,2000])

    for i in range(len(file2)):
        #f = open('water_j_distribution_s{:d}.txt'.format(i+1), 'w')
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[1]*219474.63, bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=3,label = label[i])
        count = 0

        #f.write('#Bin center    #Number of count \n')
        #for num in y:
        #    f.write('{:f},{:d} \n'.format(binEdges[count],num))
        #    count = count + 1


        #f.close()
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)
file2 = ['./result_water/result_HWW_W.s1','./result_water/result_HWW_W.s2','./result_water/result_HWW_W.s3','./result_water/result_HWW_W.s4']
save = './result_water/water_erot_distribution.eps'
label = ['No Constraint','Hard ZPE on monomer','Soft ZPE on both','Hard ZPE on both']
xtitle = 'Water monomer rotational energy (cm$^{-1}$)'
w_erot_plot(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label)



#
# file1 = 'j=4.prn'
# file2 = 'result_H_WWW.j4_j4'
# save = 'No'# 'j=4_c4.eps'
# ymax1 = ymax1
# ymax2 = 50
# xtitle='Speed (m/s)'
# ytitle1='Experiment signal intensity '
# ytitle2='Count (with soft ZPE constraint)'
# linewidth=1
# bin = 10
# ax_title='HCl channel, J(HCl)=4'
# #j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2)
#
# file1 = 'j=4.prn'
# file2 = 'result_H_WWW.j3_j4'
# save = 'j=4_c3.eps'
# ymax1 = ymax1
# ymax2 = 50
# xtitle='Speed (m/s)'
# ytitle1='Experiment signal intensity '
# ytitle2='Count (with hard ZPE constraint)'
# linewidth=1
# bin = 10
# ax_title='HCl channel, J(HCl)=4'
# #j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2)
#
#
#
#
# file1 = 'water_channel/waterj321.prn'
# file2 = 'water_channel/result_HWW_W.j4_j321'
# save = 'water_channel/j4_j321.eps'
# ymax1 = 125
# ymax2 = 20
# xtitle='Speed (m/s)'
# ytitle1='Experiment signal intensity '
# ytitle2='Count (with soft ZPE constraint)'
# linewidth=1
# bin = 10
# ax_title='Water channel, $J_{K_aK_c}=3_{21}$'
# #j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2)
#
#
# file1 = 'water_channel/waterj321.prn'
# file2 = 'water_channel/result_HWW_W.j1_j321'
# save = 'water_channel/j1_j321.eps'
# ymax1 = 125
# ymax2 = 20
# xtitle='Speed (m/s)'
# ytitle1='Experiment signal intensity '
# ytitle2='Count (with no ZPE constraint)'
# linewidth=1
# bin = 10
# ax_title='Water channel, $J_{K_aK_c}=3_{21}$'
# #j_plot(file1,file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2)
#
#
#

