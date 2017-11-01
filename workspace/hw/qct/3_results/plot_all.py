import matplotlib.pyplot as plt
import numpy as np
from pypes.plot import plot

ph = '../2_final_condition/hcl/result_H_WWW'
ph2 = '../2_final_condition/hcl/result_H_WWW'
file_hclj4 = [ph+'_s1_j4.txt',ph+'_s2_j4.txt',ph+'_s3_j4.txt',ph+'_s4_j4.txt']
file_hclj6 = [ph+'_s1_j6.txt',ph+'_s2_j6.txt',ph+'_s3_j6.txt',ph+'_s4_j6.txt']
file_hcl = [ph+'_s1.txt',ph+'_s2.txt',ph+'_s3.txt',ph+'_s4.txt']
file_expj4 = './expdata/hcl_speed_j=4.prn'
file_expj6 = './expdata/hcl_speed_j=6.prn'
ph = '../2_final_condition/water/result_HWW_W'
file_w = [ph+'_s1.txt',ph+'_s2.txt',ph+'_s3.txt',ph+'_s4.txt']

def hcl_j_speed_plot(file1,file2,save,ymax2,bin,index):
    label = ['No Constraint', 'HCl only ZPE', 'Soft ZPE', 'Hard ZPE']
    xtitle = 'Velocity (m/s)'
    ytitle1 = 'Experiment signal intensity'
    ytitle2 = 'Trajectory count'
    ymax1 = 150
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.text(note[1], note[2], note[0], fontsize=18)
    #ax.set_title(ax_title)
    data1 = np.loadtxt(file1, unpack=True)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    plt.plot(data1[0],data1[index],'k-',linewidth=1,label='Exp')
    ax.legend()
    ax2 = ax.twinx()
    #data2 = list()
    ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)
    for i in range(len(file2)):
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[2], bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=3,label = label[i])

    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)

"J selected speed distribution for HCl"

save = './result_hcl/j4speed.eps'
ymax2 = 200
bin = [20,15,15,10]
note = ['J=4',900,120]
hcl_j_speed_plot(file_expj4,file_hclj4,save,ymax2,bin,index=2)

save = './result_hcl/j6speed.eps'
ymax2 = 120
bin = [20,10,10,10]
label = ['No Constraint','HCl only ZPE','Soft ZPE','Hard ZPE']
note = ['J=6',900,120]
hcl_j_speed_plot(file_expj6,file_hclj6,save,ymax2, bin,index=1)


def hcl_j_distribution(file2,save):
    import numpy as np
    """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ymax1 = 1.2
    ymax2 = 1.2
    xtitle = 'HCl rotational energy (cm$^{-1})$'
    ytitle1 = 'Experiment relative population'
    ytitle2 = 'Theory (normalized at J=3)'
    linewidth = 1
    bin = 10
    label = ['No Constraint', 'HCl only ZPE', 'Soft ZPE', 'Hard ZPE']
    Be_HCl = 10.59341

    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    x = list()
    for num in [3, 4, 5, 6, 7, 8]:
        x.append(num*(num+1) * Be_HCl)
    #x = [3, 4, 5, 6, 7, 8] * Be_HCl
    y = [1, 0.6, 0.4, 0.09, 0.034, 0.02]
    yerr = [0.07, 0.072, 0.1, 0.08, 0.02, 0.012]

    ax.errorbar(x, y, yerr=yerr, fmt='ko', capsize=5, elinewidth=3, markeredgewidth=2, label='Exp')
    ax.legend()

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
        #print(data2[0])
        x1 = np.array(data2[0])
        unique[i], counts[i] = np.unique(x1, return_counts=True)
        count = 0
        #for j in unique[i]:
         #   if j == 3:
          #      loc.append(count)
         #   count = count + 1
        print(unique[i],counts[i])
        plt.plot(unique[i],counts[i]/counts[i][3],'-',linewidth=2,label = label[i])
    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)
    return None



hcl_j_distribution(file_hcl,save = './result_hcl/hcl_j_distribution.eps')



def w_erot_plot(file2,save):
    ytitle2 = 'Trajectory count'
    xtitle = 'Water monomer rotational energy (cm$^{-1}$)'
    plt.rc('xtick', labelsize=12)  # fontsize of the axes title
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    label = ['No Constraint', 'Hard ZPE on monomer', 'Soft ZPE on both', 'Hard ZPE on both']

    ax2.set_xlabel(xtitle)

    bin = [50,50,50,20]
    ax2.set_ylabel(ytitle2)
    ax2.set_xlim([0,2000])

    for i in range(len(file2)):
        #f = open('water_j_distribution_s{:d}.txt'.format(i+1), 'w')
        data2 = np.loadtxt(file2[i], unpack=True)

        y, binEdges = np.histogram(data2[0], bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=3,label = label[i])
        #count = 0
        #f.write('#Bin center    #Number of count \n')
        #for num in y:
        #    f.write('{:f},{:d} \n'.format(binEdges[count],num))
        #    count = count + 1
        #f.close()
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)


w_erot_plot(file_w,save = './result_water/water_erot_distribution.eps')





"Not used"


def hcl_erot_plot(file2,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label):
    plt.rc('xtick', labelsize=12)  # fontsize of the axes title
    fig = plt.figure()
    ax2 = fig.add_subplot(111)

    ax2.set_xlabel(xtitle)

    bin = [50,50,50,20]
    ax2.set_ylabel(ytitle2)
    ax2.set_xlim([0,2000])

    for i in range(len(file2)):
        #f = open('water_j_distribution_s{:d}.txt'.format(i+1), 'w')
        data2 = np.loadtxt(file2[i], unpack=True)

        y, binEdges = np.histogram(data2[0], bins=bin[i])
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
save = './result_hcl/hcl_erot_distribution.eps'
label = ['No Constraint', 'Hard ZPE on monomer', 'Soft ZPE on both', 'Hard ZPE on both']
xtitle = 'HCl monomer rotational energy (cm$^{-1}$)'
#hcl_erot_plot(file_hcl, save, ymax1, ymax2, xtitle, ytitle1, ytitle2, note, label)


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
        y, binEdges = np.histogram(data2[0], bins=bin[i])
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
#file2 = #['./result_water/result_HWW_W.s1','./result_water/result_HWW_W.s2','./result_water/result_HWW_W.s3','./result_water/result_HWW_W.s4']
save = 'water_j_distribution.eps'
label = ['No Constraint','Hard ZPE on monomer','Soft ZPE on both','Hard ZPE on both']
xtitle = 'J(water monomer)'
#j_w_plot(file_w,save,ymax1,ymax2,xtitle,ytitle1,ytitle2,note,label)


def joverallspped_plot2(file_hcl,save,xlim,xtitle,ytitle1,ytitle2,label,index1,index2):
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
#joverallspped_plot2(file1,file2,save,xlim,xtitle,ytitle1,ytitle2,label, index1=2, index2 = 4000)
xlim=[10000,20000]
bin = [100, 100, 100, 50]
save = 'WWW_internal.eps'
xtitle='(H$_2$O)$_3$ Internal Energy (cm$^{-1}$)'
#j_plot2(file1,file2,save,xlim,xtitle,ytitle1,ytitle2,label, index1=3, index2=10150)

