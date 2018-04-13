import matplotlib.pyplot as plt
import numpy as np
from pypes.plot import plot
from matplotlib import rc
from pylab import rcParams
rcParams['figure.figsize'] = 3.37, 2.8
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)
#plt.rcParams.update({'font.family': 'Times'})
plt.rc('font', family='serif')
plt.rc('font', serif='Times New Roman')

ph = '../2_final_condition/hcl/result_H_WWW'
file_hcl = [ph+'_s1.txt',ph+'_s2.txt',ph+'_s3.txt',ph+'_s4.txt']
file_hclj4 = [ph+'_s1_j4.txt',ph+'_s2_j4.txt',ph+'_s3_j4.txt',ph+'_s4_j4.txt']
file_hclj6 = [ph+'_s1_j6.txt',ph+'_s2_j6.txt',ph+'_s3_j6.txt',ph+'_s4_j6.txt']

ph2 = '../2_final_condition/water/result_HWW_W'
file_w = [ph2+'_s1.txt',ph2+'_s2.txt',ph2+'_s3.txt',ph2+'_s4.txt']
file_wj221 = [ph2+'_s1_j221.txt',ph2+'_s2_j221.txt',ph2+'_s3_j221.txt',ph2+'_s4_j221.txt']
file_wj321 = [ph2+'_s1_j321.txt',ph2+'_s2_j321.txt',ph2+'_s3_j321.txt',ph2+'_s4_j321.txt']


file_expj4 = './expdata/hcl_speed_j=4.prn'
file_expj6 = './expdata/hcl_speed_j=6.prn'

file_expj221 = './expdata/water_speed_j221.prn'
file_expj321 = './expdata/water_speed_j321.prn'


f = open('../3_results/excel/HCl_bin.txt','w')
f2 = open('../3_results/excel/water_bin.txt','w')


def hcl_j_speed_plot(file1,file2,save,ymax2,bin,index):
    from scipy.stats import gaussian_kde
    import numpy as np
    SMALL_SIZE = 8
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 10
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


    label = ['No Constraint', 'HCl only ZPE', 'Soft ZPE', 'Hard ZPE']
    xtitle = 'Speed (m/s)\n'+note[0]
    ytitle1 = 'Experimental Intensity (a.u.)'
    ytitle2 = 'Theory (Count)'# \n(normalized at '+note[2] +')'
    ymax1 = 160
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.text(note[1], note[2], note[0], fontsize=18)
    #ax.set_title(ax_title)
    data1 = np.loadtxt(file1, unpack=True)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    plt.plot(data1[0],data1[index],'k-',linewidth=1,label='Exp')
    #ax.legend()
    ax2 = ax.twinx()
    #data2 = list()
    #ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)
    ax2.plot(np.nan, 'k-', label='Exp')  # Make an agent in ax


    for i in range(len(file2)):
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[2], bins=bin[i])


        data = data2[2]
        density = gaussian_kde(data)
        xs = np.linspace(0, 2000)
        density.covariance_factor = lambda: .25
        density._compute_covariance()
        #plt.plot(xs, density(xs))
        #plt.show()


        print('The number of configs: {:d}'.format(len(data2[2])))
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        #plt.plot(xs,density(xs),'-',linewidth=1.5,label = label[i])
        #plt.plot(bincenters,y/note[1],'-',linewidth=1.5,label = label[i])
        plt.plot(bincenters, y , '-', linewidth=1.5, label=label[i])
        total_count = sum(y)
        m_hcl = 0.03646  # (kg/mol)
        weight_Et = 0
        weight_sum = 0
        deno = 0
        for j in range(len(y)):
            f.write("{:5d}   {:6.2f}\n".format(y[j],bincenters[j]))
            weight_sum = weight_sum + y[j]*bincenters[j]
            #weight_Et = weight_Et + 0.5 * m_hcl * y[j] * bincenters[j] * bincenters[j]
            weight_Et = weight_Et + 0.5 * m_hcl * y[j] * bincenters[j]
            deno = deno +   y[j] / bincenters[j]
        f.write("\n")
        print(weight_Et,deno)
        speed_avg = weight_sum / total_count #(m/s)
        #Et_avg = weight_Et / total_count /1000 * 83.593473 #(cm-1)
        Et_avg = weight_Et / deno / 1000 * 83.593473  # (cm-1)
        print('average hcl speed is: {:f} '.format(speed_avg))
        print('average hcl Et is: {:f} \n '.format(Et_avg))

    ax2.legend()
    a = plot()
    a.save(fig=fig, save=save)

"J selected speed distribution for HCl"

save = './result_hcl/FIG_S1a_j4speed.eps'
ymax2 = 1.2
bin = [50,50,50,50]
note = ['(a)',171,'311 m/s']
#hcl_j_speed_plot(file_expj4,file_hclj4,save,ymax2,bin,index=2)

save = './result_hcl/FIG_S1b_j6speed.eps'
ymax2 = 1.2
bin = [40,40,40,40]
label = ['No Constraint','HCl only ZPE','Soft ZPE','Hard ZPE']
note = ['(b)',96,'260 m/s']
#hcl_j_speed_plot(file_expj6,file_hclj6,save,ymax2, bin,index=1)


def w_speed_plot(file1,file2,save,ymax1,ymax2,bin,index):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 10
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    label = ['No Constraint', 'Hard ZPE on water', 'Soft ZPE', 'Hard ZPE']
    xtitle = 'Speed (m/s)\n'+note[0]
    ytitle1 = 'Experiment Intensity (a.u.)'
    ytitle2 = 'Trajectory count'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.text(note[1], note[2], note[0], fontsize=18)
    #ax.set_title(ax_title)
    data1 = np.loadtxt(file1, unpack=True)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    #plt.plot(data1[0],data1[1],'r-',linewidth=linewidth,label='IR ON')
    plt.plot(data1[0],data1[index],'k-',linewidth=1,label='Exp')
    #ax.legend()
    ax2 = ax.twinx()
    #data2 = list()
    ax2.plot(np.nan, 'k-', label='Exp')
    ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)
    for i in range(len(file2)):
        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[2], bins=bin[i])
        print('The number of configs: {:d}'.format(len(data2[2])))
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y,'-',linewidth=1,label = label[i])
        total_count = sum(y)
        m_hcl = 0.03646  # (kg/mol)
        weight_Et = 0
        weight_sum = 0
        for j in range(len(y)):
            f2.write("{:5d}   {:6.2f}\n".format(y[j],bincenters[j]))
            weight_sum = weight_sum + y[j]*bincenters[j]
            weight_Et = weight_Et + 0.5 * m_hcl * y[j] * bincenters[j] * bincenters[j]
        f2.write("\n")
        speed_avg = weight_sum / total_count #(m/s)
        Et_avg = weight_Et / total_count /1000 * 83.593473 #(cm-1)
        print('average hcl speed is: {:f} '.format(speed_avg))
        print('average hcl Et is: {:f} \n '.format(Et_avg))
    ax2.legend(loc=7)
    a = plot()
    a.save(fig=fig, save=save)
"J selected speed distribution for water"

save = './result_water/j221speed.eps'
ymax1 = 250
ymax2 = 40
bin = [30,30,30,30]
note = ['(a)',800,140]
#w_speed_plot(file_expj221,file_wj221,save,ymax1,ymax2,bin,index=1)

save = './result_water/j321speed.eps'
ymax1 = 140
ymax2 = 70
bin = [30,30,30,30]
note = ['(b)',800,140]
#w_speed_plot(file_expj321,file_wj321,save,ymax1, ymax2,bin,index=1)





def hcl_j_distribution(file2,save):
    import numpy as np
    """This is to plot simple scatter plot for n columns. order in the first column will be take as x, all other columns are taken as y"""
    SMALL_SIZE = 8
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 10
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ymax1 = 1.2
    ymax2 = 1.2
    xtitle = 'Rotational Energy (wavenumbers)'
    ytitle1 = 'Experimental Relative Population'
    ytitle2 = 'Theory \n(normalized at J\"=3)'
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

    ax.errorbar(x, y, yerr=yerr, fmt='ko', capsize=3,elinewidth=0.7, markeredgewidth=0.5, label='Exp',markersize=3)
   # ax.legend()

    ax2 = ax.twinx()
    ax2.errorbar(np.nan,np.nan,yerr=np.nan, fmt='ko',capsize=3,elinewidth=0.7, markeredgewidth=0.5, label='Exp',markersize=3)
    ax2.set_ylim([0, ymax2])
    ax2.set_ylabel(ytitle2)
    multi=[1,5,10,100]
    unique = [0,0,0,0]
    counts = [0,0,0,0]
    loc=list()
    markers = ['o','s','*','v']
    max = list()
    for i in range(len(file2)):
        data2 = np.loadtxt(file2[i], unpack=True)
        #print(data2[0])
        #x1 = np.array(data2[0])
        x1 = np.array(data2[1])
        unique[i], counts[i] = np.unique(x1, return_counts=True)
        count = 0


        #print(unique[i],counts[i])
        plt.plot((unique[i]*(unique[i]+1)*Be_HCl),counts[i]/counts[i][3],'-',marker=markers[i],linewidth=0.7,label = label[i],markersize=4, fillstyle='none',mew=0.5)
        print(counts[i])
    ax2.legend()
    a = plot()
    a.save(fig=fig, save=save)
    return None



#hcl_j_distribution(file_hcl,save = './result_hcl/FIG4_hcl_j_distribution.eps')



def w_erot_plot(file2,save):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 10
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    ytitle2 = 'Theory'# \n(normalized at 145 wavenumbers)'
    xtitle = 'Rotational Energy (wavenumbers)'
    plt.rc('xtick', labelsize=12)  # fontsize of the axes title
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    label = ['No Constraint', 'Hard ZPE on monomer', 'Soft ZPE on both', 'Hard ZPE on both']

    ax2.set_xlabel(xtitle)

    bin = [50,50,50,20]
    ax2.set_ylabel(ytitle2)
    ax2.set_xlim([0,3000])

    for i in range(len(file2)):
        #f = open('water_j_distribution_s{:d}.txt'.format(i+1), 'w')
        data2 = np.loadtxt(file2[i], unpack=True)

        y, binEdges = np.histogram(data2[0], bins=bin[i])
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincenters,y/2119,'-',linewidth=1.5,label = label[i])
        for j in range(len(y)):
            f2.write("{:5d}   {:6.2f}\n".format(y[j],bincenters[j]))
        f2.write("\n")

    ax2.legend()
    a = plot()
    a.save(fig=fig, save=save)


#w_erot_plot(file_w,save = './result_water/FIG5_water_erot_distribution.eps')



def hcl_j2_speed_plot(file1,file2,save,ymax2,bin,index):
    from scipy.stats import gaussian_kde
    import numpy as np
    SMALL_SIZE = 8
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 10
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


    label = ['No Constraint', 'HCl only ZPE', 'Soft ZPE', 'Hard ZPE']
    xtitle = 'Speed (m/s)'#\n'+note[0]
    ytitle1 = 'Intensity (a.u.)'
    ytitle2 = 'Theory'# \n(normalized at '+note[2] +')'
    ymax1 = 1.2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.text(note[1], note[2], note[0], fontsize=18)
    #ax.set_title(ax_title)
    data1 = np.loadtxt(file1, unpack=True)


    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle1)
    ax.set_ylim([0, ymax1])
    ax.set_xlim([0, 1400])
    plt.plot(data1[0],data1[index]/note[1],'k-',linewidth=0.5,label='Exp')

    for i in range(len(file2)):
        if i is not 2:
            continue

        data2 = np.loadtxt(file2[i], unpack=True)
        y, binEdges = np.histogram(data2[2], bins=bin[i])


        # data = data2[2]
        # density = gaussian_kde(data)
        # xs = np.linspace(0, 2000)
        # density.covariance_factor = lambda: .25
        # density._compute_covariance()
        # #plt.plot(xs, density(xs))
        # #plt.show()


        print('The number of configs: {:d}'.format(len(data2[2])))
        bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
        #plt.plot(xs,density(xs),'-',linewidth=1.5,label = label[i])
        for j in range(len(y)):
            print('binceters = {:f}, y = {:d}'.format(bincenters[j], y[j]))

        plt.plot(bincenters,y/note[2],'-',linewidth=1.5,label = label[i])
        total_count = sum(y)
        m_hcl = 0.03646  # (kg/mol)
        weight_Et = 0
        weight_sum = 0
        for j in range(len(y)):
            f.write("{:5d}   {:6.2f}\n".format(y[j],bincenters[j]))
            weight_sum = weight_sum + y[j]*bincenters[j]
            weight_Et = weight_Et + 0.5 * m_hcl * y[j] * bincenters[j] * bincenters[j]
        f.write("\n")
        speed_avg = weight_sum / total_count #(m/s)
        Et_avg = weight_Et / total_count /1000 * 83.593473 #(cm-1)
        print('average hcl speed is: {:f} '.format(speed_avg))
        print('average hcl Et is: {:f} \n '.format(Et_avg))


    ax.legend()
    a = plot()
    a.save(fig=fig, save=save)

"J selected speed distribution for HCl"

save = './result_hcl/temp_FIG3a_j4speed_norm.eps'
ymax2 = 1.2
b = 50
bin = [b, b, b, b]
note = ['(a)',150.41193, 118]
hcl_j2_speed_plot(file_expj4,file_hclj4,save,ymax2,bin,index=2)

save = './result_hcl/temp_FIG3b_j6speed_norm.eps'
ymax2 = 1.2
bin = [40,40,40,40]
label = ['No Constraint','HCl only ZPE','Soft ZPE','Hard ZPE']
note = ['(b)',167.52983,73]
hcl_j2_speed_plot(file_expj6,file_hclj6,save,ymax2, bin,index=1)





























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

