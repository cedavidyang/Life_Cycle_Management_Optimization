import numpy as np
import os
import sys
#import matplotlib.pyplot as plt
#from mpldatacursor import datacursor
#plt.rc('font', family='serif', size=12)
#plt.rc('text', usetex=True)
import matplotlib as mpl
#mpl.use("pgf")
pgf_with_custom_preamble = {
    "figure.figsize": [3.54, 2.655],
    #"figure.subplot.bottom": 0.14,
    #"figure.subplot.top": 0.93,
    "font.family": "serif", # use serif/main font for text elements
    "font.size": 9, # use font size
    "text.usetex": True,    # use inline math for ticks
    #'text.latex.unicode': True,
    #"pgf.rcfonts": False,   # don't setup fonts from rc parameters
    #"pgf.preamble": [
        #r'\usepackage{fontspec}',         # load additional packages
        #]
}
mpl.rcParams.update(pgf_with_custom_preamble)
import matplotlib.pyplot as plt
from mpldatacursor import datacursor
annotate_text = u'{label}'

def rate2suffix(icorr_mean_list):
    suffix = ''
    for icorr in icorr_mean_list:
        if icorr == 0.5:
            suffix += 'a'
        elif icorr == 1.0:
            suffix += 'b'
        else:
            print 'illegal corrosion rate'
            sys.exit(1)
    return suffix


def front_deprecated(icorr_mean_list):
    suffix = rate2suffix(icorr_mean_list)
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename = 'popdata_'+suffix+'.npz'
    datafile = os.path.join(datapath,filename)
    if os.path.isfile(datafile) is False:
        print 'no data available, execute Life_Cycle_Optimization.py with \
                icorr_mean_list={} fist'.format(icorr_mean_list)
        sys.exit(1)
    else:
        popdata = np.load(datafile)
        # all pop is the same as pop
        allpop = popdata['allpop']
        allfits = popdata['allfits']
        front = popdata['front']
        frontfits = popdata['frontfits']
        pop = popdata['pop']
        popfits = popdata['popfits']

    plt.ion()
    plt.figure()

    ##plt.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], 'bo', markeredgecolor='b')
    for ind, popfit in zip(front, frontfits):
        plt.semilogx(popfit[0], popfit[1], 'bo',
                label=u'flexure: {:d}, shear: {:d}, deck: {:d}'.format(ind[0], ind[1], ind[2]))

    plt.ylim((-1,np.max(popfits)*1.01))
    ax = plt.gca()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

    plt.xlabel(u'Failure probability (log scale)', fontsize=12)
    plt.ylabel(u'Strengthening cost (mm\\textsuperscript{3})', fontsize=12)

    annotate_text = u'P={x:.2e}, C={y:.2e}\n {{ {label} }}'
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=12,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    pause = raw_input('press any key after annotation...')

    #plt.semilogx(np.array(allfits)[:,0], np.array(allfits)[:,1], 'o', markerfacecolor='lightgrey',
            #markeredgecolor='lightgrey', alpha=0.8)
    #plt.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], 'bo', markeredgecolor='b')

def front(icorr_mean_list):
    suffix = rate2suffix(icorr_mean_list)
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename = 'popdata_'+suffix+'_MOPSO2.npz'
    filename2 = 'popdata_'+suffix+'.npz'
    datafile2 = os.path.join(datapath,filename)
    if os.path.isfile(datafile) is False:
        print 'no data available, execute Dependent_Optimization_MOPSO2.py with \
                icorr_mean_list={} fist'.format(icorr_mean_list)
        sys.exit(1)
    else:
        popdata = np.load(datafile)
        popdata2 = np.load(datafile2)
        # all pop is the same as pop
        allpop = popdata['allpop']
        allfits = popdata['allfits']
        front2 = popdata2['front']
        frontfits2 = popdata2['frontfits']
        front = popdata['front']
        frontfits = popdata['frontfits']
        pop = popdata['pop']
        popfits = popdata['popfits']

    plt.ion()
    plt.figure()

    ##plt.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], 'bo', markeredgecolor='b')
    for ind, popfit in zip(front2, frontfits2):
        ## journal version
        #plt.semilogx(popfit[0], popfit[1], 'bo',markeredgecolor='b',
                #label=u'flexure: {:d}, shear: {:d}, deck: {:d}'.format(ind[0], ind[1], ind[2]))
        # conference version
        plt.semilogx(popfit[0], popfit[1], 'b.',markeredgecolor='b',
                label=u'({:d},{:d},{:d})'.format(ind[0], ind[1], ind[2]))

    plt.ylim((-1,np.max(popfits)*1.01))
    ax = plt.gca()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

    plt.xlabel(u'Failure probability')
    plt.ylabel(u'Strengthening cost (mm\\textsuperscript{3})')

    ## journal version
    #annotate_text = u'P={x:.2e}, C={y:.2e}\n {{ {label} }}'
    # conference version
    annotate_text = u'{label}'
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    pause = raw_input('press any key after annotation...')

    plt.semilogx(np.array(allfits)[:,0], np.array(allfits)[:,1], 'o', markerfacecolor='lightgrey',
            markeredgecolor='lightgrey', alpha=0.8)
    plt.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], 'bo', markeredgecolor='b')

def compare_indicator(icorr_mean_list):
    suffix = rate2suffix(icorr_mean_list)
    filename1 = 'popdata_'+suffix+'_beta_NSGA.npz'
    filename2 = 'popdata_'+suffix+'_MOPSO2.npz'
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    datafile1 = os.path.join(datapath,filename1)
    datafile2 = os.path.join(datapath,filename2)
    datafiles = [datafile1, datafile2]
    fronts = []
    for datafile in datafiles:
        if os.path.isfile(datafile1) is False:
            print 'no data available, execute optimization first'
            sys.exit(1)
        else:
            popdata = np.load(datafile)
            front = popdata['front']
            fronts.append(front)

    distance=[]
    if len(fronts[0])>len(fronts[1]):
        front1 = fronts[1]
        front2 = fronts[0]
    else:
        front1 = fronts[0]
        front2 = fronts[1]
    for sol,solfit in zip(front1,frontfit1):
        vec = np.array(front2)-np.array(sol)
        distance.append(min(np.linalg.norm(vec, axis=1)))
    print 'distance in the parameter space: {}'.format(np.mean(distance))

    from mpl_toolkits.mplot3d import Axes3D
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for front,c,m in zip(fronts, ['b','r'], ['o','^'], ['Point-in-time', 'Cumulative-time']):
        ax.scatter(front[:,0], front[:,1], front[:,2], c=c, marker=m, label=lb)
    ax.set_xlim([-1,60])
    ax.set_ylim([-1,60])
    ax.set_zlim([-1,40])
    ax.set_xlabel('Flexure')
    ax.set_ylabel('Shear')
    ax.set_zlabel('Deck',rotation=90)
    ax.xaxis._axinfo['label']['space_factor'] = 2.8
    ax.yaxis._axinfo['label']['space_factor'] = 2.8
    ax.zaxis._axinfo['label']['space_factor'] = 2.8
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))


def compare_optmization(icorr_mean_list,indicator=None):
    suffix = rate2suffix(icorr_mean_list)
    if indicator is None or indicator.lower() == 'cumulative':
        filename1 = 'popdata_'+suffix+'_NSGA.npz'
        filename2 = 'popdata_'+suffix+'_MOPSO2.npz'
    else:
        filename1 = 'popdata_'+suffix+'beta_NSGA.npz'
        filename2 = 'popdata_'+suffix+'beta_MOPSO2.npz'
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    datafile1 = os.path.join(datapath,filename1)
    datafile2 = os.path.join(datapath,filename2)
    datafiles = [datafile1, datafile2]
    frontfits_list = []
    fronts = []
    for datafile in datafiles:
        if os.path.isfile(datafile1) is False:
            print 'no data available, execute optimization first'
            sys.exit(1)
        else:
            popdata = np.load(datafile)
            frontfits = popdata['frontfits']
            front = popdata['front']
            fronts.append(front)
            frontfits_list.append(frontfits)

    distance=[]
    fitdistance=[]
    frontfit1 = frontfits_list[0]
    frontfit2 = frontfits_list[1]
    if len(fronts[0])>len(fronts[1]):
        frontfit1 = frontfits_list[1]
        frontfit2 = frontfits_list[0]
        front1 = fronts[1]
        front2 = fronts[0]
    else:
        frontfit1 = frontfits_list[0]
        frontfit2 = frontfits_list[1]
        front1 = fronts[0]
        front2 = fronts[1]
    for sol,solfit in zip(front1,frontfit1):
        vec = np.array(front2)-np.array(sol)
        fitvec = np.array(frontfit2)-np.array(solfit)
        distance.append(min(np.linalg.norm(vec, axis=1)))
        fitdistance.append(min(np.linalg.norm(fitvec, axis=1)))
    print 'distance in the parameter space: {}'.format(np.mean(distance))
    print 'distance in the objective space: {}'.format(np.mean(fitdistance))

    from mpl_toolkits.mplot3d import Axes3D
    plt.ion()
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111, projection='3d')
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax1.set_xlabel(u'Failure probability')
    ax1.set_ylabel(u'Strengthening cost (mm\\textsuperscript{3})')
    ax2.set_xlabel(u'Flexure'); #ax2.set_xlim([-1,100])
    ax2.set_ylabel(u'Shear'); #ax2.set_ylim([-1,100])
    ax2.set_zlabel(u'Deck'); #ax2.set_zlim([-1,100])
    for front, frontfits,c,m,lb in zip(fronts, frontfits_list, ['b','r'], ['o','^'], ['NSGA-II', 'MOPSO-II']):
        ax1.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], c=c, marker=m, ls='None',label=lb)
        ax2.scatter(front[:,0], front[:,1], front[:,2], c=c, marker=m,label=lb)



def costkeeping(icorr_mean_list):
    suffix = rate2suffix(icorr_mean_list)
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename = 'costkeeping_'+suffix+'.npz'
    datafile = os.path.join(datapath,filename)
    if os.path.isfile(datafile) is False:
        print 'no data available, execute Life_Cycle_Optimization.py with \
                icorr_mean_list={} fist'.format(icorr_mean_list)
        sys.exit(1)
    else:
        costkeeping = np.load(datafile)

    plt.ion()
    plt.figure()
    plt.plot(costkeeping['flexure'][0,:], costkeeping['flexure'][1,:], 'b-',
            label='Flexure')
    plt.plot(costkeeping['shear'][0,:], costkeeping['shear'][1,:], 'r--',
            label='Shear')
    plt.plot(costkeeping['deck'][0,:], costkeeping['deck'][1,:], 'g-.',
            label='Deck')
    service_life = np.max((np.max(costkeeping['flexure'][0,:]),
            np.max(costkeeping['shear'][0,:]), np.max(costkeeping['deck'][0,:])))
    max_cost = np.max((np.max(costkeeping['flexure'][1,:]),
            np.max(costkeeping['shear'][1,:]), np.max(costkeeping['deck'][1,:])))
    plt.xlim((-1,service_life+1))
    plt.ylim((-1,max_cost*1.01))

    plt.xlabel('Strengthening time (year)')
    plt.ylabel('Cost (mm\\textsuperscript{3})')

    ax = plt.gca()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    pause = raw_input('press any key after annotation...')


def pfkeeping(icorr_mean_list):
    suffix = rate2suffix(icorr_mean_list)
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename = 'pfkeeping_'+suffix+'.npz'
    datafile = os.path.join(datapath,filename)
    if os.path.isfile(datafile) is False:
        print 'no data available, execute Life_Cycle_Optimization.py with \
                icorr_mean_list={} fist'.format(icorr_mean_list)
        sys.exit(1)
    else:
        pfkeeping = np.load(datafile)

    plt.ion()
    fig = plt.figure()
    ax0 = plt.subplot2grid((2,4), (0,0), colspan=2)
    ax1 = plt.subplot2grid((2,4), (0,2), colspan=2)
    ax2 = plt.subplot2grid((2,4), (1,1), colspan=2)
    ax0.semilogy(pfkeeping['flexure'][0,:], pfkeeping['flexure'][1,:], 'bo')
    ax1.semilogy(pfkeeping['shear'][0,:], pfkeeping['shear'][1,:], 'bo')
    ax2.semilogy(pfkeeping['deck'][0,:], pfkeeping['deck'][1,:], 'bo')

    # labels
    ax0.set_title('(a) Flexure', fontsize=12)
    ax0.set_xlabel('Strengthening time (year)', fontsize=12)
    ax0.set_ylabel('Failure prob.', fontsize=12)
    ax1.set_title('(b) Shear', fontsize=12)
    ax1.set_xlabel('Strengthening time (year)', fontsize=12)
    ax1.set_ylabel('Failure prob.', fontsize=12)
    ax2.set_title('(c) Deck', fontsize=12)
    ax2.set_xlabel('Strengthening time (year)', fontsize=12)
    ax2.set_ylabel('Failure prob.', fontsize=12)

    # change w/hspace
    plt.subplots_adjust(wspace=0.6, hspace=0.3)
    plt.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.08)


def history(icorr_mean_list, str_yr_2dlist):
    if len(str_yr_2dlist) == 1:
        # only one strengthening option is plotted, thus component pfs are
        # added

        # load data
        suffix = rate2suffix(icorr_mean_list)
        filename = 'pfhistory_str_'
        for ti in str_yr_2dlist[0]:
            filename = filename + str(int(ti)) + '_'
        datapath = os.path.join(os.path.abspath('./'), 'data')
        filename = filename+suffix+'.npz'
        datafile = os.path.join(datapath,filename)

        if os.path.isfile(datafile) is False:
            print 'no data available, execute Life_Cycle_History.py with \
                    icorr_mean_list={} and str_yr_list={} fist'.format(icorr_mean_list,
                            str_yr_list)
            sys.exit(1)
        else:
            pfhistory = np.load(datafile)
            time_array = pfhistory['time']
            pf_sys = pfhistory['system']
            pf_flex = pfhistory['flexure']
            pf_shear = pfhistory['shear']
            pf_deck = pfhistory['deck']

        plt.ion()
        plt.figure()
        plt.semilogy(time_array, pf_flex, 'b', ls='--', label='Flexure')
        plt.semilogy(time_array, pf_shear, 'r', ls='-.', label='Shear')
        plt.semilogy(time_array, pf_deck, 'g', ls='-', label='Deck and system')
        plt.semilogy(time_array, pf_sys, 'ko', ls='-', label='Deck and system')

    else:
        plt.ion()
        plt.figure()
        # multiple strengthening options, only system pfs are plotted
        ls_list = ['-', '--', '-.', ':', '  ', ' ']
        for str_yr_list,ls in zip(str_yr_2dlist, ls_list):
            # load data
            suffix = rate2suffix(icorr_mean_list)
            filename = 'pfhistory_str_'
            for ti in str_yr_list:
                filename = filename + str(int(ti)) + '_'
            datapath = os.path.join(os.path.abspath('./'), 'data')
            filename = filename+suffix+'.npz'
            datafile = os.path.join(datapath,filename)

            if os.path.isfile(datafile) is False:
                print 'no data available, execute Life_Cycle_History.py with \
                        icorr_mean_list={} and str_yr_list={} fist'.format(icorr_mean_list,
                                str_yr_list)
                sys.exit(1)
            else:
                pfhistory = np.load(datafile)
                time_array = pfhistory['time']
                pf_sys = pfhistory['system']
                ## journal
                #plt.semilogy(time_array, pf_sys, ls=ls,
                        #label=u'flexure: {:d}, shear: {:d}, deck: {:d}'.format(
                            #str_yr_list[0], str_yr_list[1], str_yr_list[2]))
                # conference
                plt.semilogy(time_array, pf_sys, ls=ls,
                        label=u'({:d},{:d},{:d})'.format(
                            str_yr_list[0], str_yr_list[1], str_yr_list[2]))

    plt.xlabel('Time (year)')
    plt.ylabel('Failure probability')

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))


def pointintimehistory(icorr_mean_list, str_yr_2dlist):
    if len(str_yr_2dlist) == 1:
        # only one strengthening option is plotted, thus component pfs are
        # added

        # load data
        suffix = rate2suffix(icorr_mean_list)
        filename = 'pointpfhistory_str_'
        for ti in str_yr_2dlist[0]:
            filename = filename + str(int(ti)) + '_'
        datapath = os.path.join(os.path.abspath('./'), 'data')
        filename = filename+suffix+'.npz'
        datafile = os.path.join(datapath,filename)

        if os.path.isfile(datafile) is False:
            print 'no data available, execute Point_Life_Cycle_History.py with \
                    icorr_mean_list={} and str_yr_list={} fist'.format(icorr_mean_list,
                            str_yr_list)
            sys.exit(1)
        else:
            pfhistory = np.load(datafile)
            time_array = pfhistory['time']
            pf_sys = pfhistory['system']
            pf_flex = pfhistory['flexure']
            pf_shear = pfhistory['shear']
            pf_deck = pfhistory['deck']

        plt.ion()
        plt.figure()
        plt.semilogy(time_array, pf_flex, 'b', ls='--', label='Flexure')
        plt.semilogy(time_array, pf_shear, 'r', ls='-.', label='Shear')
        plt.semilogy(time_array, pf_deck, 'g', ls='-', label='Deck')
        plt.semilogy(time_array, pf_sys, 'ko', ls='-', label='System')

    else:
        plt.ion()
        plt.figure()
        # multiple strengthening options, only system pfs are plotted
        ls_list = ['-', '--', '-.', ':', '  ', ' ']
        for str_yr_list,ls in zip(str_yr_2dlist, ls_list):
            # load data
            suffix = rate2suffix(icorr_mean_list)
            filename = 'pfhistory_str_'
            for ti in str_yr_list:
                filename = filename + str(int(ti)) + '_'
            datapath = os.path.join(os.path.abspath('./'), 'data')
            filename = filename+suffix+'.npz'
            datafile = os.path.join(datapath,filename)

            if os.path.isfile(datafile) is False:
                print 'no data available, execute Point_Life_Cycle_History.py with \
                        icorr_mean_list={} and str_yr_list={} fist'.format(icorr_mean_list,
                                str_yr_list)
                sys.exit(1)
            else:
                pfhistory = np.load(datafile)
                time_array = pfhistory['time']
                pf_sys = pfhistory['system']
                ## journal
                #plt.semilogy(time_array, pf_sys, ls=ls,
                        #label=u'flexure: {:d}, shear: {:d}, deck: {:d}'.format(
                            #str_yr_list[0], str_yr_list[1], str_yr_list[2]))
                # conference
                plt.semilogy(time_array, pf_sys, ls=ls,
                        label=u'({:d},{:d},{:d})'.format(
                            str_yr_list[0], str_yr_list[1], str_yr_list[2]))

    plt.xlabel('Time (year)')
    plt.ylabel('Failure probability')

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))


def lifetimefitting(icorr_mean_list, str_yr_list=[0., 0., 0.]):
    from scipy import stats
    from scipy.optimize import curve_fit

    # load data
    suffix = rate2suffix(icorr_mean_list)
    filename = 'pfhistory_str_'
    for ti in str_yr_list:
        filename = filename + str(int(ti)) + '_'
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename = filename+suffix+'.npz'
    datafile = os.path.join(datapath,filename)

    if os.path.isfile(datafile) is False:
        print 'no data available, execute Life_Cycle_History.py with \
                icorr_mean_list={} and str_yr_list={} fist'.format(icorr_mean_list,
                        str_yr_list)
        sys.exit(1)
    else:
        pfhistory = np.load(datafile)
        time_array = pfhistory['time']
        pf_sys = pfhistory['system']
        # lifetime fitting
        # candidate cdfs
        def fitexpon(xdata, *params):
            lbd = params[0]
            return stats.expon.cdf(xdata, scale=1./lbd)
        def fitweibull(xdata, *params):
            k = params[0]    # shape param k in Weibull wiki
            lbd = params[1]    # scale param lbd in Weibull wiki
            return stats.weibull_min.cdf(xdata, k, scale=lbd)
        def fitgamma(xdata, *params):
            k=params[0]
            theta = params[1]
            return stats.gamma.cdf(xdata, k, scale=theta)
        poptExpon, pcovExpon = curve_fit(fitexpon, time_array, pf_sys, p0=[1.], bounds=(0.,np.inf))
        poptWbl, pcovWbl = curve_fit(fitweibull, time_array, pf_sys, p0=[1.,5.], bounds=([0.,0.],[np.inf, np.inf]))
        poptGamma, pcovGamma = curve_fit(fitgamma, time_array, pf_sys, p0=[1.,1.], bounds=([0.,0.],[np.inf,np.inf]))

    plt.ion()
    plt.figure()
    plt.semilogy(time_array, pf_sys, 'o', label='$T_f$ data')
    plt.semilogy(time_array, fitexpon(time_array, poptExpon[0]), ls='--', label='Exponential')
    plt.semilogy(time_array, fitweibull(time_array,poptWbl[0],poptWbl[1]), ls='-', label='Weibull')
    plt.semilogy(time_array, fitgamma(time_array,poptGamma[0],poptGamma[1]), ls=':', label='Gamma')
    plt.xlabel('Time (year)')
    plt.ylabel('Failure probability')
    #plt.legend(loc='lower right', fontsize=9)

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))
