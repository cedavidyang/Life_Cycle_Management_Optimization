import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpldatacursor import datacursor
plt.rc('font', family='serif', size=12)
plt.rc('text', usetex=True)
#import matplotlib as mpl
#mpl.use("pgf")
#pgf_with_custom_preamble = {
    #"font.family": "serif", # use serif/main font for text elements
    #"font.size": 12, # use font size
    #"text.usetex": True,    # use inline math for ticks
    #"pgf.rcfonts": False,   # don't setup fonts from rc parameters
    #"pgf.preamble": [
        #"\\usepackage{siunitx}",         # load additional packages
        #]
#}
#mpl.rcParams.update(pgf_with_custom_preamble)

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


def front(icorr_mean_list):
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
        allpop = popdata['allpop']
        allfits = popdata['allfits']
        front = popdata['front']
        frontfits = popdata['frontfits']
        pop = popdata['pop']
        popfits = popdata['popfits']

    plt.ion()
    plt.figure()

    ##plt.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], 'bo', markeredgecolor='b')
    for ind, popfit in zip(pop, popfits):
        plt.semilogx(popfit[0], popfit[1], 'b.',
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

    plt.semilogx(np.array(allfits)[:,0], np.array(allfits)[:,1], 'o', markerfacecolor='lightgrey',
            markeredgecolor='lightgrey', alpha=0.8)
    plt.semilogx(np.array(popfits)[:,0], np.array(popfits)[:,1], 'bo', markeredgecolor='b')


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
    plt.plot(costkeeping['flexure'][0,:], costkeeping['flexure'][1,:], 'bo',
            label='flexure cost')
    plt.plot(costkeeping['shear'][0,:], costkeeping['shear'][1,:], 'r^',
            label='shear cost')
    plt.plot(costkeeping['deck'][0,:], costkeeping['deck'][1,:], 'gv',
            label='deck cost')
    service_life = np.max((np.max(costkeeping['flexure'][0,:]),
            np.max(costkeeping['shear'][0,:]), np.max(costkeeping['deck'][0,:])))
    max_cost = np.max((np.max(costkeeping['flexure'][1,:]),
            np.max(costkeeping['shear'][1,:]), np.max(costkeeping['deck'][1,:])))
    plt.xlim((-1,service_life+1))
    plt.ylim((-1,max_cost*1.01))

    plt.xlabel('strengthening time (year)', fontsize=12)
    plt.ylabel('cost (mm\\textsuperscript{3})', fontsize=12)

    ax = plt.gca()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=12,
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
    ax0.set_title('(a) flexure', fontsize=12)
    ax0.set_xlabel('strengthening time (year)', fontsize=12)
    ax0.set_ylabel('failure prob.', fontsize=12)
    ax1.set_title('(b) shear', fontsize=12)
    ax1.set_xlabel('strengthening time (year)', fontsize=12)
    ax1.set_ylabel('failure prob.', fontsize=12)
    ax2.set_title('(c) deck', fontsize=12)
    ax2.set_xlabel('strengthening time (year)', fontsize=12)
    ax2.set_ylabel('failure prob.', fontsize=12)

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
        plt.semilogy(time_array, pf_flex, 'b', ls='--', label='flexure')
        plt.semilogy(time_array, pf_shear, 'r', ls='-.', label='shear')
        plt.semilogy(time_array, pf_deck, 'g', ls='-', label='deck and system')
        plt.semilogy(time_array, pf_sys, 'ko', ls='-', label='deck and system')

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
                plt.semilogy(time_array, pf_sys, ls=ls,
                        label=u'flexure: {:d}, shear: {:d}, deck: {:d}'.format(
                            str_yr_list[0], str_yr_list[1], str_yr_list[2]))

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=12,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    pause = raw_input('press any key after annotation...')

    plt.xlabel('time (year)', fontsize=12)
    plt.ylabel('failure probability', fontsize=12)
