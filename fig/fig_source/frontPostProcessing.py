#=========================================================================
# plot the Pareto fronts
#=========================================================================
import numpy as np
import os
import sys
import matplotlib as mpl
mpl.use("pgf")
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "font.size": 12, # use font size
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
        "\\usepackage{siunitx}",         # load additional packages
        ]
}
mpl.rcParams.update(pgf_with_custom_preamble)
import matplotlib.pyplot as plt
from mpldatacursor import datacursor

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

def main(icorr_mean_list):
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
        popfits = popfitsdata['popfits']

    plt.ion()
    plt.close('all')
    #plt.rc('font', family='serif', size=12)
    #plt.rc('text', usetex=False)

    ##plt.semilogx(np.array(frontfits)[:,0], np.array(frontfits)[:,1], 'bo', markeredgecolor='b')
    #plt.semilogx(np.array(popfits)[:,0], np.array(popfits)[:,1], 'bo', markeredgecolor='b')
    for popfit in popfits:
        plt.semilogx(popfit[0], popfit[1], 'bo',markeredgecolor='b',
                label=u'flexure: {:d}, shear: {:d}, deck: {:d}'.format(ind[0], ind[1], ind[2]))

    plt.ylim((-1,np.max(popfits)*1.1))
    ax = plt.gca()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

    plt.xlabel(u'Failure probability (log scale)')
    plt.ylabel(u'Strengthening cost (mm\\textsuperscript{3})')

    annotate_text = u'P={x:.2e}, C={y:.2e}\n {{ {label} }}'
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=12,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    pause = raw_input('press any key after annotation...')

    plt.semilogx(np.array(allfits)[:,0], np.array(allfits)[:,1], 'o', markerfacecolor='lightgrey',
            markeredgecolor='lightgrey', alpha=0.8)

if __name__ == '__main__':
    icorr_mean_list = input('mean icorr list:\n')
    main(icorr_mean_list)
