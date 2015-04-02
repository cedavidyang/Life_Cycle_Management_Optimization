#=========================================================================
# post processing of pfkeeping
#=========================================================================
import numpy as np
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
#from mpldatacursor import datacursor
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=12)
plt.rc('text', usetex=True)
plt.ion()
plt.close('all')

if __name__ == '__main__':
    pfkeeping = np.load('pfkeeping.npz')

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

    plt.show()
