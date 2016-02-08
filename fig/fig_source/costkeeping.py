#=========================================================================
# post processing of costkeeping
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
plt.ion()
plt.close('all')

if __name__ == '__main__':
    costkeeping = np.load('costkeeping_bbb.npz')
    plt.figure()
    plt.plot(costkeeping['flexure'][0,:], costkeeping['flexure'][1,:], 'bo',
            label='flexure cost')
    plt.plot(costkeeping['shear'][0,:], costkeeping['shear'][1,:], 'r^',
            label='shear cost')
    plt.plot(costkeeping['deck'][0,:], costkeeping['deck'][1,:], 'gv',
            label='deck cost')
    plt.xlim((-1,101))
    plt.ylim((-1,8e7))

    plt.xlabel('strengthening time (year)', fontsize=12)
    plt.ylabel('cost (mm\\textsuperscript{3})', fontsize=12)

    ax = plt.gca()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

    datacursor(formatter='{label}'.format,display='multiple', draggable=True,
            bbox=None, fontsize=12,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))
