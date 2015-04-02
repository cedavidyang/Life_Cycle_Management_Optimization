# CE-based importance sampling
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
from sklearn.mixture.gmm import GMM

from constants.sampleConstants import COVT0_COV, INIT_K, K_STEPS
from CEfunctions import combine2smps
from robustness import getInitialPeaks, getInitialCovar
from mainSmp import MainSmp
from preSmp import PreSmp
from pyre.distributions import *

class CEImportanceSmp(object):
    def __init__(self, analysis_type, ncomp=None, nadp=None, nsmp=None, nmain=None):
        self.analysis_type = analysis_type
        self.ncomp = ncomp
        self.nadp = nadp
        self.nsmp = nsmp
        self.nmain = nmain
        # sampling objects
        self.mainSmp = None
        self.preSmp = None
        # sampling results
        self.pfpre = None
        self.Spf2pre = None
        self.kopt = None
        self.pfmain = None
        self.Spf2main = None
        self.pf = None
        self.Spf2 = None

    def setPreSmp(self, r_array, sl_array, rate_array, gd_array, gcov_array):
        # initial peaks
        corrTarget = np.eye(self.ncomp)
        rv_list = []
        for icomp in xrange(self.ncomp):
            cov_variable = Normal('cov', mean=0., stdv=COVT0_COV)
            rv_list.append(cov_variable.rv)
        rv_array = np.array(rv_list)
        peaks0 = getInitialPeaks(rv_array, self.ncomp, corrTarget=corrTarget, options={'disp':False, 'anneal':True})
        weights0 = 1./self.ncomp * np.ones(self.ncomp)
        covar0 = getInitialCovar(rv_array, self.ncomp)

        # initial values
        gmdistr = GMM(n_components=self.ncomp , covariance_type='full')
        gmdistr.weights_ = weights0
        gmdistr.means_ = peaks0
        gmdistr.covars_ = covar0

        self.preSmp = PreSmp(self.analysis_type, self.nadp, self.nsmp, gmdistr, rv_array, r_array,
            sl_array, rate_array, gd_array, gcov_array)

    def conductPreSmp(self, timepoint):
        if self.preSmp is None:
            print '[ERROR:] PreSmp is not initialized, MainSmp cannot be initialized before PreSmp'
            sys.exit(1)
        else:
            # preliminary sampling
            pfpre, Spf2pre = self.preSmp.adaptation(timepoint)
            k_array = np.arange(INIT_K, 1.0-K_STEPS, -K_STEPS)
            # get kopt
            kopt = self.preSmp.getKopt(k_array, timepoint)

        self.pfpre = pfpre
        self.Spf2pre = Spf2pre
        self.kopt = kopt

    def setMainSmp(self):
        if self.pfpre is None:
            print '[ERROR:] PreSmp has not been conducted, MainSmp cannot be initialized'
            sys.exit(1)
        else:
            self.mainSmp = MainSmp(preSmp=self.preSmp)
            self.mainSmp.setSmpNum(self.nmain)

    def conductMainSmp(self, timepoint):
        if self.preSmp is None:
            print '[ERROR:] MainSmp is not initialized, cannot conduct main sampling'
            sys.exit(1)
        else:
            smps = self.mainSmp.sample()
            pfmain, Spf2main = self.mainSmp.getPf(smps, timepoint)
            pf, Spf2, COVpf = combine2smps(self.pfpre, self.Spf2pre, pfmain, Spf2main)
            self.pfmain = pfmain
            self.Spf2main = Spf2main
            self.pf = pf
            self.Spf2 = Spf2

    def reset(self):
        # reset results, must be implemented after conductMainSmp
        self.preSmp = None
        self.pfpre = None
        self.Spf2pre = None
        self.pfmain = None
        self.Spf2main = None
        self.pf = None
        self.Spf2 = None
