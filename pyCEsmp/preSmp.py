# cross-entropy-based adaptive importance sampling
import numpy as np
import scipy.stats as stats

from scipy.integrate import quad
from sklearn.mixture.gmm import GMM
from constants import STD_RATIO1, STD_RATIO2, STD_RATIO3, MAX_NUM_UPDATE, KRATIO_CRITERIA
from pyre.distributions import *

from pyCEsmp.mainSmp import MainSmp

class PreSmp(MainSmp):
    def __init__(self, analysis_type, nadp, nsmp, gmdistr, rv_array, r_array, sl_array, rate_array, gd_array, gcov_array):
        self.nadp = nadp
        MainSmp.__init__(self, analysis_type, nsmp, gmdistr, rv_array, r_array, sl_array, rate_array, gd_array, gcov_array)

    def setAdpNum(self, new_nadp):
        self.nadp = new_nadp

    def getAdpNum(self):
        return self.nadp

    def updateParameter(self, smps, responsibilities, gf2hv):

        nV = self.gmdistr.means_.shape[1]
        nCompo = self.gmdistr.get_params()['n_components']

        # update weights
        if nCompo == 1:
            piV = np.array([1.0])
        else:
            if np.sum(gf2hv) == 0 or np.isnan(np.sum(gf2hv)) or np.isinf(np.sum(gf2hv)):
                raise ValueError("[updating weights]: zero division")
            piV = np.sum( np.dot(gf2hv[:, np.newaxis], np.ones((1,nCompo))) * responsibilities, axis=0 ) / np.sum(gf2hv)
            piV = piV / np.sum(piV)

        # update mean
        #if nCompo == 1:
        #    rm_mean = self.r_array[0].stats(moments='m')
        #    rv_mean = self.r_array[1].stats(moments='m')
        #    muV = np.array([[rm_mean, rv_mean]])
        #else:
        muV = np.zeros(self.gmdistr.means_.shape)
        for iCompo in range(nCompo):
            tmps = gf2hv * responsibilities[:, iCompo]
            if np.sum(tmps) == 0 or np.isnan(np.sum(tmps)) or np.isinf(np.sum(tmps)):
                raise ValueError("[updating means]: zero division")
            for iV in range(nV):
                muV[iCompo, iV] = np.dot( smps[:, iV], tmps ) / np.sum(tmps)

        # update covariance
        covarV = np.zeros(self.gmdistr.covars_.shape)
        muVold = self.gmdistr.means_
        for iCompo in range(nCompo):
            tmps = gf2hv * responsibilities[:, iCompo]
            if np.sum(tmps) == 0 or np.isnan(np.sum(tmps)) or np.isinf(np.sum(tmps)):
                raise ValueError("[updating covariance]: zero division")
            for iV in range(nV):
                for jV in range(nV):
                    smps_i = smps[:, iV]
                    smps_j = smps[:, jV]
                    covarV[iCompo, iV, jV] = np.dot( (smps_i-muVold[iCompo, iV]) * (smps_j-muVold[iCompo, jV]), tmps ) / np.sum(tmps)

        self.gmdistr = GMM(n_components = nCompo, covariance_type='full')
        self.gmdistr.weights_ = piV
        self.gmdistr.means_ = muV
        self.gmdistr.covars_ = covarV

    def adaptation(self, ti):

        kadp = np.ones(self.nadp)
        kadp[0] = STD_RATIO1
        kadp[1:3] = STD_RATIO2
        kadp[3:5] = STD_RATIO3

        sumX = 0; sumX2 = 0;
        for iadp in range(self.nadp):
            nCompo = self.gmdistr.get_params()['n_components']
            tmp_gmdistr = GMM(n_components=nCompo,covariance_type='full')
            tmp_gmdistr.weights_ = self.gmdistr.weights_
            tmp_gmdistr.means_ = self.gmdistr.means_
            tmp_gmdistr.covars_ = self.gmdistr.covars_ * kadp[iadp]**2
            subSmp = MainSmp(self.analysis_type, self.nsmp, tmp_gmdistr, self.rv_array, self.r_array, self.sl_array, self.rate_array, self.gd_array, self.gcov_array)
            i = 0
            while i<=MAX_NUM_UPDATE:
                try:
                    smps = subSmp.sample()
                    f = subSmp.priorPDF(smps)
                    g = 1 - subSmp.condAvailability(smps, ti)
                    hv, responsibilities = subSmp.score_samples(smps)
                    gf2hv = g*f/hv

                    self.updateParameter(smps, responsibilities, gf2hv)

                    i += 1
                    break
                except ValueError:
                    i += 1
                    continue

            sumX += np.sum(gf2hv)
            sumX2 += np.sum(gf2hv**2)

        npre = self.nadp * self.nsmp
        pfpre = 1./npre * sumX
        Spf2pre = 1./(npre*(npre-1)) * (sumX2 - 2*sumX*pfpre + npre * pfpre**2)

        return pfpre, Spf2pre

    def getKopt(self, k_array, ti):
        pf = np.zeros(k_array.shape)
        for k in k_array:
            nCompo = self.gmdistr.get_params()['n_components']
            tmp_gmdistr = GMM(n_components=nCompo,covariance_type='full')
            tmp_gmdistr.weights_ = self.gmdistr.weights_
            tmp_gmdistr.means_ = self.gmdistr.means_
            tmp_gmdistr.covars_ = self.gmdistr.covars_ * k**2
            subSmp = MainSmp(self.analysis_type, self.nsmp, tmp_gmdistr, self.rv_array, self.r_array, self.sl_array, self.rate_array, self.gd_array, self.gcov_array)
            smps = subSmp.sample()

            indx = np.where(k_array==k)[0][0]
            pf[indx], dummy = subSmp.getPf(smps, ti)
            if indx>0 and pf[indx] < KRATIO_CRITERIA * np.mean(pf[:indx]):
                k_indx = indx-1
                break
            else:
                k_indx = -1

        return k_array[k_indx]


if __name__ == '__main__':
    ti = 1
    nCompo = 4
    nsmp = 400
    nadp = 20
    nk = 100
    muV0 = np.array([[2768.03327821769, 1519.21878505205, 2060.98121494796, 812.166721782312],
                       [812.166721782312, 2060.98121494796, 1519.21878505205, 2768.03327821769]]).T
    muV0 = np.load('peaks0.npy')
    piV0 = np.array([0.25, 0.25, 0.25, 0.25])
    covarV0 = np.array([[[80300.0, 0.0],[0.0, 80300.0]], [[80300.0, 0.0],[0.0, 80300.0]],
                          [[80300.0, 0.0],[0.0, 80300.0]], [[80300.0, 0.0],[0.0, 80300.0]]])
    gmdistr = GMM(n_components=4, covariance_type='full')
    gmdistr.weights_ = piV0
    gmdistr.means_ = muV0
    gmdistr.covars_ = covarV0
    resistRv = Lognormal('r', mean=1790.1, stdv=283.37)
    rDistr = resistRv.rv
    #r_array = np.array([stats.norm(loc=1790.1, scale=283.37), stats.norm(loc=1790.1, scale=283.37)])
    r_array = np.array([rDistr, rDistr])
    sl_array = np.array([stats.norm(loc=301.1, scale=120.44)])
    ti_mean = 15.91
    a =  0.0135;    b = 0.8580
    gd = lambda x: (1-a*(x-ti_mean)**b) if x>ti_mean else 1.0
    gd_array = np.array([gd, gd])
    k_array = np.arange(2.0, 1.0-0.05/2, -0.05)
    
    testSmp = PreSmp(nadp, nsmp, gmdistr, r_array, sl_array, gd_array)
    pf, Spf2 = testSmp.adaptation(ti)
    piVopt = testSmp.gmdistr.weights_
    muVopt = testSmp.gmdistr.means_
    covarVopt = testSmp.gmdistr.covars_
    kopt = testSmp.getKopt(k_array, ti)
