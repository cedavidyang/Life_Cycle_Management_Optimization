# likelihood functions
import numpy as np
import scipy.stats as stats
import scipy.special as spec

from pyre.distributions import *
#from functions import *
from constants.evidenceConstants import *


def halfcellLikelihood(isIni_smp, half_cell_evidence):
    half_cell_likelihood_isIni = stats.norm.pdf(half_cell_evidence, ISINI_HALF_CELL_MEAN, ISINI_HALF_CELL_STD) / stats.norm.pdf(ISINI_HALF_CELL_MEAN, ISINI_HALF_CELL_MEAN, ISINI_HALF_CELL_STD)
    half_cell_likelihood_noIni = stats.norm.pdf(half_cell_evidence, NOINI_HALF_CELL_MEAN, NOINI_HALF_CELL_STD) / stats.norm.pdf(NOINI_HALF_CELL_MEAN, NOINI_HALF_CELL_MEAN, NOINI_HALF_CELL_STD)
    half_cell_likelihood = np.empty(isIni_smp.shape)
    half_cell_likelihood[isIni_smp] = half_cell_likelihood_isIni
    half_cell_likelihood[np.logical_not(isIni_smp)] = half_cell_likelihood_noIni

    return half_cell_likelihood


def icorrLikelihood(icorr_smp, icorr_evidence):
    icorr_likelihood = np.zeros(icorr_smp.shape)
    icorr_likelihood[np.logical_and(icorr_smp==0, icorr_evidence==0)] = 1.0
    icorr_likelihood[np.logical_and(icorr_smp==0, icorr_evidence!=0)] = 0.0
    icorr_likelihood[icorr_smp!=0] = stats.norm.pdf(icorr_evidence, icorr_smp[icorr_smp!=0], icorr_smp[icorr_smp!=0]*MEASURED_ICORR_COV) / stats.norm.pdf(0, 0, MEASURED_ICORR_COV)

    return icorr_likelihood
