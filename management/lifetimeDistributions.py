# calculate time-variant failure probability and fit the results to Weibull
# lifetime distribution
import os
import sys
import numpy as np
import scipy.stats as stats
from scipy.special import gammainc, erf
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import scipy.stats as stats
from sklearn.mixture.gmm import GMM

import time
import datetime

from pyre.distributions import *
from constants.beamConstants import SERVICE_LIFE, RELIABILITY_DT
from constants.beamConstants import M_LLIM_MEAN, M_LLIM_COV
from constants.beamConstants import V_LLIM_MEAN, V_LLIM_COV
from constants.beamConstants import M_LLIM_DECK_MEAN, M_LLIM_DECK_COV
from constants.beamConstants import LL_ARRIVAL_RATE, LL_ARRIVAL_RATE_DECK
from constants.sampleConstants import *
from pyCEsmp import *

def getPreSmpObject(analysis_type, new_service_time, frp_mean_history, frp_cov_history):

    # time array
    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)

    # resistance array and live load array
    flexure_mean_history = frp_mean_history[0]
    flexure_cov_history = frp_cov_history[0]
    shear_mean_history = frp_mean_history[1]
    shear_cov_history = frp_cov_history[1]
    rm_mean = flexure_mean_history[0]
    rv_mean = shear_mean_history[0]
    r_array = np.array([rm_mean, rv_mean])

    slmRv = Normal('slm', mean=M_LLIM_MEAN, stdv=M_LLIM_MEAN*M_LLIM_COV)
    slvRv = Normal('slv', mean=V_LLIM_MEAN, stdv=V_LLIM_MEAN*V_LLIM_COV)
    sldRv = Normal('sld', mean=M_LLIM_DECK_MEAN, stdv=M_LLIM_DECK_MEAN*M_LLIM_DECK_COV)
    slmDistr = slmRv.rv
    slvDistr = slvRv.rv
    sldDistr = sldRv.rv
    sl_array = np.array([slmDistr, slvDistr, sldDistr])

    rate_m = LL_ARRIVAL_RATE
    rate_v = LL_ARRIVAL_RATE
    rate_d = LL_ARRIVAL_RATE_DECK
    rate_array = np.array([rate_m, rate_v, rate_d])

    # function gd(t) and gcov(t)
    #def gd_func(x, a, b):
    #    return a+b*x
    def gd_func(x, a, b, c, d, e):
        return a+b*x**1+c*x**2+d*x**3+e*x**4
    gt_flex = flexure_mean_history/flexure_mean_history[0] / np.sqrt(1+flexure_cov_history**2)
    gt_shear = shear_mean_history/shear_mean_history[0] / np.sqrt(1+shear_cov_history**2)
    gd_flex_popt, pcov = curve_fit(gd_func, new_service_time, gt_flex)
    gd_shear_popt,pcov = curve_fit(gd_func, new_service_time, gt_shear)
    gd_flexure = lambda x: gd_func(x, gd_flex_popt[0], gd_flex_popt[1], gd_flex_popt[2], gd_flex_popt[3], gd_flex_popt[4])
    gd_shear = lambda x: gd_func(x, gd_shear_popt[0], gd_shear_popt[1], gd_shear_popt[2], gd_shear_popt[3], gd_shear_popt[4])
    #gd_flexure = lambda x: gd_func(x, gd_flex_popt[0], gd_flex_popt[1])
    #gd_shear = lambda x: gd_func(x, gd_shear_popt[0], gd_shear_popt[1])

    def gcov_func(x, a, b):    # dummy function, not used
        return a+b*x
    at_flex = np.sqrt(np.log(flexure_cov_history**2+1)/COVT0_COV**2)
    at_shear = np.sqrt(np.log(shear_cov_history**2+1)/COVT0_COV**2)
    gcov_flex_popt, pcov = curve_fit(gcov_func, new_service_time, at_flex)
    gcov_shear_popt,pcov = curve_fit(gcov_func, new_service_time, at_shear)
    #gcov_flexure = lambda x: gcov_func(x, gcov_flex_popt[0], gcov_flex_popt[1], gcov_flex_popt[2], gcov_flex_popt[3], gcov_flex_popt[4])
    #gcov_shear = lambda x: gcov_func(x, gcov_shear_popt[0], gcov_shear_popt[1], gcov_shear_popt[2], gcov_shear_popt[3], gcov_shear_popt[4])
    gcov_flexure = lambda x: at_flex[0] + 0*x
    gcov_shear = lambda x: at_shear[0] + 0*x

    gd_array = np.array([gd_flexure, gd_shear])
    gcov_array = np.array([gcov_flexure, gcov_shear])

    #plt.close('all')
    #plt.figure()
    #plt.plot(new_service_time, gt_flex, 'bo',
    #         new_service_time, gt_shear, 'rs',
    #         new_service_time, gd_flexure(new_service_time), 'b-',
    #         new_service_time, gd_shear(new_service_time), 'r-')
    #plt.figure()
    #plt.plot(new_service_time, at_flex, 'bo',
    #         new_service_time, at_shear, 'rs',
    #         new_service_time, gcov_flexure(new_service_time), 'b-',
    #         new_service_time, gcov_shear(new_service_time), 'r-')
    #plt.show()

    # initial peaks
    corrTarget = np.eye(2)
    cov_variable1 = Normal('cov1', mean=0., stdv=COVT0_COV)
    cov_rv1 = cov_variable1.rv
    cov_variable2 = Normal('cov2', mean=0., stdv=COVT0_COV)
    cov_rv2 = cov_variable2.rv
    rv_array = np.array([cov_rv1, cov_rv2])
    peaks0 = getInitialPeaks(rv_array, NUM_COMPONENT, corrTarget=corrTarget, options={'disp':False, 'anneal':True})
    weights0 = 1./NUM_COMPONENT * np.ones(NUM_COMPONENT)
    covar0 = getInitialCovar(rv_array, NUM_COMPONENT)

    # initial values
    gmdistr = GMM(n_components=NUM_COMPONENT , covariance_type='full')
    gmdistr.weights_ = weights0
    gmdistr.means_ = peaks0
    gmdistr.covars_ = covar0

    # CE-based smpling
    nTimePoint = time_array.size
    #pfpre = np.zeros(nTimePoint)
    #pfmain = np.zeros(nTimePoint)
    Spf2pre = np.zeros(nTimePoint)
    Spf2main = np.zeros(nTimePoint)
    kopt = np.zeros(nTimePoint)
    #pf = np.zeros(nTimePoint)
    #Spf2 = np.zeros(nTimePoint)
    # preliminary sampling
    preSmp = PreSmp(analysis_type, NUM_ADAPTATION, NUM_PRE_SMP, gmdistr, rv_array, r_array, sl_array, rate_array, gd_array, gcov_array)

    return preSmp

#def timeVariantReliability(preSmp, timepoint):
#
#    # preliminary sampling
#    pfpre, Spf2pre = preSmp.adaptation(timepoint)
#    k_array = np.arange(INIT_K, 1.0-K_STEPS, -K_STEPS)
#    # get kopt
#    kopt = preSmp.getKopt(k_array, timepoint)
#    # main sampling
#    mainSmp = MainSmp(preSmp=preSmp)
#    mainSmp.setSmpNum(NUM_MAIN_SMP)
#    smps = mainSmp.sample()
#    pfmain, Spf2main = mainSmp.getPf(smps, timepoint)
#    # combine preliminary and main sampling
#    pf, Spf2, COVpf = combine2smps(pfpre, Spf2pre, pfmain, Spf2main)
#
#    return {'pfpre': pfpre, 'Spf2pre': Spf2pre, 'pfmain': pfmain,
#            'Spf2main': Spf2main, 'pf': pf, 'Spf2': Spf2}

def getLifetimeDistributionParameter(lifetime_name, time_array, pf_array, p0=None):

    if lifetime_name.lower() == 'weibull':
        def weibull_cdf(x, lbda, k):
            return 1.0 - np.exp(-(x/lbda)**k)
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(weibull_cdf, time_array, pf_array, p0 = [200.0, 0.5] if p0 is None else p0, sigma = weights)
        #lifetime_parameter, pcov = curve_fit(weibull_cdf, time_array, pf_array, p0 = [157.9, 159], sigma = weights)
        #lifetime_parameter = np.array([157.9, 159])    # initial value for const resistance
        lifetime_cdf = lambda x: weibull_cdf(x, lifetime_parameter[0], lifetime_parameter[1])
    elif lifetime_name.lower() == 'exponential':
        def exponential_cdf(x, lbda):
            return 1.0 - np.exp(-(x*lbda))
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(exponential_cdf, time_array, pf_array, p0=[0.01] if p0 is None else p0, sigma = weights)
        lifetime_cdf = lambda x: exponential_cdf(x, lifetime_parameter[0])
    elif lifetime_name.lower() == 'exponential_power':
        def power_cdf(x, lbda, k):
            return 1.0 - np.exp(1.0-np.exp(lbda*x**k))
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(power_cdf, time_array, pf_array, p0=[0.01, 1.0] if p0 is None else p0, sigma = weights)
        lifetime_cdf = lambda x: power_cdf(x, lifetime_parameter[0], lifetime_parameter[1])
    elif lifetime_name.lower() == 'gamma':
        def gamma_cdf(x, k, theta):
            return gammainc(k, x/theta)
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(gamma_cdf, time_array, pf_array, p0=p0, sigma=weights)
        lifetime_cdf = lambda x: gamma_cdf(x, lifetime_parameter[0], lifetime_parameter[1])
    elif 'loglog_poly' in lifetime_name.lower():
        deg = int(lifetime_name[-1])
        loglogSt = np.log(-np.log(1.0-pf_array))
        lifetime_parameter = np.polyfit(time_array, loglogSt, deg)
        lifetime_cdf = lambda x: 1.0 - np.exp(-np.exp(np.poly1d(lifetime_parameter)(x)))
    elif  lifetime_name.lower() == 'normal':
        def normal_cdf(x, mu, s):
            return 0.5*(1+erf((x-mu)/(s*np.sqrt(2))))
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(normal_cdf, time_array, pf_array, p0=[150.0, 50.0] if p0 is None else p0, sigma=weights)
        lifetime_cdf = lambda x: normal_cdf(x, lifetime_parameter[0], lifetime_parameter[1])
    elif lifetime_name.lower() == 'weibull_mixture2':
        def weibull_cdf(x, lbda, k1, k2, p1):
            return 1.0 - p1*np.exp(-(x/lbda)**k1) - (1-p1)*np.exp(-(x/lbda)**k2)
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(weibull_cdf, time_array, pf_array, p0 = [200.0, 1, 10, 1./3.] if p0 is None else p0, sigma = weights)
        lifetime_cdf = lambda x: weibull_cdf(x, lifetime_parameter[0], lifetime_parameter[1],
            lifetime_parameter[2], lifetime_parameter[3])
    elif lifetime_name.lower() == 'weibull_mixture3':
        def weibull_cdf(x, lbda, k1, k2, k3, p1, p2):
            return 1.0 - p1*np.exp(-(x/lbda)**k1) - p2*np.exp(-(x/lbda)**k2) - (1-p1-p2)*np.exp(-(x/lbda)**k3)
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(weibull_cdf, time_array, pf_array, p0 = [200.0, 1, 5, 10, 1./8., 1./3.] if p0 is None else p0, sigma = weights)
        lifetime_cdf = lambda x: weibull_cdf(x, lifetime_parameter[0], lifetime_parameter[1],
            lifetime_parameter[2], lifetime_parameter[3], lifetime_parameter[4], lifetime_parameter[5])
    elif lifetime_name.lower() == 'weibull_mixture4':
        def weibull_cdf(x, lbda, k1, k2, k3, k4, p1, p2, p3):
            return 1.0 - p1*np.exp(-(x/lbda)**k1) - p2*np.exp(-(x/lbda)**k2) - \
                p3*np.exp(-(x/lbda)**k3) - (1-p1-p2-p3)*np.exp(-(x/lbda)**k4)
        weights = np.copy(pf_array)
        #weights[0] = 1.0
        lifetime_parameter, pcov = curve_fit(weibull_cdf, time_array, pf_array, p0 = [200.0, 1, 5, 15, 20, 2./3., 1./6., 1./9.] if p0 is None else p0, sigma = weights)
        lifetime_cdf = lambda x: weibull_cdf(x, lifetime_parameter[0], lifetime_parameter[1],
            lifetime_parameter[2], lifetime_parameter[3], lifetime_parameter[4],
            lifetime_parameter[5], lifetime_parameter[6], lifetime_parameter[7])
    else:
        print '[ERROR:] {} illegal lifetime name'.format(lifetime_name.lower())
        sys.exit(1)


    return lifetime_parameter, lifetime_cdf

if __name__ == '__main__':
    from scipy.special import gamma
    import matplotlib.pyplot as plt
    plt.ion()
    ## lower tail data
    #time_array = np.arange(10.,110,10.)
    #time_array = np.insert(time_array[time_array>0], 0, 1.)
    #res = np.load(os.path.join(os.path.abspath('./'), 'data', 'frp_cs', 'shear', 'reliability_parallel_stryr_52.npz'))
    #pf_array = res['pf']

    pf_init = np.array([2.87855688e-11,   2.50411980e-09,   4.62755761e-09,
         6.65924872e-09,   8.47060055e-09,   1.02859066e-08,
         1.20380322e-08,   1.35961518e-08,   1.52916619e-08,
         1.66467175e-08])
    pf_array = np.array([  1.83821390e-08,   1.54197926e-07,   5.14594048e-07,\
        1.68850680e-06,   6.02364972e-06,   2.01626483e-05,\
        6.53350320e-05,   1.85676755e-04,   4.94043808e-04,\
        1.16988724e-03,   2.66841549e-03,   5.93252252e-03,\
        1.31120241e-02,   3.20902966e-02,   8.39448756e-02,\
        2.35869513e-01,   6.13201146e-01,   9.69313898e-01,\
        9.99030612e-01,   9.98048137e-01,   9.99030622e-01])
    time_array = np.array([   1.,   10.,   20.,   30.,   40.,   50.,   60.,   70.,   80.,\
        90.,  100.,  110.,  120.,  130.,  140.,  150.,  160.,  170.,\
        180.,  190.,  200.])
    pf_array = pf_array[:10]
    time_array = time_array[:10]

    ## with initial <1yr data
    #pf_array = np.hstack((pf_init, pf_array[:18]))
    #time_array = np.hstack((np.linspace(0.001, 0.9, 10), time_array[:18]))

    ## data with constant resistance
    #pf_array = np.array([1.59650071e-13,  1.02051700e-12,  8.12105938e-12,
    #6.49772458e-11,  5.54705282e-10,  4.85687046e-09,  4.34790219e-08,
    #3.93959657e-07,  3.57256801e-06,  3.20011154e-05,  2.78885682e-04,
    #2.32281600e-03,  1.80321113e-02,  1.22653748e-01,  5.70597288e-01,
    #9.91690383e-01])
    #time_array = np.arange(130, 161, 2.)

    plt.close('all')
    plt.rc('font', family='serif', size=12)

    #lifetime_name_list = ['weibull', 'exponential', 'exponential_power', 'gamma', 'loglog_poly1', 'loglog_poly5', 'normal', 'weibull_mixture2', 'weibull_mixture3']
    #lifetime_name_list = ['loglog_poly1', 'loglog_poly4', 'normal', 'weibull_mixture2', 'weibull_mixture3', 'weibull_mixture4']
    lifetime_name_list = ['weibull_mixture2', 'weibull_mixture3']
    for lifetime_name in lifetime_name_list:
        lifetime_parameter, lifetime_cdf = getLifetimeDistributionParameter(lifetime_name, time_array, pf_array)
        if lifetime_name == 'weibull_mixture2':
            lbda = lifetime_parameter[0]
            k1 = lifetime_parameter[1]
            k2 = lifetime_parameter[2]
            p1 = lifetime_parameter[-1]
            p2 = 1.0-p1
            mean_life = lbda*(p1*gamma(1+1/k1)+p2*gamma(1+1/k2))
            print '{} parameters:'.format(lifetime_name)
            print lifetime_parameter
            print mean_life
        elif lifetime_name == 'weibull_mixture3':
            lbda = lifetime_parameter[0]
            k1 = lifetime_parameter[1]
            k2 = lifetime_parameter[2]
            k3 = lifetime_parameter[3]
            p1 = lifetime_parameter[-2]
            p2 = lifetime_parameter[-1]
            p3 = 1.0-p1-p2
            p2 = 1.0-p1
            mean_life = lbda*(p1*gamma(1+1/k1)+p2*gamma(1+1/k2)+p3*gamma(1+1/k3))
            print '{} parameters:'.format(lifetime_name)
            print lifetime_parameter
            print mean_life


        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        ax1.semilogy(time_array, pf_array, 'o')
        ax1.semilogy(time_array, lifetime_cdf(time_array), 'r-')
        relative_error = (lifetime_cdf(time_array) - pf_array) / pf_array
        ax2.stem(time_array, relative_error, linefmt='b-', markerfmt='bo', basefmt='r-')
        ax2.set_xlabel(r'service time (yr)')
        ax1.set_ylabel('CDF in log')
        ax2.set_ylabel('relative error')

        plt.show()
