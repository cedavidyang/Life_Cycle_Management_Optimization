# Life-Cycle reliability
import os
import sys
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
plt.ion()
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from sklearn.mixture.gmm import GMM

from multiprocessing import Pool

import time
import datetime

from constants import *
from pyre.distributions import *
from pyCEsmp import *

type_of_component = raw_input('type of component (system):')

# reliability output folder
print 'type of component: '+type_of_component
if 's' in type_of_component.lower():
    type_of_component = 'shear'
elif 'f' in type_of_component.lower():
    type_of_component = 'flexure'
elif 'd' in type_of_component.lower():
    type_of_component = 'deck'
else:
    print 'ERROR: illegal type of component, must be flexure, shear or deck'
    sys.exit(1)
if 'd' in type_of_component.lower():
    FLEXTURE_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'deck')
else:
    FLEXTURE_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'flexure')
SHEAR_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'shear')
RE_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', type_of_component)

# read degradation data
flexure_datafile = os.path.join(FLEXTURE_DATAFILE_PATH, 'LWS_results.txt')
shear_datafile = os.path.join(SHEAR_DATAFILE_PATH, 'LWS_results.txt')
service_time = np.loadtxt(flexure_datafile)[0,:]
#end_indx = np.where(service_time == SERVICE_LIFE)[0][0]+1
end_indx = service_time.size
service_time = service_time[:end_indx]
# flexural strength rc
flexure_mean_history = np.loadtxt(flexure_datafile)[20,:][:end_indx]
flexure_std_history = np.loadtxt(flexure_datafile)[21,:][:end_indx]
flexure_cov_history = flexure_std_history/flexure_mean_history
# shear strength rc
shear_mean_history = np.loadtxt(shear_datafile)[22,:][:end_indx]
shear_std_history = np.loadtxt(shear_datafile)[23,:][:end_indx]
shear_cov_history = shear_std_history/shear_mean_history
# assume flexural and shear strength follow lognormal
rm_mean = flexure_mean_history[0]
rm_stdv = flexure_std_history[0]
rv_mean = shear_mean_history[0]
rv_stdv = shear_std_history[0]

# resistance array and live load array
r_array = np.array([rm_mean, rv_mean])
slmRv = Normal('slm', mean=M_LLIM_MEAN, stdv=M_LLIM_MEAN*M_LLIM_COV)
slvRv = Normal('slv', mean=V_LLIM_MEAN, stdv=V_LLIM_MEAN*V_LLIM_COV)
sldRv = Normal('sld', mean=M_LLIM_DECK_MEAN, stdv=M_LLIM_DECK_MEAN*M_LLIM_DECK_COV)
slmDistr = slmRv.rv
slvDistr = slvRv.rv
sldDistr = sldRv.rv
sl_array = np.array([slmDistr, slvDistr, sldDistr])

# live load arrival rate
rate_m = LL_ARRIVAL_RATE
rate_v = LL_ARRIVAL_RATE
rate_d = LL_ARRIVAL_RATE_DECK
rate_array = np.array([rate_m, rate_v, rate_d])

# mean degradation and cov evolution
def gd_func(x, a, b, c, d, e):
    return a +b*x**1+c*x**2+d*x**3+e*x**4
gt_flex = flexure_mean_history/flexure_mean_history[0] / np.sqrt(1+flexure_cov_history**2)
gt_shear = shear_mean_history/shear_mean_history[0] / np.sqrt(1+shear_cov_history**2)
gd_flex_popt, pcov = curve_fit(gd_func, service_time, gt_flex)
gd_shear_popt,pcov = curve_fit(gd_func, service_time, gt_shear)
gd_flexure = lambda x: gd_func(x, gd_flex_popt[0], gd_flex_popt[1], gd_flex_popt[2], gd_flex_popt[3], gd_flex_popt[4])
gd_shear = lambda x: gd_func(x, gd_shear_popt[0], gd_shear_popt[1], gd_shear_popt[2], gd_shear_popt[3], gd_shear_popt[4])

def gcov_func(x, a, b, c, d, e):
    return a +b*x+c*x**2+d*x**3+e*x**4
at_flex = np.sqrt(np.log(flexure_cov_history**2+1)/COVT0_COV**2)
at_shear = np.sqrt(np.log(shear_cov_history**2+1)/COVT0_COV**2)
gcov_flex_popt, pcov = curve_fit(gcov_func, service_time, at_flex)
gcov_shear_popt,pcov = curve_fit(gcov_func, service_time, at_shear)
gcov_flexure = lambda x: gcov_func(x, gcov_flex_popt[0], gcov_flex_popt[1], gcov_flex_popt[2], gcov_flex_popt[3], gcov_flex_popt[4])
gcov_shear = lambda x: gcov_func(x, gcov_shear_popt[0], gcov_shear_popt[1], gcov_shear_popt[2], gcov_shear_popt[3], gcov_shear_popt[4])

gd_array = np.array([gd_flexure, gd_shear])
gcov_array = np.array([gcov_flexure, gcov_shear])

# compare
nb = raw_input('if continue?')
if nb == 'y' or nb == 'Y':
    print 'continue analysis'
else:
    plt.plot(service_time, gt_flex, 'bo',
            service_time, gt_shear, 'rs',
            service_time, gd_flexure(service_time), 'b-',
            service_time, gd_shear(service_time), 'r-')
    plt.figure()
    plt.plot(service_time, at_flex, 'bo',
            service_time, at_shear, 'rs',
            service_time, gcov_flexure(service_time), 'b-',
            service_time, gcov_shear(service_time), 'r-')
    plt.show()
    sys.exit(1)

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

#preSmp = PreSmp(NUM_ADAPTATION, 1, gmdistr, rv_array, r_array, sl_array, gd_array, gcov_array)
#time_array = np.linspace(0, 200, 101)
#pf_with_const_resistance = []
#for ti in time_array:
#    smps = np.array([[0.0, 0.0]])
#    pf = 1.0 - preSmp.condAvailability(smps, ti)
#    pf_with_const_resistance.append(pf[0])
#print 'result assuming constant shear resistance:'
#print np.array(pf_with_const_resistance)

def timeVariantReliability(ti):
    # CE-based smpling
    # preliminary sampling
    preSmp = PreSmp(type_of_component, NUM_ADAPTATION, NUM_PRE_SMP, gmdistr, rv_array, r_array, sl_array, rate_array, gd_array, gcov_array)
    pfpre, Spf2pre = preSmp.adaptation(ti)
    k_array = np.arange(INIT_K, 1.0-K_STEPS, -K_STEPS)
    kopt = preSmp.getKopt(k_array, ti)
    # main sampling
    mainSmp = MainSmp(preSmp=preSmp)
    mainSmp.setSmpNum(NUM_MAIN_SMP)
    smps = mainSmp.sample()
    pfmain, Spf2main = mainSmp.getPf(smps, ti)

    ## save sampling functions
    #optimalSmpFunc(rv_array, mainSmp, ti)

    pf, Spf2, COVpf = combine2smps(pfpre, Spf2pre, pfmain, Spf2main)

    return {'pfpre': pfpre, 'Spf2pre': Spf2pre, 'pfmain': pfmain,
            'Spf2main': Spf2main, 'pf': pf, 'Spf2': Spf2}
    # serial version return
    #return {'pfpre': pfpre, 'Spf2pre': Spf2pre, 'pfmain': pfmain,
    #        'Spf2main': Spf2main, 'pf': pf, 'Spf2': Spf2,
    #        'preSmp': preSmp, 'mainSmp': mainSmp}

tmp = timeVariantReliability(10)

def main():

    # time array
    #time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.arange(RELIABILITY_DT,END_AGE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.insert(time_array, 0, 1.)
    nTimePoint = time_array.shape[0]

    # initial values
    pfpre = np.zeros(nTimePoint)
    pfmain = np.zeros(nTimePoint)
    Spf2pre = np.zeros(nTimePoint)
    Spf2main = np.zeros(nTimePoint)
    kopt = np.zeros(nTimePoint)
    pf = np.zeros(nTimePoint)
    Spf2 = np.zeros(nTimePoint)
    smps_func_array = np.empty(shape=pf.shape, dtype=object)

    #print 'CALC: Serial version'
    #start_delta_time = time.time()
    ## serial version
    #for ti in time_array:
    #    print('calculating time year {} ...'.format(ti))
    #    # CE-based smpling
    #    indx = np.where(time_array==ti)[0][0]
    #    #pfpre[indx], Spf2pre[indx], pfmain[indx], Spf2main[indx], pf[indx], Spf2[indx] = timeVariantReliability(ti)
    #    res = timeVariantReliability(ti)
    #    smps_func_array[indx] = res['mainSmp']
    #    pfpre[indx] = res['pfpre']
    #    Spf2pre[indx] = res['Spf2pre']
    #    pfmain[indx] = res['pfmain']
    #    Spf2main[indx] = res['Spf2main']
    #    pf[indx] = res['pf']
    #    Spf2[indx] = res['Spf2']
    #np.savez('results_serial.npz', pfpre=pfpre, Spf2pre=Spf2pre,pfmain=pfmain, Spf2main=Spf2main, pf=pf, Spf2=Spf2)
    ##np.savez('sample_func_serial.npz', smps_func = smps_func_array)
    #result_dict = np.load('results_serial.npz')
    #printResults(time_array, result_dict)
    #delta_time = time.time() - start_delta_time
    #print 'DONE: Serial version',str(datetime.timedelta(seconds=delta_time))
    ### display the shape of optimal sampling function at year 80
    ##optimalSmpFunc(rv_array, smps_func_array[0], np.array([80.0]))

    print 'CALC: Parallel version'
    try:
        start_delta_time = time.time()
        pool = Pool(processes=3)
        res = pool.map_async(timeVariantReliability, time_array).get(0xFFFF)
        pool.close()
        pool.join()
        for ti in time_array:
            indx = np.where(time_array==ti)[0][0]
            pfpre[indx] = res[indx]['pfpre']
            Spf2pre[indx] = res[indx]['Spf2pre']
            pfmain[indx] = res[indx]['pfmain']
            Spf2main[indx] = res[indx]['Spf2main']
            pf[indx] = res[indx]['pf']
            Spf2[indx] = res[indx]['Spf2']
        save_file = os.path.join(RE_DATAFILE_PATH, 'reliability_results_parallel.npz')
        np.savez(save_file, pfpre=pfpre, Spf2pre=Spf2pre,pfmain=pfmain, Spf2main=Spf2main, pf=pf, Spf2=Spf2)
        result_dict = np.load(save_file)
        printResults(time_array, result_dict)
        delta_time = time.time() - start_delta_time
        print 'DONE: Parallel version',str(datetime.timedelta(seconds=delta_time))
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()

if __name__ == '__main__':
    main()
