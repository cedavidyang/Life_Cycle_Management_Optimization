# Life-Cylce Management of FRP strengthening
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import scipy.stats as stats
from multiprocessing import Pool

import time
import datetime

from constants import *
from pyre.distributions import *
from management import *
from pyCEsmp import *

type_of_component = raw_input('type of component (system):')

RC_CS_FLEX_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'flexure')
flex_datafile = os.path.join(RC_CS_FLEX_PATH, 'LWS_results.txt')
RC_CS_DECK_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'deck')
deck_datafile = os.path.join(RC_CS_DECK_PATH, 'LWS_results.txt')
RC_CS_SHEAR_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'shear')
shear_datafile = os.path.join(RC_CS_SHEAR_PATH, 'LWS_results.txt')

# characteristic strengthening time
strengthen_time = np.arange(FRP_DESIGN_YR, SERVICE_LIFE, 2.)
new_service_time_list = getServiceTimeFromFile(shear_datafile, strengthen_time)

if 'f' in type_of_component.lower():
    # get corrosion loss
    flex_loss_dict = getCorrosionLossFromFile(flex_datafile, strengthen_time)
    # get ratio lists
    flex_mean_strength_array, tfrp_array = strengtheningFromCorrosionLoss(flex_loss_dict, 'flexure')
    # get gt and at histories
    flex_mean_history_array, flex_cov_history_array = getDeteriorationFromFile(flex_datafile, strengthen_time, flex_mean_strength_array, 'flexure')
    shear_mean_history_array = np.copy(flex_mean_history_array)
    shear_cov_history_array = np.copy(flex_cov_history_array)
elif 's' in type_of_component.lower():
    # get corrosion loss
    shear_loss_dict = getCorrosionLossFromFile(shear_datafile, strengthen_time)
    # get ratio lists
    shear_mean_strength_array, tfrp_array = strengtheningFromCorrosionLoss(shear_loss_dict, 'shear')
    # get gt and at histories
    shear_mean_history_array, shear_cov_history_array = getDeteriorationFromFile(shear_datafile, strengthen_time, shear_mean_strength_array, 'shear')
    flex_mean_history_array = np.copy(shear_mean_history_array)
    flex_cov_history_array = np.copy(shear_cov_history_array)
elif 'd' in type_of_component.lower():
    # get corrosion loss
    deck_loss_dict = getCorrosionLossFromFile(deck_datafile, strengthen_time)
    # get ratio lists
    deck_mean_strength_array, tfrp_array = strengtheningFromCorrosionLoss(deck_loss_dict, 'deck')
    # get gt and at histories
    flex_mean_history_array, flex_cov_history_array = getDeteriorationFromFile(deck_datafile, strengthen_time, deck_mean_strength_array, 'deck')
    shear_mean_history_array = np.copy(flex_mean_history_array)
    shear_cov_history_array = np.copy(flex_cov_history_array)

# get time-variant failure probability
time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
time_array = np.insert(time_array[time_array>0], 0, 1.)
nTimePoint = time_array.shape[0]

pfpre = np.zeros(nTimePoint)
pfmain = np.zeros(nTimePoint)
Spf2pre = np.zeros(nTimePoint)
Spf2main = np.zeros(nTimePoint)
kopt = np.zeros(nTimePoint)
pf = np.zeros(nTimePoint)
Spf2 = np.zeros(nTimePoint)

if 'f' in type_of_component.lower():
    subfolder = 'flexure'
elif 's' in type_of_component.lower():
    subfolder = 'shear'
elif 'd' in type_of_component.lower():
    subfolder = 'deck'
else:
    print '[Error:] illegal type of component, must be flexure, shear or deck'
    sys.exit(1)
datapath = os.path.join(os.path.abspath('./'), 'data', 'frp_cs', subfolder)

for str_yr, new_service_time,\
    flex_mean_history, shear_mean_history, flex_cov_history, shear_cov_history\
    in zip(strengthen_time, new_service_time_list, flex_mean_history_array, shear_mean_history_array,
        flex_cov_history_array, shear_cov_history_array):

    filename = 'reliability_parallel_stryr_'+str(int(str_yr))+'.npz'
    datafile = os.path.join(datapath, filename)

    frp_mean_history = np.array([flex_mean_history, shear_mean_history])
    frp_cov_history = np.array([flex_cov_history, shear_cov_history])

    preSmp = getPreSmpObject(type_of_component, new_service_time, frp_mean_history, frp_cov_history)

    if os.path.isfile(datafile):
        print 'Reliability for strengthening at year '+str(int(str_yr))
        result_dict = np.load(datafile)
        printResults(time_array, result_dict)
        continue

    def timeVariantReliability(timepoint):
        # preliminary sampling
        pfpre, Spf2pre = preSmp.adaptation(timepoint)
        k_array = np.arange(INIT_K, 1.0-K_STEPS, -K_STEPS)
        # get kopt
        kopt = preSmp.getKopt(k_array, timepoint)
        # main sampling
        mainSmp = MainSmp(preSmp=preSmp)
        mainSmp.setSmpNum(NUM_MAIN_SMP)
        smps = mainSmp.sample()
        pfmain, Spf2main = mainSmp.getPf(smps, timepoint)
        # combine preliminary and main sampling
        pf, Spf2, COVpf = combine2smps(pfpre, Spf2pre, pfmain, Spf2main)
        return {'pfpre': pfpre, 'Spf2pre': Spf2pre, 'pfmain': pfmain,
                'Spf2main': Spf2main, 'pf': pf, 'Spf2': Spf2}

    print 'CALC: reliability (parallel version) for strengthening at year '+str(int(str_yr))
    try:
        start_delta_time = time.time()
        pool = Pool(processes=3)
        res = pool.map_async(timeVariantReliability, time_array).get(0xFFFF)
        pool.close()
        pool.join()
        for indx, ti in enumerate(time_array):
            pfpre[indx] = res[indx]['pfpre']
            Spf2pre[indx] = res[indx]['Spf2pre']
            pfmain[indx] = res[indx]['pfmain']
            Spf2main[indx] = res[indx]['Spf2main']
            pf[indx] = res[indx]['pf']
            Spf2[indx] = res[indx]['Spf2']
        delta_time = time.time() - start_delta_time
        print 'DONE: Parallel version',str(datetime.timedelta(seconds=delta_time))
        np.savez(datafile, pfpre=pfpre, Spf2pre=Spf2pre,pfmain=pfmain, Spf2main=Spf2main, pf=pf, Spf2=Spf2)
        result_dict = np.load(datafile)
        printResults(time_array, result_dict)
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()

# fit distribution
lifetime_parameter_list = []
plt.close('all')
plt.rc('font', family='serif', size=12)

for str_yr, new_service_time,\
    flex_mean_history, shear_mean_history, flex_cov_history, shear_cov_history\
    in zip(strengthen_time, new_service_time_list, flex_mean_history_array, shear_mean_history_array,
        flex_cov_history_array, shear_cov_history_array):

    filename = 'reliability_parallel_stryr_'+str(int(str_yr))+'.npz'
    datafile = os.path.join(datapath, filename)

    if os.path.isfile(datafile):
        res = np.load(datafile)
        pf_array = res['pf']

        #lifetime_name_list = ['weibull', 'gamma', 'loglog_poly5', 'normal', 'weibull_mixture3']
        #lifetime_name_list = ['weibull_mixture2', 'weibull_mixture3']
        lifetime_name_list = ['weibull_mixture3']
        for lifetime_name in lifetime_name_list:
            lifetime_parameter, lifetime_cdf = getLifetimeDistributionParameter(lifetime_name, time_array, pf_array)
            print '{} fitting: lifetime distribution parameters:\n {}'.format(lifetime_name, lifetime_parameter)

            f, (ax1, ax2) = plt.subplots(2, sharex=True)
            ax1.semilogy(time_array, pf_array, 'o')
            ax1.semilogy(time_array, lifetime_cdf(time_array), 'r-')
            relative_error = (lifetime_cdf(time_array) - pf_array) / pf_array
            ax2.stem(time_array, relative_error, linefmt='b-', markerfmt='bo', basefmt='r-')

            plt.show()
        wait = raw_input('strengthen year {}: press to continue'.format(str_yr))


if __name__ == '__main__':
    pass
