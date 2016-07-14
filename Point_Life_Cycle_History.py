# get the pf history for all the components and the system
import os
import sys
import numpy as np
from multiprocessing import Pool
import itertools

import time
import datetime

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from management.performanceFuncs import pointintimeHistory

str_yr = [0, 0, 0]
icorr_mean_list=[1,1,1]
nprocess = 2

def getPf(param_list):
    pf_list = pointintimeHistory(*param_list)
    return pf_list

def main():
    np.random.seed(64)

    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.append(time_array, str_yr)
    time_array = np.insert(time_array[time_array>0], 0, 1.)
    time_array = np.sort(time_array)
    time_array[time_array==0] = 1.
    time_array = np.unique(time_array)
    nTimePoint = time_array.shape[0]

    iter_param = itertools.izip(itertools.repeat(str_yr), itertools.repeat(icorr_mean_list),
            time_array, itertools.repeat(False))

    print "parallel computing begins"
    start_delta_time = time.time()

    # calculate pfs
    pool = Pool(processes=nprocess)
    res = pool.map_async(getPf, iter_param).get(0xFFFF)
    pool.close()
    pool.join()

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))

    return time_array, res


if __name__ == '__main__':
    time_array, pfs = main()
    pfs = np.array(pfs)

    # save data
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
    suffix = rate2suffix(icorr_mean_list)
    filename = 'pointpfhistory_str_'
    for ti in str_yr:
        filename = filename + str(int(ti)) + '_'
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename = filename+suffix+'.npz'
    datafile = os.path.join(datapath,filename)

    np.savez(datafile, time=time_array, system=pfs[:,0], flexure=pfs[:,1],
            shear=pfs[:,2], deck=pfs[:,3])
