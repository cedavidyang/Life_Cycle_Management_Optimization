# postprocessing after optimization
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from multiprocessing import Pool
import itertools

import time
import datetime

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from management.performanceFuncs import performanceHistory

def getPf(param_list):
    pf_list = performanceHistory(*param_list)
    return pf_list

def main():
    np.random.seed(64)

    str_yr = [0, 0, 0]
    icorr_mean_list=[1,1,1]

    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.insert(time_array[time_array>0], 0, 1.)
    nTimePoint = time_array.shape[0]

    iter_param = itertools.izip(itertools.repeat(str_yr), itertools.repeat(icorr_mean_list),
            time_array, itertools.repeat(False))

    print "parallel computing begins"
    start_delta_time = time.time()

    # calculate pfs
    pool = Pool(processes=3)
    res = pool.map_async(getPf, iter_param).get(0xFFFF)
    pool.close()
    pool.join()

    #res = map(getPf,iter_param)

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))

    #save data
    datapath = os.path.join(os.path.abspath('./'), 'data', 'simpleCorrosion')
    filename = 'pfhistory_str_'
    for ti in str_yr:
        filename = filename + str(int(ti))
    datafile = os.path.join(datapath, filename)
    np.save(datafile, res)

    return res


if __name__ == '__main__':
    pfs = main()
