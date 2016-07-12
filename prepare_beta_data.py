import itertools
import math
import operator
import random
import numpy as np
from multiprocessing import Pool

import time
import datetime
import os
import sys

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
from management.performanceFuncs import pointintimeFunc, generateCompData
from management.component import Component
from management.system import System

def parallelbetaflex(timecombo):
    stryr = timecombo[0]
    year = timecombo[1]
    pf, cost = generateCompData('flex', stryr, year)
    return pf, cost

def parallelbetashear(timecombo):
    stryr = timecombo[0]
    year = timecombo[1]
    pf, cost = generateCompData('shear', stryr, year)
    return pf, cost

def parallelbetadeck(timecombo):
    stryr = timecombo[0]
    year = timecombo[1]
    pf, cost = generateCompData('deck', stryr, year)
    return pf, cost

if __name__ == '__main__':

    icorr_mean_list = [1.,1.,1.]
    life = 100
    num_processes = 20


    flexdata = np.zeros((life,life+2))
    sheardata = np.zeros((life,life+2))
    deckdata = np.zeros((life,life+2))

    lifetime = np.arange(0, life+1)
    strtime = np.arange(0, life)

    print "test"
    start_delta_time = time.time()

    # parallel edition
    timecombo = np.transpose(np.array(np.meshgrid(strtime, lifetime)))
    timecombo = np.reshape(timecombo, ((life)*(life+1),2))
    pool = Pool(processes=num_processes)
    resflex = pool.map_async(parallelbetaflex, timecombo).get(0xFFFF)
    resshear = pool.map_async(parallelbetashear, timecombo).get(0xFFFF)
    resdeck = pool.map_async(parallelbetadeck, timecombo).get(0xFFFF)
    pool.close()
    pool.join()
    for data, res in zip([flexdata, sheardata, deckdata], [resflex, resshear, resdeck]):
        pfdata = np.array(res)[:,0]
        costdata = np.array(res)[:,1]
        pfdata = np.reshape(pfdata, (life, life+1))
        costdata = np.reshape(costdata, (life, life+1))
        data[:,:-1] = pfdata
        data[:,-1] = costdata[:,-1]

    # # non-parallel edition
    # for i, stryr in enumerate(strtime):
        # # girder flexure
        # pf1, cost1 = generateCompData('flex', stryr, 0)
        # flexdata[i,-1] = cost1
        # # girder shear
        # pf2, cost2 = generateCompData('shear', stryr, 0)
        # sheardata[i,-1] = cost2
        # # deck
        # pf3, cost3 = generateCompData('deck', stryr, 0)
        # deckdata[i,-1] = cost3
        # for j, year in enumerate(lifetime):
            # # girder flexure
            # pf1, cost1 = generateCompData('flex', stryr, year)
            # flexdata[i,j] = pf1
            # # girder shear
            # pf2, cost2 = generateCompData('shear', stryr, year)
            # sheardata[i,j] = pf2
            # # deck
            # pf3, cost3 = generateCompData('deck', stryr, year)
            # deckdata[i,j] = pf3

    # save data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    np.savez(os.path.join(datapath,'beta_data'), flex=flexdata, shear=sheardata, deck=deckdata)

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))
