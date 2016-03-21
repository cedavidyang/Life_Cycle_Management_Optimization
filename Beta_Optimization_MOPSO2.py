import itertools
import math
import operator
import random
import numpy as np
from multiprocessing import Pool, Manager

import time
import datetime
import os
import sys

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
from management.performanceFuncs import pointintimeFunc
from management.component import Component
from management.system import System

from deap import algorithms
from deap import base
from deap import creator
from deap import tools


NPARAM = 3
INRTLMAX = 0.8
INRTLMIN = 0.8
PCONF = 1.5
SCONF = 1.5
PMIN = 0.0
PMAX = 100.0
SMIN = -0.2
SMAX = 0.2
NREP = 500
NDIV = 50
NPOP = 500
# stop criteria
NGEN = 10
NMAX = 100
TOL1 = 0.001
TOL2 = 0.001
NCR = 10

#icorr_mean_list = np.array(input('corrosoin rate:')).astype('double')
#year = np.array(input('expected life:')).astype('double')
#num_processes = np.array(input('number of processes:')).astype('int')
icorr_mean_list = [1.,1.,1.]
year = 100
num_processes = 40

creator.create("FitnessMulti", base.Fitness, weights=(-1.0,-1.0))
creator.create("Particle", list, fitness=creator.FitnessMulti, speed=list,
        pmin=None, pmax=None, smin=None, smax=None, geoinfo=None,
        best=None)
creator.create("Swarm", list, gbest=None, gbestfit=list)

def init_particle(pcls, size, pmin, pmax, smin, smax, geoinfo):
    part = pcls(random.randint(pmin, pmax) for _ in xrange(size))
    spacing = pmax-pmin
    part.speed = [random.uniform(smin*spacing, smax*spacing) for _ in xrange(size)]
    part.pmin = pmin
    part.pmax = pmax
    part.smin = smin
    part.smax = smax
    part.geoinfo = geoinfo
    return part

def set_geoinfo(part, indx, subindx):
    part.geoinfo = (indx,subindx)

def update_particle(part, best, w, phi1, phi2):
    tmpu1 = random.uniform(0,phi1)
    tmpu2 = random.uniform(0,phi2)
    u1 = [tmpu1 for _ in range(len(part))]
    u2 = [tmpu2 for _ in range(len(part))]
    w = [w for _ in range(len(part))]
    v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
    v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
    part.speed = list(map(operator.add, map(operator.mul, w, part.speed), map(operator.add, v_u1, v_u2)))
    newpos = map(operator.add, part, part.speed)
    newposndarray = np.round(newpos)
    newposndarray = newposndarray.astype(int)
    newpos = np.ndarray.tolist(newposndarray)
    for i, pos in enumerate(newpos):
        # stay in the domain and reverse velocity if hit the boundaries
        if pos < part.pmin:
            newpos[i] = part.pmin
            #part.speed[i] = -part.speed[i]
        elif pos > part.pmax:
            newpos[i] = part.pmax
            #part.speed[i] = -part.speed[i]
    for i, speed in enumerate(part.speed):
        if speed < part.smin*(part.pmax-part.pmin):
            part.speed[i] = part.smin*(part.pmax-part.pmin)
        elif speed > part.smax*(part.pmax-part.pmin):
            part.speed[i] = part.smax*(part.pmax-part.pmin)
    #part[:] = list(map(operator.add, part, part.speed))
    part[:] = newpos

def select_leader(gbest):
    crowdingDist = np.array([part.fitness.crowding_dist for part in gbest])
    if all(np.isinf(crowdingDist)):
        crowdingDist = np.ones(crowdingDist.shape)
    else:
        maxdist = np.max(crowdingDist[np.logical_not(np.isinf(crowdingDist))])
        crowdingDist[np.isinf(crowdingDist)] = maxdist
    rouletteWheel = np.cumsum(crowdingDist, axis = 0)
    finalIdx = np.where(((rouletteWheel.max() - rouletteWheel.min())*random.random()) <= rouletteWheel)[0]
    h = finalIdx[0]

    return gbest[h]

def unique_rows(A, return_index, return_inverse):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    Returns I if return_index is True
    Returns J if return_inverse is True
    """
    A = np.require(A, requirements='C')
    assert A.ndim == 2, "array must be 2-dim'l"

    B = np.unique(A.view([('', A.dtype)]*A.shape[1]),
            return_index=return_index,
            return_inverse=return_inverse)

    if return_index or return_inverse:
        return (B[0].view(A.dtype).reshape((-1, A.shape[1]), order='C'),) \
            + B[1:]
    else:
        return B.view(A.dtype).reshape((-1, A.shape[1]), order='C')

toolbox = base.Toolbox()
toolbox.register("particle", init_particle, creator.Particle, size=NPARAM, pmin=PMIN, pmax=PMAX,
        smin=SMIN, smax=SMAX, geoinfo=None)
toolbox.register("swarm", tools.initRepeat, creator.Swarm, toolbox.particle)
toolbox.register("update", update_particle, phi1=PCONF, phi2=SCONF)
toolbox.register("evaluate", pointintimeFunc, icorr_mean_list=icorr_mean_list, year=year)
toolbox.register("sort", tools.sortNondominated)
toolbox.register("assign_crowding_dist", tools.emo.assignCrowdingDist)
toolbox.register("select_leader", select_leader)
toolbox.register("unique_rows", unique_rows, return_index=False, return_inverse=False)


def main():
    # reset bookkeeping
    System.resetBookKeeping()

    ## use existing bookkeeping data
    #bookkeeping = np.load('bookkeeping.npz')

    manager = Manager()
    System.bookkeeping = manager.dict(System.bookkeeping)

    pool = Pool(processes=num_processes)
    toolbox.register("map", pool.map)

    #System.bookkeeping = dict(System.bookkeeping)
    #toolbox.register("map", map)

    print "MULTIOBJECTIVE OPTIMIZATION: parallel version"
    start_delta_time = time.time()

    # optimization
    random.seed(64)

    logbook = tools.Logbook()
    logbook.header = ["gen", "evals", "nfront", "mean1", "mean2", "tol1", "tol2", "time"]

    # initialization: first generation
    swarm = toolbox.swarm(n=NPOP)
    fitnessvalues = toolbox.map(toolbox.evaluate,swarm)
    for i,part in enumerate(swarm):
        part.fitness.values = fitnessvalues[i]
        part.best = toolbox.clone(part)
    firstfront = toolbox.sort(swarm, NPOP, first_front_only=True)[0]
    swarm.gbest = []
    swarm.gbestfit = []
    for i,part in enumerate(firstfront):
        swarm.gbest.append(toolbox.clone(part))
        swarm.gbestfit.append(part.fitness.values)
    delta_time = time.time() - start_delta_time
    logbook.record(gen=1, evals=NPOP,nfront=len(swarm.gbest),
            mean1='nan', mean2='nan', tol1=TOL1, tol2=TOL2,time=delta_time)
    print(logbook.stream)

    gbestfitlast = swarm.gbestfit
    evolStop = False
    # start the evolution
    g = 2
    distances1 = []
    distances2 = []
    frontfitlast = np.array([[1., 0.]])
    evolStop = False
    while not evolStop:
        toolbox.assign_crowding_dist(swarm.gbest)
        for part in swarm:
            leader = toolbox.select_leader(swarm.gbest)
            w = INRTLMAX - ((INRTLMAX - INRTLMIN)/NMAX)*g
            toolbox.update(part, leader, w)
        fitnessvalues = toolbox.map(toolbox.evaluate,swarm)
        for i,part in enumerate(swarm):
            part.fitness.values = fitnessvalues[i]
            # update personal best
            if part.fitness.dominates(part.best.fitness):
                part.best = toolbox.clone(part)
            elif part.best.fitness.dominates(part.fitness):
                pass
            else:
                if (random.randint(0,1) == 0):
                    part.best = toolbox.clone(part)

        # update external archive
        firstfront = toolbox.sort(swarm, NPOP, first_front_only=True)[0]
        firstfront = toolbox.sort(swarm.gbest+firstfront, NREP+NPOP, first_front_only=True)[0]
        RepX = np.array([part[:] for part in firstfront])
        uniqueIdx = np.sort(toolbox.unique_rows(RepX, return_index=True)[1])
        swarm.gbest = []
        for idx in uniqueIdx:
            swarm.gbest.append(firstfront[idx])
        swarm.gbestfit = [part.fitness.values for part in swarm.gbest]

        #:::: If the external population has reached its maximum allowable ::::#
        if (len(swarm.gbest) > NREP):
            toolbox.assign_crowding_dist(swarm.gbest)
            crowding_dist = np.array([part.fitness.crowding_dist for part in swarm.gbest])
            kept_indx = np.argsort(crowding_dist)[-NREP:]
            gbest = []
            for indx in kept_index:
                gbest.append(swarm.gbest[idx])
            swarm.gbest = gbest

        # check if stop evolution
        distance1=[]
        distance2=[]
        frontfit = np.array(swarm.gbestfit)
        frontfit[frontfit[:,1] == 0,1] = 1.
        frontfitlast[frontfitlast[:,1] == 0,1] = 1.
        for obj in frontfit:
            distance1.append(min(np.abs(np.log10(frontfitlast[:,0])-np.log10(obj[0]))))
            #distance2.append(min(np.abs(frontfitlast[:,1]-obj[1])))
            distance2.append(min(np.abs(np.log10(frontfitlast[:,1])-np.log10(obj[1]))))
        distances1.append(np.mean(distance1))
        distances2.append(np.mean(distance2))
        longest1 = 0.
        longest2 = 0.
        for point1 in frontfit:
            for point2 in frontfit:
                dist1 = np.abs(np.log10(point1[0])-np.log10(point2[0]))
                #dist2 = np.abs(point1[1]-point2[1])
                dist2 = np.abs(np.log10(point1[1])-np.log10(point2[1]))
                if dist1 > longest1:
                    longest1 = dist1
                if dist2 > longest2:
                    longest2 = dist2
        tol1 = np.maximum(longest1/NPOP,TOL1)
        tol2 = np.maximum(longest2/NPOP,TOL2)
        evolStop = (len(distances1)>NGEN and all(np.array(distances1[-NGEN:])<tol1)
            and all(np.array(distances2[-NGEN:])<tol2)) or g>NMAX
        frontfitlast = frontfit

        # Gather all the fitnesses in one list and print the stats
        #record = stats.compile(pop)
        delta_time = time.time() - start_delta_time
        logbook.record(gen=g, evals=NPOP,nfront=len(swarm.gbest),
                mean1=distances1[-1], mean2=distances2[-1],
                tol1=tol1, tol2=tol2, time=delta_time)
        print(logbook.stream)

        g+=1

    pool.close()
    pool.join()

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))


    return swarm, logbook

if __name__ == "__main__":

    swarm, logbook = main()
    allfits = [part.fitness.values for part in swarm]

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
    # load data
    datapath = os.path.join(os.path.abspath('./'), 'data')
    filename_list = ['logbook_'+suffix+'_beta_MOPSO2.npz', 'bookkeeping_'+suffix+'.npz',
            'popdata_'+suffix+'_beta_MOPSO2.npz']
    datafiles = []
    for filename in filename_list:
        datafile = os.path.join(datapath,filename)
        datafiles.append(datafile)

    np.savez(datafiles[0], logbook=logbook)
    np.savez(datafiles[1], bookkeeping=System.bookkeeping)
    np.savez(datafiles[-1], allpop=swarm, allfits=allfits, front=swarm.gbest,
            frontfits=swarm.gbestfit, pop=swarm, popfits=allfits)
