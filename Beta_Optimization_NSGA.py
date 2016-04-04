import os
import sys
import itertools
import math
import operator
import random
import numpy as np
from multiprocessing import Pool, Manager

import time
import datetime

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
from management.performanceFuncs import pointintimeFunc
from management.component import Component
from management.system import System

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

#icorr_mean_list = np.array(input('corrosoin rate:')).astype('double')
#year = np.array(input('expected life:')).astype('double')
#num_processes = np.array(input('number of processes:')).astype('int')
icorr_mean_list = [1.,1.,0.5]
year = 100
num_processes = 40

creator.create("FitnessMulti", base.Fitness, weights=(-1.0,-1.0))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_plan", random.randint, FRP_DESIGN_YR-TIME_INTERVAL,
        SERVICE_LIFE-TIME_INTERVAL)

# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_plan, 3)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", pointintimeFunc, icorr_mean_list=icorr_mean_list, year=year)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutUniformInt,
        low=FRP_DESIGN_YR-TIME_INTERVAL, up=SERVICE_LIFE-TIME_INTERVAL, indpb=0.33)
toolbox.register("select", tools.selNSGA2)
toolbox.register("sort", tools.sortNondominated)

NPOP = 500
# stop criteria
NGEN = 10
NMAX = 100
TOL1 = 0.001
TOL2 = 0.001
NCR = 10


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

    pop = toolbox.population(n=NPOP)
    fits = toolbox.map(toolbox.evaluate, pop)
    for fit,ind in zip(fits, pop):
        ind.fitness.values = fit

    nevals = NPOP

    g = 1
    distances1 = []
    distances2 = []
    frontfitlast = np.array([[1., 0.]])
    nevalsum = 0
    evolStop = False
    halloffame = tools.ParetoFront()
    while not evolStop:
        offspring = algorithms.varOr(pop, toolbox, lambda_=NPOP, cxpb=0.9, mutpb=0.1)
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        nevals = len(invalid_ind)
        nevalsum += nevals
        fits = toolbox.map(toolbox.evaluate, invalid_ind)
        for fit,ind in zip(fits, invalid_ind):
            ind.fitness.values = fit
        pop = toolbox.select(offspring+pop, k=NPOP)
        front = toolbox.sort(offspring+pop, k=NPOP, first_front_only=True)[0]
        halloffame.update(front)

        ## check if stop evolution
        #distance=[]
        #frontfit = [ind.fitness.values for ind in halloffame]
        #for obj in frontfit:
            #vector = np.array(frontfitlast)-np.array(obj)
            #distance.append(min(np.linalg.norm(vector, axis=1)))
        #distances.append(np.mean(distance))
        #longest = 0.
        #distances.append(np.mean(distance))
        #vertlongest = 0.
        #horilongest = 0.
        #for point1 in frontfit:
            #for point2 in frontfit:
                ##dist = np.linalg.norm(np.array(point1)-np.array(point2))
                #if dist > longest:
                    #longest = dist
        #tol = longest/NPOP
        #tol = np.maximum(tol, TOL)
        #evolStop = (len(distances)>NGEN and all(np.array(distances[-NGEN:])<tol)) or g>NMAX
        #frontfitlast = frontfit
        # check if stop evolution
        distance1=[]
        distance2=[]
        frontfit = np.array([ind.fitness.values for ind in halloffame])
        frontfit[frontfit[:,1] == 0,1] = 1.
        frontfitlast[frontfitlast[:,1] == 0,1] = 1.
        for obj in frontfit:
            #vector = np.array(frontfitlast)-np.array(obj)
            #distance.append(min(np.linalg.norm(vector, axis=1)))
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
        delta_time = time.time() - start_delta_time
        logbook.record(gen=g, evals=nevals,nfront=len(halloffame),
                mean1=distances1[-1], mean2=distances2[-1],
                tol1=tol1, tol2=tol2, time=delta_time)
        print(logbook.stream)

        g+=1

    pool.close()
    pool.join()

    System.bookkeeping = System.bookkeeping.copy()

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))


    return pop, logbook, halloffame, nevalsum

if __name__ == "__main__":

    allpop, logbook, halloffame, nevalsum = main()
    front_parallel = halloffame

    allfits = [ind.fitness.values for ind in allpop]
    frontfits = [ind.fitness.values for ind in front_parallel]
    pop = allpop[-NPOP:]
    popfits = [ind.fitness.values for ind in pop]

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
    filename_list = ['logbook_'+suffix+'_beta_NSGA.npz','bookkeeping_'+suffix+'.npz',
            'popdata_'+suffix+'_beta_NSGA.npz']
    datafiles = []
    for filename in filename_list:
        datafile = os.path.join(datapath,filename)
        datafiles.append(datafile)

    np.savez(datafiles[0], logbook=logbook)
    np.savez(datafiles[1], bookkeeping=System.bookkeeping)
    np.savez(datafiles[-1], allpop=allpop, allfits=allfits, front=front_parallel,
            frontfits=frontfits, pop=pop, popfits=popfits)
