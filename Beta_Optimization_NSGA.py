import os
import sys
import numpy as np
from multiprocessing import Pool, Manager

import time
import datetime

import random

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
from management.performanceFuncs import pointintimeFunc
from management.component import Component

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

#icorr_mean_list = np.array(input('corrosoin rate:')).astype('double')
#year = np.array(input('expected life:')).astype('double')
#num_processes = np.array(input('number of processes:')).astype('int')
icorr_mean_list = [1.,1.,1.]
year = 100
num_processes = 10

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

NPOP = 300
NGEN = 30


def main():
    # reset bookkeeping
    Component.resetCostKeeping()
    Component.resetPfKeeping()
    Component.resetRiskKeeping()

    ## use existing pf data
    #pfkeeping = np.load('pfkeeping.npz')
    #Component.pfkeeping['flexure'] = pfkeeping['flexure']
    #Component.pfkeeping['shear'] = pfkeeping['shear']
    #Component.pfkeeping['deck'] = pfkeeping['deck']
    ## use existing cost data
    #costkeeping = np.load('costkeeping.npz')
    #Component.costkeeping['flexure'] = costkeeping['flexure']
    #Component.costkeeping['shear'] = costkeeping['shear']
    #Component.costkeeping['deck'] = costkeeping['deck']

    manager = Manager()
    Component.pfkeeping = manager.dict(Component.pfkeeping)
    Component.costkeeping = manager.dict(Component.costkeeping)
    Component.riskkeeping = manager.dict(Component.riskkeeping)

    pool = Pool(processes=num_processes)
    toolbox.register("map", pool.map)

    print "MULTIOBJECTIVE OPTIMIZATION: parallel version"
    start_delta_time = time.time()

    # optimization
    random.seed(64)

    logbook = tools.Logbook()
    logbook.header = ["gen", "evals", "nfront", "mean", "tol"]

    pop = toolbox.population(n=NPOP)
    fits = toolbox.map(toolbox.evaluate, pop)
    for fit,ind in zip(fits, pop):
        ind.fitness.values = fit

    nevals = NPOP

    g = 1
    distances = []
    frontfitlast = np.zeros((1,2))
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

        # check if stop evolution
        distance=[]
        frontfit = [ind.fitness.values for ind in halloffame]
        for obj in frontfit:
            vector = np.array(frontfitlast)-np.array(obj)
            distance.append(min(np.linalg.norm(vector, axis=1)))
        distances.append(np.mean(distance))
        longest = 0.
        for point1 in frontfit:
            for point2 in frontfit:
                dist = np.linalg.norm(np.array(point1)-np.array(point2))
                if dist > longest:
                    longest = dist
        tol = longest/NPOP
        tol = np.maximum(tol, TOL)
        evolStop = (len(distances)>NGEN and all(np.array(distances[-NGEN:])<tol)) or g>NMAX
        frontfitlast = frontfit

        # Gather all the fitnesses in one list and print the stats
        record = stats.compile(pop)
        logbook.record(gen=g, evals=nevals,nfront=len(halloffame),
                mean=distances[-1], tol=tol, **stats.compile(pop))
        print(logbook.stream)

        g+=1

    pool.close()
    pool.join()

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))


    return pop, logbook, halloffame, nevalsum

if __name__ == "__main__":

    allpop, log, halloffame, nevalsum = main()
    front_parallel = halloffame

    # sort pfbooking
    pf_flex = Component.pfkeeping['flexure']
    indx = np.argsort(pf_flex[0])
    pf_flex = pf_flex[:,indx]
    pf_shear = Component.pfkeeping['shear']
    indx = np.argsort(pf_shear[0])
    pf_shear = pf_shear[:,indx]
    pf_deck = Component.pfkeeping['deck']
    indx = np.argsort(pf_deck[0])
    pf_deck = pf_deck[:,indx]
    # save costbooking
    cost_flex = Component.costkeeping['flexure']
    indx = np.argsort(cost_flex[0])
    cost_flex = cost_flex[:,indx]
    cost_shear = Component.costkeeping['shear']
    indx = np.argsort(cost_shear[0])
    cost_shear = cost_shear[:,indx]
    cost_deck = Component.costkeeping['deck']
    indx = np.argsort(cost_deck[0])
    cost_deck = cost_deck[:,indx]

    allfits = [ind.fitness.values for ind in allpop]
    frontfits = [ind.fitness.values for ind in front_parallel[0]]
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
    filename_list = ['pfkeeping_'+suffix+'.npz', 'costkeeping_'+suffix+'.npz',
            'popdata_'+suffix+'.npz']
    datafiles = []
    for filename in filename_list:
        datafile = os.path.join(datapath,filename)
        datafiles.append(datafile)

    np.savez(datafiles[0], flexure=pf_flex, shear=pf_shear, deck=pf_deck)
    np.savez(datafiles[1], flexure=cost_flex, shear=cost_shear, deck=cost_deck)
    np.savez(datafiles[2], allpop=allpop, allfits=allfits, front=front_parallel,
            frontfits=frontfits, pop=pop, popfits=popfits)
