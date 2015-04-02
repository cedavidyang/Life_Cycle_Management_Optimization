import os
import sys
import numpy as np
import array
from multiprocessing import Pool, Manager
import matplotlib.pyplot as plt
plt.ion()

from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
from management.component import Component
from management.system import System
from management.performanceFuncs import evalFitness

import time
import datetime

import random

from deap import algorithms
from deap import base
from deap import creator
from deap import tools


icorr_mean_list = [1,1,1]
year = 100

creator.create("FitnessMulti", base.Fitness, weights=(-1.0,-1.0))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_plan", random.randint, FRP_DESIGN_YR, SERVICE_LIFE-TIME_INTERVAL)

# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_plan, 3)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", evalFitness)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutUniformInt,
        low=FRP_DESIGN_YR, up=SERVICE_LIFE-TIME_INTERVAL, indpb=0.05)
toolbox.register("select", tools.selNSGA2)
toolbox.register("sort", tools.sortNondominated)

#class MyManager(BaseManager):
#    pass
#MyManager.register('pfkeeping', Component, exposed=('pfkeeping'))
#MyManager.register('costkeeping', Component, exposed=('costkeeping'))

def main_series():
    Component.resetPfKeeping()
    Component.resetCostKeeping()

    print "MULTIOBJECTIVE OPTIMIZATION: series version"
    start_delta_time = time.time()

    # optimization
    random.seed(64)

    npop = 1000
    ngen = 200

    stats = tools.Statistics(key=lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    logbook = tools.Logbook()
    logbook.header = "gen", "evals", "avg", "std", "min", "max"

    pop = toolbox.population(n=npop)
    fits = toolbox.map(toolbox.evaluate, pop)
    for fit,ind in zip(fits, pop):
        ind.fitness.values = fit

    nevals = npop
    allpop = []
    for gen in range(ngen):
        allpop = allpop+pop
        record = stats.compile(pop)
        logbook.record(gen=gen, evals=nevals, **record)
        print(logbook.stream)

        offspring = algorithms.varOr(pop, toolbox, lambda_=npop, cxpb=0.5, mutpb=0.1)
        #invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #invalid_ind = offspring
        nevals = len(offspring)
        fits = toolbox.map(toolbox.evaluate, offspring)
        for fit,ind in zip(fits, offspring):
            ind.fitness.values = fit
        pop = toolbox.select(offspring+pop, k=npop)

    front = toolbox.sort(allpop, k=int(ngen*npop), first_front_only=True)

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))

    return pop, logbook, front

def main_parallel():
    Component.resetPfKeeping()
    Component.resetCostKeeping()

    manager = Manager()
    Component.pfkeeping = manager.dict(Component.pfkeeping)
    Component.costkeeping = manager.dict(Component.costkeeping)

    pool = Pool(processes=3)
    toolbox.register("map", pool.map)

    print "MULTIOBJECTIVE OPTIMIZATION: parallel version"
    start_delta_time = time.time()

    # optimization
    random.seed(64)

    npop = 100
    ngen = 50

    stats = tools.Statistics(key=lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    logbook = tools.Logbook()
    logbook.header = "gen", "evals", "avg", "std", "min", "max"

    pop = toolbox.population(n=npop)
    fits = toolbox.map(toolbox.evaluate, pop)
    for fit,ind in zip(fits, pop):
        ind.fitness.values = fit

    nevals = npop
    allpop = []
    for gen in range(ngen):
        allpop = allpop+pop
        record = stats.compile(pop)
        logbook.record(gen=gen, evals=nevals, **record)
        print(logbook.stream)

        offspring = algorithms.varOr(pop, toolbox, lambda_=npop, cxpb=0.5, mutpb=0.1)
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        nevals = len(invalid_ind)
        fits = toolbox.map(toolbox.evaluate, invalid_ind)
        for fit,ind in zip(fits, invalid_ind):
            ind.fitness.values = fit
        pop = toolbox.select(offspring+pop, k=npop)

    front = toolbox.sort(allpop, k=int(ngen*npop), first_front_only=True)

    pool.close()
    pool.join()

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))

    return allpop, logbook, front

if __name__ == "__main__":
    #pop,log,front_series = main_series()
    #flex = Component.pfkeeping['flexure']
    #indx = np.argsort(flex[0])
    #flex = flex[:,indx]
    #np.save('flex_series1.npy', flex)
    #pop,log,front_series = main_series()
    #flex = Component.pfkeeping['flexure']
    #indx = np.argsort(flex[0])
    #flex = flex[:,indx]
    #np.save('flex_series2.npy', flex)
    allpop,log,front_parallel = main_series()
    flex = Component.pfkeeping['flexure']
    indx = np.argsort(flex[0])
    flex = flex[:,indx]
    #np.save('flex_parallel.npy', flex)

    allfits = [ind.fitness.values for ind in allpop]
    frontfits = [ind.fitness.values for ind in front_parallel[0]]

    #toolbox.register("map", map)
    #pop_res = toolbox.map(toolbox.evaluate, pop)
    #pop_res = np.array(pop_res)
    #front_series_res = toolbox.map(toolbox.evaluate, front_series[0])
    #front_parallel_res = toolbox.map(toolbox.evaluate, front_parallel[0])
    #front_series_res = np.array(front_series_res)
    #front_parallel_res = np.array(front_parallel_res)

    plt.close('all')
    plt.rc('font', family='serif', size=12)

    #plt.scatter(pop_res[:,0], pop_res[:,1], facecolors='none')
    #plt.scatter(front_series_res[:,0], front_series_res[:,1],
    #        marker='^', facecolors='b', edgecolors='b', label='series w/ bookkeeping')
    #plt.scatter(front_parallel_res[:,0], front_parallel_res[:,1],
    #        marker='o', facecolors='r', edgecolors='r', label='parallel w/ bookkeeping and shared memory')
    plt.scatter(np.array(allfits)[:,0], np.array(allfits)[:,1], facecolors='none')
    plt.scatter(np.array(frontfits)[:,0], np.array(frontfits)[:,1],
            marker='o', facecolors='r', edgecolors='r', label='parallel w/ bookkeeping and shared memory')

    plt.xlabel('Surrogate failure probability')
    plt.ylabel('Surrogate strengthening cost')
    #plt.legend(loc='upper right', prop={'size':12})

    for front in front_parallel[0]:
        if front[0] == front[1] and front[1] == front[2]:
            print '{}'.format(front)
        else:
            continue
