import os
import sys
import numpy as np
import array
from multiprocessing import Pool, Manager
import matplotlib.pyplot as plt
plt.ion()

import time
import datetime

import random
from management.performanceFuncs import evalFitnessMatlab

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

def mutUniformUser(individual, low, high, indpb):
    for indx, chrom in enumerate(individual):
        if random.random()<indpb:
            individual[indx] = np.random.uniform(low, high)

    return individual,

creator.create("FitnessMulti", base.Fitness, weights=(-1.0,-1.0))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_plan", np.random.uniform, low=-5, high=5)
#toolbox.register("attr_plan", random.randint, -5, 5)

# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual,
        toolbox.attr_plan, 3)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", evalFitnessMatlab)
#toolbox.register("mate", tools.cxUniform, indpb=0.9)
toolbox.register("mate", tools.cxOnePoint)
#toolbox.register("mutate", tools.mutUniformInt, low=-5, up=5, indpb=0.05)
toolbox.register("mutate", mutUniformUser, low=-5, high=5, indpb=0.05)
#toolbox.register("mutate", tools.mutPolynomialBounded, eta=1,
        #low=-5, up=5, indpb=0.33)
#toolbox.register("mutate", tools.mutUniformInt, low=-5, up=5, indpb=0.05)
toolbox.register("select", tools.selNSGA2)
toolbox.register("sort", tools.sortNondominated)

def main():
    print "MULTIOBJECTIVE OPTIMIZATION: series version"
    start_delta_time = time.time()

    # optimization
    #random.seed(64)

    npop = 500
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

        offspring = algorithms.varOr(pop, toolbox, lambda_=npop, cxpb=0.9, mutpb=0.1)
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

if __name__ == "__main__":
    pop,log,front = main()

    popfits = [ind.fitness.values for ind in pop]
    frontfits = [ind.fitness.values for ind in front[0]]

    plt.close('all')
    plt.rc('font', family='serif', size=12)

    plt.scatter(np.array(popfits)[:,0], np.array(popfits)[:,1],
            marker='v', facecolors='b', edgecolors='b', label='parallel w/ bookkeeping and shared memory')
    plt.scatter(np.array(frontfits)[:,0], np.array(frontfits)[:,1],
            marker='^', facecolors='r', edgecolors='r', label='parallel w/ bookkeeping and shared memory')

    plt.xlabel('Objective 1')
    plt.ylabel('Objective 2')
