import numpy as np
import time
import datetime

from constants.sampleConstants import NUM_COMPONENT, NUM_ADAPTATION,\
    NUM_MAIN_SMP, NUM_KOPT_SMP, NUM_PRE_SMP
from constants import END_AGE, RELIABILITY_DT, SERVICE_LIFE, FRP_DESIGN_YR
from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
from corrosion.simpleCorrosion import simpleCorrosionLHS
from management.component import Component
from management.system import System
from pyre.distributions import *

def evalFitness(individual):
    comp_type_list = ['flexure', 'shear', 'deck']
    pf_list = []
    cost_list = []
    for str_yr,comp_type in zip(individual, comp_type_list):
        component = Component(comp_type, str_yr=str_yr)
        if str_yr in Component.pfkeeping[comp_type][0,:]:
            indx = np.where(Component.pfkeeping[comp_type][0,:]==str_yr)[0][0]
            pf = Component.pfkeeping[comp_type][1,:][indx]
        else:
            time.sleep(0.01)
            pf = 1./str_yr
            tmp = np.array([[str_yr], [pf]])
            Component.pfkeeping[comp_type] = np.hstack((Component.pfkeeping[comp_type], tmp))
        pf_list.append(pf)

        if str_yr in Component.costkeeping[comp_type][0,:]:
            indx = np.where(Component.costkeeping[comp_type][0,:]==str_yr)[0][0]
            cost = Component.costkeeping[comp_type][1,:][indx]
        else:
            time.sleep(0.01)
            cost = 10.*str_yr
            tmp = np.array([[str_yr], [cost]])
            Component.costkeeping[comp_type] = np.hstack((Component.costkeeping[comp_type], tmp))
        cost_list.append(cost)

    return max(pf_list),sum(cost_list)

def evalFitnessMatlab(individual):
    y = np.zeros(2)
    for i in xrange(2):
        y[0] = y[0] - 10*np.exp(-0.2*np.sqrt(individual[i]**2+ individual[i+1]**2))

    for i in xrange(3):
        y[1] = y[1] + np.abs(individual[i])**0.8 + 5*np.sin(individual[i]**3)

    return y[0], y[1]

def performanceFunc(str_yr_list, icorr_mean_list=[1,1,1], year=100, register=True):
    # Structural age
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    comp_type_list = ['flexure', 'shear', 'deck']
    # construct virgin component
    no_str_list = [0, 0, 0]
    virgin_component_list = []
    for comp_type,str_yr,icorr_mean,str_yr2 in zip(comp_type_list, no_str_list, icorr_mean_list, str_yr_list):
        virgin_component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
        if str_yr in Component.costkeeping[comp_type][0,:] and\
            str_yr in Component.pfkeeping[comp_type][0,:] and\
                str_yr2 in Component.costkeeping[comp_type][0,:] and\
                   str_yr2 in Component.pfkeeping[comp_type][0,:]:
                    pass
        else:
            resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
            #if str_yr not in Component.costkeeping[comp_type][0,:]:
            #    # register the cost
            #    Component.registerCost(comp_type, str_yr, cost)
            # construc the list
            virgin_component.setServiceTime(service_time)
            virgin_component.setResistanceMean(resistance_mean)
            virgin_component.setResistanceCov(resistance_cov)
            virgin_component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
        virgin_component_list.append(virgin_component)

    if str_yr_list == no_str_list:
        # construct system and return if no strengthening plan
        system = System(virgin_component_list)
        pf, pf_dict = system.accumulativePf(year, register)
        total_cost = 0.0
        for comp_type in comp_type_list:
            if str_yr not in Component.costkeeping[comp_type][0,:]:
                # register the cost
                Component.registerCost(comp_type, 0, 0)
    else:
        # update pf and total_cost if having strengthening
        # strengthened component
        cost_list = []
        component_list = []
        for comp_type,str_yr,icorr_mean in zip(comp_type_list, str_yr_list, icorr_mean_list):
            component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
            # cost and pf
            if str_yr in Component.costkeeping[comp_type][0,:] and\
                    str_yr in Component.pfkeeping[comp_type][0,:]:
                # no need to get mean and cov of resistance since pf is in
                # bookkeeping database
                indx = np.where(Component.costkeeping[comp_type][0,:]==str_yr)[0][0]
                cost = Component.costkeeping[comp_type][1,indx]
            elif str_yr not in Component.costkeeping[comp_type][0,:] and\
                    str_yr in Component.pfkeeping[comp_type][0,:]:
                # no need to set service_time and resistance since pf in
                # bookkeeping database
                # calculate resistance and cost
                resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
                # register the cost
                Component.registerCost(comp_type, str_yr, cost)
            else:
                resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
                # register the cost
                Component.registerCost(comp_type, str_yr, cost)
                # construc the list
                component.setServiceTime(service_time)
                component.setResistanceMean(resistance_mean)
                component.setResistanceCov(resistance_cov)
                component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
            cost_list.append(cost)
            component_list.append(component)
        # construct system
        system = System(virgin_component_list, component_list)
        pf,pf_dict = system.accumulativePf(year, register)
        total_cost = sum(cost_list)

    return pf,total_cost

def pointintimeFunc(str_yr_list, icorr_mean_list=[1,1,1], year=100, register=True):
    key = np.array2string(np.array(str_yr_list, dtype=int))
    if key in iter(System.bookkeeping.keys()):
        pf = System.bookkeeping[key][0]
        total_cost = System.bookkeeping[key][1]
    else:
        # Structural age
        service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
        comp_type_list = ['flexure', 'shear', 'deck']
        # construct virgin component
        no_str_list = [0, 0, 0]
        virgin_component_list = []
        for comp_type,str_yr,icorr_mean in zip(comp_type_list, no_str_list, icorr_mean_list):
            virgin_component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
            resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
            virgin_component.setServiceTime(service_time)
            virgin_component.setResistanceMean(resistance_mean)
            virgin_component.setResistanceCov(resistance_cov)
            virgin_component_list.append(virgin_component)

        if str_yr_list == no_str_list:
            # construct system and return if no strengthening plan
            system = System(virgin_component_list)
            pf, pf_dict = system.pointintimePf(year, register)
            total_cost = 0.0
            for comp_type in comp_type_list:
                if str_yr not in Component.costkeeping[comp_type][0,:]:
                    # register the cost
                    Component.registerCost(comp_type, 0, 0)
        else:
            # update pf and total_cost if having strengthening
            # strengthened component
            total_cost = 0.0
            component_list = []
            for comp_type,str_yr,icorr_mean in zip(comp_type_list, str_yr_list, icorr_mean_list):
                component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
                resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
                # register the cost
                Component.registerCost(comp_type, str_yr, cost)
                # construc the list
                component.setServiceTime(service_time)
                component.setResistanceMean(resistance_mean)
                component.setResistanceCov(resistance_cov)
                total_cost += cost
                component_list.append(component)
            # construct system
            system = System(virgin_component_list, component_list)
            pf,pf_dict = system.pointintimePf(year, register)
        if register:
            System.bookkeeping[key] = [pf,total_cost]

    return pf,total_cost

def performanceHistory(str_yr_list, icorr_mean_list, year, register=False):
    # Structural age
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    comp_type_list = ['flexure', 'shear', 'deck']
    # construct virgin component
    no_str_list = [0, 0, 0]
    virgin_component_list = []
    for comp_type,str_yr,icorr_mean in zip(comp_type_list, no_str_list, icorr_mean_list):
        virgin_component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
        resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
        virgin_component.setServiceTime(service_time)
        virgin_component.setResistanceMean(resistance_mean)
        virgin_component.setResistanceCov(resistance_cov)
        virgin_component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
        virgin_component_list.append(virgin_component)

    if str_yr_list == no_str_list:
        # construct system and return if no strengthening plan
        system = System(virgin_component_list)
        pf_sys,pf_dict = system.accumulativePf(year, register)
        total_cost = 0.0
    else:
        # update pf and total_cost if having strengthening
        # strengthened component
        cost_list = []
        component_list = []
        for comp_type,str_yr,icorr_mean in zip(comp_type_list, str_yr_list, icorr_mean_list):
            component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
            # cost and pf
            if str_yr in Component.costkeeping[comp_type][0,:] and\
                    str_yr in Component.pfkeeping[comp_type][0,:]:
                indx = np.where(Component.costkeeping[comp_type][0,:]==str_yr)[0][0]
                cost = Component.costkeeping[comp_type][1,indx]
            elif str_yr not in Component.costkeeping[comp_type][0,:] and\
                    str_yr in Component.pfkeeping[comp_type][0,:]:
                resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
            else:
                resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
                # construc the list
                component.setServiceTime(service_time)
                component.setResistanceMean(resistance_mean)
                component.setResistanceCov(resistance_cov)
                component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
            component_list.append(component)
        # construct system
        system = System(virgin_component_list, component_list)
        pf_sys, pf_dict = system.accumulativePf(year, register)
        total_cost = sum(cost_list)

    pf_list = [pf_sys, pf_dict['flexure'], pf_dict['shear'], pf_dict['deck']]

    return pf_list

def pointintimeHistory(str_yr_list, icorr_mean_list=[1,1,1], year=100, register=True):
    key = np.array2string(np.array(str_yr_list, dtype=int))
    # Structural age
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    comp_type_list = ['flexure', 'shear', 'deck']
    # construct virgin component
    no_str_list = [0, 0, 0]
    virgin_component_list = []
    for comp_type,str_yr,icorr_mean in zip(comp_type_list, no_str_list, icorr_mean_list):
        virgin_component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
        resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
        virgin_component.setServiceTime(service_time)
        virgin_component.setResistanceMean(resistance_mean)
        virgin_component.setResistanceCov(resistance_cov)
        virgin_component_list.append(virgin_component)

    if str_yr_list == no_str_list:
        # construct system and return if no strengthening plan
        system = System(virgin_component_list)
        pf, pf_dict = system.pointintimePf(year, register)
    else:
        # update pf and total_cost if having strengthening
        # strengthened component
        component_list = []
        for comp_type,str_yr,icorr_mean in zip(comp_type_list, str_yr_list, icorr_mean_list):
            component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
            resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
            # construc the list
            component.setServiceTime(service_time)
            component.setResistanceMean(resistance_mean)
            component.setResistanceCov(resistance_cov)
            component_list.append(component)
        # construct system
        system = System(virgin_component_list, component_list)
        pf,pf_dict = system.pointintimePf(year, register)

    pf_sys = np.max(pf_dict['system'])
    pf_flex = np.max(pf_dict['flexure'])
    pf_shear = np.max(pf_dict['shear'])
    pf_deck = np.max(pf_dict['deck'])
    pf_list = [pf_sys, pf_flex, pf_shear, pf_deck]

    return pf_list

def generateCompData(comp_type, str_yr, year, icorr_mean=1., life=100):
    # Structural age
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    comp_type_list = ['flexure', 'shear', 'deck']
    # construct virgin component
    virgin_component = Component(comp_type, maintain_tag=False, str_yr=0)
    resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
    virgin_component.setServiceTime(service_time)
    virgin_component.setResistanceMean(resistance_mean)
    virgin_component.setResistanceCov(resistance_cov)
    pf = virgin_component.pointintimePf(year, register=False)

    #strengthened component
    if str_yr != 0:
        if year<str_yr:
            pf = virgin_component.pointintimePf(year, register=False)
        else:
            component = Component(comp_type, str_yr=str_yr)
            resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
            component.setServiceTime(service_time)
            component.setResistanceMean(resistance_mean)
            component.setResistanceCov(resistance_cov)
            pf = component.pointintimePf(year, register=False)

    return pf,cost

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    print "debugging of performance function"
    start_delta_time = time.time()

    year=100
    pf_list = []
    cost_list = []
    str_yr_list = np.arange(0, 81, 20)
    str_schemes = [[0, 0, 0], [0, 0, 10], [0, 23, 29], [50, 23, 29]]
    for str_yr_list in str_schemes:
        #pf_sys, total_cost = performanceFunc([i, i, i], year=year)
        #pf_sys, total_cost = pointintimeFunc([i, i, i], year=year)
        pf_sys, total_cost = pointintimeFunc(str_yr_list, year=year)
        #pf = Component.pfkeeping['flexure'][1,1]
        pf_list.append(pf_sys)
        cost_list.append(total_cost)
    #print 'failure probability at year {:d}: {:.4e}'.format(year, pf)
    #print 'total cost during year {:d}: {:.4e} mm^3'.format(year, total_cost)

    print 'mean {}'.format(np.mean(pf_list))
    print 'std {}'.format(np.std(pf_list))

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))

    plt.plot(pf_list, 'bo')
