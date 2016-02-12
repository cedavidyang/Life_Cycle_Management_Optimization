# system class for bridge superstructure
import os
import sys
import numpy as np

from component import Component

class System(object):

    def __init__(self, virgin_component_list, component_list=None):
        self.virgin_component_list = virgin_component_list
        if component_list is None:
            self.component_list = virgin_component_list
        else:
            self.component_list = component_list

    def addComponent(self, component):
        self.component_list.append(component)

    def accumulativePf(self, timepoint, register=True):
        pf_dict = {}
        for component,virgin_component in zip(self.component_list, self.virgin_component_list):
            str_yr = component.str_yr
            comp_type = component.comp_type
            if str_yr in Component.pfkeeping[comp_type][0,:]:
                indx = np.where(Component.pfkeeping[comp_type][0,:]==str_yr)[0][0]
                pf_dict[component.comp_type] = Component.pfkeeping[comp_type][1,indx]
            elif component.maintain_tag == False or str_yr>=timepoint:
                pf_dict[comp_type] = virgin_component.accumulativePf(timepoint, register=False)
                if register == True:
                    tmp = np.array([[str_yr], [pf_dict[comp_type]]])
                    Component.pfkeeping[comp_type] = np.hstack((Component.pfkeeping[comp_type], tmp))
            else:
                frp_At = 1 - component.accumulativePf(timepoint-str_yr, register=False)
                rc_St = 1 - virgin_component.accumulativePf(str_yr, register=False)
                St = rc_St*frp_At
                pf_dict[component.comp_type] = 1-St
                # write survival to bookkeeping
                if register == True:
                    tmp = np.array([[str_yr], [1-St]])
                    Component.pfkeeping[comp_type] = np.hstack((Component.pfkeeping[comp_type], tmp))

        flex_survivor = 1-pf_dict['flexure']
        shear_survivor = 1-pf_dict['shear']
        beam_survivor = flex_survivor * shear_survivor
        beam_fail = 1 - beam_survivor
        deck_survivor = 1-pf_dict['deck']

        sys_survivor = deck_survivor * (beam_survivor**3 + beam_survivor**3*beam_fail
            + beam_survivor**2*beam_fail + beam_survivor**3*beam_fail + beam_survivor**3*beam_fail**2)
        sys_pf = 1-sys_survivor

        return sys_pf, pf_dict

    def pointintimePf(self, timepoint, register=True):
        pf_dict = {}
        #str_yr_array = np.zeros((len(self.component_list),1), dtype=int)
        str_yr_list = []
        pf_dummy = 0.
        for i,component in enumerate(self.component_list):
            str_yr_list.append(component.str_yr)
        #str_yr_list.append(timepoint)
        str_yr_array = np.array(str_yr_list)

        nstr = str_yr_array.size
        syspf_list = []
        for i,t in enumerate(np.append(str_yr_array,timepoint)):
            frpindx = str_yr_array<t
            components = np.array(self.component_list)[frpindx]
            virgin_components = np.array(self.virgin_component_list)[np.logical_not(frpindx)]
            for component in components:
                comp_type = component.comp_type
                if t in Component.pfkeeping[comp_type][0,:]:
                    indx = np.where(Component.pfkeeping[comp_type][0,:]==t)[0][0]
                    pf_dict[component.comp_type] = Component.pfkeeping[comp_type][1,indx]
                else:
                    if component.service_time is None:
                        import ipdb; ipdb.set_trace() # BREAKPOINT
                        pass
                    pf_dict[component.comp_type] = component.pointintimePf(t-component.str_yr)
                    # write survival to bookkeeping
                    if register == True:
                        tmp = np.array([[t], [pf_dict[component.comp_type]]])
                        Component.pfkeeping[comp_type] = np.hstack((Component.pfkeeping[comp_type], tmp))
            for component in virgin_components:
                comp_type = component.comp_type
                if t in Component.pfkeeping['virgin'+comp_type][0,:]:
                    indx = np.where(Component.pfkeeping['virgin'+comp_type][0,:]==t)[0][0]
                    pf_dict[component.comp_type] = Component.pfkeeping['virgin'+comp_type][1,indx]
                else:
                    if component.service_time is None:
                        import ipdb; ipdb.set_trace() # BREAKPOINT
                        pass
                    pf_dict[component.comp_type] = component.pointintimePf(t)
                    # write survival to bookkeeping
                    if register == True:
                        tmp = np.array([[t], [pf_dict[component.comp_type]]])
                        Component.pfkeeping['virgin'+comp_type] = np.hstack((Component.pfkeeping['virgin'+comp_type], tmp))

            flex_survivor = 1-pf_dict['flexure']
            shear_survivor = 1-pf_dict['shear']
            beam_survivor = flex_survivor * shear_survivor
            beam_fail = 1 - beam_survivor
            deck_survivor = 1-pf_dict['deck']

            sys_survivor = deck_survivor * (beam_survivor**3 + beam_survivor**3*beam_fail
                + beam_survivor**2*beam_fail + beam_survivor**3*beam_fail + beam_survivor**3*beam_fail**2)
            sys_pf = 1-sys_survivor
            syspf_list.append(sys_pf)

        sys_pf_max = np.max(syspf_list)

        return sys_pf_max, pf_dict


if __name__ == '__main__':
    import time
    import datetime
    from constants.sampleConstants import NUM_COMPONENT, NUM_ADAPTATION,\
        NUM_MAIN_SMP, NUM_KOPT_SMP, NUM_PRE_SMP
    from constants import END_AGE, RELIABILITY_DT
    from constants.simpleCorrosionConstants import START_AGE, TIME_INTERVAL, END_AGE
    from corrosion.simpleCorrosion import simpleCorrosionLHS
    from pyre.distributions import *

    print "debugging of system class"
    start_delta_time = time.time()

    # Structural age
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    comp_type_list = ['flexure', 'shear', 'deck']
    icorr_mean_list = [1, 1, 1]

    # construct virgin component
    str_yr_list = [0, 0, 0]
    cost_list = []
    virgin_component_list = []
    for comp_type,str_yr,icorr_mean in zip(comp_type_list, str_yr_list, icorr_mean_list):
        virgin_component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
        resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
        virgin_component.setServiceTime(service_time)
        virgin_component.setResistanceMean(resistance_mean)
        virgin_component.setResistanceCov(resistance_cov)
        virgin_component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
        virgin_component_list.append(virgin_component)
        cost_list.append(cost)

    # strengthened component
    str_yr_list = [0, 50, 0]
    cost_list = []
    component_list = []
    for comp_type,str_yr,icorr_mean in zip(comp_type_list, str_yr_list, icorr_mean_list):
        component = Component(comp_type, maintain_tag=False, str_yr=str_yr)
        resistance_mean,resistance_cov,cost = simpleCorrosionLHS(comp_type, service_time, icorr_mean, str_yr)
        component.setServiceTime(service_time)
        component.setResistanceMean(resistance_mean)
        component.setResistanceCov(resistance_cov)
        component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
        component_list.append(component)
        cost_list.append(cost)

    year = 100
    system = System(virgin_component_list, component_list)
    pf = system.accumulativePf(year)
    total_cost = sum(cost_list)
    print 'failure probability at year {:d}: {:.4e}'.format(year, pf)
    print 'total cost during year {:d}: {:.4e} mm^3'.format(year, total_cost)

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))
