# life-cycle indices of component and system, in preparation for optimization
import os
import sys
import numpy as np

from scipy.interpolate import interp1d

from constants import FRP_DESIGN_YR, SERVICE_LIFE, END_AGE
from constants import RELIABILITY_DT
from constants.beamConstants import BFRP_NOM, SPAN_NOM, PLATE2SUPPORT_NOM
from constants.beamConstants import BEAM_HEIGHT_NOM, WEB_WIDTH_NOM, WFRPV_NOM, SFRPV_NOM
from constants.beamConstants import FRP_DECK_RATIO, SPAN_DECK_NOM, PLATE2SUPPORT_DECK_NOM
from constants.manageConstants import FRP_UNIT_PRICE, GIRDER_NUM
from preprocessing import getCorrosionLossFromFile, strengtheningFromCorrosionLoss


def loadRcSurvivor(type_of_component):
    if 'f' in type_of_component.lower():
        subfolder = 'flexure'
    elif 's' in type_of_component.lower():
        subfolder = 'shear'
    elif 'd' in type_of_component.lower():
        subfolder = 'deck'
    else:
        print "[ERROR:] illegal type_of_component, must be flexure, shear or deck"
        sys.exit(1)
    datapath = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', subfolder)
    filename = 'reliability_results_parallel.npz'
    datafile = os.path.join(datapath, filename)

    if os.path.isfile(datafile):
        res = np.load(datafile)
        pf = res['pf']
        rc_St = 1.-pf
    else:
        print "[ERROR:] {} does not exist for component {}".format(datafile, subfolder)
        print "run Life_Cycle_Reliability_RC.py first with parameter {}".format(type_of_component)
        sys.exit(1)

    return rc_St

def loadComponentAvailability(type_of_component, str_yr):
    if 'f' in type_of_component.lower():
        subfolder = 'flexure'
    elif 's' in type_of_component.lower():
        subfolder = 'shear'
    elif 'd' in type_of_component.lower():
        subfolder = 'deck'
    else:
        print "[ERROR:] illegal type_of_component, must be flexure, shear or deck"
        sys.exit(1)
    datapath = os.path.join(os.path.abspath('./'), 'data', 'frp_cs', subfolder)
    filename = 'reliability_parallel_stryr_' + str(int(str_yr)) + '.npz'
    datafile = os.path.join(datapath, filename)

    if os.path.isfile(datafile):
        res = np.load(datafile)
        pf = res['pf']
        At = 1.-pf
    else:
        print "[ERROR:] {} does not exist for component {}".format(datafile, subfolder)
        print "run Life_Cycle_Management_Availability.py first with parameter {}".format(type_of_component)
        sys.exit(1)

    return At


def componentSurvivor(type_of_component, str_yr, time_array):
    # time arrays
    rc_time_array = np.arange(RELIABILITY_DT,END_AGE+RELIABILITY_DT,RELIABILITY_DT)
    rc_time_array = np.insert(rc_time_array, 0, (0.,1.))
    frp_time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    frp_time_array = np.insert(frp_time_array, 0, (0.,1.))

    rc_St = loadRcSurvivor(type_of_component)
    rc_St = np.insert(rc_St, 0, 1)
    rc_St_interp = interp1d(rc_time_array, rc_St)
    if str_yr >= FRP_DESIGN_YR:
        frp_At = loadComponentAvailability(type_of_component, str_yr)
        frp_At = np.insert(frp_At, 0, 1)
        frp_At = 1-(1-frp_At)*(1+frp_time_array/10)
        frp_At_interp = interp1d(frp_time_array, frp_At)

        St = rc_St_interp(time_array)
        mask = time_array>str_yr
        St[mask] = rc_St_interp(str_yr) * frp_At_interp(time_array[mask]-str_yr)
    else:
        St = rc_St_interp(time_array)


    return St

def systemSurvivor(str_yr_dict, time_array):
    survivor_dict = {}
    for type_of_component, str_yr in str_yr_dict.iteritems():
        survivor_dict[type_of_component] = componentSurvivor(type_of_component, str_yr, time_array)

    flex_survivor = survivor_dict['f']
    shear_survivor = survivor_dict['s']
    beam_survivor = flex_survivor * shear_survivor
    beam_fail = 1 - beam_survivor
    deck_survivor = survivor_dict['d']

    sys_survivor = deck_survivor * (beam_survivor**3 + beam_survivor**3*beam_fail
        + beam_survivor**2*beam_fail + beam_survivor**3*beam_fail + beam_survivor**3*beam_fail**2)
    #sub_survivor = 1-beam_fail**2
    #sys_survivor = deck_survivor * sub_survivor**4

    return sys_survivor

def strengtheningCost(str_yr_dict):
    cost_dict = {}
    for type_of_component, str_yr in str_yr_dict.iteritems():
        if 'f' in type_of_component.lower():
            RC_CS_FLEX_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'flexure')
            datafile = os.path.join(RC_CS_FLEX_PATH, 'LWS_results.txt')
            compute_cost = lambda tfrp: GIRDER_NUM*tfrp*BFRP_NOM*(SPAN_NOM-2*PLATE2SUPPORT_NOM)*FRP_UNIT_PRICE
        elif 's' in type_of_component.lower():
            RC_CS_SHEAR_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'shear')
            datafile = os.path.join(RC_CS_SHEAR_PATH, 'LWS_results.txt')
            compute_cost = lambda tfrp: GIRDER_NUM*tfrp*WFRPV_NOM*(BEAM_HEIGHT_NOM*2+WEB_WIDTH_NOM)*(SPAN_NOM/3.*2/SFRPV_NOM)*FRP_UNIT_PRICE
        elif 'd' in type_of_component.lower():
            RC_CS_DECK_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'deck')
            datafile = os.path.join(RC_CS_DECK_PATH, 'LWS_results.txt')
            compute_cost = lambda tfrp: (GIRDER_NUM-1)*tfrp*(SPAN_NOM/FRP_DECK_RATIO)*(SPAN_DECK_NOM-2*PLATE2SUPPORT_DECK_NOM)*FRP_UNIT_PRICE
        if str_yr >= FRP_DESIGN_YR:
            # get corrosion loss
            loss_dict = getCorrosionLossFromFile(datafile, str_yr)
            # get tfrp
            dummy, tfrp_array = strengtheningFromCorrosionLoss(loss_dict, type_of_component)
            # compute cost
            cost = compute_cost(tfrp_array[0])
        else:
            cost = 0
        cost_dict[type_of_component] = cost

    total_cost = sum(cost_dict.values())

    return total_cost


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()

    # debugging of componentSurvivor()
    design_service_life = FRP_DESIGN_YR + SERVICE_LIFE
    time_array0 = np.arange(RELIABILITY_DT,design_service_life+RELIABILITY_DT,RELIABILITY_DT)
    str_yr_array = np.arange(FRP_DESIGN_YR, SERVICE_LIFE, 2)
    str_yr_array = np.insert(str_yr_array, 0, 0)

    plt.close('all')
    plt.rc('font', family='serif', size=12)

    component_list = ['f', 's', 'd']
    for type_of_component in component_list:
        plt.figure()
        plt.xlabel('service time (year)')
        plt.ylabel('CDF (log scale)')
        for str_yr in str_yr_array:
            time_array = time_array0 if str_yr==0 else np.sort(np.append(time_array0, str_yr))
            St = componentSurvivor(type_of_component, str_yr, time_array)
            plt.semilogy(time_array, 1-St)
        from matplotlib.patches import Ellipse, Circle
        from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea,AnchoredAuxTransformBox
        if type_of_component == 'f':
            #plt.figure()
            #fig=plt.figure(1, figsize=(3,3))
            #ax = plt.subplot(111)
            #ada = AnchoredDrawingArea(40, 20, 0, 0,
            #                        loc=1, pad=0., frameon=False)
            #p1 = Circle((10, 10), 10)
            #ada.drawing_area.add_artist(p1)
            #p2 = Circle((30, 10), 5, fc="r")
            #ada.drawing_area.add_artist(p2)
            #ax.add_artist(ada)
            ax = plt.gca()
            xy = (20,-35)
            ada = AnchoredDrawingArea(40,40,0,0,loc=1, pad=0., frameon=False)
            el = Ellipse(xy=xy, width=10, height=80, fc='none', ec='k', ls='dashed')
            ada.drawing_area.add_artist(el)
            ax.add_artist(ada)

    # debugging of systemSurvivor()
    #str_yr_dict = {'f':60, 's':0, 'd':56}
    str_yr_dict = {'f':0, 's':0, 'd':0}
    #str_yr_dict = {'f':0, 's':70, 'd':60}
    print 'Total strengthening cost: {} mm^3'.format(strengtheningCost(str_yr_dict))
    plt.figure()
    flex_St = componentSurvivor('f', str_yr_dict['f'], time_array0)
    shear_St = componentSurvivor('s', str_yr_dict['s'], time_array0)
    deck_St = componentSurvivor('d', str_yr_dict['d'], time_array0)
    sys_St = systemSurvivor(str_yr_dict, time_array0)
    print 'Failure probability during lifetime: {:.4e}'.format(1-sys_St[-1])
    plt.semilogy(time_array0, 1-flex_St, 'b--', label='flexure')
    plt.semilogy(time_array0, 1-shear_St, 'r--', label='shear')
    plt.semilogy(time_array0, 1-deck_St, 'g--', label='deck')
    plt.semilogy(time_array0, 1-sys_St, 'k-', label='system')
    plt.legend(loc=2)
    plt.xlabel('service time (year)')
    plt.ylabel('CDF (log scale)')
