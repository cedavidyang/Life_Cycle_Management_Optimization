# preprocessing data: preparation for optimization and management
import os
import sys
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
import warnings

from pyre.distributions import *
from constants import N_LHS_SMP, SERVICE_LIFE
from constants import FRP_FLEX_COV, FRP_SHEAR_COV, COVT0_COV
from constants import FRP_DECK_RATIO, FRP_DECK_COV
from frpDesignRatio import *

def getServiceTimeFromFile(datafile, strengthen_time):
    time_span = np.loadtxt(datafile)[0,:]
    new_service_time_list = []
    for timepoint in strengthen_time:
        start_indx = np.where(time_span==timepoint)[0][0]
        end_indx = np.where(time_span==timepoint+SERVICE_LIFE)[0][0]+1
        new_service_time = time_span[start_indx:end_indx] - timepoint
        new_service_time_list.append(new_service_time)

    return np.array(new_service_time_list)

def getCorrosionLossFromFile(datafile, strengthen_time):
    """ get corrosion loss from datafile
    """
    #from corrosion.initialVariables import region1_no, region2_no, ds_region1_mean, ds_region2_mean
    rebar_type = os.path.basename(os.path.dirname(datafile))
    if 'f' in rebar_type.lower():
        ds_region1_mean = DS_REGION1_MEAN
        ds_region2_mean = DS_REGION2_MEAN
        region1_no = REGION1_NO
        region2_no = REGION2_NO
    elif 's' in rebar_type.lower():
        ds_region1_mean = DS_REGION1_SHEAR_MEAN
        ds_region2_mean = DS_REGION2_SHEAR_MEAN
        region1_no = REGION1_NO_SHEAR
        region2_no = REGION2_NO_SHEAR
    elif 'd' in rebar_type.lower():
        ds_region1_mean = DS_REGION1_DECK_MEAN
        ds_region2_mean = DS_REGION2_DECK_MEAN
        region1_no = REGION1_NO_DECK
        region2_no = REGION2_NO_DECK
    else:
        print '[ERROR:] Illegal input for rebar type, must be flexure, shear or deck'
    time_span = np.loadtxt(datafile)[0,:]
    corrosion_prob1 = np.loadtxt(datafile)[3,:]
    corrosion_prob2 = np.loadtxt(datafile)[4,:]
    ds1_mean = np.loadtxt(datafile)[7,:]
    ds1_std = np.loadtxt(datafile)[8,:]
    ds2_mean = np.loadtxt(datafile)[9,:]
    ds2_std = np.loadtxt(datafile)[10,:]
    flexure_mean = np.loadtxt(datafile)[20,:]
    flexure_std = np.loadtxt(datafile)[21,:]
    shear_mean = np.loadtxt(datafile)[22,:]
    shear_std = np.loadtxt(datafile)[23,:]

    p1t_array = corrosion_prob1[np.in1d(time_span, strengthen_time)]
    p2t_array = corrosion_prob2[np.in1d(time_span, strengthen_time)]
    ds1t_mean_array = ds1_mean[np.in1d(time_span, strengthen_time)]
    ds1t_std_array = ds1_std[np.in1d(time_span, strengthen_time)]
    ds2t_mean_array = ds2_mean[np.in1d(time_span, strengthen_time)]
    ds2t_std_array = ds2_std[np.in1d(time_span, strengthen_time)]

    # nominal and sample values of dst and Ast
    ds1t_nom_array = ds1t_mean_array
    ds2t_nom_array = ds2t_mean_array
    ds1t_smp_array = []
    ds2t_smp_array = []
    for p1t, p2t, ds1t_mean, ds1t_std, ds2t_mean, ds2t_std\
        in zip(p1t_array, p2t_array, ds1t_mean_array, ds1t_std_array,
            ds2t_mean_array, ds2t_std_array):
        # region 1
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                corroded_ds1_mean = 1./p1t * (ds1t_mean - (1.-p1t)*ds_region1_mean)
                corroded_ds1_std = np.sqrt(ds1t_std**2/p1t**2)
                corroded_ds1 = Lognormal('corroded_ds1', corroded_ds1_mean, corroded_ds1_std)
                corroded_ds1_smp = corroded_ds1.rv.rvs(size=N_LHS_SMP)
                intact_no = int(N_LHS_SMP*(1-p1t))
                corroded_ds1_smp[:intact_no] = ds_region1_mean
                ds1t_smp = corroded_ds1_smp
                ds1t_smp_array.append(ds1t_smp)
            except (ValueError, Warning):
                ds1t_smp = np.ones(N_LHS_SMP)*ds_region1_mean
                ds1t_smp_array.append(ds1t_smp)
            # region 2
            try:
                corroded_ds2_mean = 1./p2t * (ds2t_mean - (1.-p2t)*ds_region2_mean)
                corroded_ds2_std = np.sqrt(ds2t_std**2/p2t**2)
                corroded_ds2 = Lognormal('corroded_ds2', corroded_ds2_mean, corroded_ds2_std)
                corroded_ds2_smp = corroded_ds2.rv.rvs(size=N_LHS_SMP)
                intact_no = int(N_LHS_SMP*(1-p2t))
                corroded_ds2_smp[:intact_no] = ds_region2_mean
                ds2t_smp = corroded_ds2_smp
                ds2t_smp_array.append(ds2t_smp)
            except (ValueError, Warning):
                ds2t_smp = np.ones(N_LHS_SMP)*ds_region2_mean
                ds2t_smp_array.append(ds2t_smp)
    ds1t_smp_array = np.array(ds1t_smp_array)
    ds2t_smp_array = np.array(ds2t_smp_array)
    Ast_nom_array = np.pi/4*ds1t_nom_array**2*region1_no + np.pi/4*ds2t_nom_array**2*region2_no
    Ast_nom_array = Ast_nom_array[:, np.newaxis]
    Ast_smp_array = np.pi/4*ds1t_smp_array**2*region1_no + np.pi/4*ds2t_smp_array**2*region2_no

    corrosion_loss = {'ds1t_nom_array':ds1t_nom_array, 'ds2t_nom_array':ds2t_nom_array,
        'ds1t_smp_array':ds1t_smp_array, 'ds2t_smp_array':ds2t_smp_array,
        'Ast_nom_array':Ast_nom_array, 'Ast_smp_array':Ast_smp_array}

    return corrosion_loss

def strengtheningFromCorrosionLoss(corrosion_loss_dict, type_of_component):

    frpBeamDesign = generateFrpBeamForDesign(type_of_component)
    frpBeamMc = generateFrpBeamForMC(type_of_component)

    ds1t_nom_array = corrosion_loss_dict['ds1t_nom_array']
    ds2t_nom_array = corrosion_loss_dict['ds2t_nom_array']
    ds1t_smp_array = corrosion_loss_dict['ds1t_smp_array']
    ds2t_smp_array = corrosion_loss_dict['ds2t_smp_array']
    Ast_nom_array = corrosion_loss_dict['Ast_nom_array']
    Ast_smp_array = corrosion_loss_dict['Ast_smp_array']

    # get initial design capacity
    if 'f' in type_of_component.lower():
        #frpBeamDesign.tfrp = np.zeros(frpBeamDesign.tfrp.shape)
        #initial_design_capacity = frpBeamDesign.flexCapacityDesign()
        initial_design_capacity = 1.25*M_DCDW_NOM + 1.75*(M_LLIM_NOM+M_DESIGN_LANE)*LDF_FLEX
        # get mean/design ratio array
        tfrp_list = []
        ier_list = []
        mean_strength_list = []
        print 'print tfrp for {} strengthening'.format(type_of_component)
        for ds1t_nom,ds2t_nom,ds1t_smp,ds2t_smp,Ast_nom,Ast_smp \
            in zip(ds1t_nom_array,ds2t_nom_array,ds1t_smp_array,ds2t_smp_array,
                Ast_nom_array,Ast_smp_array):
            frpBeamDesign.ds = ds1t_nom
            frpBeamDesign.Ast = Ast_nom
            frpBeamMc.ds = ds1t_smp
            frpBeamMc.Ast = Ast_smp
            # design FRP strengthening
            design_frp = lambda tfrp: designFlexStrength(frpBeamDesign,
                tfrp) - initial_design_capacity
            tfrp, infodict, ier, mesg = fsolve(design_frp, 0.11, full_output=True)
            print '{} mm'.format(tfrp[0])
            #tfrp = 0.0 if tfrp[0]<0 else tfrp[0]
            if tfrp[0]<0:
                tfrp = 0.0
            elif tfrp[0]>MAX_LAYER * 0.166:
                tfrp = MAX_LAYER * 0.166
            else:
                tfrp = tfrp[0]
            tfrp_list.append(tfrp)
            ier_list.append(ier)
            # set tfrpv for mc beam object
            frpBeamDesign.tfrp = np.array([tfrp])
            frpBeamMc.tfrp = tfrp * np.ones(N_LHS_SMP)
            # get mean
            mufrp_smp = frpBeamMc.flexCapacity()
            mean_strength = np.mean(mufrp_smp)
            mean_strength_list.append(mean_strength)
    elif 's' in type_of_component.lower():
        frpBeamDesign.tfrpv = np.zeros(frpBeamDesign.tfrpv.shape)
        #initial_design_capacity = frpBeamDesign.shearCapacityDesign()
        initial_design_capacity = 1.25*V_DCDW_NOM + 1.75*(V_LLIM_NOM+V_DESIGN_LANE)*LDF_SHEAR
        # get mean/design ratio array
        tfrp_list = []
        ier_list = []
        mean_strength_list = []
        print 'print tfrpv for {} strengthening'.format(type_of_component)
        for ds1t_nom,ds2t_nom,ds1t_smp,ds2t_smp,Ast_nom,Ast_smp \
            in zip(ds1t_nom_array,ds2t_nom_array,ds1t_smp_array,ds2t_smp_array,
                Ast_nom_array,Ast_smp_array):
            frpBeamDesign.dsv = ds1t_nom
            frpBeamDesign.Asvt = Ast_nom
            frpBeamMc.dsv = ds1t_smp
            frpBeamMc.Asvt = Ast_smp
            # design FRP strengthening
            design_frp = lambda tfrpv: designShearStrength(frpBeamDesign,
                tfrpv) - initial_design_capacity
            tfrpv, infodict, ier, mesg = fsolve(design_frp, 0.001, full_output=True)
            print '{} mm'.format(tfrpv[0])
            #tfrpv = 0.0 if tfrpv[0]<0 else tfrpv[0]
            if tfrpv[0]<0:
                tfrpv = 0.0
            elif tfrpv[0]>MAX_LAYER * 0.166:
                tfrpv = MAX_LAYER * 0.166
            else:
                tfrpv = tfrpv[0]
            tfrp_list.append(tfrpv)
            ier_list.append(ier)
            # set tfrpv for mc beam object
            frpBeamDesign.tfrpv = np.array([tfrpv])
            frpBeamMc.tfrpv = tfrpv * np.ones(N_LHS_SMP)
            # get mean
            vufrp_smp = frpBeamMc.shearCapacity()
            mean_strength = np.mean(vufrp_smp)
            mean_strength_list.append(mean_strength)
    elif 'd' in type_of_component.lower():
        #frpBeamDesign.tfrp = np.zeros(frpBeamDesign.tfrp.shape)
        #initial_design_capacity = frpBeamDesign.flexCapacityDesign()
        initial_design_capacity = 1.25*M_DCDW_DECK_NOM + 1.75*(1.2*M_LLIM_DECK_NOM)
        # get mean/design ratio array
        tfrp_list = []
        ier_list = []
        mean_strength_list = []
        frpBeamDesign.bfrp = frpBeamDesign.b/FRP_DECK_RATIO
        frpBeamMc.bfrp = frpBeamMc.b/FRP_DECK_RATIO
        print 'print tfrp for {} strengthening'.format(type_of_component)
        for ds1t_nom,ds2t_nom,ds1t_smp,ds2t_smp,Ast_nom,Ast_smp \
            in zip(ds1t_nom_array,ds2t_nom_array,ds1t_smp_array,ds2t_smp_array,
                Ast_nom_array,Ast_smp_array):
            frpBeamDesign.ds = ds1t_nom
            frpBeamDesign.Ast = Ast_nom
            frpBeamMc.ds = ds1t_smp
            frpBeamMc.Ast = Ast_smp
            # design FRP strengthening
            design_frp = lambda tfrp: designFlexStrength(frpBeamDesign,
                tfrp) - initial_design_capacity
            tfrp, infodict, ier, mesg = fsolve(design_frp, 0.11, full_output=True)
            print '{} mm'.format(tfrp[0])
            #tfrp = 0.0 if tfrp[0]<0 else tfrp[0]
            if tfrp[0]<0:
                tfrp = 0.0
            elif tfrp[0]>MAX_LAYER * 0.166:
                tfrp = MAX_LAYER * 0.166
            else:
                tfrp = tfrp[0]
            tfrp_list.append(tfrp)
            ier_list.append(ier)
            # set tfrpv for mc beam object
            frpBeamDesign.tfrp = np.array([tfrp])
            frpBeamMc.tfrp = tfrp * np.ones(N_LHS_SMP)
            # get mean
            mufrp_smp = frpBeamMc.flexCapacity()
            mean_strength = np.mean(mufrp_smp)
            mean_strength_list.append(mean_strength)
    else:
        print '[ERROR:] illegal type of component, must be flexure, shear or deck'
        sys.exit(1)

    return np.array(mean_strength_list), np.array(tfrp_list)

def getDeteriorationFromFile(datafile, strengthen_time, mean_strength_array, type_of_component):

    time_span = np.loadtxt(datafile)[0,:]
    frp_mean_history_list = []
    frp_cov_history_list = []
    for timepoint,mean_strength in zip(strengthen_time,mean_strength_array):
        start_indx = np.where(time_span==timepoint)[0][0]
        end_indx = np.where(time_span==timepoint+SERVICE_LIFE)[0][0]+1
        if 'f' in type_of_component.lower():
            rc_mean_history = np.loadtxt(datafile)[20,start_indx:end_indx]
            rc_std_history = np.loadtxt(datafile)[21,start_indx:end_indx]
            frp_contribution = mean_strength - rc_mean_history[0]
            frp_mean_history = rc_mean_history + frp_contribution
            frp_cov_history = np.ones(frp_mean_history.shape) * FRP_FLEX_COV
        elif 's' in type_of_component.lower():
            rc_mean_history = np.loadtxt(datafile)[22,start_indx:end_indx]
            rc_std_history = np.loadtxt(datafile)[23,start_indx:end_indx]
            frp_contribution = mean_strength - rc_mean_history[0]
            frp_mean_history = rc_mean_history + frp_contribution
            frp_cov_history = np.ones(frp_mean_history.shape) * FRP_SHEAR_COV
        elif 'd' in type_of_component.lower():
            rc_mean_history = np.loadtxt(datafile)[20,start_indx:end_indx]
            rc_std_history = np.loadtxt(datafile)[21,start_indx:end_indx]
            frp_contribution = mean_strength - rc_mean_history[0]
            frp_mean_history = rc_mean_history + frp_contribution
            frp_cov_history = np.ones(frp_mean_history.shape) * FRP_DECK_COV
        else:
            print '[ERROR:] illegal type of component, must be flexure, shear or deck'
            sys.exit(1)


        frp_mean_history_list.append(frp_mean_history)
        frp_cov_history_list.append(frp_cov_history)

    return np.array(frp_mean_history_list), np.array(frp_cov_history_list)
