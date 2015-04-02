# simple corrosion according to Enright 1998
import os
import sys
import numpy as np
import scipy.special as spc
from scipy.optimize import fsolve

from pyre.distributions import *
from management.frpDesignRatio import generateFrpBeamForMC, generateFrpBeamForDesign,\
    generateRcBeamForMC, designFlexStrength, designShearStrength
#from constants.beamConstants import MAX_LAYER, FRP_DESIGN_YR
#from constants.beamConstants import M_DCDW_NOM, M_LLIM_NOM, M_DESIGN_LANE, LDF_FLEX
#from constants.beamConstants import V_DCDW_NOM, V_LLIM_NOM, V_DESIGN_LANE, LDF_SHEAR
#from constants.beamConstants import M_DCDW_DECK_NOM, M_LLIM_DECK_NOM
#from constants.beamConstants import BFRP_NOM, SPAN_NOM, PLATE2SUPPORT_NOM
#from constants.beamConstants import WFRPV_NOM, BEAM_HEIGHT_NOM, WEB_WIDTH_NOM, SFRPV_NOM
#from constants.beamConstants import FRP_DECK_RATIO, SPAN_DECK_NOM, PLATE2SUPPORT_DECK_NOM
from constants.beamConstants import *
from constants.manageConstants import GIRDER_NUM, FRP_UNIT_PRICE
from constants.corrosionConstants import N_LHS_SMP
from constants.simpleCorrosionConstants import *

def simpleCorrosionLHS(component_type, service_time, icorr_mean, str_yr=None,
        fix_cov=True):

    if str_yr is None or str_yr < FRP_DESIGN_YR:
        str_yr = 0

    # beam object
    frpBeam = generateFrpBeamForMC(component_type)
    frpBeamDesign = generateFrpBeamForDesign(component_type)
    # get initial_design_capacity
    if 'f' in component_type.lower() or 'd' in component_type.lower():
        initial_design_capacity = designFlexStrength(frpBeamDesign, np.array([0]))
    else:
        initial_design_capacity = designShearStrength(frpBeamDesign, np.array([0]))

    # initial variables
    if 'f' in component_type.lower():
        cover_mean = FLEX_COVER_MEAN
        cover_cov = FLEX_COVER_COV
        distance12 = FLEX_DISTANCE12
        region1_no = FLEX_REGION1_NO
        region2_no = FLEX_REGION2_NO
        ds1_mean = FLEX_DS1_MEAN
        ds2_mean = FLEX_DS2_MEAN
        ds1_cov = FLEX_DS1_COV
        ds2_cov = FLEX_DS2_COV
    elif 's' in component_type.lower():
        cover_mean = SHEAR_COVER_MEAN
        cover_cov = SHEAR_COVER_COV
        distance12 = SHEAR_DISTANCE12
        region1_no = SHEAR_REGION1_NO
        region2_no = SHEAR_REGION2_NO
        ds1_mean = SHEAR_DS1_MEAN
        ds2_mean = SHEAR_DS2_MEAN
        ds1_cov = SHEAR_DS1_COV
        ds2_cov = SHEAR_DS2_COV
    elif 'd' in component_type.lower():
        cover_mean = DECK_COVER_MEAN
        cover_cov = DECK_COVER_COV
        distance12 = DECK_DISTANCE12
        region1_no = DECK_REGION1_NO
        region2_no = DECK_REGION2_NO
        ds1_mean = DECK_DS1_MEAN
        ds2_mean = DECK_DS2_MEAN
        ds1_cov = DECK_DS1_COV
        ds2_cov = DECK_DS2_COV

    ## MC with LHS of structural deterioration
    # percentage array for LHS
    cdf_lhs = np.linspace(0.5/N_LHS_SMP, 1.-0.5/N_LHS_SMP, num=N_LHS_SMP)
    # concrete cover
    concrete_cover = Lognormal('Xc', cover_mean, cover_mean*cover_cov)
    xc1_smp = concrete_cover.rv.ppf(cdf_lhs)
    xc1_smp = np.random.permutation(xc1_smp)
    xc2_smp = xc1_smp + distance12
    # diffusion coefficient
    diffusion_coefficient = Lognormal('Dc', DC_MEAN, DC_MEAN*DC_COV)
    dc_smp = diffusion_coefficient.rv.ppf(cdf_lhs)
    dc_smp = np.random.permutation(dc_smp)
    # surface chloride
    surface_chloride = Lognormal('C0', C0_MEAM, C0_MEAM*C0_COV)
    c0_smp = surface_chloride.rv.ppf(cdf_lhs)
    c0_smp = np.random.permutation(c0_smp)
    # critical chloride
    critical_chloride = Lognormal('Ccr', CCR_MEAN, CCR_MEAN*CCR_COV)
    ccr_smp = critical_chloride.rv.ppf(cdf_lhs)
    ccr_smp = np.random.permutation(ccr_smp)
    # corrosion initiation time
    ti1_smp = xc1_smp**2 / (4.*dc_smp) * (spc.erfinv((c0_smp-ccr_smp)/c0_smp))**(-2)
    ti2_smp = xc2_smp**2 / (4.*dc_smp) * (spc.erfinv((c0_smp-ccr_smp)/c0_smp))**(-2)
    # initial reinforcement
    rebar_region1 = Lognormal('ds1', ds1_mean, ds1_cov)
    ds1_smp = rebar_region1.rv.ppf(cdf_lhs)
    ds1_smp = np.random.permutation(ds1_smp)
    rebar_region2 = Lognormal('ds2', ds2_mean, ds2_cov)
    ds2_smp = rebar_region2.rv.ppf(cdf_lhs)
    ds2_smp = np.random.permutation(ds2_smp)
    resistance_mean = []
    resistance_cov = []
    # corrosion rate
    corrosion_current = Lognormal('icorr', icorr_mean, icorr_mean*ICORR_COV)
    icorr_smp = corrosion_current.rv.ppf(cdf_lhs)
    icorr_smp = np.random.permutation(icorr_smp)
    corrosion_rate = icorr_smp * 0.0232
    for indx,t in enumerate(service_time+str_yr):
        # at the time of strengthening
        str_corrosion1 = str_yr - ti1_smp
        str_corrosion1[str_corrosion1<0] = 0
        str_corrosion2 = str_yr - ti2_smp
        str_corrosion2[str_corrosion2<0] = 0
        # residual reinforcement
        ds1str_smp = ds1_smp - str_corrosion1*corrosion_rate
        ds1str_smp[ds1str_smp<0] = 0
        ds2str_smp = ds2_smp - str_corrosion2*corrosion_rate
        ds2str_smp[ds2str_smp<0] = 0
        As_str_smp = np.pi/4 * (ds1str_smp**2*region1_no + ds2str_smp**2*region2_no)

        # time after corrosion
        t_corrosion1 = t+str_yr - ti1_smp
        t_corrosion1[t_corrosion1<0] = 0
        t_corrosion2 = t+str_yr - ti2_smp
        t_corrosion2[t_corrosion2<0] = 0
        # residual reinforcement
        ds1t_smp = ds1_smp - t_corrosion1*corrosion_rate
        ds1t_smp[ds1t_smp<0] = 0
        ds2t_smp = ds2_smp - t_corrosion2*corrosion_rate
        ds2t_smp[ds2t_smp<0] = 0
        Ast_smp = np.pi/4 * (ds1t_smp**2*region1_no + ds2t_smp**2*region2_no)
        if str_yr == 0:
            ME_flex = Lognormal('ME_flex', ME_FLEX_RC_MEAN, ME_FLEX_RC_MEAN*ME_FLEX_RC_COV)
            ME_flex_smp = ME_flex.rv.ppf(cdf_lhs)
            ME_flex_smp = np.random.permutation(ME_flex_smp)
            ME_shear = Lognormal('ME_shear', ME_SHEAR_RC_MEAN, ME_SHEAR_RC_MEAN*ME_SHEAR_RC_COV)
            ME_shear_smp = ME_shear.rv.ppf(cdf_lhs)
            ME_shear_smp = np.random.permutation(ME_shear_smp)
            frpBeam.tfrp = np.zeros(N_LHS_SMP)
            frpBeam.tfrpv = np.zeros(N_LHS_SMP)
            frpBeam.ME_flex = ME_flex_smp
            frpBeam.ME_shear = ME_shear_smp
            if 'f' in component_type.lower() or 'd' in component_type.lower():
                frpBeam.ds = ds1t_smp
                frpBeam.Ast = Ast_smp
                rsmp = frpBeam.flexCapacity()
                if fix_cov == True:
                    rcov = RC_FLEX_COV
            else:
                frpBeam.dsv = ds1t_smp
                frpBeam.Asvt = Ast_smp
                rsmp = frpBeam.shearCapacity()
                if fix_cov == True:
                    rcov = RC_SHEAR_COV
            cost = 0
        elif 'f' in component_type.lower():
            #initial_design_capacity = 1.25*M_DCDW_NOM + 1.75*(M_LLIM_NOM+M_DESIGN_LANE)*LDF_FLEX
            if indx==0:    # at the first iteration, calculate tfrp
                frpBeamDesign.ds = np.array([np.mean(ds1str_smp)])
                frpBeamDesign.Ast = np.array([np.mean(As_str_smp)])
                frpBeam.ds = ds1str_smp
                frpBeam.Ast = As_str_smp
                # design FRP strengthening
                design_frp = lambda tfrp: designFlexStrength(frpBeamDesign,
                    tfrp) - initial_design_capacity
                tfrp, infodict, ier, mesg = fsolve(design_frp, 0.11, full_output=True)
                if tfrp[0]<0:
                    tfrp = 0.0
                elif tfrp[0]>MAX_LAYER * 0.166:
                    tfrp = MAX_LAYER * 0.166
                else:
                    tfrp = tfrp[0]
                cost = GIRDER_NUM*tfrp*BFRP_NOM*(SPAN_NOM-2*PLATE2SUPPORT_NOM)*FRP_UNIT_PRICE
                frpBeam.tfrp = tfrp * np.ones(N_LHS_SMP)
                # from tfrp deduce model error
                if tfrp < FRP_MIN_THICK:
                    mu = (ME_FLEX_FRPIC_MEAN-ME_FLEX_RC_MEAN)/FRP_MIN_THICK*tfrp + ME_FLEX_RC_MEAN
                    cov = (ME_FLEX_FRPIC_COV-ME_FLEX_RC_COV)/FRP_MIN_THICK*tfrp + ME_FLEX_RC_COV
                    ME_flex = Lognormal('ME_flex', mu, mu*cov)
                    ME_flex_smp = ME_flex.rv.ppf(cdf_lhs)
                    ME_flex_smp = np.random.permutation(ME_flex_smp)
                    frpBeam.ME_flex = ME_flex_smp
            frpBeamDesign.ds = np.array([np.mean(ds1t_smp)])
            frpBeamDesign.Ast = np.array([np.mean(Ast_smp)])
            frpBeam.ds = ds1t_smp
            frpBeam.Ast = Ast_smp
            rsmp = frpBeam.flexCapacity()
            # use fixed cov
            if fix_cov == True:
                rcov = (FRP_FLEX_COV-RC_FLEX_COV)/FRP_MIN_THICK*tfrp + RC_FLEX_COV
                rcov = np.min((rcov,FRP_FLEX_COV))
        elif 's' in component_type.lower():
            #initial_design_capacity = 1.25*V_DCDW_NOM + 1.75*(V_LLIM_NOM+V_DESIGN_LANE)*LDF_SHEAR
            if indx==0:    # at the first iteration, calculate tfrp
                frpBeamDesign.dsv = np.array([np.mean(ds1str_smp)])
                frpBeamDesign.Asvt = np.array([np.mean(As_str_smp)])
                frpBeam.dsv = ds1str_smp
                frpBeam.Asvt = As_str_smp
                # design FRP strengthening
                design_frp = lambda tfrpv: designShearStrength(frpBeamDesign,
                    tfrpv) - initial_design_capacity
                tfrpv, infodict, ier, mesg = fsolve(design_frp, 0.001, full_output=True)
                if tfrpv[0]<0:
                    tfrpv = 0.0
                elif tfrpv[0]>MAX_LAYER * 0.166:
                    tfrpv = MAX_LAYER * 0.166
                else:
                    tfrpv = tfrpv[0]
                cost = GIRDER_NUM*tfrpv*WFRPV_NOM*(BEAM_HEIGHT_NOM*2+WEB_WIDTH_NOM)*\
                    (SPAN_NOM/3.*2/SFRPV_NOM)*FRP_UNIT_PRICE
                frpBeam.tfrpv = tfrpv * np.ones(N_LHS_SMP)
                # from tfrpv deduce model error
                if tfrpv < FRP_MIN_THICK:
                    mu = (ME_SHEAR_FRPU_MEAN-ME_SHEAR_RC_MEAN)/FRP_MIN_THICK*tfrpv + ME_SHEAR_RC_MEAN
                    cov = (ME_SHEAR_FRPU_COV-ME_SHEAR_RC_COV)/FRP_MIN_THICK*tfrpv + ME_SHEAR_RC_COV
                    ME_shear = Lognormal('ME_shear', mu, mu*cov)
                    ME_shear_smp = ME_shear.rv.ppf(cdf_lhs)
                    ME_shear_smp = np.random.permutation(ME_shear_smp)
                    frpBeam.ME_shear = ME_shear_smp
            frpBeamDesign.dsv = np.array([np.mean(ds1t_smp)])
            frpBeamDesign.Asvt = np.array([np.mean(Ast_smp)])
            frpBeam.dsv = ds1t_smp
            frpBeam.Asvt = Ast_smp
            rsmp = frpBeam.shearCapacity()
            # use fixed cov
            if fix_cov == True:
                rcov = (FRP_SHEAR_COV-RC_SHEAR_COV)/FRP_MIN_THICK*tfrpv + RC_SHEAR_COV
                rcov = np.min((rcov,FRP_SHEAR_COV))
        elif 'd' in component_type.lower():
            #initial_design_capacity = 1.25*M_DCDW_DECK_NOM + 1.75*(1.2*M_LLIM_DECK_NOM)
            frpBeamDesign.bfrp = frpBeamDesign.b/FRP_DECK_RATIO
            frpBeam.bfrp = frpBeam.b/FRP_DECK_RATIO
            if indx==0:    # at the first iteration, calculate tfrp
                frpBeamDesign.ds = np.array([np.mean(ds1str_smp)])
                frpBeamDesign.Ast = np.array([np.mean(As_str_smp)])
                frpBeam.ds = ds1str_smp
                frpBeam.Ast = As_str_smp
                # design FRP strengthening
                design_frp = lambda tfrp: designFlexStrength(frpBeamDesign,
                    tfrp) - initial_design_capacity
                tfrp, infodict, ier, mesg = fsolve(design_frp, 0.11, full_output=True)
                if tfrp[0]<0:
                    tfrp = 0.0
                elif tfrp[0]>MAX_LAYER * 0.166:
                    tfrp = MAX_LAYER * 0.166
                else:
                    tfrp = tfrp[0]
                cost = (GIRDER_NUM-1)*tfrp*(SPAN_NOM/FRP_DECK_RATIO)*\
                    (SPAN_DECK_NOM-2*PLATE2SUPPORT_DECK_NOM)*FRP_UNIT_PRICE
                frpBeam.tfrp = tfrp * np.ones(N_LHS_SMP)
                # from tfrp deduce model error
                if tfrp < FRP_MIN_THICK:
                    # for beam b/bfrp=2, for deck b/bfrp=10, so times 5
                    frp_min_thick = DECK_DEPTH_NOM/(BEAM_DEPTH_NOM/FRP_MIN_THICK)*5
                    mu = (ME_FLEX_FRPIC_MEAN-ME_FLEX_RC_MEAN)/frp_min_thick*tfrp + ME_FLEX_RC_MEAN
                    cov = (ME_FLEX_FRPIC_COV-ME_FLEX_RC_COV)/frp_min_thick*tfrp + ME_FLEX_RC_COV
                    ME_flex = Lognormal('ME_flex', mu, mu*cov)
                    ME_flex_smp = ME_flex.rv.ppf(cdf_lhs)
                    ME_flex_smp = np.random.permutation(ME_flex_smp)
                    frpBeam.ME_flex = ME_flex_smp
            frpBeamDesign.ds = np.array([np.mean(ds1t_smp)])
            frpBeamDesign.Ast = np.array([np.mean(Ast_smp)])
            frpBeam.ds = ds1t_smp
            frpBeam.Ast = Ast_smp
            rsmp = frpBeam.flexCapacity()
            # use fixed cov
            if fix_cov == True:
                # for beam b/bfrp=2, for deck b/bfrp=10, so times 5
                frp_min_thick = DECK_DEPTH_NOM/(BEAM_DEPTH_NOM/FRP_MIN_THICK)*5
                rcov = (FRP_DECK_COV-RC_DECK_COV)/frp_min_thick*tfrp + RC_DECK_COV
                rcov = np.min((rcov,FRP_DECK_COV))

        rmean = np.mean(rsmp)
        if fix_cov == True:
            resistance_mean.append(rmean)
            resistance_cov.append(rcov)
        else:
            rstd = np.std(rsmp)
            rcov = rstd/rmean
            resistance_mean.append(rmean)
            resistance_cov.append(rcov)

    return np.array(resistance_mean), np.array(resistance_cov), cost


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from constants import COVT0_COV
    import time
    import datetime

    print "debugging simple corrosion"
    start_delta_time = time.time()

    #np.random.seed(64)
    # Structural age
    str_yr = input('strengthening time=?')
    comp_type = raw_input('comp_type?')
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    resistance_mean, resistance_cov, cost = simpleCorrosionLHS(comp_type, service_time, 1, str_yr)

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))
    if comp_type == 'flexure':
        tfrp = cost/(GIRDER_NUM*BFRP_NOM*(SPAN_NOM-2*PLATE2SUPPORT_NOM)*FRP_UNIT_PRICE)
    elif comp_type == 'shear':
        tfrp = cost/(GIRDER_NUM*WFRPV_NOM*(BEAM_HEIGHT_NOM*2+WEB_WIDTH_NOM)*
                (SPAN_NOM/3.*2/SFRPV_NOM)*FRP_UNIT_PRICE)
    elif comp_type == 'deck':
        tfrp = cost/((GIRDER_NUM-1)*(SPAN_NOM/FRP_DECK_RATIO)*
                (SPAN_DECK_NOM-2*PLATE2SUPPORT_DECK_NOM)*FRP_UNIT_PRICE)
    print 'strengthenging frp thickness: {:6f} mm'.format(tfrp)

    gd_data = resistance_mean/resistance_mean[0] / np.sqrt(1.+resistance_cov**2)
    gd_param = np.polyfit(service_time, gd_data, 4)
    gd = lambda x: np.poly1d(gd_param)(x)

    at_data = np.sqrt(np.log(resistance_cov**2+1.)/COVT0_COV**2)
    gcov_param = np.polyfit(service_time, at_data, 4)
    gcov = lambda x: np.poly1d(gcov_param)(x)

    plt.plot(service_time, gd_data, 'b+',
            service_time, gd(service_time), 'b-')
    #plt.plot(service_time, at_data, 'b+',
            #service_time, gcov(service_time), 'b-')
    plt.show()
