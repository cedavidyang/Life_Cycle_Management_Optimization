# get relation between diameter loss and the ratio of mean/design
import os
import sys
import numpy as np
import math
import scipy.io
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

from resistance import *
from constants.beamConstants import *
from constants.corrosionConstants import N_SMP, N_LHS_SMP
#from corrosion.initialVariables import rebar_type
from corrosion.corrosionFuncs import concStrengthVariable, concEffEcVariable, concTensileVariable

def generateFrpBeamForDesign(rebar_type):
    """ prepare for the beam design, frp strengthening in beamConstants.py is
        temperarily used
    """
    ## FRP deterioration function
    #mat = scipy.io.loadmat( os.path.join(os.path.abspath('./'),
        #'matlab_data', 'frp_degradation.mat') )
    #tdegrade = mat['timeEnvYR'].flatten()
    #sfrpdegrade = mat['degradeFrp2Env'].flatten()
    #frp_interp1d = interp1d(tdegrade, sfrpdegrade, bounds_error=False,
        #fill_value=sfrpdegrade[-1])
    #frp_degrade = lambda t: frp_interp1d(t)
    # overwrite frp_degrade, degradation not considered
    frp_degrade = lambda t: t*0 + 1.0

    ## bonding deterioratio function
    #bond_degrade = lambda t: BOND_P1*t**4 + BOND_P2*t**3 + BOND_P3*t**2 + BOND_P4*t + BOND_P5
    # overwrite bond_degrade, degradation not considered
    bond_degrade = lambda t: t*0 + 1.0

    nsmp = 1
    # concrete strength
    fc_smp = np.array([FC_NOM])
    ## concrete elastic modulus
    #Ec_smp = 4733.0 * np.sqrt(FC_NOM)
    # concrete tensile strength
    ft_smp = np.array([0.6227 * np.sqrt(FC_NOM)])
    # model error of flexural strength of RC beams
    ME_flex_smp = np.array([1.0])
    # yielding strength of steel
    fy_smp = np.array([FSY_NOM])
    if 'd' in rebar_type.lower():
        # beam width and depth
        b_smp = np.array([DECK_WIDTH_NOM])
        d_smp = np.array([DECK_DEPTH_NOM])
        # flange width and depth
        bf_smp = np.array([DECK_WIDTH_NOM])
        hf_smp = np.array([DECK_DEPTH_NOM])
    else:
        # beam width and depth
        b_smp = np.array([WEB_WIDTH_NOM])
        d_smp = np.array([BEAM_DEPTH_NOM])
        # flange width and depth
        bf_smp = np.array([FLANGE_WIDTH_NOM])
        hf_smp = np.array([FLANGE_THICK_NOM])
    # model error of shear strength of RC beams
    ME_shear_smp = np.array([1.0])
    # yielding strength of shear reinforcement
    fyv_smp = np.array([FSYV_NOM])
    # shear depth
    dv_smp = np.array([BEAM_DEPTH_NOM])
    # shear interval
    sv_smp = np.array([SHEAR_INTERVAL_NOM])
    ## FRP-related random variables
    # model error of flexural strengthening
    ME_flex_frp_smp = np.array([1.0])
    # model error of shear strengthening
    ME_shear_frp_smp = np.array([1.0])
    # material properties
    ecu_smp = np.array([ECU_NOM])
    Es_smp = np.array([ES_NOM])
    fyc_smp = np.array([FYC_NOM])
    # flexural strengthening
    Efrp_smp = np.array([EFRP_NOM])
    ffrp_smp = np.array([FFRP_NOM])
    tfrp_smp = np.array([TFRP_MEAN])
    bfrp_smp = np.array([BFRP_MEAN])

    if EM_FLEX_FORM.lower() == 'ic':
        fanchor_smp = np.array([0.0])
    else:
        fanchor_smp = np.array([FANCHOR_NOM])

    eini_smp = np.array([INI_FRP_STRAIN])

    # shear strengthening
    Efrpv_smp = np.array([EFRPV_NOM])
    ffrpv_smp = np.array([FFRPV_NOM])
    tfrpv_smp = np.array([TFRPV_NOM])
    # geometric properties
    if 'd' in rebar_type.lower():
        h_smp = d_smp +  DECK_HEIGHT_NOM - DECK_DEPTH_NOM
        a_smp = np.array([PLATE2SUPPORT_DECK_NOM])
        ss_smp = np.array([SHEAR_SPAN_DECK_NOM])
        l_smp = np.array([SPAN_DECK_NOM])
    else:
        h_smp = d_smp +  BEAM_HEIGHT_NOM - BEAM_DEPTH_NOM
        a_smp = np.array([PLATE2SUPPORT_NOM])
        ss_smp = np.array([SHEAR_SPAN_NOM])
        l_smp = np.array([SPAN_NOM])
    dsc_smp = np.array([DSC_NOM])
    Asc_smp = np.pi/4*dsc_smp**2 * DSC_NO
    dfrp_smp = h_smp
    dfrpt_smp = hf_smp
    wfrpv_smp = np.array([WFRPV_NOM])
    sfrpv_smp = np.array([SFRPV_NOM])
    frpIncline_smp = np.array([FRP_INCLINE_NOM]) * np.pi/180.

    if 'd' in rebar_type.lower():
        ds_smp = np.array([DS_REGION1_DECK_NOM])
        ds_smp2 = np.array([DS_REGION2_DECK_NOM])
        ds_prev = np.array([DS_REGION1_DECK_NOM])
        ds_prev2 = np.array([DS_REGION2_DECK_NOM])
    else:
        ds_smp = np.array([DS_REGION1_NOM])
        ds_smp2 = np.array([DS_REGION2_NOM])
        ds_prev = np.array([DS_REGION1_NOM])
        ds_prev2 = np.array([DS_REGION2_NOM])
    dsv_smp = np.array([DS_REGION1_SHEAR_NOM])
    dsv_smp2 = np.array([DS_REGION2_SHEAR_NOM])
    dsv_prev = np.array([DS_REGION1_SHEAR_NOM])
    dsv_prev2 = np.array([DS_REGION2_SHEAR_NOM])

    tAfterStr = 0.0
    if EM_SHEAR_FORM.lower() == 'w':
        # FRP degradation
        ffrpvt_smp = frp_degrade(tAfterStr) * ffrpv_smp
        env_smp = np.ones(nsmp)
    else:
        # bonding degradation
        env_smp = bond_degrade(tAfterStr) * np.ones(nsmp)
        ffrpvt_smp = np.copy(ffrpv_smp)

    # get resistance
    # shear
    Ast_smp = np.pi/4 * (ds_smp**2 * REGION1_NO + ds_smp2**2 * REGION2_NO)
    Asvt_smp = np.pi/4 * (dsv_smp**2 * REGION1_NO_SHEAR +
        dsv_smp2**2 * REGION2_NO_SHEAR)
    ME_dict, material_dict, geo_dict = assembleBeamDict(
        ME_flex_smp, ME_shear_smp,
        fc_smp, LAMBDA_FC, fy_smp, fyv_smp,
        Ast_smp, Asvt_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
    ME_dict, material_dict, geo_dict = assembleFrpBeamDict(
        ME_dict, material_dict, geo_dict, ME_flex_frp_smp, ME_shear_frp_smp,
        ft_smp, FC_NOM/0.8, ecu_smp, Es_smp, fyc_smp,
        Efrp_smp, ffrp_smp, tfrp_smp, bfrp_smp, fanchor_smp, eini_smp,
        Efrpv_smp, tfrpv_smp, ffrpvt_smp, EM_SHEAR_FORM,
        h_smp, a_smp, ss_smp, l_smp, dsc_smp, Asc_smp, dfrp_smp, dfrpt_smp,
        BAR_TYPE, dsv_smp, wfrpv_smp, sfrpv_smp, frpIncline_smp, env_smp)

    frpBeamDesign= FRPBeam('frp_beam_design', nsmp, ME_dict, material_dict, geo_dict)

    return frpBeamDesign

def generateFrpBeamForMC(rebar_type):
    """generate FRP-strengthened beam objects for MC simulation
    """
    ## FRP deterioration function
    #mat = scipy.io.loadmat( os.path.join(os.path.abspath('./'),
        #'matlab_data', 'frp_degradation.mat') )
    #tdegrade = mat['timeEnvYR'].flatten()
    #sfrpdegrade = mat['degradeFrp2Env'].flatten()
    #frp_interp1d = interp1d(tdegrade, sfrpdegrade, bounds_error=False,
        #fill_value=sfrpdegrade[-1])
    #frp_degrade = lambda t: frp_interp1d(t)
    # overwrite frp_degrade, degradation not considered
    frp_degrade = lambda t: t*0 + 1.0

    ## bonding deterioratio function
    #bond_degrade = lambda t: BOND_P1*t**4 + BOND_P2*t**3 + BOND_P3*t**2 + BOND_P4*t + BOND_P5
    # overwrite bond_degrade, degradation not considered
    bond_degrade = lambda t: t*0 + 1.0

    # percentage array for LHS
    cdf_lhs = np.linspace(0.5/N_LHS_SMP, 1.-0.5/N_LHS_SMP, num=N_LHS_SMP)

    # concrete strength
    compressive_strength = concStrengthVariable()
    #np.random.seed(CONCRETE_STRENGTH_SEED)
    #fc_smp = compressive_strength.rv.rvs(size = N_SMP)
    fc_smp = compressive_strength.rv.ppf(cdf_lhs)
    fc_smp = np.random.permutation(fc_smp)

    # effective concrete elastic modulus
    elastic_modulus = concEffEcVariable()
    #np.random.seed(CONCRETE_MODULUS_SEED)
    #Ec_smp = elastic_modulus.rv.rvs(size = N_SMP)
    Ec_smp = elastic_modulus.rv.ppf(cdf_lhs)
    Ec_smp = np.random.permutation(Ec_smp)

    # concrete tensile strength
    tensile_strength = concTensileVariable()
    #np.random.seed(CONCRETE_TENSION_SEED)
    #ft_smp = tensile_strength.rv.rvs(size = N_SMP)
    ft_smp = tensile_strength.rv.ppf(cdf_lhs)
    ft_smp = np.random.permutation(ft_smp)

    # model error of flexural strength of RC beams
    ME_flex = modelErrorRCFlexVariable()
    #np.random.seed(ME_FLEX_RC_SEED)
    #ME_flex_smp = ME_flex.rv.rvs(size=N_SMP)
    ME_flex_smp = ME_flex.rv.ppf(cdf_lhs)
    ME_flex_smp = np.random.permutation(ME_flex_smp)

    # yielding strength of steel
    fy_steel = steelYieldingVariable()
    #np.random.seed(FSY_SEED)
    #fy_smp = fy_steel.rv.rvs(size=N_SMP)
    fy_smp = fy_steel.rv.ppf(cdf_lhs)
    fy_smp = np.random.permutation(fy_smp)

    # beam width and depth
    beam_width = beamWidthVariable(rebar_type)
    #np.random.seed(BEAM_WIDTH_SEED)
    #b_smp = beam_width.rv.rvs(size=N_SMP)
    b_smp = beam_width.rv.ppf(cdf_lhs)
    b_smp = np.random.permutation(b_smp)
    beam_depth = beamDepthVariable(rebar_type)
    #np.random.seed(BEAM_DEPTH_SEED)
    #d_smp = beam_depth.rv.rvs(size=N_SMP)
    d_smp = beam_depth.rv.ppf(cdf_lhs)
    d_smp = np.random.permutation(d_smp)

    # flange width and depth
    flange_width = flangeWidthVariable(rebar_type)
    #np.random.seed(FLANGE_WIDTH_SEED)
    #bf_smp = flange_width.rv.rvs(size=N_SMP)
    bf_smp = flange_width.rv.ppf(cdf_lhs)
    bf_smp = np.random.permutation(bf_smp)
    flange_depth = flangeDepthVariable(rebar_type)
    #np.random.seed(FLANGE_DEPTH_SEED)
    #hf_smp = flange_depth.rv.rvs(size=N_SMP)
    hf_smp = flange_depth.rv.ppf(cdf_lhs)
    hf_smp = np.random.permutation(hf_smp)

    # model error of shear strength of RC beams
    ME_shear = modelErrorRCShearVariable()
    #np.random.seed(ME_SHEAR_RC_SEED)
    #ME_shear_smp = ME_shear.rv.rvs(size=N_SMP)
    ME_shear_smp = ME_shear.rv.ppf(cdf_lhs)
    ME_shear_smp = np.random.permutation(ME_shear_smp)

    # yielding strength of shear reinforcement
    fyv_steel = shearYieldingVariable()
    #np.random.seed(FSYV_SEED)
    #fyv_smp = fyv_steel.rv.rvs(size=N_SMP)
    fyv_smp = fyv_steel.rv.ppf(cdf_lhs)
    fyv_smp = np.random.permutation(fyv_smp)

    # shear depth
    dv = shearDepthVariable()
    #np.random.seed(SHEAR_DEPTH_SEED)
    #dv_smp = dv.rv.rvs(size=N_SMP)
    dv_smp = dv.rv.ppf(cdf_lhs)
    dv_smp = np.random.permutation(dv_smp)

    # shear interval
    sv = shearIntervalVariable()
    #np.random.seed(SHEAR_INTERVAL_SEED)
    #sv_smp = sv.rv.rvs(size=N_SMP)
    sv_smp = sv.rv.ppf(cdf_lhs)
    sv_smp = np.random.permutation(sv_smp)

    ## FRP-related random variables
    # model error of flexural strengthening
    ME_flex_frp = modelErrorFRPFlexVariable()
    #np.random.seed(ME_FLEX_FRP_SEED)
    #ME_flex_frp_smp = ME_flex_frp.rv.rvs(size=N_SMP)
    ME_flex_frp_smp = ME_flex_frp.rv.ppf(cdf_lhs)
    ME_flex_frp_smp = np.random.permutation(ME_flex_frp_smp)
    # model error of shear strengthening
    ME_shear_frp = modelErrorFRPShearVariable()
    #np.random.seed(ME_SHEAR_FRP_SEED)
    #ME_shear_frp_smp = ME_shear_frp.rv.rvs(size=N_SMP)
    ME_shear_frp_smp = ME_shear_frp.rv.ppf(cdf_lhs)
    ME_shear_frp_smp = np.random.permutation(ME_shear_frp_smp)
    # material properties
    ecu_smp = ECU_MEAN * np.ones(N_LHS_SMP)
    Es_smp = ES_MEAN * np.ones(N_LHS_SMP)
    fyc_smp = FYC_MEAN * np.ones(N_LHS_SMP)
    # flexural strengthening
    Efrp = EfrpVariable()
    #np.random.seed(EFRP_SEED)
    #Efrp_smp = Efrp.rv.rvs(size=N_SMP)
    Efrp_smp = Efrp.rv.ppf(cdf_lhs)
    Efrp_smp = np.random.permutation(Efrp_smp)

    ffrp = ffrpVariable()
    #np.random.seed(FFRP_SEED)
    #ffrp_smp = ffrp.rv.rvs(size=N_SMP)
    ffrp_smp = ffrp.rv.ppf(cdf_lhs)
    ffrp_smp = np.random.permutation(ffrp_smp)

    tfrp_smp = TFRP_MEAN * np.ones(N_LHS_SMP)
    bfrp_smp = BFRP_MEAN * np.ones(N_LHS_SMP)

    if EM_FLEX_FORM.lower() == 'ic':
        fanchor_smp = np.zeros(N_LHS_SMP)
    else:
        fanchor_smp = FANCHOR_MEAN * np.ones(N_LHS_SMP)

    eini_smp = INI_FRP_STRAIN * np.ones(N_LHS_SMP)

    # shear strengthening
    Efrpv = EfrpShearVariable()
    #np.random.seed(EFRPV_SEED)
    #Efrpv_smp = Efrpv.rv.rvs(size=N_SMP)
    Efrpv_smp = Efrpv.rv.ppf(cdf_lhs)
    Efrpv_smp = np.random.permutation(Efrpv_smp)

    ffrpv = ffrpShearVariable()
    #np.random.seed(FFRPV_SEED)
    #ffrpv_smp = ffrpv.rv.rvs(size=N_SMP)
    ffrpv_smp = ffrpv.rv.ppf(cdf_lhs)
    ffrpv_smp = np.random.permutation(ffrpv_smp)

    tfrpv_smp = TFRPV_MEAN * np.ones(N_LHS_SMP)

    # geometric properties
    if 'd' in rebar_type.lower():
        h_smp = d_smp +  BEAM_HEIGHT_MEAN - BEAM_DEPTH_MEAN
        a_smp = PLATE2SUPPORT_DECK_MEAN * np.ones(N_LHS_SMP)
        ss_smp = SHEAR_SPAN_DECK_MEAN * np.ones(N_LHS_SMP)
        l_smp = SPAN_DECK_MEAN * np.ones(N_LHS_SMP)
    else:
        h_smp = d_smp +  BEAM_HEIGHT_MEAN - BEAM_DEPTH_MEAN
        a_smp = PLATE2SUPPORT_MEAN * np.ones(N_LHS_SMP)
        ss_smp = SHEAR_SPAN_MEAN * np.ones(N_LHS_SMP)
        l_smp = SPAN_MEAN * np.ones(N_LHS_SMP)
    dsc_smp = DSC_MEAN * np.ones(N_LHS_SMP)
    Asc_smp = np.pi/4*dsc_smp**2 * DSC_NO
    dfrp_smp = np.copy(h_smp)
    dfrpt_smp = hf_smp
    wfrpv_smp = WFRPV_MEAN * np.ones(N_LHS_SMP)
    sfrpv_smp = SFRPV_MEAN * np.ones(N_LHS_SMP)
    frpIncline_smp = FRP_INCLINE_MEAN * np.ones(N_LHS_SMP) * np.pi/180.

    if 'd' in rebar_type.lower():
        ds_smp = np.ones(N_LHS_SMP) * DS_REGION1_DECK_MEAN
        ds_smp2 = np.ones(N_LHS_SMP) * DS_REGION2_DECK_MEAN
        ds_prev = np.ones(N_LHS_SMP) * DS_REGION1_DECK_MEAN
        ds_prev1 = np.ones(N_LHS_SMP) * DS_REGION2_DECK_MEAN
    else:
        ds_smp = np.ones(N_LHS_SMP) * DS_REGION1_MEAN
        ds_smp2 = np.ones(N_LHS_SMP) * DS_REGION2_MEAN
        ds_prev = np.ones(N_LHS_SMP) * DS_REGION1_MEAN
        ds_prev2 = np.ones(N_LHS_SMP) * DS_REGION2_MEAN
    dsv_smp = np.ones(N_LHS_SMP) * DS_REGION1_SHEAR_MEAN
    dsv_smp2 = np.ones(N_LHS_SMP) * DS_REGION2_SHEAR_MEAN
    dsv_prev = np.ones(N_LHS_SMP) * DS_REGION1_SHEAR_MEAN
    dsv_prev2 = np.ones(N_LHS_SMP) * DS_REGION2_SHEAR_MEAN

    tAfterStr = 0.0
    if EM_SHEAR_FORM.lower() == 'w':
        # FRP degradation
        ffrpvt_smp = frp_degrade(tAfterStr) * ffrpv_smp
        env_smp = np.ones(N_LHS_SMP)
    else:
        # bonding degradation
        env_smp = bond_degrade(tAfterStr) * np.ones(N_LHS_SMP)
        ffrpvt_smp = np.copy(ffrpv_smp)

    # get beam object
    Ast_smp = np.pi/4 * (ds_smp**2 * REGION1_NO + ds_smp2**2 * REGION2_NO)
    Asvt_smp = np.pi/4 * (ds_smp**2 * REGION1_NO_SHEAR +
        ds_smp2**2 * REGION2_NO_SHEAR)
    ME_dict, material_dict, geo_dict = assembleBeamDict(
        ME_flex_smp, ME_shear_smp,
        fc_smp, LAMBDA_FC, fy_smp, fyv_smp,
        Ast_smp, Asvt_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
    ME_dict, material_dict, geo_dict = assembleFrpBeamDict(
        ME_dict, material_dict, geo_dict, ME_flex_frp_smp, ME_shear_frp_smp,
        ft_smp, FC_NOM/0.8, ecu_smp, Es_smp, fyc_smp,
        Efrp_smp, ffrp_smp, tfrp_smp, bfrp_smp, fanchor_smp, eini_smp,
        Efrpv_smp, tfrpv_smp, ffrpvt_smp, EM_SHEAR_FORM,
        h_smp, a_smp, ss_smp, l_smp, dsc_smp, Asc_smp, dfrp_smp, dfrpt_smp,
        BAR_TYPE, dsv_smp, wfrpv_smp, sfrpv_smp, frpIncline_smp, env_smp)
    frpBeamMc= FRPBeam('frp_beam_mc', N_LHS_SMP, ME_dict, material_dict, geo_dict)

    return frpBeamMc

def designShearStrength(frpBeamDesign, tfrpv=None, wfrpv=None, sfrpv=None):
    """get design (factored) shear strengthe of FRP-strengthened beams
       (HK model).
       :Args:
       frpBeamDesignShear --- frpBeam object
       tfrpv              --- thickness of FRP for shear
       wfrpv              --- width of FRP for shear
       sfrpv              --- interval of FRP for shear (center to center)
    """
    # change parameters to provided values, more secure methods recommended
    if tfrpv is not None:
        frpBeamDesign.tfrpv = tfrpv
    if wfrpv is not None:
        frpBeamDesign.wfrpv = wfprv
    if sfrpv is not None:
        frpBeamDesign.sfrpv = sfrpv

    vufrp_design = frpBeamDesign.shearCapacityDesign()

    return vufrp_design

def sampleShearStrength(frpBeamMc, tfrpv_smp=None, wfrpv_smp=None, sfrpv_smp=None):
    """get sample values of shear strength of FRP-strengthened beams
       (HK model).
       :Args:
       frpBeamMc      --- frpBeam object
       tfrpv_smp      --- thickness of FRP for shear
       wfrpv_smp      --- width of FRP for shear
       sfrpv_smp      --- interval of FRP for shear (center to center)
    """

    # change parameters to provided values, more secure methods recommended
    if tfrpv_smp is not None:
        frpBeamMc.tfrpv = tfrpv_smp
    if wfrpv_smp is not None:
        frpBeamMc.wfrpv = wfprv_smp
    if sfrpv_smp is not None:
        frpBeamMc.sfrpv = sfrpv_smp

    vufrp_smp = frpBeamMc.shearCapacity()

    return vufrp_smp

def designFlexStrength(frpBeamDesign, tfrp=None):
    """get design (factored) flexural strengthe of FRP-strengthened beams
       (HK model).
       :Args:
       frpBeamDesign --- frpBeam object
       tfrp          --- thickness of FRP for flexure
    """
    # change parameters to provided values, more secure methods recommended
    if tfrp is not None:
        frpBeamDesign.tfrp = tfrp

    mufrp_design = frpBeamDesign.flexCapacityDesign()

    return mufrp_design

def sampleFlexStrength(frpBeamMc, tfrpv_smp=None):
    """get sample values of flexural strength of FRP-strengthened beams
       (HK model).
       :Args:
       frpBeamMc --- frpBeam object
       tfrpv     --- thickness of FRP for flexure
    """
    # change parameters to provided values, more secure methods recommended
    if tfrp_smp is not None:
        frpBeamMc.tfrp = tfrp_smp

    mufrp_smp = frpBeamMc.flexCapacity()

    return mufrp_smp

def generateRcBeamForMC(rebar_type):
    """generate RC beam objects for MC simulation
    """
    # percentage array for LHS
    cdf_lhs = np.linspace(0.5/N_LHS_SMP, 1.-0.5/N_LHS_SMP, num=N_LHS_SMP)

    # concrete strength
    compressive_strength = concStrengthVariable()
    fc_smp = compressive_strength.rv.ppf(cdf_lhs)
    fc_smp = np.random.permutation(fc_smp)

    # model error of flexural strength of RC beams
    ME_flex = modelErrorRCFlexVariable()
    ME_flex_smp = ME_flex.rv.ppf(cdf_lhs)
    ME_flex_smp = np.random.permutation(ME_flex_smp)

    # yielding strength of steel
    fy_steel = steelYieldingVariable()
    fy_smp = fy_steel.rv.ppf(cdf_lhs)
    fy_smp = np.random.permutation(fy_smp)

    # beam width and depth
    beam_width = beamWidthVariable(rebar_type)
    b_smp = beam_width.rv.ppf(cdf_lhs)
    b_smp = np.random.permutation(b_smp)
    beam_depth = beamDepthVariable(rebar_type)
    d_smp = beam_depth.rv.ppf(cdf_lhs)
    d_smp = np.random.permutation(d_smp)

    # flange width and depth
    flange_width = flangeWidthVariable(rebar_type)
    bf_smp = flange_width.rv.ppf(cdf_lhs)
    bf_smp = np.random.permutation(bf_smp)
    flange_depth = flangeDepthVariable(rebar_type)
    hf_smp = flange_depth.rv.ppf(cdf_lhs)
    hf_smp = np.random.permutation(hf_smp)

    # model error of shear strength of RC beams
    ME_shear = modelErrorRCShearVariable()
    ME_shear_smp = ME_shear.rv.ppf(cdf_lhs)
    ME_shear_smp = np.random.permutation(ME_shear_smp)

    # yielding strength of shear reinforcement
    fyv_steel = shearYieldingVariable()
    fyv_smp = fyv_steel.rv.ppf(cdf_lhs)
    fyv_smp = np.random.permutation(fyv_smp)

    # shear depth
    dv = shearDepthVariable()
    dv_smp = dv.rv.ppf(cdf_lhs)
    dv_smp = np.random.permutation(dv_smp)

    # shear interval
    sv = shearIntervalVariable()
    sv_smp = sv.rv.ppf(cdf_lhs)
    sv_smp = np.random.permutation(sv_smp)

    if 'd' in rebar_type.lower():
        ds_smp = np.ones(N_LHS_SMP) * DS_REGION1_DECK_MEAN
        ds_smp2 = np.ones(N_LHS_SMP) * DS_REGION2_DECK_MEAN
    else:
        ds_smp = np.ones(N_LHS_SMP) * DS_REGION1_MEAN
        ds_smp2 = np.ones(N_LHS_SMP) * DS_REGION2_MEAN
    dsv_smp = np.ones(N_LHS_SMP) * DS_REGION1_SHEAR_MEAN
    dsv_smp2 = np.ones(N_LHS_SMP) * DS_REGION2_SHEAR_MEAN

    # get beam object
    Ast_smp = np.pi/4 * (ds_smp**2 * REGION1_NO + ds_smp2**2 * REGION2_NO)
    Asvt_smp = np.pi/4 * (ds_smp**2 * REGION1_NO_SHEAR +
        ds_smp2**2 * REGION2_NO_SHEAR)
    ME_dict, material_dict, geo_dict = assembleBeamDict(
        ME_flex_smp, ME_shear_smp,
        fc_smp, LAMBDA_FC, fy_smp, fyv_smp,
        Ast_smp, Asvt_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)

    rcBeam = RCBeam('rc', ME_dict, material_dict, geo_dict)

    return rcBeam


def main():
    """find mean to design ration of FRP-strengthened girders
    """
    data_path = os.path.join(os.path.abspath('./'), 'data', 'frp_design_ratio')
    datafile = os.path.join(data_path , 'frp_shear_design_ratio.npz')

    frpBeamDesignShear = generateFrpBeamForDesign()
    frpBeamMcShear = generateFrpBeamForMC()

    # check design models with hand calculation
    if os.path.isfile(datafile):
        tfrpv = np.load(datafile)['tfrpv'][18]
        dsvt = np.load(datafile)['dsvt'][18]
        # set design and mc beamobject
        Asvt = np.pi/4 * (dsvt**2 * REGION1_NO_SHEAR + dsvt**2 * REGION2_NO_SHEAR)
        frpBeamDesignShear.dsv = np.array([dsvt])
        frpBeamDesignShear.Asvt = np.array([Asvt])
        frpBeamDesignShear.tfrpv = np.array([tfrpv])

        frpBeamMcShear.dsv = dsvt * np.ones(N_SMP)
        frpBeamMcShear.Asvt = Asvt * np.ones(N_SMP)
        frpBeamMcShear.tfrpv = tfrpv * np.ones(N_SMP)

        vufrp_design = frpBeamDesignShear.shearCapacityDesign()
        vufrp_smp = frpBeamMcShear.shearCapacity()
        ratio = np.mean(vufrp_smp)/vufrp_design
        print 'residual shear reinforcement: {:10.4f}'.format(dsvt)
        print 'tfrpv_array: {:10.4f}'.format(tfrpv)
        print 'vufrp_design: {:10.4e}'.format(vufrp_design[0])
        print 'vufrp_mean: {:10.4e}'.format(np.mean(vufrp_smp))
        print 'ratio: {:6.4f}'.format(ratio[0])
        sys.exit(1)

    # design tfrpv
    FACTORED_SHEAR_LOAD = 624    #kN
    design_tfrpv = lambda tfrpv: design_shear_strength(frpBeamDesignShear,
            tfrpv) - FACTORED_SHEAR_LOAD

    # get ratio for different diameter loss
    residual_diameter_array = np.linspace(DS_REGION1_SHEAR_NOM, 0.5*DS_REGION1_SHEAR_NOM, 50)
    tfrpv_list = []
    ier_list = []
    ratio_list = []
    for dsvt in residual_diameter_array:
        # set dsv for design and mc beam object
        Asvt = np.pi/4 * (dsvt**2 * REGION1_NO_SHEAR + dsvt**2 * REGION2_NO_SHEAR)
        frpBeamDesignShear.dsv = np.array([dsvt])
        frpBeamDesignShear.Asvt = np.array([Asvt])
        frpBeamMcShear.dsv = dsvt * np.ones(N_SMP)
        frpBeamMcShear.Asvt = Asvt * np.ones(N_SMP)
        # design FRP strengthening
        tfrpv, infodict, ier, mesg = fsolve(design_tfrpv, 0.001, full_output=True)
        #tfrpv = 0.0 if tfrpv[0]<0 else tfrpv[0]
        if tfrpv[0]<0:
            tfrpv = 0.0
        elif tfrpv[0]>MAX_LAYER * 0.166:
            tfrpv = MAX_LAYER * 0.166
        else:
            tfrpv = tfrpv[0]
        tfrpv_list.append(tfrpv)
        ier_list.append(ier)
        # set tfrpv for design and mc beam object
        frpBeamDesignShear.tfrpv = np.array([tfrpv])
        frpBeamMcShear.tfrpv = tfrpv * np.ones(N_SMP)
        # get ratio of mean/design
        vufrp_design = frpBeamDesignShear.shearCapacityDesign()
        vufrp_smp = frpBeamMcShear.shearCapacity()
        ratio = np.mean(vufrp_smp)/vufrp_design
        ratio_list.append(ratio[0])

    tfrpv_array = np.array(tfrpv_list)
    ier_array = np.array(ier_list)
    ratio_array = np.array(ratio_list)

    #save to data
    frp_shear_design_ratio = np.savez(datafile, dsvt=residual_diameter_array,
            tfrpv=tfrpv_array, ier=ier_array, ratio=ratio_array)

if __name__ == '__main__':
    main()
