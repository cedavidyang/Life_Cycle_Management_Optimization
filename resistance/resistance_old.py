# flexura and shear resistance of RC and FRP-strengthened RC beams
import sys
import numpy as np

from constants import *
from pyre import *


def modelErrorRCFlexVariable():
    mu = ME_FLEX_RC_MEAN
    sigma = mu * ME_FLEX_RC_COV
    ME_flex = Lognormal('ME_flex', mu, sigma)
    
    return ME_flex
    
    
def steelYieldingVariable():
    mu = FSY_MEAN
    sigma = mu * FSY_COV
    if FSY_DISTR == 'deterministic':
        fy_steel = Deterministic('fy', mu, sigma)
    elif FSY_DISTR == 'lognormal':
        fy_steel = Lognormal('fy', mu, sigma)
    
    return fy_steel
    
    
def beamWidthVariable():
    mu = WEB_WIDTH_MEAN
    sigma = mu * WEB_WIDTH_COV 
    if WEB_WIDTH_DISTR == 'deterministic':
        beam_width = Deterministic('beam_width', mu, sigma)
    elif WEB_WIDTH_DISTR == 'lognormal':
        beam_width = Lognormal('beam_width', mu, sigma)
    
    return beam_width
    
    
def beamDepthVariable():
    mu = BEAM_DEPTH_MEAN
    sigma = mu * BEAM_DEPTH_COV
    if BEAM_DEPTH_DISTR == 'deterministic':
        beam_depth = Deterministic('beam_depth', mu, sigma)
    elif BEAM_DEPTH_DISTR == 'normal':
        beam_depth = Normal('beam_depth', mu, sigma)
    
    
    return beam_depth
    
    
def flangeWidthVariable():
    mu = FLANGE_WIDTH_MEAN
    sigma = mu * FLANGE_WIDTH_COV 
    if WEB_WIDTH_DISTR == 'deterministic':
        flange_width = Deterministic('flange_width', mu, sigma)
    elif WEB_WIDTH_DISTR == 'lognormal':
        flange_width = Lognormal('flange_width', mu, sigma)
    
    return flange_width
    
    
def flangeDepthVariable():
    mu = FLANGE_THICK_MEAN
    sigma = mu * FLANGE_THICK_COV
    if BEAM_DEPTH_DISTR == 'deterministic':
        flange_depth = Deterministic('flange_depth', mu, sigma)
    elif BEAM_DEPTH_DISTR == 'normal':
        flange_depth = Normal('flange_depth', mu, sigma)
    
    
    return flange_depth
    
    
def residualSteelArea(section_loss_smp, section_loss_smp2):
    As0_region1_onebar = np.pi/4 * DS_REGION1_MEAN**2    #[mm2]
    As0_region2_onebar = np.pi/4 * DS_REGION2_MEAN**2
    Ast_smp = REGION1_NO * (As0_region1_onebar - section_loss_smp) + REGION2_NO *(As0_region2_onebar - section_loss_smp2)
    
    return Ast_smp
    
    
def modelErrorRCShearVariable():
    mu = ME_SHEAR_RC_MEAN
    sigma = mu * ME_SHEAR_RC_COV
    ME_flex = Lognormal('ME_shear', mu, sigma)
    
    return ME_flex
    
    
def shearYieldingVariable():
    mu = FSYV_MEAN
    sigma = mu*FSYV_COV
    fsyv_steel = Lognormal('fsyv_steel', mu, sigma)
    
    return fsyv_steel
    
    
def shearDepthVariable():
    mu = SHEAR_DEPTH_MEAN
    sigma = mu*SHEAR_DEPTH_COV
    dv = Normal('dv', mu, sigma)
    
    return dv 
    
       
def shearIntervalVariable():
    mu = SHEAR_INTERVAL_MEAN
    sigma = mu*SHEAR_INTERVAL_COV
    if SHEAR_INTERVAL_DISTR == 'normal':
        sv = Normal('sv', mu, sigma)
    elif SHEAR_INTERVAL_DISTR == 'deterministic':
        sv = Deterministic('sv', mu, sigma)
    
    return sv     


def modelErrorFRPShearVariable():
    if EM_SHEAR_FORM.lower() == 'side' or EM_SHEAR_FORM .lower() == 's':
        mu = ME_SHEAR_FRPS_MEAN
        sigma = mu*ME_SHEAR_FRPS_COV
        if ME_SHEAR_FRPS_DISTR == 'normal':
            me = Normal('me_shear_frp', mu, sigma)
        elif ME_SHEAR_FRPS_DISTR == 'lognormal':
            me = Lognormal('me_shear_frp', mu, sigma)
    elif EM_SHEAR_FORM.lower() == 'u':
        mu = ME_SHEAR_FRPU_MEAN
        sigma = mu*ME_SHEAR_FRPU_COV
        if ME_SHEAR_FRPU_DISTR == 'normal':
            me = Normal('me_shear_frp', mu, sigma)
        elif ME_SHEAR_FRPU_DISTR == 'lognormal':
            me = Lognormal('me_shear_frp', mu, sigma)
    elif EM_SHEAR_FORM.lower() == 'w':
        mu = ME_SHEAR_FRPW_MEAN
        sigma = mu*ME_SHEAR_FRPW_COV
        if ME_SHEAR_FRPW_DISTR == 'normal':
            me = Normal('me_shear_frp', mu, sigma)
        elif ME_SHEAR_FRPW_DISTR == 'lognormal':
            me = Lognormal('me_shear_frp', mu, sigma)
    else:
        print '[Error:] unknown shear strengthening form'
        sys.exit(1)
        
    return me
    
    
def modelErrorFRPFlexVariable():
    if EM_FLEX_FORM.lower() == 'ic':
        mu = ME_FLEX_FRPIC_MEAN
        sigma = mu*ME_FLEX_FRPIC_COV
        if ME_FLEX_FRPIC_DISTR == 'normal':
            me = Normal('me_flex_frp', mu, sigma)
        elif ME_FLEX_FRPIC_DISTR == 'lognormal':
            me = Lognormal('me_flex_frp', mu, sigma)
    elif EM_SHEAR_FORM.lower() == 'anchored' or EM_SHEAR_FORM.lower() == 'rupture':
        mu = ME_FLEX_FRPRT_MEAN
        sigma = mu*ME_FLEX_FRPRT_COV
        if ME_FLEX_FRPRT_DISTR == 'normal':
            me = Normal('me_flex_frp', mu, sigma)
        elif ME_FLEX_FRPRT_DISTR == 'lognormal':
            me = Lognormal('me_flex_frp', mu, sigma)
    else:
        print '[Error:] unknown shear strengthening form'
        sys.exit(1)
        
    return me
    

def EfrpVariable():
    mu = EFRP_MEAN
    sigma = mu*EFRP_COV
    if EFRP_DISTR == 'lognormal':
        Efrp = Lognormal('Efrp', mu, sigma)
    elif EFRP_DISTR == 'deterministic':
        Efrp = Deterministic('Efrp', mu, sigma)
    
    return Efrp
    
      
def ffrpVariable():
    mu = FFRP_MEAN
    sigma = mu*FFRP_COV
    if FFRP_DISTR == 'lognormal':
        ffrp = Lognormal('ffrp', mu, sigma)
    elif FFRP_DISTR == 'weibull':
        ffrp = Weibull('ffrp', mu, sigma)
    elif FFRP_DISTR == 'deterministic':
        ffrp = Deterministic('ffrp', mu, sigma)
    
    return ffrp


def EfrpShearVariable():
    mu = EFRPV_MEAN
    sigma = mu*EFRPV_COV
    if EFRP_DISTR == 'lognormal':
        Efrpv = Lognormal('Efrpv', mu, sigma)
    elif EFRP_DISTR == 'deterministic':
        Efrpv = Deterministic('Efrpv', mu, sigma)
    
    return Efrpv
    
      
def ffrpShearVariable():
    mu = FFRPV_MEAN
    sigma = mu*FFRPV_COV
    if FFRP_DISTR == 'lognormal':
        ffrp = Lognormal('ffrpvv', mu, sigma)
    elif FFRP_DISTR == 'weibull':
        ffrp = Weibull('ffrpv', mu, sigma)
    elif FFRP_DISTR == 'deterministic':
        ffrp = Deterministic('ffrpv', mu, sigma)
    
    return ffrp  
        
        
    
def assembleBeamDict(ME_flex, ME_shear, fc, fy, fyv, Ast, Asvt, b, d, bf, hf, dv, sv):
    ME_dict = {'ME_flex': ME_flex, 'ME_shear':ME_shear}
    material_dict = {'fc':fc, 'fy':fy, 'fyv':fyv}
    geo_dict = {'b':b, 'd':d, 'bf':bf, 'hf':hf, 'dv':dv, 'sv':sv, 'Ast':Ast, 'Asvt':Asvt}
    
    return ME_dict, material_dict, geo_dict
    
    
def assembleFrpBeamDict(ME_dict, material_dict, geo_dict, ME_flex, ME_shear,
                        ft, ecu, Es, fyc, Efrp, ffrp, tfrp, bfrp, fanchor, eini,
                        Efrpv, tfrpv, ffrpv, strform,
                        h, a, ss, l, dc, Asc, dfrp, dfrpt, barType, dsv, wfrpv, sfrpv, frpIncline):
        #self.Efrpv  = np.array(material_dict['Efrpv'])
        #self.tfrpv  = np.array(material_dict['tfrpv'])
        #self.ffrpv  = np.array(material_dict['ffrpv'])
        #self.strform = np.array(material_dict['strfrom'])
    # re-assign model errors
    ME_dict['ME_flex'] = ME_flex
    ME_dict['ME_shear'] = ME_shear
    # add material properties
    new_material_dict = {'ft':ft, 'ecu':ecu, 'Es':Es, 'fyc':fyc, 'Efrp':Efrp,
                         'ffrp':ffrp, 'tfrp':tfrp, 'bfrp':bfrp, 'fanchor': fanchor, 'eini':eini,    # for FRP flexural strengthening
                         'Efrpv':Efrpv, 'tfrpv':tfrpv, 'ffrpv':ffrpv, 'strform': strform}
    material_dict.update(new_material_dict)
    # add geometric properties
    new_ME_dict = {'h':h, 'a':a, 'ss':ss, 'l':l, 'dc':dc, 'Asc':Asc,     # for FRP flexural strengthening;
                'dfrp':dfrp, 'dfrpt':dfrpt, 'barType':barType, 'dsv':dsv,
                'wfrpv':wfrpv, 'sfrpv':sfrpv, 'frpIncline':frpIncline}    # for FRP shear strengthening  
    ME_dict.update(new_ME_dict)
    
    return ME_dict, material_dict, geo_dict