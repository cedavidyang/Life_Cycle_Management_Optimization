# functions used in corrosion models
import numpy as np
import scipy.stats as stats
import scipy.special as spec

#from corrosion import cover_mean, cover_cov, cover_distr
#from corrosion import ds_region1_nom, ds_region1_mean, ds_region1_cov, ds_region1_distr
#from corrosion import ds_region2_nom, ds_region2_mean, ds_region2_cov, ds_region2_distr
#from corrosion import region1_no, region2_no, distance_12
from constants.beamConstants import *
from constants.corrosionConstants import *
from pyre.distributions import *

def water2cementRatio(fc_cylinder):
    """return w2c ratio (Bolomey's formula) Vu and Stewart (2000)
    :Args:
        fc_cylinder: mean cylinder strength [MPa]
    :Returns:
        w2c: water-to-cement ratio [-]
    """
    w2c = 27./(fc_cylinder+13.5)

    return w2c


def chlorideContent(age, temp_smp, Drcm_smp, dc_smp, Cs_smp):
    # get parameters
    ke = np.exp(B_E*(1./TEMP_REF - 1./temp_smp))
    at = (TIME_REF/age)**AGE_EXP
    kt = TRANSFER_VARIABLE
    Dapp_smp = ke * Drcm_smp * kt * at

    # redefine parameters
    z = dc_smp
    t = age
    Dapp = Dapp_smp
    Cs = Cs_smp
    C0 = 0

    # calculation
    chloride_concentration = C0 + (Cs - C0) * ( 1 - spec.erf(z / (2*np.sqrt(Dapp*t)) ) )

    return chloride_concentration


def radialPressure(mass_loss_smp, Ec_smp, gamma_smp, delta0_smp, dc_smp):
    from initialVariables import ds_region1_mean

    Ml = mass_loss_smp
    Eef = Ec_smp     # Ec_smp are smps of effective concrete modulus, creep coefficient has been considered
    gamma = gamma_smp
    delta0 = delta0_smp
    d0 = ds_region1_mean
    nu = CONCRETE_POISSON
    pi = np.pi
    x = dc_smp

    psi = (d0+2*delta0)**2 / (2*x*(x+(d0+2*delta0)))

    pcor_smp = 2*(Ml/(RHO_S*KG_TO_G)/(CORROSION_MM_TO_M)**2)*Eef*(gamma-1.) / (pi*d0*(1+psi+nu)*(d0+2*delta0)) - 2*delta0*Eef / ((1+psi+nu)*(d0+2*delta0))
    pcor_smp = np.maximum(0,pcor_smp)

    return pcor_smp


def refDiffusionCoefVariable():
    """Return referrence diffusion coefficient Drcm (Papadakis et al. 1996)
    :Args:
    """
    fc_cylinder = FC_MEAN
    wc = water2cementRatio(fc_cylinder)
    ac = A_2_C
    mu_d0 = D_H2O * 0.15 * (1 + RHO_C * wc) / (1 + RHO_C * wc + RHO_C/RHO_A * ac) * ( (RHO_C*wc - 0.85) / (1 + RHO_C*wc) )**3
    sigma_d0 = mu_d0 * D_COV
    ref_diffusion_coefficient = Lognormal('Drcm',mu_d0,sigma_d0)

    return ref_diffusion_coefficient


def concCoverVariable():
    from initialVariables import cover_mean, cover_cov, cover_distr
    mu = cover_mean
    sigma = mu * cover_cov
    if cover_distr == 'lognormal':
        concrete_cover = Lognormal('dc', mu, sigma) 
    elif cover_distr == 'normal':
        concrete_cover = Normal('dc', mu, sigma) 

    return concrete_cover


def surfaceClVariable():
    d_coast = DISTANCE_COAST
    if 0<=d_coast<=0.1:
        muC = 2.95
    elif 0.1<d_coast<=2.84:
        muC = 1.15 - 1.80 * np.log10(d_coast)
    elif d_coast>2.84:
        muC = 1.15 - 1.80 * np.log10(2.84)
    else:
        sys.exit("Error: 'Cs' value is not in the list")
    sigmaC = muC * CS_COV
    surface_chloride = Lognormal('Cs',muC,sigmaC)

    return surface_chloride


def criticalClVariable():
    mu = C_CRIT_MEAN
    sigma = mu * C_CRIT_COV
    if C_CRIT_DISTR == 'normal':
        critical_chloride = Normal('Ccr',mu,sigma)
    elif C_CRIT_DISTR == 'uniform':
        lb = C_CRIT_MEAN - np.sqrt(3)*C_CRIT_MEAN*C_CRIT_COV
        ub = C_CRIT_MEAN + np.sqrt(3)*C_CRIT_MEAN*C_CRIT_COV
        critical_chloride = Uniform('Ccr',lb,ub, input_type='LB UB')

    return critical_chloride


def concStrengthVariable():
    mu = FC_MEAN
    sigma = FC_COV * mu
    compressive_strength = Normal('fc',mu,sigma)

    return compressive_strength


def concEffEcVariable():
    mu = 4733.0 * np.sqrt(FC_NOM) / (1+CREEP_COEFFICIENT)
    sigma = mu * EC_COV
    if EC_DISTR == 'normal':
        elastic_modulus = Normal('Ec', mu, sigma)
    elif EC_DISTR == 'lognormal':
        elastic_modulus = Lognormal('Ec', mu, sigma)

    return elastic_modulus


def concTensileVariable():
    mu = 0.6227 * np.sqrt(FC_NOM) * (FC_MEAN/FC_NOM)
    sigma = mu*FT_COV
    if FT_DISTR == 'normal':
        tensile_strength = Normal('ft', mu, sigma)
    elif FT_DISTR == 'lognormal':
        tensile_strength = Lognormal('ft', mu, sigma)

    return tensile_strength


def criticalPressure(dc_smp, ft_smp):
    from initialVariables import ds_region1_mean
    pcr_smp = 2*dc_smp*ft_smp / ds_region1_mean

    return pcr_smp


def porousLayerThickVariable():
    # porous layer thickness: delta0
    porous_zone = Uniform('porous_zone', POROUS_LB, POROUS_UB, input_type='LB UB')

    return porous_zone


def modelErrorRcVariable():
    # model error of Rc
    mu = 0
    sigma = 0.1203
    log_rc_var = Normal('delta_Rc', mu, sigma)

    return log_rc_var


def modelErrorIcorrVariable():
    mu = 0
    sigma = 0.3312
    log_icorr_var = Normal('delta_icorr', mu, sigma)

    return log_icorr_var


def modelErrorCrackVariable():
    crk_var = Lognormal('crk_var', mean=1, stdv=CRACK_COV)

    return crk_var


def modelErrorDiffusionVariable():
    Dck_var = Lognormal('Dck_var', mean=1, stdv=1/np.sqrt(DIFFUSION_K))

    return Dck_var


def tempVariable():
    temperature = Normal('Temperature', TEMP_MEAN, TEMP_STDV)

    return temperature


def volRatioVariable():
    volume_ratio = Beta('gamma_val',GAMMA_VOL_ALPHA, GAMMA_VOL_BETA, a=GAMMA_VOL_LB, b=GAMMA_VOL_UB, input_type = 'alpha beta')

    return volume_ratio


def logRcSmp(Cl_smp, log_rc_var_smp):
    log_rc_smp = 8.03 - 0.549 * np.log(1+1.69*Cl_smp) + log_rc_var_smp

    return log_rc_smp


def logIcorrSmp(Cl_smp, temp_smp, Rc_smp, tcorr_smp, log_icorr_var_smp):
    log_icorr_smp = 7.89 + 0.7771 * np.log(1.69 * Cl_smp) - 3006.0/temp_smp - 0.000116*Rc_smp + 2.24*tcorr_smp**(-0.215) + log_icorr_var_smp

    return log_icorr_smp


def corrosionLossSmp(ds_smp, imean_smp):
    from initialVariables import ds_region1_mean
    dm_loss_smp = CORROSION_CONST1 * np.pi * ds_smp * imean_smp * CORROSION_CONST2 * TIME_INTERVAL * CORROSION_NACM2_TO_AM2 * CORROSION_MM_TO_M    # g/m
    dVolume_loss_smp = dm_loss_smp / (RHO_S * KG_TO_G)  / CORROSION_MM_TO_M**2    # in [mm^2], is also delta sectional loss
    ds_curr_sq = (ds_smp)**2 - 4.0*dVolume_loss_smp/np.pi
    ds_smp[ds_curr_sq<0] = 0.0
    ds_smp[ds_curr_sq>=0] =  np.sqrt( ds_curr_sq[ds_curr_sq>=0] )
    mass_loss_smp = (RHO_S*KG_TO_G) * np.pi / 4 * ( (ds_region1_mean*CORROSION_MM_TO_M)**2 - (ds_smp*CORROSION_MM_TO_M)**2)    #[g/m]
    section_loss_smp = np.pi / 4 * ( ds_region1_mean**2 - ds_smp**2)    # [mm2]

    return mass_loss_smp, section_loss_smp, ds_smp


def crackWidthSmp(ds_crack_smp, ds_smp, crk_var_smp):
    from initialVariables import ds_region1_mean
    As0 = np.pi/4 * (ds_region1_mean**2 - ds_crack_smp**2)
    As = np.pi/4 * (ds_region1_mean**2 - ds_smp**2)
    deltaA = np.maximum(0,As-As0)
    wc_det = CRACK_K * deltaA
    wc_smp = wc_det * crk_var_smp

    return wc_smp


def diffusionRatioSmp(Dck_var_smp, wc_smp):
    fw_smp = (DIFFUSION_CONST1 * wc_smp**2 + DIFFUSION_CONST2 * wc_smp)*Dck_var_smp + 1
    fw_smp[fw_smp>20.0] = 20.0

    return fw_smp


def residualSteelArea(section_loss_smp, section_loss_smp2):
    from initialVariables import ds_region1_mean, ds_region2_mean
    from initialVariables import region1_no, region2_no
    As0_region1_onebar = np.pi/4 * ds_region1_mean**2    #[mm2]
    As0_region2_onebar = np.pi/4 * ds_region2_mean**2
    Ast_smp = region1_no * (As0_region1_onebar - section_loss_smp) + region2_no *(As0_region2_onebar - section_loss_smp2)

    return Ast_smp
