# DBN with LWS v0_1 (with definitive condition state)
import os
import numpy as np
import scipy.stats as stats
import scipy.special as spec

from constants import *
from evidence import *
from corrosion import *

from resistance import *
from functions import *
from pyre.distributions import *

import time
import datetime

def corrosionMC():
    #analysis_type = raw_input('analysis type: flexure, shear or deck?')
    analysis_type = rebar_type
    if 'f' in analysis_type.lower():
        datafile_path = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'flexure')
    elif 's' in analysis_type.lower():
        datafile_path = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'shear')
    elif 'd' in analysis_type.lower():
        datafile_path = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'deck')
    else:
        print '[ERROR:] illegal analysis type, must be flexure, shear or deck'
        sys.exit(1)

    # Structural age
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)

    weight_sum = 0
    n_iter = 1
    seed_indx = np.arange(1,45,2)

    chloride_sums = np.zeros(1)
    corrosion_state1_sums = np.zeros(1)
    corrosion_state2_sums = np.zeros(1)
    corrosion_rate_sums = np.zeros(1)
    residual_diameter1_sums = np.zeros(1)
    residual_diameter2_sums = np.zeros(1)
    radial_pressure_sums = np.zeros(1)
    crack_prob_sums = np.zeros(1)
    ds_crack_sums = np.zeros(1)
    crack_width_sums = np.zeros(1)
    diffusion_crack_sums = np.zeros(1)
    rc_flexure_sums = np.zeros(1)
    rc_shear_sums = np.zeros(1)
    #rc_deck_sums = np.zeros(1)

    #chloride_data = np.array([]).reshape(service_time.size,0)
    #corrosion_state_data = np.array([]).reshape(service_time.size,0)
    #corrosion_rate_data = np.array([]).reshape(service_time.size,0)
    #mean_corrosion_rate_data = np.array([]).reshape(service_time.size,0)
    #residual_diameter_data = np.array([]).reshape(service_time.size,0)
    #radial_pressure_data = np.array([]).reshape(service_time.size,0)
    #crack_initiation_data = np.array([]).reshape(service_time.size,0)
    #crack_width_data = np.array([]).reshape(service_time.size,0)
    #diffusion_crack_data = np.array([]).reshape(service_time.size,0)
    #rc_flexure_data = np.array([]).reshape(service_time.size, 0)
    #rc_shear_data = np.array([]).reshape(service_time.size, 0)
    #likelihood_weighting_data = np.array([])

    while weight_sum <= SUM_WEIGHT and n_iter<=MAX_ITER_LW:
        ## seeds for sample generation
        DIFFUSION_REF_SEED = seed_indx[0]**2
        CONCRETE_COVER_SEED = seed_indx[1]**2
        SURFACE_CL_SEED = seed_indx[2]**2
        CRITICAL_CL_SEED = seed_indx[3]**2
        CONCRETE_RESISTANCE_SEED = seed_indx[4]**2
        CORROSION_RATE_SEED = seed_indx[5]**2
        CONCRETE_STRENGTH_SEED = seed_indx[6]**2
        CONCRETE_MODULUS_SEED = seed_indx[7]**2
        CONCRETE_TENSION_SEED = seed_indx[8]**2
        POROUS_LAYER_SEED = seed_indx[9]**2
        CRACK_WIDTH_SEED = seed_indx[10]**2
        CRACK_DIFFUSION_SEED = seed_indx[11]**2
        # flexural resistance
        ME_FLEX_RC_SEED = seed_indx[12]**2
        FSY_SEED = seed_indx[13]**2
        BEAM_WIDTH_SEED = seed_indx[14]**2
        BEAM_DEPTH_SEED = seed_indx[15]**2
        FLANGE_WIDTH_SEED = seed_indx[16]**2
        FLANGE_DEPTH_SEED = seed_indx[17]**2
        # shear resistance
        ME_SHEAR_RC_SEED = seed_indx[18]**2
        SHEAR_DEPTH_SEED = seed_indx[19]**2
        FSYV_SEED = seed_indx[20]**2
        SHEAR_INTERVAL_SEED = seed_indx[21]**2

        seed_indx = seed_indx + 45

        # history of interest
        chloride_history = np.zeros((service_time.size, N_SMP))
        corrosion_state1_history = np.zeros((service_time.size, N_SMP))
        corrosion_state2_history = np.zeros((service_time.size, N_SMP))
        corrosion_rate_history = np.zeros((service_time.size, N_SMP))
        residual_diameter1_history = np.zeros((service_time.size, N_SMP))
        residual_diameter2_history = np.zeros((service_time.size, N_SMP))
        radial_pressure_history = np.zeros((service_time.size, N_SMP))
        crack_initiation_history = np.zeros((service_time.size, N_SMP))
        ds_crack_history = np.zeros((service_time.size, N_SMP))
        crack_width_history = np.zeros((service_time.size, N_SMP))
        diffusion_crack_history = np.zeros((service_time.size, N_SMP))
        rc_flexure_history = np.zeros((service_time.size, N_SMP))
        rc_shear_history = np.zeros((service_time.size, N_SMP))
        #rc_deck_history = np.zeros((service_time.size, N_SMP))
        # initial likelihood weighting
        likelihood_weighting = np.ones(N_SMP)

        ## initial samples
        # reference diffusion coefficient
        ref_diffusion_coefficient = refDiffusionCoefVariable()
        np.random.seed(DIFFUSION_REF_SEED)
        Drcm_smp = ref_diffusion_coefficient.rv.rvs(size=N_SMP)    # [mm^2/year]

        # concrete cover
        concrete_cover = concCoverVariable()
        np.random.seed(CONCRETE_COVER_SEED)
        dc_smp = concrete_cover.rv.rvs(size = N_SMP)    # [mm]

        # surface chloride
        surface_chloride = surfaceClVariable()
        np.random.seed(SURFACE_CL_SEED)
        Cs_smp = surface_chloride.rv.rvs(size = N_SMP)     # [kg/m3]

        # critical chloride
        critical_chloride = criticalClVariable()
        np.random.seed(CRITICAL_CL_SEED)
        Ccr_smp = critical_chloride.rv.rvs(size = N_SMP)

        # concrete strength
        compressive_strength = concStrengthVariable()
        np.random.seed(CONCRETE_STRENGTH_SEED)
        fc_smp = compressive_strength.rv.rvs(size = N_SMP)

        # effective concrete elastic modulus
        elastic_modulus = concEffEcVariable()
        np.random.seed(CONCRETE_MODULUS_SEED)
        Ec_smp = elastic_modulus.rv.rvs(size = N_SMP)

        # concrete tensile strength
        tensile_strength = concTensileVariable()
        np.random.seed(CONCRETE_TENSION_SEED)
        ft_smp = tensile_strength.rv.rvs(size = N_SMP)

        # critical radial pressure
        pcr_smp = criticalPressure(dc_smp, ft_smp)

        # porous layer thickness: delta0
        porous_zone = porousLayerThickVariable()
        np.random.seed(POROUS_LAYER_SEED)
        delta0_smp = porous_zone.rv.rvs(size = N_SMP)

        # model error of Rc
        log_rc_var = modelErrorRcVariable()
        np.random.seed(CONCRETE_RESISTANCE_SEED)
        log_rc_var_smp = log_rc_var.rv.rvs(size=N_SMP)

        # model error of icorr
        log_icorr_var = modelErrorIcorrVariable()
        np.random.seed(CORROSION_RATE_SEED)
        log_icorr_var_smp = log_icorr_var.rv.rvs(size=N_SMP) - np.log(1.08)

        # model error of CRACK_K
        crk_var = modelErrorCrackVariable()
        np.random.seed(CRACK_WIDTH_SEED)
        crk_var_smp = crk_var.rv.rvs(size=N_SMP)

        # model error of fw
        Dck_var = modelErrorDiffusionVariable()
        np.random.seed(CRACK_DIFFUSION_SEED)
        Dck_var_smp = Dck_var.rv.rvs(size=N_SMP)

        # model error of flexural strength of RC beams
        ME_flex = modelErrorRCFlexVariable()
        np.random.seed(ME_FLEX_RC_SEED)
        ME_flex_smp = ME_flex.rv.rvs(size=N_SMP)

        ## model error of deck strength of RC beams
        #ME_deck = modelErrorRCFlexVariable()
        #ME_deck_smp = ME_flex.rv.rvs(size=N_SMP)

        # yielding strength of steel
        fy_steel = steelYieldingVariable()
        np.random.seed(FSY_SEED)
        fy_smp = fy_steel.rv.rvs(size=N_SMP)

        # beam width and depth
        beam_width = beamWidthVariable()
        np.random.seed(BEAM_WIDTH_SEED)
        b_smp = beam_width.rv.rvs(size=N_SMP)
        beam_depth = beamDepthVariable()
        np.random.seed(BEAM_DEPTH_SEED)
        d_smp = beam_depth.rv.rvs(size=N_SMP)

        # flange width and depth
        flange_width = flangeWidthVariable()
        np.random.seed(FLANGE_WIDTH_SEED)
        bf_smp = flange_width.rv.rvs(size=N_SMP)
        flange_depth = flangeDepthVariable()
        np.random.seed(FLANGE_DEPTH_SEED)
        hf_smp = flange_depth.rv.rvs(size=N_SMP)

        # model error of shear strength of RC beams
        ME_shear = modelErrorRCShearVariable()
        np.random.seed(ME_SHEAR_RC_SEED)
        ME_shear_smp = ME_shear.rv.rvs(size=N_SMP)

        # yielding strength of shear reinforcement
        fyv_steel = shearYieldingVariable()
        np.random.seed(FSYV_SEED)
        fyv_smp = fyv_steel.rv.rvs(size=N_SMP)

        # shear depth
        dv = shearDepthVariable()
        np.random.seed(SHEAR_DEPTH_SEED)
        dv_smp = dv.rv.rvs(size=N_SMP)

        # shear interval
        sv = shearIntervalVariable()
        np.random.seed(SHEAR_INTERVAL_SEED)
        sv_smp = sv.rv.rvs(size=N_SMP)


        ## temperary nodes
        # temperature
        temp = tempVariable()
        # volume ratio: gamma
        volume_ratio = volRatioVariable()


        ## start the process
        Cl_prev = np.zeros(N_SMP)
        Cl_prev2 = np.zeros(N_SMP)
        isIni_prev = np.zeros(N_SMP).astype(bool)
        isIni_prev2 = np.zeros(N_SMP).astype(bool)
        iniYr_smp = np.ones(N_SMP) * (END_AGE*2)
        iniYr_smp2 = np.ones(N_SMP) * (END_AGE*2)
        icorr_prev = np.zeros(N_SMP)
        icorr_prev2 = np.zeros(N_SMP)
        ds0_smp = np.ones(N_SMP) * ds_region1_mean
        ds_smp = np.ones(N_SMP) * ds_region1_mean
        ds_smp2 = np.ones(N_SMP) * ds_region2_mean
        ds_prev = np.ones(N_SMP) * ds_region1_mean
        ds_prev2 = np.ones(N_SMP) * ds_region2_mean
        pcor_prev = np.ones(N_SMP)
        isCrack_prev = np.zeros(N_SMP).astype(bool)
        ds_crack_smp = np.zeros(N_SMP)
        wc_prev = np.zeros(N_SMP)
        fw_prev = np.ones(N_SMP)
        gamma_pre = np.zeros(N_SMP)
        # temperature smps
        #np.random.seed(int(age)*n_iter)
        temp_smp = temp.rv.rvs(size=N_SMP)
        #np.random.seed(int(age)*n_iter)
        gamma_smp = volume_ratio.rv.rvs(size = N_SMP)
        #gamma_smp = np.maximum(gamma_pre, gamma_smp)
        #gamma_pre = gamma_smp

        for age in service_time:
            ## temperature smps
            #np.random.seed(int(age)*n_iter)
            #temp_smp = temp.rv.rvs(size=N_SMP)

            # get Cl smps
            D_smp = np.copy(Drcm_smp)
            # region 1
            D_smp[isCrack_prev] = fw_prev[isCrack_prev] * Drcm_smp[isCrack_prev]
            Cl_smp = chlorideContent(age, temp_smp, D_smp, dc_smp, Cs_smp)
            Cl_smp = np.maximum(Cl_smp, Cl_prev)
            chloride_history[service_time==age, :] = Cl_smp
            Cl_prev = Cl_smp
            # region 2
            D_smp2 = np.copy(Drcm_smp)
            Cl_smp2 = chlorideContent(age, temp_smp, D_smp2, dc_smp+DISTANCE_12, Cs_smp)
            Cl_smp2 = np.maximum(Cl_smp2, Cl_prev2)
            Cl_prev2 = Cl_smp2

            # get corrosion state
            # region 1
            isIni_smp = Cl_smp>=Ccr_smp
            isIni_smp = np.logical_or(isIni_smp, isIni_prev)
            corrosion_state1_history[service_time==age, :] = isIni_smp
            # check evidence
            if not np.isnan(evidence_dict['iniStat'][service_time==age]):    # with evidence
                isIni_smp_priori = np.copy(isIni_smp)
                if evidence_dict['iniStat'][service_time==age].astype(bool):
                    isIni_smp = np.ones(N_SMP, dtype=bool)
                    condition_prob = isIni_smp_priori
                else:
                    isIni_smp = np.zeros(N_SMP, dtype=bool)
                    condition_prob = np.logical_not(isIni_smp_priori)
                likelihood_weighting = likelihood_weighting * condition_prob
            elif not np.isnan(evidence_dict['halfCell'][service_time==age]):
                half_cell_evidence = evidence_dict['halfCell'][service_time==age]
                condition_prob = halfcellLikelihood(isIni_smp, half_cell_evidence)
                likelihood_weighting = likelihood_weighting * condition_prob
            # region 2
            isIni_smp2 = Cl_smp2>=Ccr_smp
            isIni_smp2 = np.logical_or(isIni_smp2, isIni_prev2)
            corrosion_state2_history[service_time==age, :] = isIni_smp2

            # calculate tcorr
            # region 1
            iniYr_smp[np.logical_and(np.logical_not(isIni_prev), isIni_smp)] = age
            tcorr_smp = age - iniYr_smp + TIME_INTERVAL
            tcorr_smp[tcorr_smp<0] = 0.0
            isIni_prev = isIni_smp
            # region 2
            iniYr_smp2[np.logical_and(np.logical_not(isIni_prev2), isIni_smp2)] = age
            tcorr_smp2 = age - iniYr_smp2 + TIME_INTERVAL
            tcorr_smp2[tcorr_smp2<0] = 0.0
            isIni_prev2 = isIni_smp2

            # get Rc smp
            # region 1
            log_smp = logRcSmp(Cl_smp, log_rc_var_smp)
            Rc_smp = np.exp(log_smp)
            # region 2
            log_smp2 = logRcSmp(Cl_smp2, log_rc_var_smp)
            Rc_smp2 = np.exp(log_smp2)

            # get corrosion current density (corrosion rate) smps
            # region 1
            log_smp = logIcorrSmp(Cl_smp, temp_smp, Rc_smp, tcorr_smp, log_icorr_var_smp)
            icorr_smp = np.exp(log_smp)
            icorr_smp[tcorr_smp==0] = 0
            corrosion_rate_history[service_time==age, :] = icorr_smp
            if not np.isnan(evidence_dict['icorr'][service_time==age]):
                icorr_evidence = evidence_dict['icorr'][service_time==age]
                condition_prob = icorrLikelihood(icorr_smp, icorr_evidence)
                likelihood_weighting = likelihood_weighting * condition_prob
            #region 2
            log_smp2 = logIcorrSmp(Cl_smp2, temp_smp, Rc_smp2, tcorr_smp2, log_icorr_var_smp)
            icorr_smp2 = np.exp(log_smp2)
            icorr_smp2[tcorr_smp2==0] = 0

            # compute imean
            # region 1
            #imean_smp = 0.5*(icorr_prev + icorr_smp)
            imean_smp = icorr_smp
            icorr_prev = icorr_smp
            # region 2
            imean_smp2 = icorr_smp2

            # get mass loss and residual diameter
            # region 1
            mass_loss_smp, section_loss_smp, ds_smp = corrosionLossSmp(ds_smp, imean_smp)
            ds_smp = np.minimum(ds_prev, ds_smp)
            residual_diameter1_history[service_time==age, :] = ds_smp
            ds_prev = ds_smp
            # region 2
            mass_loss_smp2, section_loss_smp2, ds_smp2 = corrosionLossSmp(ds_smp2, imean_smp2)
            ds_smp2 = np.minimum(ds_prev2, ds_smp2)
            residual_diameter2_history[service_time==age, :] = ds_smp2
            ds_prev2 = ds_smp2

            # get radial pressure
            pcor_smp = radialPressure(mass_loss_smp, Ec_smp, gamma_smp, delta0_smp, dc_smp)
            pcor_smp[pcor_smp<pcor_prev] = pcor_prev[pcor_smp<pcor_prev]
            radial_pressure_history[service_time==age, :] = pcor_smp
            pcor_prev = pcor_smp

            # get crack state and ds at crack
            isCrack_smp = pcor_smp>pcr_smp
            isCrack_smp = np.logical_or(isCrack_smp, isCrack_prev)
            crack_initiation_history[service_time==age, :] = isCrack_smp
            # check evidence
            if not np.isnan(evidence_dict['crkStat'][service_time==age]):    # with evidence
                isCrack_smp_priori = np.copy(isCrack_smp)
                if evidence_dict['crkStat'][service_time==age].astype(bool):
                    isCrack_smp = np.ones(N_SMP, dtype=bool)
                    condition_prob = isCrack_smp_priori
                else:
                    isCrack_smp = np.zeros(N_SMP, dtype=bool)
                    condition_prob = np.logical_not(isCrack_smp_priori)
                likelihood_weighting = likelihood_weighting * condition_prob
            ds_crack_smp[np.logical_and(np.logical_not(isCrack_prev), isCrack_smp)] = ds_smp[np.logical_and(np.logical_not(isCrack_prev), isCrack_smp)] 
            ds_crack_history[service_time==age, :] = ds_crack_smp
            isCrack_prev = isCrack_smp

            # get crack width
            wc_smp = crackWidthSmp(ds_crack_smp, ds_smp, crk_var_smp)
            wc_smp = np.maximum(wc_smp, wc_prev)
            crack_width_history[service_time == age, :] = wc_smp
            wc_prev = wc_smp
            ## check evidence
            if not np.isnan(evidence_dict['conditionState'][service_time==age]):    # with evidence
                crk_std = np.maximum(CRACK_MEASUREMENT_ABSERROR, wc_smp*CRACK_MEASUREMENT_RELERROR)
                #deltaAsloss = np.pi/4*ds_crack_smp**2-np.pi/4*ds_smp**2
                #deltaAsloss = np.maximum(deltaAsloss, 0)
                #crk_std = wc_smp * 0.4*np.exp(0.5*CRACK_GAMMA1+0.5*CRACK_GAMMA2*(deltaAsloss))
                #crk_std = np.maximum(CRACK_MEASUREMENT_ABSERROR, crk_std)
                if evidence_dict['conditionState'][service_time==age] == 1:
                    wc_smp_to_cdf_indx = (CS1_UB - wc_smp) / crk_std
                    condition_prob = stats.norm.cdf(wc_smp_to_cdf_indx)
                elif evidence_dict['conditionState'][service_time==age] == 2:
                    wc_smp_to_cdf_indx1 = (CS2_LB - wc_smp) / crk_std
                    wc_smp_to_cdf_indx2 = (CS2_UB - wc_smp) / crk_std
                    condition_prob = stats.norm.cdf(wc_smp_to_cdf_indx2) - stats.norm.cdf(wc_smp_to_cdf_indx1)
                elif evidence_dict['conditionState'][service_time==age] == 3:
                    wc_smp_to_cdf_indx1 = (CS3_LB - wc_smp) / crk_std
                    wc_smp_to_cdf_indx2 = (CS3_UB - wc_smp) / crk_std
                    condition_prob = stats.norm.cdf(wc_smp_to_cdf_indx2) - stats.norm.cdf(wc_smp_to_cdf_indx1)
                else:    # cnndition state 4
                    wc_smp_to_cdf_indx = (CS4_LB - wc_smp) / crk_std
                    condition_prob = 1 - stats.norm.cdf(wc_smp_to_cdf_indx)
                likelihood_weighting = likelihood_weighting * condition_prob

            # get diffusion coefficient of cracked concrete
            fw_smp = diffusionRatioSmp(Dck_var_smp, wc_smp)
            fw_smp[np.logical_not(isCrack_smp)] = 1.0
            fw_smp = np.maximum(fw_prev, fw_smp)
            diffusion_crack_history[service_time == age, :] = fw_smp
            fw_prev = fw_smp

            # get resistance
            # name, ME_flex, fc, fy, Ast, b, d, bf, hf
            # flexural
            Ast_smp = residualSteelArea(section_loss_smp, section_loss_smp2)
            ME_dict, material_dict, geo_dict = assembleBeamDict(ME_flex_smp, ME_shear_smp,\
                                                                fc_smp, LAMBDA_FC, fy_smp, fyv_smp,\
                                                                Ast_smp, Ast_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
            rcBeam = RCBeam('rc_beam', ME_dict, material_dict, geo_dict)
            mu_smp = rcBeam.flexCapacity()
            rc_flexure_history[service_time == age, :] = mu_smp
            # shear
            Asvt_smp = residualSteelArea(section_loss_smp, section_loss_smp2)
            ME_dict, material_dict, geo_dict = assembleBeamDict(ME_flex_smp, ME_shear_smp,\
                                                                fc_smp, LAMBDA_FC, fy_smp, fyv_smp,\
                                                                Asvt_smp, Asvt_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
            rcBeam = RCBeam('rc_beam', ME_dict, material_dict, geo_dict)
            vu_smp = rcBeam.shearCapacity()
            rc_shear_history[service_time == age, :] = vu_smp
            ## deck
            #Ast_smp = residualSteelArea(section_loss_smp, section_loss_smp2)
            #ME_dict, material_dict, geo_dict = assembleBeamDict(ME_deck_smp, ME_shear_smp,\
            #                                                    fc_smp, LAMBDA_FC, fy_smp, fyv_smp,\
            #                                                    Ast_smp, Ast_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
            #rcBeam = RCBeam('rc_beam', ME_dict, material_dict, geo_dict)
            #mu_smp = rcBeam.flexCapacity()
            #rc_deck_history[service_time == age, :] = mu_smp

        # calculate sums
        sums = weightedSum(chloride_history, likelihood_weighting, axis_data=-1)
        chloride_sums = chloride_sums + sums

        sums = weightedSum(corrosion_state1_history, likelihood_weighting, axis_data=-1)
        corrosion_state1_sums = corrosion_state1_sums + sums

        sums = weightedSum(corrosion_state2_history, likelihood_weighting, axis_data=-1)
        corrosion_state2_sums = corrosion_state2_sums + sums

        sums = weightedSum(corrosion_rate_history, likelihood_weighting, axis_data=-1)
        corrosion_rate_sums = corrosion_rate_sums + sums

        sums = weightedSum(residual_diameter1_history, likelihood_weighting, axis_data=-1)
        residual_diameter1_sums = residual_diameter1_sums + sums

        sums = weightedSum(residual_diameter2_history, likelihood_weighting, axis_data=-1)
        residual_diameter2_sums = residual_diameter2_sums + sums

        sums = weightedSum(radial_pressure_history, likelihood_weighting, axis_data=-1)
        radial_pressure_sums = radial_pressure_sums + sums

        sums = weightedSum(crack_initiation_history, likelihood_weighting, axis_data=-1)
        crack_prob_sums = crack_prob_sums + sums

        sums = weightedSum(ds_crack_history, likelihood_weighting, axis_data=-1)
        ds_crack_sums = ds_crack_sums + sums

        sums = weightedSum(crack_width_history, likelihood_weighting, axis_data=-1)
        crack_width_sums = crack_width_sums + sums

        sums = weightedSum(diffusion_crack_history, likelihood_weighting, axis_data=-1)
        diffusion_crack_sums = diffusion_crack_sums + sums

        sums = weightedSum(rc_flexure_history, likelihood_weighting, axis_data=-1)
        rc_flexure_sums = rc_flexure_sums + sums

        sums = weightedSum(rc_shear_history, likelihood_weighting, axis_data=-1)
        rc_shear_sums = rc_shear_sums + sums

        #sums = weightedSum(rc_deck_history, likelihood_weighting, axis_data=-1)
        #rc_deck_sums = rc_deck_sums + sums

        ## save data
        ##save_percentile = np.sum(likelihood_weighting) / SUM_WEIGHT * 100.0
        ##accept_weight = np.percentile(likelihood_weighting, 100 - save_percentile)
        #chloride_data = np.hstack((chloride_data, chloride_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #corrosion_state_data = np.hstack((corrosion_state_data, corrosion_state_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #corrosion_rate_data = np.hstack((corrosion_rate_data, corrosion_rate_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #mean_corrosion_rate_data = np.hstack((mean_corrosion_rate_data, mean_corrosion_rate_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #residual_diameter_data = np.hstack((residual_diameter_data, residual_diameter_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))        
        #radial_pressure_data = np.hstack((radial_pressure_data, radial_pressure_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #crack_initiation_data = np.hstack((crack_initiation_data, crack_initiation_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #crack_width_data = np.hstack((crack_width_data, crack_width_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #diffusion_crack_data = np.hstack((diffusion_crack_data, diffusion_crack_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #rc_flexure_data = np.hstack((rc_flexure_data, rc_flexure_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #rc_shear_data = np.hstack((rc_shear_data, rc_shear_history[:, likelihood_weighting>=ACCEPT_WEIGHT]))
        #likelihood_weighting_data = np.hstack((likelihood_weighting_data, likelihood_weighting[likelihood_weighting>=ACCEPT_WEIGHT]))

        weight_sum = weight_sum + np.sum(likelihood_weighting)
        n_iter += 1

    if n_iter > MAX_ITER:
        print 'Warning: accumulated weight sum is smaller than ' + str(SUM_WEIGHT) 

    ## save data to binary files
    #np.save(os.path.join(DATAFILE_PATH,'chloride_history.npy'), chloride_data)
    #np.save(os.path.join(DATAFILE_PATH,'corrosion_state_history.npy'), corrosion_state_data)
    #np.save(os.path.join(DATAFILE_PATH,'corrosion_rate_history.npy'), corrosion_rate_data)
    #np.save(os.path.join(DATAFILE_PATH,'mean_corrosion_rate_history.npy'), mean_corrosion_rate_data)
    #np.save(os.path.join(DATAFILE_PATH,'residual_diameter_history.npy'), residual_diameter_data)
    #np.save(os.path.join(DATAFILE_PATH,'radial_pressure_history.npy'), radial_pressure_data)
    #np.save(os.path.join(DATAFILE_PATH,'crack_initiation_history.npy'), crack_initiation_data)
    #np.save(os.path.join(DATAFILE_PATH,'crack_width_history.npy'), crack_width_data)
    #np.save(os.path.join(DATAFILE_PATH,'diffusion_crack_history.npy'), diffusion_crack_data)
    #np.save(os.path.join(DATAFILE_PATH,'rc_flexure_history.npy'), rc_flexure_data)
    #np.save(os.path.join(DATAFILE_PATH,'rc_shear_history.npy'), rc_shear_data)
    #np.save(os.path.join(DATAFILE_PATH,'likelihood_weighting.npy'), likelihood_weighting_data)

    # calculate mean and stds
    chloride_mean_history, chloride_std_history = weightedAvgAndStdFromSum(chloride_sums)
    corrosion_prob1_history, dummy = weightedAvgAndStdFromSum(corrosion_state1_sums)
    corrosion_prob2_history, dummy = weightedAvgAndStdFromSum(corrosion_state2_sums)
    corrosion_rate_mean_history, corrosion_rate_std_history = weightedAvgAndStdFromSum(corrosion_rate_sums)
    residual_diameter1_mean_history, residual_diameter1_std_history = weightedAvgAndStdFromSum(residual_diameter1_sums)
    residual_diameter2_mean_history, residual_diameter2_std_history = weightedAvgAndStdFromSum(residual_diameter2_sums)
    radial_pressure_mean_history, radial_pressure_std_history = weightedAvgAndStdFromSum(radial_pressure_sums)
    crack_prob_history, dummy = weightedAvgAndStdFromSum(crack_prob_sums)
    ds_crack_mean_history, ds_crack_std_history = weightedAvgAndStdFromSum(ds_crack_sums)
    crack_width_mean_history, crack_width_std_history = weightedAvgAndStdFromSum(crack_width_sums)
    diffusion_crack_mean_history, diffusion_crack_std_history = weightedAvgAndStdFromSum(diffusion_crack_sums)
    rc_flexure_mean_history, rc_flexure_std_history = weightedAvgAndStdFromSum(rc_flexure_sums)
    rc_shear_mean_history, rc_shear_std_history = weightedAvgAndStdFromSum(rc_shear_sums)
    #rc_deck_mean_history, rc_deck_std_history = weightedAvgAndStdFromSum(rc_deck_sums)

    # save data to text file
    datafile = os.path.join(datafile_path, 'LWS_results.txt')
    with open(datafile, 'w') as f_handle:
        np.savetxt(f_handle, np.array(['# service life']), fmt="%s")
        np.savetxt(f_handle, service_time.reshape(1, service_time.size), fmt='%d')
        np.savetxt(f_handle, np.array(['# chloride history (mean)']), fmt='%s')
        np.savetxt(f_handle, chloride_mean_history.reshape(1, chloride_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# chloride history (std)']), fmt='%s')
        np.savetxt(f_handle, chloride_std_history.reshape(1, chloride_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# corrosion prob 1 history']), fmt='%s')
        np.savetxt(f_handle, corrosion_prob1_history.reshape(1, corrosion_prob1_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# corrosion prob 2 history']), fmt='%s')
        np.savetxt(f_handle, corrosion_prob2_history.reshape(1, corrosion_prob2_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# corrosion rate history (mean)']), fmt='%s')
        np.savetxt(f_handle, corrosion_rate_mean_history.reshape(1, corrosion_rate_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# corrosion rate history (std)']), fmt='%s')
        np.savetxt(f_handle, corrosion_rate_std_history.reshape(1, corrosion_rate_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# residual diameter 1 history (mean)']), fmt='%s')
        np.savetxt(f_handle, residual_diameter1_mean_history.reshape(1, residual_diameter1_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# residual diameter 1 history (std)']), fmt='%s')
        np.savetxt(f_handle, residual_diameter1_std_history.reshape(1, residual_diameter1_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# residual diameter 2 history (mean)']), fmt='%s')
        np.savetxt(f_handle, residual_diameter2_mean_history.reshape(1, residual_diameter2_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# residual diameter 2 history (std)']), fmt='%s')
        np.savetxt(f_handle, residual_diameter2_std_history.reshape(1, residual_diameter2_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# radial pressure history (mean)']), fmt='%s')
        np.savetxt(f_handle, radial_pressure_mean_history.reshape(1, radial_pressure_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# radial pressure history (std)']), fmt='%s')
        np.savetxt(f_handle, radial_pressure_std_history.reshape(1, radial_pressure_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# crack prob history']), fmt='%s')
        np.savetxt(f_handle, crack_prob_history.reshape(1, crack_prob_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# ds at crack history (mean)']), fmt='%s')
        np.savetxt(f_handle, ds_crack_mean_history.reshape(1, ds_crack_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# ds at crack history (std)']), fmt='%s')
        np.savetxt(f_handle, ds_crack_std_history.reshape(1, ds_crack_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# crack width history (mean)']), fmt='%s')
        np.savetxt(f_handle, crack_width_mean_history.reshape(1, crack_width_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# crack width history (std)']), fmt='%s')
        np.savetxt(f_handle, crack_width_std_history.reshape(1, crack_width_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# cracked concrete diffusion history (fw, mean)']), fmt='%s')
        np.savetxt(f_handle, diffusion_crack_mean_history.reshape(1, diffusion_crack_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# cracked concrete diffusion history (fw, std)']), fmt='%s')
        np.savetxt(f_handle, diffusion_crack_std_history.reshape(1, diffusion_crack_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# flexural resistance history (mean, kNm)']), fmt='%s')
        np.savetxt(f_handle, rc_flexure_mean_history.reshape(1, rc_flexure_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# flexural resistance history (std, kNm)']), fmt='%s')
        np.savetxt(f_handle, rc_flexure_std_history.reshape(1, rc_flexure_std_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# shear resistance history (mean, kN)']), fmt='%s')
        np.savetxt(f_handle, rc_shear_mean_history.reshape(1, rc_shear_mean_history.size), fmt='%.4e')
        np.savetxt(f_handle, np.array(['# shear resistance history (std, kN)']), fmt='%s')
        np.savetxt(f_handle, rc_shear_std_history.reshape(1, rc_shear_std_history.size), fmt='%.4e')
        #np.savetxt(f_handle, np.array(['# deck resistance history (mean, kNm)']), fmt='%s')
        #np.savetxt(f_handle, rc_deck_mean_history.reshape(1, rc_deck_mean_history.size), fmt='%.4e')
        #np.savetxt(f_handle, np.array(['# deck resistance history (std, kNm)']), fmt='%s')
        #np.savetxt(f_handle, rcstd_history.reshape(1, rcstd_history.size), fmt='%.4e')


if __name__ == '__main__':
    print 'CALC: begin'
    start_delta_time = time.time()
    corrosionMC()
    delta_time = time.time() - start_delta_time
    print 'DONE: ',str(datetime.timedelta(seconds=delta_time))
