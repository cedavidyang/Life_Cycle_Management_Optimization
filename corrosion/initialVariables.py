import os
import sys
import numpy as np
from constants.beamConstants import *
from constants.corrosionConstants import START_AGE, END_AGE, TIME_INTERVAL
from constants.evidenceConstants import *

def importUserRebar():
    """import user defined variables and refine constants for corrosion analysis
    """
    # define rebar under consideration: flexure or shear
    rebar_type = raw_input('define rebar under consideration: flexure, shear or deck?')
    # corrosion related variables
    if 'f' in rebar_type.lower():
        cover_nom = COVER_NOM;            cover_mean = COVER_MEAN;            cover_cov = COVER_COV;            cover_distr = COVER_DISTR;
        ds_region1_nom = DS_REGION1_NOM;  ds_region1_mean = DS_REGION1_MEAN;  ds_region1_cov = DS_REGION1_COV;  ds_region1_distr = DS_REGION1_DISTR;
        ds_region2_nom = DS_REGION2_NOM;  ds_region2_mean = DS_REGION2_MEAN;  ds_region2_cov = DS_REGION2_COV;  ds_region2_distr = DS_REGION2_DISTR;
        region1_no = REGION1_NO
        region2_no = REGION2_NO
        distance_12 = DISTANCE_12
    elif 's' in rebar_type.lower():
        cover_nom = COVER_SHEAR_NOM;            cover_mean = COVER_SHEAR_MEAN;            cover_cov = COVER_SHEAR_COV;            cover_distr = COVER_SHEAR_DISTR;
        ds_region1_nom = DS_REGION1_SHEAR_NOM;  ds_region1_mean = DS_REGION1_SHEAR_MEAN;  ds_region1_cov = DS_REGION1_SHEAR_COV;  ds_region1_distr = DS_REGION1_SHEAR_DISTR;
        ds_region2_nom = DS_REGION2_SHEAR_NOM;  ds_region2_mean = DS_REGION2_SHEAR_MEAN;  ds_region2_cov = DS_REGION2_SHEAR_COV;  ds_region2_distr = DS_REGION2_SHEAR_DISTR;
        region1_no = REGION1_NO_SHEAR
        region2_no = REGION2_NO_SHEAR
        distance_12 = DISTANCE_12_SHEAR
    elif 'd' in rebar_type.lower():
        cover_nom = COVER_DECK_NOM;            cover_mean = COVER_DECK_MEAN;            cover_cov = COVER_DECK_COV;            cover_distr = COVER_DECK_DISTR;
        ds_region1_nom = DS_REGION1_DECK_NOM;  ds_region1_mean = DS_REGION1_DECK_MEAN;  ds_region1_cov = DS_REGION1_DECK_COV;  ds_region1_distr = DS_REGION1_DECK_DISTR;
        ds_region2_nom = DS_REGION2_DECK_NOM;  ds_region2_mean = DS_REGION2_DECK_MEAN;  ds_region2_cov = DS_REGION2_DECK_COV;  ds_region2_distr = DS_REGION2_DECK_DISTR;
        region1_no = REGION1_NO_DECK
        region2_no = REGION2_NO_DECK
        distance_12 = DISTANCE_12_DECK
    else:
        print '[ERROR:] Illegal input for rebar type, must be flexure, shear or deck'
        sys.exit(1)
    ## resistance related variables
    #if 'd' in rebar_type.lower():
    #    beam_height_nom = DECK_HEIGHT_NOM;          beam_height_mean = DECK_HEIGHT_MEAN;          beam_height_cov = DECK_HEIGHT_COV;         beam_height_distr = DECK_HEIGHT_DISTR;
    #    beam_depth_nom = DECK_DEPTH_NOM;          beam_depth_mean = DECK_DEPTH_MEAN;          beam_depth_cov = DECK_DEPTH_COV;        beam_depth_distr = DECK_DEPTH_DISTR;
    #    flange_width_nom = DECK_WIDTH_NOM;        flange_width_mean = DECK_WIDTH_MEAN;        flange_width_cov = DECK_WIDTH_COV;        flange_width_distr = DECK_WIDTH_DISTR;
    #    flange_thick_nom = DECK_DEPTH_NOM;         flange_thick_mean = DECK_DEPTH_MEAN;         flange_thick_cov = DECK_DEPTH_COV;        flange_thick_distr = DECK_DEPTH_DISTR; 
    #    web_width_nom = DECK_WIDTH_NOM;            web_width_mean = DECK_WIDTH_MEAN;            web_width_cov = DECK_WIDTH_COV;           web_width_distr = DECK_WIDTH_DISTR;
    #else:
    #    beam_height_nom = BEAM_HEIGHT_NOM;          beam_height_mean = BEAM_HEIGHT_MEAN;          beam_height_cov = BEAM_HEIGHT_COV;         beam_height_distr = BEAM_HEIGHT_DISTR;
    #    beam_depth_nom = BEAM_DEPTH_NOM;          beam_depth_mean = BEAM_DEPTH_MEAN;          beam_depth_cov = BEAM_DEPTH_COV;        beam_depth_distr = BEAM_DEPTH_DISTR;
    #    flange_width_nom = FLANGE_WIDTH_NOM;        flange_width_mean = FLANGE_WIDTH_MEAN;        flange_width_cov = FLANGE_WIDTH_COV;        flange_width_distr = FLANGE_WIDTH_DISTR;
    #    flange_thick_nom = FLANGE_THICK_NOM;         flange_thick_mean = FLANGE_THICK_MEAN;         flange_thick_cov = FLANGE_THICK_COV;        flange_thick_distr = FLANGE_THICK_DISTR;
    #    web_width_nom = WEB_WIDTH_NOM;            web_width_mean = WEB_WIDTH_MEAN;            web_width_cov = WEB_WIDTH_COV;           web_width_distr = WEB_WIDTH_DISTR;

    return {'rebar_type': rebar_type,
            'cover_nom':cover_nom, 'cover_mean':cover_mean, 'cover_cov':cover_cov, 'cover_distr':cover_distr,
            'ds_region1_nom':ds_region1_nom, 'ds_region1_mean':ds_region1_mean, 'ds_region1_cov':ds_region1_cov, 'ds_region1_distr':ds_region1_distr,
            'ds_region2_nom':ds_region2_nom, 'ds_region2_mean':ds_region2_mean, 'ds_region2_cov':ds_region2_cov, 'ds_region2_distr':ds_region2_distr,
            'region1_no':region1_no, 'region2_no':region2_no, 'distance_12':distance_12}

# import user defined rebar: for flexure, shear or deck
res = importUserRebar()
rebar_type = res['rebar_type']
cover_nom = res['cover_nom']
cover_mean = res['cover_mean']
cover_cov = res['cover_cov']
cover_distr = res['cover_distr']

ds_region1_nom = res['ds_region1_nom']
ds_region1_mean = res['ds_region1_mean']
ds_region1_cov = res['ds_region1_cov']
ds_region1_distr = res['ds_region1_distr']

ds_region2_nom = res['ds_region2_nom']
ds_region2_mean = res['ds_region2_mean']
ds_region2_cov = res['ds_region2_cov']
ds_region2_distr = res['ds_region2_distr']

region1_no = res['region1_no']
region2_no = res['region2_no']
distance_12 = res['distance_12']

def importUserEvidence():
    """import user defined variables and refine constants for evidence constants
       define evidence type: (1) corrosion initiation state
                             (2) half-cell potential
                             (3) corrosion rate icorr
                             (4) crack initiation state
                             (5) condition state
    """
    #from corrosion import rebar_type
    # create evidence dict
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    evidence_list = np.empty(service_time.shape)
    evidence_list[:] = np.nan
    evidence_dict = {'iniStat': np.copy(evidence_list),
                    'halfCell': np.copy(evidence_list),
                    'icorr': np.copy(evidence_list),
                    'crkStat': np.copy(evidence_list),
                    'conditionState': np.copy(evidence_list)}

    print 'define type of evidence:'
    print '(1) corrosion initiation state'
    print '(2) half-cell potential'
    print '(3) corrosion rate icorr'
    print '(4) crack initiation state'
    print '(5) condition state'
    evidence_type = input()

    if evidence_type is 1:
        evidence_dict['iniStat'][service_time==NO_INI_YR] = False
        evidence_dict['iniStat'][service_time==IS_INI_YR] = True
    elif evidence_type is 2:
        evidence_dict['halfCell'][service_time==HC_DATA1['year']] = HC_DATA1['potential']
        evidence_dict['halfCell'][service_time==HC_DATA2['year']] = HC_DATA2['potential']
    elif evidence_type is 3:
        evidence_dict['icorr'][service_time==CORROSION_RATE_DATA1['year']] = CORROSION_RATE_DATA1['current']
        evidence_dict['icorr'][service_time==CORROSION_RATE_DATA2['year']] = CORROSION_RATE_DATA2['current']
    elif evidence_type is 4:
        evidence_dict['crkStat'][service_time==NO_CRK_YR] = False
        evidence_dict['crkStat'][service_time==IS_CRK_YR] = True
    elif evidence_type is 5:
        #rebar_type = raw_input('rebar under consideration: flexure, shear or deck?')
        if 'f' in rebar_type.lower():
            cs1_to_cs2 = CS1_TO_CS2_FLEX
            cs2_to_cs3 = CS2_TO_CS3_FLEX
            cs3_to_cs4 = CS3_TO_CS4_FLEX
            evidence_end_yr = EVIDENCE_END_YR_FLEX
        elif 's' in rebar_type.lower():
            cs1_to_cs2 = CS1_TO_CS2_SHEAR
            cs2_to_cs3 = CS2_TO_CS3_SHEAR
            cs3_to_cs4 = CS3_TO_CS4_SHEAR
            evidence_end_yr = EVIDENCE_END_YR_SHEAR
        elif 'd' in rebar_type.lower():
            cs1_to_cs2 = CS1_TO_CS2_DECK
            cs2_to_cs3 = CS2_TO_CS3_DECK
            cs3_to_cs4 = CS3_TO_CS4_DECK
            evidence_end_yr = EVIDENCE_END_YR_DECK
        else:
             print '[ERROR:] illegal type of rebar, must be flexure, shear or deck'
             sys.exit(1)
        evidence_dict['conditionState'][service_time<=cs1_to_cs2] = 1
        evidence_dict['conditionState'][np.logical_and(service_time>cs1_to_cs2, service_time<=cs2_to_cs3)] = 2
        evidence_dict['conditionState'][np.logical_and(service_time>cs2_to_cs3, service_time<=cs3_to_cs4)] = 3
        evidence_dict['conditionState'][np.logical_and(service_time>cs3_to_cs4, service_time<=evidence_end_yr)] = 4

    return evidence_dict

# import user defined evidence dictionary
evidence_dict = importUserEvidence()

if __name__ == '__main__':
    res = importUserRebar()
