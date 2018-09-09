"""constants used in the life-cycle FRP design"""


"""model information
# BEAM_HEIGHT: beam height [mm], Enright and Frangopol (1999)

# BEAM_DEPTH: effective depth [mm], Enright and Frangopol (1999)
# SHEAR_DEPTH_NOM: effective shear depth [mm], Enright and Frangopol (1999)

# FLANGE_WIDTH: effective flange width [mm], Enright and Frangopol (1999) and AASHTO 4.6.2.6.1
# FLANGE_THICK: effective flange thickness [mm], Enright and Frangopol (1999)
# WEB_WIDTH: web width of T-beam (beam width of rectangular beam) [mm], Enright and Frangopol (1999)
# COVER: cover thickness [mm], mean: Enright and Frangopol (1999), cov and distribution: Enright and Frangopol (1997)
# DS_REGION1: diameter of reinforcing steel at region 1 (outer layer) [mm], Enright and Frangopol (1999)
# DS_REGION2: diameter of reinforcing steel at region 2 (inner layer) [mm], Enright and Frangopol (1999)
# REGION1_NO: number of rebars in region 1, Enright and Frangopol (1999)
# REGION2_NO: number of rebars in region 2, Enright and Frangopol (1999)
# DISTANCE_12: distance between 1 and 2, Enright and Frangopol (1999)
# FC: concrete strength at the 28th day [MPa], Enright and Frangopol (1997), mean value has been altered according to nominal value, COV is also different.
      the original values are based on Ellingwood et al. (1980), but 0.85 in 0.85fc' is used mistakenly twice in Enright and Frangopol (1997)
# FCT: concrete tensile strength [MPa], Papakonstantinou and Shinozuka (2013)
# ES: steel modulus [MPa], assupmtion
# FSY: yielding strength of steel [MPa], Enright and Frangopol (1997) and Ellingwood et al. (1980)

# FSYV_NOM, yielding strneght of shear reinforement [MPa]

# ME_FLEX_RC: model error of flexural strength of conventional RC beams, ACI model [-], Ellingwood et al. (1980), distribution is assumed
# ME_SHEAR_RC: model error of shear strength of conventional RC beams, ACI model [-], Ellingwood et al. (1980), distribution is assumed

# M_DCDW: moment due to DC and DW [kNm], Enright and Frangopol (1999), evenly divided by 5 girders
# M_LLIM: moment due to truck load and IM [kNm], Enright and Frangopol (1999), cov has been altered to match V_LLIM
# V_DCDW: shear due to DC and DW [kN], Enright and Frangopol (1999), evenly divided by 5 girders
# V_LLIM: shear due to truck load and IM [kN], Enright and Frangopol (1999)
# LL_ARRIVAL_RATE: arrival rate of maximum live load [times/year]
"""

# geometric properteis
BEAM_HEIGHT_NOM = 790.;          BEAM_HEIGHT_MEAN = 790.;          BEAM_HEIGHT_COV = 0.;         BEAM_HEIGHT_DISTR = 'deterministic';
BEAM_DEPTH_NOM = 663.0;          BEAM_DEPTH_MEAN = 663.0;          BEAM_DEPTH_COV = 0.03;        BEAM_DEPTH_DISTR = 'normal';
SHEAR_DEPTH_NOM = 663.0;         SHEAR_DEPTH_MEAN = 663.0;         SHEAR_DEPTH_COV = 0.03;       SHEAR_DEPTH_DISTR = 'normal';
FLANGE_WIDTH_NOM = 2600.;        FLANGE_WIDTH_MEAN = 2600.;        FLANGE_WIDTH_COV = 0.;        FLANGE_WIDTH_DISTR = 'deterministic';
FLANGE_THICK_NOM = 190.;         FLANGE_THICK_MEAN = 190.;         FLANGE_THICK_COV = 0.;        FLANGE_THICK_DISTR = 'deterministic';
WEB_WIDTH_NOM = 400.;            WEB_WIDTH_MEAN = 400.;            WEB_WIDTH_COV = 0.;           WEB_WIDTH_DISTR = 'deterministic';
# general information about reinforcement
BAR_TYPE = 'deformed'
DSC_NO = 0
DSC_NOM = 0.0;                   DSC_MEAN = 0.0;                   DSC_COV = 0.;                 DSC_DISTR = 'deterministic';
SHEAR_INTERVAL_NOM = 95.;        SHEAR_INTERVAL_MEAN = 95.;        SHEAR_INTERVAL_COV = 0;       SHEAR_INTERVAL_DISTR = 'deterministic';
# data for flexure limit state
COVER_NOM = 51.0;                COVER_MEAN = 51.0;                COVER_COV = 0.20;             COVER_DISTR = 'normal';            
DS_REGION1_NOM = 35.8;           DS_REGION1_MEAN = 35.8;           DS_REGION1_COV = 0.02;        DS_REGION1_DISTR = 'normal';
DS_REGION2_NOM = 35.8;           DS_REGION2_MEAN = 35.8;           DS_REGION2_COV = 0.02;        DS_REGION2_DISTR = 'normal';
REGION1_NO = 6
REGION2_NO = 2
DISTANCE_12 = 87.
# data for shear limit state
COVER_SHEAR_NOM = 38.0;          COVER_SHEAR_MEAN = 38.0;          COVER_SHEAR_COV = 0.20;       COVER_SHEAR_DISTR = 'normal';                       
DS_REGION1_SHEAR_NOM = 12.7;     DS_REGION1_SHEAR_MEAN = 12.7;     DS_REGION1_SHEAR_COV = 0.02;  DS_REGION1_SHEAR_DISTR = 'normal';
DS_REGION2_SHEAR_NOM = 12.7;     DS_REGION2_SHEAR_MEAN = 12.7;     DS_REGION2_SHEAR_COV = 0.02;  DS_REGION2_SHEAR_DISTR = 'normal';
REGION1_NO_SHEAR = 2
REGION2_NO_SHEAR = 0
DISTANCE_12_SHEAR = 0
# data for deck limit state
SHEAR_SPAN_DECK_NOM = 1250.;        SHEAR_SPAN_DECK_MEAN = 1250.;   SHEAR_SPAN_DECK_COV = 0.;    SHEAR_SPAN_DECK_DISTR = 'deterministic';
SPAN_DECK_NOM = 2500.;              SPAN_DECK_MEAN = 2500.;         SPAN_DECK_COV = 0.;          SPAN_DECK_DISTR = 'deterministic';
PLATE2SUPPORT_DECK_NOM = 0.;        PLATE2SUPPORT_DECK_MEAN = 0.;   PLATE2SUPPORT_DECK_COV = 0.; PLATE2SUPPORT_DECK_DISTR = 'deterministic';
DECK_HEIGHT_NOM = 190.;          DECK_HEIGHT_MEAN = 190.;          DECK_HEIGHT_COV = 0.;         DECK_HEIGHT_DISTR = 'deterministic';
DECK_DEPTH_NOM = 157.2;          DECK_DEPTH_MEAN = 157.2;          DECK_DEPTH_COV = 0.03;        DECK_DEPTH_DISTR = 'normal';
DECK_WIDTH_NOM = 1000.;          DECK_WIDTH_MEAN = 1000.;          DECK_WIDTH_COV = 0.;          DECK_WIDTH_DISTR = 'deterministic';
COVER_DECK_NOM = 25.4;           COVER_DECK_MEAN = 25.4;           COVER_DECK_COV = 0.20;        COVER_DECK_DISTR = 'normal';
DS_REGION1_DECK_NOM = 15.9;      DS_REGION1_DECK_MEAN = 15.9;      DS_REGION1_DECK_COV = 0.02;   DS_REGION1_DECK_DISTR = 'normal';
DS_REGION2_DECK_NOM = 15.9;      DS_REGION2_DECK_MEAN = 15.9;      DS_REGION2_DECK_COV = 0.02;   DS_REGION2_DECK_DISTR = 'normal';
REGION1_NO_DECK = 8
REGION2_NO_DECK = 0
DISTANCE_12_DECK = 0
FRP_DECK_RATIO = 10.
# material properties
LAMBDA_FC = 1    # factor for light weight concrete, for normal concrete LAMBDA_FC = 1,
ECU_NOM = 0.0035;                ECU_MEAN = 0.0035;                ECU_COV = 0.0035;             ECU_DISTR = 'deterministic';
FC_NOM = 20.7;                   FC_MEAN = 25.9;                   FC_COV = 0.15;                FC_DISTR = 'normal';\
                                                                   FT_COV = 0.20;                FT_DISTR = 'normal';\
                                                                   EC_COV = 0.12;                EC_DISTR = 'normal';
ES_NOM = 206e3;                  ES_MEAN = 206e3;                  ES_COV = 0;                   ES_DISTR = 'deterministic';
FSY_NOM = 276;                   FSY_MEAN = 310.5;                 FSY_COV = 0.12;               FSY_DISTR = 'lognormal';
FYC_NOM = 276;                   FYC_MEAN = 310.5;                 FYC_COV = 0.12;               FYC_DISTR = 'lognormal';
FSYV_NOM = 276;                  FSYV_MEAN = 310.5;                FSYV_COV = 0.12;              FSYV_DISTR = 'lognormal';
# model errors
#ME_FLEX_RC_NOM = 1.0;            ME_FLEX_RC_MEAN = 1.01;           ME_FLEX_RC_COV = 0.043;       ME_FLEX_RC_DISTR = 'lognormal';
#ME_SHEAR_RC_NOM = 1.0;           ME_SHEAR_RC_MEAN = 1.09;          ME_SHEAR_RC_COV = 0.115;      ME_SHEAR_RC_DISTR = 'lognormal';
# from Nowak 1999, ME_SHEAR_RC_COV is changed from 0.100 to 0.125 to make cov
# of shear resistance 0.155 (Nowak 1999), since cov of shear resistance is
# small due to the direct usage of sqrt(fc_smp) in beam.py
ME_FLEX_RC_NOM = 1.0;            ME_FLEX_RC_MEAN = 1.02;           ME_FLEX_RC_COV = 0.060;       ME_FLEX_RC_DISTR = 'lognormal';
ME_SHEAR_RC_NOM = 1.0;           ME_SHEAR_RC_MEAN = 1.075;         ME_SHEAR_RC_COV = 0.125;      ME_SHEAR_RC_DISTR = 'lognormal';
# load variables
M_DCDW_NOM = 217.1;              M_DCDW_MEAN = 228.0;              M_DCDW_COV = 0.00;            M_DCDW_DISTR = 'deterministic';         
V_DCDW_NOM = 94.98;              V_DCDW_MEAN = 99.72;              V_DCDW_COV = 0.00;            V_DCDW_DISTR = 'deterministic'; 
#M_LLIM_NOM = 686.2;              M_LLIM_MEAN = 365.7;              M_LLIM_COV = 0.20;            M_LLIM_DISTR = 'normal';
#V_LLIM_NOM = 336.2;              V_LLIM_MEAN = 158.8;              V_LLIM_COV = 0.20;            V_LLIM_DISTR = 'normal';
# for flexure and shear, from Nowak 1999
M_LLIM_NOM = 430.0*1.33;              M_LLIM_MEAN = 430.0*1.15*0.74;              M_LLIM_COV = 0.23;            M_LLIM_DISTR = 'normal';
V_LLIM_NOM = 220.6*1.33;              V_LLIM_MEAN = 220.6*1.15*0.68;              V_LLIM_COV = 0.23;            V_LLIM_DISTR = 'normal';
M_DESIGN_LANE = 95.0
V_DESIGN_LANE = 42.3
LL_ARRIVAL_RATE = 1000
# for deck, from Nowak 1999, single-lane governed, mean value divided by m=1.2
M_DCDW_DECK_NOM = 1.571;              M_DCDW_DECK_MEAN = 1.650;                   M_DCDW_DECK_COV = 0.00;       M_DCDW_DECK_DISTR = 'deterministic';
M_LLIM_DECK_NOM = 20.59*1.33;         M_LLIM_DECK_MEAN = 17.16*1.15*0.74;         M_LLIM_DECK_COV = 0.23;       M_LLIM_DECK_DISTR = 'normal';
LL_ARRIVAL_RATE_DECK = 220000
# LDF_FLEX = 0.78
# LDF_SHEAR = 0.86
LDF_FLEX = 0.82
LDF_SHEAR = 0.83
LL_MEAN_INCREASE_RATE = 0.00
LL_STD_INCREASE_RATE = 0.00


"""preventive maintenance options
# SERVICE_LIFE: service life [yr]
PM_FLEX_YR: preventive maintenance year for flexural rebars
"""
SERVICE_LIFE = 100.
RELIABILITY_DT = 10.
PM_FLEX_YR = 50.


""" FRP strengthening
# FRP_LAYOUT: 'none' --- no FRP strengthening
#             'flexural' --- flexural strengthening: w/o anchor (IC debonding or FRP rupture)
#             'flexural+anchor' --- flexural strengthening: w/anchor (FRP rupture)
#             'shear+side' --- shear strengthening side bonding
#             'shear+u' --- shear strengthening u-jacketing
#             'shear+w' --- shear strengthening complete wrapping
# DESIGN_MODEL:   'aci'--- ACI 440.2R
#                 'hk' --- Hong Kong 2013
# EFRP: modulus of FRP [MPa]
# FFRP: initial strength of FRP [MPa]
# TFRP: thickness of FRP [mm]
# bfrp: width of FRP, mm
# flex_ACI: model error of flexural strengthening (ACI)
# FRP_yr: time to strengthen
"""
# geometric and material
FRP_LAYOUT = 'none';
DESIGN_MODEL = 'hk';
FRP_DESIGN_YR = 1;
MAX_LAYER = 40;
FRP_MIN_THICK = 0.11;
TFRP_NEGLECT = 0.01;
FRP_FLEX_COV = 0.141;
FRP_SHEAR_COV = 0.212;
FRP_DECK_COV = 0.145;
RC_FLEX_COV = 0.133;
RC_SHEAR_COV = 0.155;
RC_DECK_COV = 0.133;
# flexural
INI_FRP_STRAIN = 0.0;
EFRP_NOM = 230e3;              EFRP_MEAN = 230e3;         EFRP_COV = 0.1;         EFRP_DISTR = 'lognormal';
FFRP_NOM = 3450.;              FFRP_MEAN = 3900.;         FFRP_COV = 0.10;        FFRP_DISTR = 'weibull';
TFRP_NOM = 0.110;              TFRP_MEAN = 0.110;         TFRP_COV = 0.;          TFRP_DISTR = 'deterministic';
BFRP_NOM = 200;                BFRP_MEAN = 200;           BFRP_COV = 0.;          BFRP_DISTR = 'deterministic';
FANCHOR_NOM = 3450.;           FANCHOR_MEAN = 3450.;      FANCHOR_COV = 0.;       FANCHOR_DISTR = 'deterministic';
PLATE2SUPPORT_NOM = 0.;        PLATE2SUPPORT_MEAN = 0.;   PLATE2SUPPORT_COV = 0.; PLATE2SUPPORT_DISTR = 'deterministic';
#model errors
ME_FLEX_FRPIC_NOM = 1.0;         ME_FLEX_FRPIC_MEAN = 1.13;        ME_FLEX_FRPIC_COV = 0.080;    ME_FLEX_FRPIC_DISTR = 'lognormal';
ME_FLEX_FRPRT_NOM = 1.0;         ME_FLEX_FRPRT_MEAN = 1.00;        ME_FLEX_FRPRT_COV = 0.120;    ME_FLEX_FRPRT_DISTR = 'lognormal';
ME_SHEAR_FRPS_NOM = 1.0;         ME_SHEAR_FRPS_MEAN = 1.54;        ME_SHEAR_FRPS_COV = 0.228;    ME_SHEAR_FRPS_DISTR = 'lognormal';
ME_SHEAR_FRPU_NOM = 1.0;         ME_SHEAR_FRPU_MEAN = 1.59;        ME_SHEAR_FRPU_COV = 0.192;    ME_SHEAR_FRPU_DISTR = 'lognormal';
ME_SHEAR_FRPW_NOM = 1.0;         ME_SHEAR_FRPW_MEAN = 1.44;        ME_SHEAR_FRPW_COV = 0.188;    ME_SHEAR_FRPW_DISTR = 'lognormal';
## shear: year 52
#EFRPV_NOM = 230e3;             EFRPV_MEAN = 230e3;        EFRPV_COV = 0.1;        EFRPV_DISTR = 'lognormal';
#FFRPV_NOM = 3450.;             FFRPV_MEAN = 3900.;        FFRPV_COV = 0.10;       FFRPV_DISTR = 'weibull';
#TFRPV_NOM = 0.166;             TFRPV_MEAN = 0.166;        TFRPV_COV = 0.;         TFRPV_DISTR = 'deterministic';
#SHEAR_SPAN_NOM = 1219.;        SHEAR_SPAN_MEAN = 1219.;   SHEAR_SPAN_COV = 0.;    SHEAR_SPAN_DISTR = 'deterministic';
#SPAN_NOM = 9100.;              SPAN_MEAN = 9100.;         SPAN_COV = 0.;          SPAN_DISTR = 'deterministic';
#FRP_INCLINE_NOM = 90.;         FRP_INCLINE_MEAN = 90.;    FRP_INCLINE_COV = 0.;   FRP_INCLINE_DISTR = 'deterministic';
#WFRPV_NOM = 50.;               WFRPV_MEAN = 50.;          WFRPV_COV = 0.;         WFRPV_DISTR = 'deterministic';
#SFRPV_NOM = 150.;              SFRPV_MEAN = 150.;         SFRPV_COV = 0.;         SFRPV_DISTR = 'deterministic';
#EM_FLEX_YR = 100.
#EM_SHEAR_YR = 52.
#EM_FLEX_FORM = 'IC'
#EM_SHEAR_FORM = 'W'
## U-jacket without anchors
#TFRPV_NOM = 0.166;             TFRPV_MEAN = 0.166;        TFRPV_COV = 0.;         TFRPV_DISTR = 'deterministic';
#WFRPV_NOM = 50.;               WFRPV_MEAN = 50.;          WFRPV_COV = 0.;         WFRPV_DISTR = 'deterministic';
#SFRPV_NOM = 50.;               SFRPV_MEAN = 50.;          SFRPV_COV = 0.;         SFRPV_DISTR = 'deterministic';
#EM_SHEAR_FORM = 'U'

# shear: year 78
EFRPV_NOM = 230e3;             EFRPV_MEAN = 230e3;        EFRPV_COV = 0.1;        EFRPV_DISTR = 'lognormal';
FFRPV_NOM = 3450.;             FFRPV_MEAN = 3900.;        FFRPV_COV = 0.10;       FFRPV_DISTR = 'weibull';
TFRPV_NOM = 0.166;             TFRPV_MEAN = 0.166;        TFRPV_COV = 0.;         TFRPV_DISTR = 'deterministic';
SHEAR_SPAN_NOM = 1219.;        SHEAR_SPAN_MEAN = 1219.;   SHEAR_SPAN_COV = 0.;    SHEAR_SPAN_DISTR = 'deterministic';
SPAN_NOM = 9100.;              SPAN_MEAN = 9100.;         SPAN_COV = 0.;          SPAN_DISTR = 'deterministic';
FRP_INCLINE_NOM = 90.;         FRP_INCLINE_MEAN = 90.;    FRP_INCLINE_COV = 0.;   FRP_INCLINE_DISTR = 'deterministic';
WFRPV_NOM = 50.;               WFRPV_MEAN = 50.;          WFRPV_COV = 0.;         WFRPV_DISTR = 'deterministic';
SFRPV_NOM = 100.;              SFRPV_MEAN = 100.;         SFRPV_COV = 0.;         SFRPV_DISTR = 'deterministic';
EM_FLEX_YR = 100.
EM_SHEAR_YR = 72.
EM_FLEX_FORM = 'IC'
EM_SHEAR_FORM = 'W'
# U-jacket without anchors
TFRPV_NOM = 0.498;             TFRPV_MEAN = 0.498;        TFRPV_COV = 0.;         TFRPV_DISTR = 'deterministic';
WFRPV_NOM = 50.;               WFRPV_MEAN = 50.;          WFRPV_COV = 0.;         WFRPV_DISTR = 'deterministic';
SFRPV_NOM = 100.;              SFRPV_MEAN = 100.;         SFRPV_COV = 0.;         SFRPV_DISTR = 'deterministic';
EM_SHEAR_FORM = 'U'

"""FRP and bonding deterioration
# FRP_degrade: 0 --- no deterioration
#              1 --- Karbhari's model
#              2 --- Davalos's model
#              3 --- Interpolation of field deterioration
# bond_degrade: 0 --- no deterioration
# bonding_deterioration: bonding deterioration data obtained from FEM
"""
FRP_DEGRADE = 0;
FRPdeg1_par_1 = -3.366; FRPdeg1_par_2 = 100;
FRPdeg2_par_1 = 0.8; FRPdeg2_par_2 = 1277.5;
# load FRP_deterioration_env.mat
BOND_DEGRADE = 0;
# load bonding_deterioration.mat
# partial safety factors
GAMMA_CONC_FLEX = 1.50
GAMMA_CONC_SHEAR = 1.25
GAMMA_STEEL = 1./0.87
GAMMA_BOND = 1.25
GAMMA_FRP = 1.40
