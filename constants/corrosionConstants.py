# Constants used in the analysis

## A_2_C: aggregate-to-cement ratio, based on the figure on the following website: http://www.cement.org/cement-concrete-basics/how-concrete-is-made
#A_2_C = (0.41 + 0.26) / 0.11
# A_2_C: aggregate-to-cement ratio, based on Jiekai Zhou's study
A_2_C = 3.0
## RHO_C, [kg/m3]: density of cement, a density of 1000-1600 kg/m3 - See more at: http://teaching.ust.hk/~civl111/CHAPTER4.pdf
# RHO_C, [-]: specific gravity of cement, 3.16 (Papadakis et al. 1996)
RHO_C = 3.16
## RHO_A, [kg/m3]: density of aggregate, Normal weight aggregate: The aggregate has unit weight of 1520-1680 kg/m3, See more at: http://teaching.ust.hk/~civl111/CHAPTER3.pdf
# RHO_A, [-]: specific gravity of aggregates, 2.6 (Papadakis et al. 1996)
RHO_A = 2.6
# AGE_EXP: aging exponent (fib 2006), for Portland cement concrete Beta(m=0.3, s=0.12, a=0.0, b=1.0), only mean value is used
AGE_EXP = 0.20
## C_CL_0: existing chloride content by the end of the previous time interval, in [kg/m3] (unit should be the same as the other chloride contents)
#C_CL_0 = 0.0

# D_H2O, [mm2/yr]: chloride diffusion coefficient in an infinite solution (Vu and Stewart 2000)
D_H2O = 1.6e-5 * 1.0e2 * (365.*24.*3600)
# D_COV: COV of apparent diffusion coefficient (Vu and Stewart 2000)
D_COV = 0.30
## D0_COV: COV of reference diffusion coefficient (Papakonstantinou and Shinozuka 2013)
#D0_COV = 0.2
# CL_INI: initial chloride content in concrete (assumed to be 0.0), in [kg/m3] (unit should be the same as the other chloride contents) 
CL_INI = 0.0
 
# TEMP_REF, [K]: reference temperature Tref = 20 C = 293.15 K, fib (2006)
TEMP_REF = 293.15
# TEMP_MEAN, [K]: mean temperature in Hong Kong
TEMP_MEAN = 293.15
# TEMP_STDV, [K]: stdv of annual temperature in [K], see more at: http://en.wikipedia.org/wiki/Hong_Kong#Geography_and_climate
TEMP_STDV = (36.1 - 0.0) / 6
# B_E, regression variable be of environmental transfer variable ke, fib(2006), only mean value is used
# B_E: be, a parameter influencing environmental vairable ke, [K] (fib 2006)
B_E = 4800.
# B_E = 35000.0/8.314 # Papakonstantinou and Shinozuka 2013
# TRANSFER_VARIABLE, transfer variable (fib 2006), kt=1
TRANSFER_VARIABLE = 1.0

# TIME_REF, [yr]: reference time Tref = 28 days = 0.0767 yr
TIME_REF = 0.0767

# DISTANCE_COAST, [km]: distance to coast
DISTANCE_COAST = 0.1

# CS_COV: COV of surface chloride content, Vu and Stewart (2000)
CS_COV = 0.5

## C_CRIT_MEAN: mean of critial chloride content, [kg/m3], Papakonstantinou and Shinozuka (2013)
#C_CRIT_MEAN = 1.2
## C_CRIT_COV: COV of critial chloride content, [kg/m3], Papakonstantinou and Shinozuka (2013)
#C_CRIT_COV = 0.24
# C_CRIT_MEAN: mean of critial chloride content, [kg/m3], Papakonstantinou and Shinozuka (2013)
C_CRIT_MEAN = 0.9
# C_CRIT_COV: COV of critial chloride content, [kg/m3], Papakonstantinou and Shinozuka (2013)
C_CRIT_COV = 0.19
C_CRIT_DISTR = 'uniform'

# COVER_COV: COV of concrete cover, Enright and Frangopol (1997)
COVER_COV = 0.20

# RC_VARIATION: std of temperary random varible to preserve variability of Rc, N(0, 0.1203) in Eq. 12 (Papakonstantinou and Shinozuka 2013)
RC_VARIATION = 0.1203
# ICORR_VARIATION: std of temperary random variable to preserve variability of icorr, N(0, 0.3312) in Eq. 11 (Papakonstantinou and Shinozuka 2013)
ICORR_VARIATION = 0.3312

# CORROSION_CONST1, CORROSION_CONST2: constants in corrosion loss equation (Papakonstantinou and Shinozuka 2013)
CORROSION_CONST1 = 2.315E-4
CORROSION_CONST2 = 31536.E3
# NACM2_TO_AM2: icorr in module propagation is in nA/cm2, change it to A/m2
# MM_TO_M: reinforcement diameter, from [mm] to [m]
CORROSION_NACM2_TO_AM2 = 1.E-2
CORROSION_MM_TO_M = 1.E-3
KG_TO_G = 1.E3

# RHO_S, [kg/m3], density of steel reinforcement, PAY ATTENTION TO THE UNIT
RHO_S = 7.85E3

# CREEP_COEFFICIENT: creep coefficient of concrete modulus (Papakonstantinou and Shinozuka 2013)
CREEP_COEFFICIENT = 2.0
# CONCRETE_POISSON: Poisson's ratio of concrete, 0.18 (Papakonstantinou and Shinozuka 2013)
CONCRETE_POISSON = 0.18

# GAMMA_VOL_MEAN: Alpha of relative volume ratio, 1.6 (Papakonstantinou and Shinozuka 2013)
GAMMA_VOL_ALPHA = 1.6
# GAMMA_VOL_STD: Beta of relative volume ration, 4.0 (Papakonstantinou and Shinozuka 2013)
GAMMA_VOL_BETA = 4.0
# GAMMA_VOL_LB: lower bound of relative volume ratio, 1.695 (Papakonstantinou and Shinozuka 2013)
GAMMA_VOL_LB = 1.695
# GAMMA_VOL_UB: upper bound of relative volume ratio, 6.3 (Papakonstantinou and Shinozuka 2013)
GAMMA_VOL_UB = 6.3
# POROUS_LB: lower bound of thickness of porous zone, 5.e-3 in [mm] (Papakonstantinou and Shinozuka 2013)
POROUS_LB = 5.e-3
# POROUS_UB: upper bound of thickness of porous zone, 12.e-2 in [mm] (Papakonstantinou and Shinozuka 2013)
POROUS_UB = 12.e-2

# CRACK PROPAGATION CONSTANTS: Papakonstantinou and Shinozuka 2013
CRACK_K = 0.1007
CRACK_GAMMA1 = -0.4388
CRACK_GAMMA2 = -0.0399
CRACK_COV = 0.54
# DIFFUSION PROPAGATION CONSTANTS: Papakonstantinou and Shinozuka 2013
DIFFUSION_CONST1 = 16.3894
DIFFUSION_CONST2 = 8.3186
DIFFUSION_K = 22.8745
DIFFUSION_COV = 0.2091
DIFFUSION_UB = 20.0

# N_SMP: number of samples for LWS
N_SMP = int(1E6)
# N_LHS_SMP = 1000, fluctuation of resistance cov: (.140-.130)/.135=0.07
# N_LHS_SMP = 2000, fluctuation of resistance cov: (.137-.129)/.133=0.06
# N_LHS_SMP = 3000, fluctuation of resistance cov: (.136-.131)/.133=0.04
# N_LHS_SMP = 4000, fluctuation of resistance cov: (.136-.131)/.133=0.04
# N_LHS_SMP = 5000, fluctuation of resistance cov: (.136-.132)/.134=0.03
# N_LHS_SMP = 10000, fluctuation of resistance cov: (.136-.132)/.134=0.03
N_LHS_SMP = 200
ACCEPT_WEIGHT = 0.0005
MAX_ITER_LW = 100
SUM_WEIGHT = 1000.0
# TIME
START_AGE = 0.
END_AGE = 200.0
TIME_INTERVAL = 2.0

## constants for concrete and steel properties
#CONCRETE_GRADE = 20.7
#CONCRETE_STRENGTH = 25.9
#STEEL_DIAMETER = 35.81
#CONCRETE_COVER = 30

## FRP deterioration


## bonding deterioration
# one layer continuous
# poly4 fitting from MATLAB: f(x) = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5
BOND_P1 = 3.393e-9
BOND_P2 = -3.889e-7
BOND_P3 = 1.504e-5
BOND_P4 = -5.223e-4
BOND_P5 = 0.9992
