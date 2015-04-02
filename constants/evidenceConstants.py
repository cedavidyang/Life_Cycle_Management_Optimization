# evidence constants

# condition state constants
CS1_UB = 0.06
CS2_LB = 0.06
CS2_UB = 1.6
CS3_LB = 1.6
CS3_UB = 3.2
CS4_LB = 3.2
CRACK_MEASUREMENT_ABSERROR = 0.02    #mm
CRACK_MEASUREMENT_RELERROR = 0.00    #cov

# half cell potential
ISINI_HALF_CELL_MEAN = -400
ISINI_HALF_CELL_STD = 100
NOINI_HALF_CELL_MEAN = -150
NOINI_HALF_CELL_STD = 50

# icorr tests
MEASURED_ICORR_COV = 1.00

## evidence dict constants
NO_INI_YR = 18
IS_INI_YR = 20

HC_DATA1 = {'year': 10, 'potential':-100}
HC_DATA2 = {'year': 30, 'potential':-400}

CORROSION_RATE_DATA1 = {'year': 30, 'current': 2.0}
CORROSION_RATE_DATA2 = {'year': 50, 'current': 1.5}

NO_CRK_YR = 48
IS_CRK_YR = 50

# for shear
CS1_TO_CS2_SHEAR = 8
CS2_TO_CS3_SHEAR = 20
CS3_TO_CS4_SHEAR = 34
EVIDENCE_END_YR_SHEAR = 50

# for flexure
CS1_TO_CS2_FLEX = 10
CS2_TO_CS3_FLEX = 26
CS3_TO_CS4_FLEX = 42
EVIDENCE_END_YR_FLEX = 50

# for deck
CS1_TO_CS2_DECK = 8
CS2_TO_CS3_DECK = 18
CS3_TO_CS4_DECK = 28
EVIDENCE_END_YR_DECK = 50
