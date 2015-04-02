# -*- coding: utf-8 -*-
# evidence
import numpy as np
from constants import *


## create evidence dict
service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)

evidence_list = np.empty(service_time.shape)
evidence_list[:] = np.nan
evidence_dict = {'iniStat': np.copy(evidence_list),
                 'halfCell': np.copy(evidence_list),
                 'icorr': np.copy(evidence_list),
                 'crkStat': np.copy(evidence_list),
                 'conditionState': np.copy(evidence_list)}

# node iniStat
#evidence_dict['iniStat'][service_time==NO_INI_YR] = False
#evidence_dict['iniStat'][service_time==IS_INI_YR] = True

# half cell potential evidence
#evidence_dict['halfCell'][service_time==HC_DATA1['year']] = HC_DATA1['potential']
#evidence_dict['halfCell'][service_time==HC_DATA2['year']] = HC_DATA2['potential']

# corrosion rate evidence
#evidence_dict['icorr'][service_time==CORROSION_RATE_DATA1['year']] = CORROSION_RATE_DATA1['current']
#evidence_dict['icorr'][service_time==CORROSION_RATE_DATA2['year']] = CORROSION_RATE_DATA2['current']

# node crkStat
#evidence_dict['crkStat'][service_time==NO_CRK_YR] = False
#evidence_dict['crkStat'][service_time==IS_CRK_YR] = True

# condition state evidence
evidence_dict['conditionState'][service_time<=CS1_TO_CS2] = 1
evidence_dict['conditionState'][np.logical_and(service_time>CS1_TO_CS2, service_time<=CS2_TO_CS3)] = 2
evidence_dict['conditionState'][np.logical_and(service_time>CS2_TO_CS3, service_time<=CS3_TO_CS4)] = 3
evidence_dict['conditionState'][np.logical_and(service_time>CS3_TO_CS4, service_time<=EVIDENCE_END_YR)] = 4

#evidence_dict['conditionState'][service_time==CS1_TO_CS2] = 2
#evidence_dict['conditionState'][service_time==CS2_TO_CS3] = 3
#evidence_dict['conditionState'][service_time==CS3_TO_CS4] = 4
