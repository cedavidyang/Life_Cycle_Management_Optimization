#import os
#strepath = os.path.join(os.path.expanduser('~'), 'FERUM')
#import sys
#sys.path.append(strepath)
from pyStRe import *

# component class
import os
import sys
import numpy as np
from scipy import stats

from constants.beamConstants import FRP_DESIGN_YR, SERVICE_LIFE
from constants.beamConstants import M_DCDW_MEAN, M_DCDW_COV, M_LLIM_MEAN, M_LLIM_COV
from constants.beamConstants import V_DCDW_MEAN, V_DCDW_COV, V_LLIM_MEAN, V_LLIM_COV
from constants.beamConstants import M_DCDW_DECK_MEAN, M_DCDW_DECK_COV,\
    M_LLIM_DECK_MEAN, M_LLIM_DECK_COV
from constants.beamConstants import LL_ARRIVAL_RATE, LL_ARRIVAL_RATE_DECK
from constants.corrosionConstants import TIME_INTERVAL
from constants.sampleConstants import COVT0_COV

from pyCEsmp import *
from pyre.distributions import *
from preprocessing import getCorrosionLossFromFile, getDeteriorationFromFile,\
    strengtheningFromCorrosionLoss

class Component(object):

    pfkeeping = {'flexure':-1.*np.ones((2,1)), 'shear':-1.*np.ones((2,1)),
        'deck':-1.*np.ones((2,1)), 'virginflexure':-1.*np.ones((2,1)),
        'virginshear':-1.*np.ones((2,1)), 'virgindeck':-1.*np.ones((2,1))}
    riskkeeping = {'flexure':-1.*np.ones((2,1)), 'shear':-1.*np.ones((2,1)),
        'deck':-1.*np.ones((2,1)), 'virginflexure':-1.*np.ones((2,1)),
        'virginshear':-1.*np.ones((2,1)), 'virgindeck':-1.*np.ones((2,1))}
    costkeeping = {'flexure':-1.*np.ones((2,1)), 'shear':-1.*np.ones((2,1)),
        'deck':-1.*np.ones((2,1)), 'virginflexure':-1.*np.ones((2,1)),
        'virginshear':-1.*np.ones((2,1)), 'virgindeck':-1.*np.ones((2,1))}

    @classmethod
    def resetPfKeeping(cls):
        cls.pfkeeping = {'flexure':-1.*np.ones((2,1)), 'shear':-1.*np.ones((2,1)),
                'deck':-1.*np.ones((2,1)), 'virginflexure':-1.*np.ones((2,1)),
                'virginshear':-1.*np.ones((2,1)), 'virgindeck':-1.*np.ones((2,1))}

    @classmethod
    def resetRiskKeeping(cls):
        cls.riskkeeping = {'flexure':-1.*np.ones((2,1)), 'shear':-1.*np.ones((2,1)),
                'deck':-1.*np.ones((2,1)), 'virginflexure':-1.*np.ones((2,1)),
                'virginshear':-1.*np.ones((2,1)), 'virgindeck':-1.*np.ones((2,1))}

    @classmethod
    def registerCost(cls, comp_type, str_yr, cost):
        tmp = np.array([[str_yr], [cost]])
        cls.costkeeping[comp_type] = \
                np.hstack((cls.costkeeping[comp_type], tmp))

    @classmethod
    def resetCostKeeping(cls):
        cls.costkeeping = {'flexure':-1.*np.ones((2,1)), 'shear':-1.*np.ones((2,1)),
                'deck':-1.*np.ones((2,1)), 'virginflexure':-1.*np.ones((2,1)),
                'virginshear':-1.*np.ones((2,1)), 'virgindeck':-1.*np.ones((2,1))}

    def __init__(self, comp_type, maintain_tag=None, str_yr=None):
        self.comp_type = comp_type
        self.maintain_tag = maintain_tag
        self.str_yr = str_yr
        self.service_time = None
        self.resistance_mean = None
        self.resistance_cov = None
        self.ceSmp = None
        # check input data
        if str_yr is not None:
            if str_yr<FRP_DESIGN_YR:
                self.str_yr = 0
                self.maintain_tag = False
            else:
                self.maintain_tag = True
        #if (maintain_tag is not None) and (str_yr is not None):
        #    if str_yr<FRP_DESIGN_YR and maintain_tag == True:
        #        self.maintain_tag = False
        #        print '[WARNING:] str_yr is conflict with maintain_tag, tag is changed'
        #    elif str_yr>=FRP_DESIGN_YR and maintain_tag == False:
        #        self.maintain_tag = True
        #        print '[WARNING:] str_yr is conflict with maintain_tag, tag is changed'

    def setServiceTime(self, time_array):
        self.service_time = time_array

    def setResistanceMean(self, resistance_mean):
        self.resistance_mean = resistance_mean

    def setResistanceCov(self, resistance_cov):
        self.resistance_cov = resistance_cov

    def loadDeterioration(self, datafile):
        if self.maintain_tag == True:
            loss_dict = getCorrosionLossFromFile(datafile, self.str_yr)
            mean_strength, tfrp = strengtheningFromCorrosionLoss(loss_dict, self.comp_type)
            mean_history, cov_history = getDeteriorationFromFile(datafile, self.str_yr,
                mean_strength_array, self.comp_type)
            self.service_time = np.arange(TIME_INTERVAL, SERVICE_LIFE+TIME_INTERVAL, TIME_INTERVAL)
            self.resistance_mean = mean_history
            self.resistance_cov = cov_history
        elif 'f' in self.comp_type.lower() or 'd' in self.comp_type.lower():
            # flexural strength rc
            mean_history = np.loadtxt(datafile)[20,:]
            std_history = np.loadtxt(datafile)[21,:]
            cov_history = std_history/mean_history
            self.service_time = np.loadtxt(datafile)[0,:]
            self.resistance_mean = mean_history
            self.resistance_cov = cov_history
        elif 's' in self.comp_type.lower():
            # shear strength rc
            mean_history = np.loadtxt(datafile)[22,:]
            std_history = np.loadtxt(datafile)[23,:]
            cov_history = std_history/mean_history
            self.service_time = np.loadtxt(datafile)[0,:]
            self.resistance_mean = mean_history
            self.resistance_cov = cov_history

    def setCESampling(self, ncomp, nadp, nsmp, nmain):
        # create smp object
        self.ceSmp = CEImportanceSmp(self.comp_type, ncomp, nadp, nsmp, nmain)

    def accumulativePf(self, timepoint, register=True):
        #if self.str_yr in Component.pfkeeping[self.comp_type][0,:]:
        #    indx = np.where(Component.pfkeeping[self.comp_type][0,:]==self.str_yr)[0][0]
        #    pf = Component.pfkeeping[self.comp_type][1,indx]
        #    return pf
        #else:
            # resistance array
            r_array = np.array([self.resistance_mean[0]])
            # live load array
            if 'f' in self.comp_type.lower():
                sl = Normal('slm', mean=M_LLIM_MEAN, stdv=M_LLIM_MEAN*M_LLIM_COV)
            elif 's' in self.comp_type.lower():
                sl = Normal('slv', mean=V_LLIM_MEAN, stdv=V_LLIM_MEAN*V_LLIM_COV)
            elif 'd' in self.comp_type.lower():
                sl = Normal('sld', mean=M_LLIM_DECK_MEAN, stdv=M_LLIM_DECK_MEAN*M_LLIM_DECK_COV)
            slDistr = sl.rv
            sl_array = np.array([slDistr])
            # rate array
            if 'f' in self.comp_type.lower():
                rate = LL_ARRIVAL_RATE
            elif 's' in self.comp_type.lower():
                rate = LL_ARRIVAL_RATE
            elif 'd' in self.comp_type.lower():
                rate = LL_ARRIVAL_RATE_DECK
            rate_array = np.array([rate])
            # gd_array
            gd_data = self.resistance_mean/self.resistance_mean[0] / \
                np.sqrt(1.+self.resistance_cov**2)
            gd_param = np.polyfit(self.service_time, gd_data, 4)
            gd = lambda x: np.poly1d(gd_param)(x)
            gd_array = np.array([gd])
            # gcov_array
            at_data = np.sqrt(np.log(self.resistance_cov**2+1.)/COVT0_COV**2)
            gcov_param = np.polyfit(self.service_time, at_data, 4)
            gcov = lambda x: np.poly1d(gcov_param)(x)
            gcov_array = np.array([gcov])

            ## create smp object
            #self.ceSmp = CEImportanceSmp(self.comp_type, ncomp, nadp, nsmp, nmain)
            self.ceSmp.setPreSmp(r_array, sl_array, rate_array, gd_array, gcov_array)

    #def accumulativePf(self, timepoint, register=True):
    #    if self.str_yr in Component.pfkeeping[self.comp_type][0,:]:
    #        indx = np.where(Component.pfkeeping[self.comp_type][0,:]==self.str_yr)[0][0]
    #        pf = Component.pfkeeping[self.comp_type][1,indx]
    #        return pf
    #    else:
            self.ceSmp.conductPreSmp(timepoint)
            self.ceSmp.setMainSmp()
            self.ceSmp.conductMainSmp(timepoint)
            ## store in pfkeeping
            #if register == True:
            #    tmp = np.array([[self.str_yr], [self.ceSmp.pf]])
            #    Component.pfkeeping[self.comp_type] = \
            #        np.hstack((Component.pfkeeping[self.comp_type], tmp))

            return self.ceSmp.pf

    def pointintimePf(self, timepoint, register=True):
        # resistance array
        if timepoint == 0:
            rmean = self.resistance_mean[0]
            rcov = self.resistance_cov[0]
        else:
            rmean = self.resistance_mean[self.service_time==timepoint][0]
            rcov = self.resistance_cov[self.service_time==timepoint][0]
        if 'f' in self.comp_type:
            sd = M_DCDW_MEAN
        elif 's' in self.comp_type:
            sd = V_DCDW_MEAN
        elif 'd' in self.comp_type:
            sd = M_DCDW_DECK_MEAN
        rstd = rmean*rcov
        rmean = rmean-sd
        R = stats.lognorm(np.sqrt(np.log(1+rcov**2)), scale=rmean/np.sqrt(1+rcov**2))
        # live load array
        if 'f' in self.comp_type.lower():
            mean=M_LLIM_MEAN/0.74
            stdv=M_LLIM_MEAN/0.74*M_LLIM_COV
        elif 's' in self.comp_type.lower():
            mean=V_LLIM_MEAN/0.68
            stdv=V_LLIM_MEAN/0.68*V_LLIM_COV
        elif 'd' in self.comp_type.lower():
            mean= M_LLIM_DECK_MEAN/0.74
            stdv=M_LLIM_DECK_MEAN/0.74*M_LLIM_DECK_COV
        loc=mean-np.sqrt(6)*stdv/np.pi*np.euler_gamma,
        scale=np.sqrt(6)*stdv/np.pi
        SL = stats.gumbel_r(loc=loc,scale=scale)
        rvs = [R, SL]
        corr = np.eye(2)
        probdata = ProbData(names=['R','SL'], rvs=rvs, corr=corr, startpoint=[rmean, mean], nataf=False)
        rvs = [R,SL]

        def gf1(x, param=None):
            return x[0]-x[1]
        def dgdq1(x, param=None):
            dgd1 = 1.
            dgd2 = -1.
            return [dgd1,dgd2]

        analysisopt = AnalysisOpt(gradflag='DDM', recordu=False, recordx=False, flagsens=False, verbose=False)
        gfunc = Gfunc(gf1, dgdq1)

        formBeta = CompReliab(probdata, gfunc, analysisopt)
        formresults = formBeta.form_result()
        return formresults.pf1


if __name__ == '__main__':
    import time
    import datetime
    from constants.sampleConstants import NUM_COMPONENT, NUM_ADAPTATION,\
        NUM_MAIN_SMP, NUM_KOPT_SMP, NUM_PRE_SMP
    from constants import END_AGE, RELIABILITY_DT

    print "debugging of component class"
    start_delta_time = time.time()

    test_component = Component('shear', maintain_tag=False, str_yr=0)
    RC_CS_SHEAR_PATH = os.path.join(os.path.abspath('./'), 'data', 'rc_cs', 'shear')
    shear_datafile = os.path.join(RC_CS_SHEAR_PATH, 'LWS_results.txt')
    test_component.loadDeterioration(shear_datafile)
    test_component.setCESampling(NUM_COMPONENT, NUM_ADAPTATION, NUM_PRE_SMP, NUM_MAIN_SMP)
    year = 100
    #pf = test_component.accumulativePf(year)
    pf = test_component.pointintimePf(year)
    print 'failure probability at year {:d}: {:.4e}'.format(year, pf)

    delta_time = time.time() - start_delta_time
    print 'DONE: {} s'.format(str(datetime.timedelta(seconds=delta_time)))

    #reliability_data = os.path.join(RC_CS_SHEAR_PATH, 'reliability_results_parallel.npz')
    #res = np.load(reliability_data)
    ## time array
    #time_array = np.arange(RELIABILITY_DT,END_AGE+RELIABILITY_DT,RELIABILITY_DT)
    #time_array = np.insert(time_array, 0, 1)
    #pf_cmp = res['pf'][time_array==year]
    #print 'failure probability at year {:d}: {:.4e}'.format(year, pf_cmp[0])
