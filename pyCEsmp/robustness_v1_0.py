# robustness related functions
import numpy as np
#from scipy import optimize
import scipy.stats as stats
from pyre import *


# latin hypercube sampling for initial peaks

def latinHypercubeSmp(nV, nCompo):
    lhsArrange = np.zeros((nCompo, nV))
    for iV in range(nV):
        lhsArrange[:, iV] = np.random.permutation(nCompo)
        
    return lhsArrange
    
def arrange2smps(r_array, arrange):
    nCompo, nV = arrange.shape
    smps = np.zeros(arrange.shape)
    for iV in range(nV):
        smps[:, iV] = [r_array[iV].ppf((arrange[i, iV] + 0.5)/nCompo) for i in range(nCompo)]
        
    return smps

# simulated annealing
def anneal_LHS(r_array, arrange0, corrTarget, params=None, options=None):
    if options == None:
        options = {
        'maxfev'       : None,  # Default, formerly `maxeval`.
        'maxiter'      : 5000,   # Non-default value.
        'maxaccept'    : None,  # Default value.
        'ftol'         : 1e-6,  # Default, formerly `feps`.
        'T0'           : initial_temperature(corrTarget),  # Default value.
        'Tf'           : 1e-12, # Default value.
        'boltzmann'    : 1.0,   # Default value.
        'learn_rate'   : 0.5,   # Default value.
        'quench'       : 1.0,   # Default value.
        'm'            : 1.0,   # Default value.
        'n'            : 1.0,   # Default value.
        'lower'        : -100,   # Non-default value.
        'upper'        : +100,   # Non-default value.
        'dwell'        : 250,   # Non-default value.
        'disp'         : True   # Default value.
        }

    nCompo, nV = arrange0.shape
    smps0 = arrange2smps(r_array, arrange0)
    corrLhs = np.corrcoef(smps0.T)
    
    normBest = np.linalg.norm( (corrLhs-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
    normParent = np.linalg.norm( (corrLhs-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
    tempCurrent = np.copy(temp0)
    arrangeBest = np.copy(arrange0)
    
    iFuncEval= 1
    iIter= 1
    nTotalLoop = 1
    temp0 = options['T0']
    tempMin = options['Tf']
    normMin = options['ftol']
    nDwell = options['dwell']
    while tempCurrent>tempMin and normBest>normMin:
        iLoop = 1
        while iLoop <= nLoop:
            [corrOffspring, arrangeOffspring] = set_change(r_array, nCompo)
            normOffspring = np.linalg.norm( (corrOffspring-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
            iFuncEval += 1
            
            if normBest>normOffspring:
                normBest = np.copy(normOffspring)
                arrangeBest = np.copy(arrangeOffspring)
        
            if normParent>normOffspring:
                normParent = np.copy(normOffspring)
                iLoop = iLoop+1
                nTotalLoop = nTotalLoop+1
            else:
                dE = normOffspring - normParent
                u = np.random.rand(1)
                if u > exp( -dE/tempCurrent ):
                    iIter = iIter + 1
                    if iIter>nMax:
                        break
                    else:
                        continue
                else:
                    normParent = np.copy(normOffspring)
                    iLoop = iLoop+1
                    nTotalLoop = nTotalLoop+1
        
            iIter = iIter + 1

        tempCurrent = tempCurrent*annealTemp
        if iIter>nMax:
            break
            
    smpsBest = arrange2smps(arrangeBest) 
    corrBest = np.corrcoef(smpsBest.T)
    normBest = np.linalg.norm( (corrBest-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
    
    return smpsBest, arrangeBest, corrBest, normBest
    

def initial_temperature(corrTarget):
    nV = corrTarget.shape[0]

    isPos = corrTarget>=0
    isNeg = corrTarget<=0

    corrRemote = np.copy(corrTarget)
    corrRemote( isPos ) = -1
    corrRemote( isNeg ) = 1
    corrRemote( np.eye(nV, dtype='bool') ) = 1

    tempIni = np.linalg.norm( corrRemote-corrTarget, 'fro' ) / (nV*(nV-1.)/2.)
    
    return tempIni
    
def set_change(r_array, nCompo):
    nV = r_array.shape[0]
    arrangeOffspring = latinHypercubeSmp(nV, nCompo)
    smps = arrange2smps(r_array, arrangeOffspring)
    corrOffspring = np.corrcoef(smps.T)

    return smps, corrOffspring
    
if __name__ == '__main__':
    myopts = {
        'maxfev'       : None,  # Default, formerly `maxeval`.
        'maxiter'      : 5000,   # Non-default value.
        'maxaccept'    : None,  # Default value.
        'ftol'         : 1e-6,  # Default, formerly `feps`.
        'T0'           : None,  # Default value.
        'Tf'           : 1e-12, # Default value.
        'boltzmann'    : 1.0,   # Default value.
        'learn_rate'   : 0.5,   # Default value.
        'quench'       : 1.0,   # Default value.
        'm'            : 1.0,   # Default value.
        'n'            : 1.0,   # Default value.
        'lower'        : -100,   # Non-default value.
        'upper'        : +100,   # Non-default value.
        'dwell'        : 250,   # Non-default value.
        'disp'         : True   # Default value.
        }
    resistRv = Lognormal('r', mean=1790.1, stdv=283.37)
    rDistr = resistRv.rv
    r_array = np.array([rDistr, rDistr]) 
    variables = arrange2variable(r_array, np.array([[0,1,2,3], [0,1,2,3]]).T)