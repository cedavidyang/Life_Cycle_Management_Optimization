# robustness related functions
import numpy as np
from scipy import optimize


# latin hypercube sampling for initial peaks

def latinHypercubeSmp(seed, nV, nCompo):
    lhsArrange = np.zeros((nCompo, nV))
    for iV in range(nV):
        np.random.seed( int((iV+1)*np.abs(seed)) )
        lhsArrange[:, iV] = np.random.permutation(nCompo)
        
    return lhsArrange
    

# simulated annealing
def initialPeaks(r_array, nCompo, targetCorr, x0):
    nV = r_array.shape[0]
    def min_func(x):
        arrange = latinHypercubeSmp(int(x), 2, 4)
        sampleCorr = np.corrcoef(arrange.T)
        norm = np.linalg.norm( sampleCorr - targetCorr, 'fro' )
    
        return norm
    # bestArrange = optimize.anneal(min_func, x0, full_output='True')
    # bestArrange = optimize.anneal(min_func, x0)
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
    bestArrange = optimize.minimize(min_func, x0, method='Anneal', options=myopts)
    
    return bestArrange
    
    
if __name__ == '__main__':
    test = initialPeaks(np.array([1,1]), 4, np.array([[1.0, -1.], [-1., 1.0]]), 0.0)
    best = latinHypercubeSmp(test['x'], 2, 4)