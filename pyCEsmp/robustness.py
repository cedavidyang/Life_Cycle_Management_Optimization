# robustness related functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#from scipy import optimize
import scipy.stats as stats
from pyre.distributions import *


## latin hypercube sampling for initial peaks

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
        #norm_rv = stats.norm(loc=r_array[iV].stats(moments='m'), scale=np.sqrt(r_array[iV].stats(moments='v')))
        #smps[:, iV] = [norm_rv.ppf((arrange[i, iV] + 0.5)/nCompo) for i in range(nCompo)]

    return smps

## simulated annealing
def anneal_LHS(r_array, nCompo, corrTarget, params=None, options=None):
    """ yields a generator
    """
    user_options = options
    options = {
        'maxfev'       : 1.e5,  # Default, formerly `maxeval`.
        'maxiter'      : 5000,   # Non-default value.
        'maxaccept'    : 1000,  # Default value.
        'ftol'         : 1e-2,  # Default, formerly `feps`.
        'T0'           : initial_temperature(corrTarget),  # Default value.
        'Tf'           : 1e-12, # Default value.

        'k'            : 0.5,   # Default value.
        'learn_rate'   : 0.5,   # Default value.
        'quench'       : 1.0,   # Default value.
        'm'            : 1.0,   # Default value.
        'n'            : 1.0,   # Default value.

        'dwell'        : 50,   # Non-default value.
        'disp'         : False   # Default value.
        }
    if user_options != None:
        options.update(user_options)

    nV = r_array.shape[0]
    arrange0 = latinHypercubeSmp(nV, nCompo)
    smps0 = arrange2smps(r_array, arrange0)
    corrLhs = np.corrcoef(smps0.T)
    temp0 = options['T0']
    normBest = np.linalg.norm( (corrLhs-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
    normParent = np.linalg.norm( (corrLhs-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
    tempCurrent = np.copy(temp0)
    arrangeBest = np.copy(arrange0)

    iFuncEval= 1
    iTotalIter = 1
    iAccept = 1
    tempRatio = 1. / (1.+options['k'])
    while not (tempCurrent <= options['Tf']      or
               normBest <= options['ftol']       or
               iFuncEval > options['maxfev']   or
               iTotalIter > options['maxiter']  or
               iAccept > options['maxaccept']
               ):

        for iDwell in range(options['dwell']):
            if options['disp'] is True:
                smpsBest = arrange2smps(r_array, arrangeBest)
                corrBest = np.corrcoef(smpsBest.T)
                yield smpsBest, arrangeBest, corrBest, normBest, normParent, tempCurrent, np.array(iTotalIter)
            arrangeOffspring, corrOffspring = set_change(r_array, nCompo)
            normOffspring = np.linalg.norm( (corrOffspring-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)
            iFuncEval += 1
            # save (update) the best
            if normBest>normOffspring:
                normBest = np.copy(normOffspring)
                arrangeBest = np.copy(arrangeOffspring)

            # check if accept the child
            if normParent>normOffspring:
                normParent = np.copy(normOffspring)
                iAccept += 1
            else:
                dE = normOffspring - normParent
                u = np.random.rand(1)
                if u > np.exp( -dE/tempCurrent ):
                    iTotalIter += 1
                    continue
                else:    # accept the child with larger norm
                    iAccept += 1
                    normParent = np.copy(normOffspring)

            iTotalIter += 1


        # lower the temperature
        tempCurrent = tempCurrent*tempRatio


    smpsBest = arrange2smps(r_array, arrangeBest)
    corrBest = np.corrcoef(smpsBest.T)
    normBest = np.linalg.norm( (corrBest-corrTarget), 'fro' ) / (nV*(nV-1.)/2.)

    if options['disp'] is False:
        yield smpsBest, arrangeBest, corrBest, normBest

def initial_temperature(corrTarget):
    nV = corrTarget.shape[0]

    isPos = corrTarget>=0
    isNeg = corrTarget<=0

    corrRemote = np.copy(corrTarget)
    corrRemote[isPos] = -1.0
    corrRemote[isNeg]= 1.0
    corrRemote[np.eye(nV, dtype='bool')] = 1.0

    tempIni = np.linalg.norm( corrRemote-corrTarget, 'fro' ) / (nV*(nV-1.)/2.)

    return tempIni

def set_change(r_array, nCompo):
    nV = r_array.shape[0]
    arrangeOffspring = latinHypercubeSmp(nV, nCompo)
    smps = arrange2smps(r_array, arrangeOffspring)
    corrOffspring = np.corrcoef(smps.T)

    return arrangeOffspring, corrOffspring

## initial peaks
def getInitialPeaks(r_array, nCompo, corrTarget=None, options=None):
    user_options = options
    options = {'anneal': False,
               'disp': False,}
    if user_options is not None:
        options.update(user_options)

    nV = r_array.shape[0]
    if options['anneal'] is True:
        if corrTarget is None:
            corrTarget = np.eye(nV)
        if options['disp'] is True:
            ## postprocessing
            data_gen = lambda: anneal_LHS(r_array, nCompo, corrTarget, options={'disp':True})
            fig, ax = plt.subplots()
            lines = ax.plot([], [], 'b', [], [], 'g', [], [], 'r')
            t0 = initial_temperature(corrTarget)
            ax.set_ylim(-0., t0*1.1)
            ax.set_xlim(0, 2)
            ax.grid()
            xdata, y1data, y2data, y3data = [], [], [], []
            def run(data):
                # update the data
                smpsBest, arrangeBest, corrBest, normBest, normParent, tempCurrent, iTotalIter = data
                xdata.append(iTotalIter)
                y1data.append(normParent)
                y2data.append(normBest)
                y3data.append(tempCurrent)
                xmin, xmax = ax.get_xlim()

                if iTotalIter >= xmax:
                    ax.set_xlim(xmin, 1.2*xmax)
                    ax.figure.canvas.draw()

                lines[0].set_data(xdata, y1data)
                lines[1].set_data(xdata, y2data)
                lines[2].set_data(xdata, y3data)

                return lines
            ani = animation.FuncAnimation(fig, run, data_gen, repeat=False)
            plt.show()
        res_generator = anneal_LHS(r_array, nCompo, corrTarget, options={'disp':False})
        for value in res_generator:
            res = value
        return res[0]
    else:
        arrange = latinHypercubeSmp(nV, nCompo)
        initialPeaks = arrange2smps(r_array, arrange)
        return initialPeaks

def getInitialCovar(r_array, nCompo):
    nV = r_array.shape[0]
    var = np.array([rv.stats(moments='v') for rv in r_array])
    covarV0 = np.zeros((nCompo, nV, nV))
    for iCompo in range(nCompo):
        covarV0[iCompo, :, :] = np.diag(var)

    return covarV0



if __name__ == '__main__':

    resistRv = Lognormal('r', mean=1790.1, stdv=283.37)
    rDistr = resistRv.rv
    r_array = np.array([])

    for i in range(2):
        r_array = np.append(r_array, rDistr)
    nCompo = 4
    corrTarget = np.array([[1., -1.], [-1., 1.]])
    peaks = getInitialPeaks(r_array, nCompo, corrTarget=corrTarget, options={'anneal':True, 'disp':False}) 
    covarV0 = getInitialCovar(r_array, nCompo)

    #for i in range(5):
    #    r_array = np.append(r_array, rDistr)
    #nCompo = 10
    #corrTarget = np.eye(5)
    #peaks = getInitialPeaks(r_array, nCompo, corrTarget=corrTarget, options={'anneal':True, 'disp':True})
    #print peaks

    #for i in range(10):
    #    r_array = np.append(r_array, rDistr)
    #nCompo = 20
    #corrTarget = np.eye(10)
    #peaks = getInitialPeaks(r_array, nCompo, corrTarget=corrTarget, options={'anneal':False, 'disp':False})

