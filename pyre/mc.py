#!/usr/bin/python -tt
# -*- coding: utf-8 -*-

import os
import numpy as np
import math
import scipy.optimize as opt
import scipy.special as spec
import matplotlib.pyplot as plt

from model import *
from correlation import *
from distributions import *
from cholesky import *
from limitstate import *
from stepsize import *
from form import *

class MonteCarlo(object):
  """Monte Carlo Simulation

  The preceding sections describe some methods for determining the reliability
  index :math:`\\beta` for some common forms of the limit state
  function. However, it is sometimes extremely difficult or impossible to find
  :math:`\\beta`. [Nowak2000]_

  In this case, the probability of failure :math:`p_f` may also be estimated
  by numerical simulation methods. A large variety of simulation techniques
  can be found in the literature, indeed, the most commonly used method is the
  Monte Carlo method. [Faber2009]_

  The principle of simulation methods is to carry out random sampling in the
  physical (or standardized) space. For each of the samples the limit state
  function is evaluated to figure out, whether the configuration is desired or
  undesired. The probability of failure :math:`p_f` is estimated by the number
  of undesired configurations, respected to the total numbers of
  samples. [Lemaire2010]_

  :Attributes:
    - analysis_option (AnalysisOption): Option for the structural analysis
    - limit_state (LimitState): Information about the limit state
    - stochastic_model (StochasticModel): Information about the model
  """

  def __init__(self,analysis_options=None,limit_state=None,stochastic_model=None):

    # The stochastic model
    if stochastic_model == None:
      self.model = StochasticModel()
    else:
      self.model = stochastic_model

    # Options for the calculation
    if analysis_options == None:
      self.options = AnalysisOptions()
    else:
      self.options = analysis_options

    # The limit state function
    if limit_state == None:
      self.limitstate = LimitState()
    else:
      self.limitstate = limit_state

    self.nrv = self.model.getLenMarginalDistributions()
    self.point = None
    self.covariance = None
    self.cholesky_covariance = None
    self.inverse_covariance = None
    self.sum_q = None
    self.sum_q2 = None
    self.q_bar = None
    self.cov_q_bar = None
    self.factors = None
    self.k = None
    self.done = None
    self.block_size = None
    self.u = None
    self.x = None
    self.beta = None
    self.Pf = None
    self.G = None
    self.I = None
    self.q = None
    self.Pf = None
    self.beta = None
    self.all_X = None
    self.all_G = None
    self.bins = None

    # Computation of modified correlation matrix R0
    computeModifiedCorrelationMatrix(self)

    # Cholesky decomposition
    computeCholeskyDecomposition(self)

  def setPoint(self,point=None):
    """Set design point"""
    if point == None:
      self.point = np.zeros((self.nrv,1))
    else:
      self.point = point

  def computeRandomNumbers(self):
    """Compute random numbers"""
    if self.options.getRandomGenerator() == 0:
      self.u = np.dot(self.point,[np.ones(self.block_size)])+np.dot(self.cholesky_covariance,np.random.randn(self.nrv, self.block_size))
    elif self.options.getRandomGenerator() == 1:
      print 'Error: function not yet implemented'

  def computeTransformation(self):
    """Compute transformation from u to x space

    .. note::

       TODO: this method takes a lot of time, find some better solution.

    """
    self.x = np.zeros((self.nrv,self.block_size))

    for i in range(self.block_size):
      self.x[:,i] = u_to_x(self.u[:,i],self.model)

  def computeLimitState(self):
    """Evaluate limit-state function"""
    G, gradient = evaluateLimitState(self.x,self.model,self.options,self.limitstate,'no')
    self.G = G

  def computeResults(self):
    """Collect result of sampling"""
    self.I = np.zeros(self.block_size)
    indx=np.where(self.G[0] < 0)
    self.I[indx] = 1;

  def computeSumUpdate(self):
    """Update summation"""
    part1 = np.zeros(self.block_size)
    for i in range(self.block_size):
      part1[i] = np.vdot(self.u[:,i],self.u[:,i])
    value1 = self.u-np.dot(self.point,[np.ones(self.block_size)])
    value2 = np.dot(self.inverse_covariance,(self.u-np.dot(self.point,[np.ones(self.block_size)])))
    part2 = np.zeros(self.block_size)
    for i in range(self.block_size):
      part2[i] = np.vdot(value1[:,i],value2[:,i])

    self.q = self.I * self.factors * np.exp(-0.5*part1+0.5*part2)

    self.sum_q += np.sum(self.q)
    self.sum_q2 += np.sum(self.q**2)

  def computeCoefficientOfVariation(self):
    """Compute Coefficient of Variation"""
    n = self.k-1
    if self.sum_q > 0:
      self.q_bar[n] = 1*self.k**(-1) * self.sum_q
      variance_q_bar = 1*self.k**(-1) * ( 1*self.k**(-1) * self.sum_q2 - (1*self.k**(-1)*self.sum_q)**2 )
      self.cov_q_bar[n] = np.sqrt(variance_q_bar) * self.q_bar[n]**(-1)
      if self.cov_q_bar[n] == 0:
        self.cov_of_q_bar[n] = 1.0
    else:
      self.q_bar[n] = 0
      self.cov_q_bar[n] = 1.0

  def computePercentDone(self):
    """Compute percent done"""
    if np.floor( self.k * self.options.getSamples()**(-1)*20) > self.done:
      self.done = np.floor( self.k * self.options.getSamples()**(-1)*20)
      if self.options.printOutput():
          print self.done*5,'% complete'

  def computeFailureProbability(self):
    """Compute probability of failure"""
    if self.sum_q > 0:
      self.Pf = self.q_bar[self.k-1]
    else:
      self.Pf = 0

  def computeBeta(self):
    """Compute beta value"""
    if self.sum_q > 0:
      self.beta = -Normal.inv_cdf(self.Pf)
    else:
      self.beta = 0

  def getBeta(self):
    """Returns the beta value

    :Returns:
      - beta (float): Returns the beta value
    """
    return self.beta

  def getFailure(self):
    """Returns the probability of failure

    :Returns:
      - Pf (float): Returns the probability of failure
    """
    return self.Pf

  def getDistributionData(self):
    """Returns data for the failure

    :Returns:
      - all_G (float): Returns data from the distribution
    """
    return self.all_G

  def getBins(self):
    """Returns the amount on bins

    :Returns:
      - bins (int): Returns the amount on bins (histogram)
    """
    return self.bins


class CrudeMonteCarlo(MonteCarlo):
  """Crude Monte Carlo simulation (CMC)

  The Crude Monte Carlo simulation (CMC) is the most simple form and
  corresponds to a direct application of Equation (24). A large number
  :math:`n` of samples are simulated for the set of random variables
  :math:`{\\bf X}`. All samples that lead to a failure are counted :math:`n_f`
  and after all simulations the probability of failure :math:`p_f` may be
  estimated by [Faber2009]_

  .. math::
     :label: eq:2_93

             \\tilde{p}_f = \\frac{n_f}{n}

  Theoretically, an infinite number of simulations will provide an exact
  probability of failure. However, time and the power of computers are
  limited; therefore, a suitable amount of simulations :math:`n` are required
  to achieve an acceptable level of accuracy.

  :Attributes:
    - analysis_option (AnalysisOption): Option for the structural analysis
    - limit_state (LimitState): Information about the limit state
    - stochastic_model (StochasticModel): Information about the model
    - point (vec): Design point for the simulation
  """

  def __init__(self,analysis_options=None,limit_state=None,stochastic_model=None,point=None):

    MonteCarlo.__init__(self,analysis_options=analysis_options,limit_state=limit_state,stochastic_model=stochastic_model)

    # Set point for crude Monte Carlo / importance sampling
    self.setPoint(point)

    # Initialize variables
    self.initializeVariables()

    self.k = 0
    while self.k < self.options.getSamples():
      self.block_size = min(self.options.getBlockSize(), self.options.getSamples() - self.k)
      self.k += self.block_size

      # Computation of the random numbers
      self.computeRandomNumbers()

      # Compute transformation from u to x space
      self.computeTransformation()

      # Evaluate limit-state function
      self.computeLimitState()

      # Collect result of sampling: if g < 0 , I = 1 , else I = 0
      self.computeResults()

      # Update sums
      self.computeSumUpdate()

      # Compute coefficient of variation (of pf)
      self.computeCoefficientOfVariation()

      # Coumpute percent done
      self.computePercentDone()

      # Check convergence
      if self.cov_q_bar[self.k-1] <= self.options.getSimulationCov():
        break

    # Compute failure probability
    self.computeFailureProbability()

    # Compute beta value
    self.computeBeta()

    # Show Results
    if self.options.printOutput():
      self.showResults()

  def initializeVariables(self):
    """Initialization of the simulation variables"""
    stdv = self.options.getSimulationStdv()
    samples = self.options.getSamples()
    # Establish covariance matrix, its Cholesky decomposition, and its inverse
    self.covariance = stdv**2 * np.eye(self.nrv)
    self.cholesky_covariance = stdv * np.eye(self.nrv);  # chol_covariance = chol(covariance);
    self.inverse_covariance = 1*(stdv**2)**(-1) * np.eye(self.nrv); # inv_covariance = inv(covariance);

    # Initializations
    self.sum_q = 0
    self.sum_q2 = 0
    self.q_bar = np.zeros(samples)
    self.cov_q_bar = np.empty(samples)
    self.cov_q_bar[:] = np.nan

    # Pre-compute some factors to minimize computations inside simulation loop
    self.factors = stdv**self.nrv
    self.cov_q_bar[0] = 1.0
    self.done = 0

  def showResults(self):
    """Show results and plots"""
    print ''
    print '=================================================='
    print ''
    print ' RESULTS FROM RUNNING CRUDE MONTE CARLO SIMULATION'
    print ''
    print ' Reliability index beta:       ',self.beta
    print ' Failure probability:          ',self.Pf
    print ' Coefficient of variation of Pf',self.cov_q_bar[self.k-1]
    print ' Number of simulations:        ',self.k
    print ''
    print '=================================================='
    print ''

    npts = 200
    x = np.arange(0, self.k, self.block_size, dtype=int)
    x[0] = 1
    idx = x-1
    if self.k*self.block_size**(-1) > npts:
      s = np.round(len(x)*200**(-1))
      idx = np.arange(0, len(x), s, dtype=int)
      x = x[idx]
      idx = x-1

    # Plot how the coefficient of variation varies with number of simulations
    plt.clf()
    plt.plot(x, self.cov_q_bar[idx],'ro')
    if self.cov_q_bar[self.block_size-1]>0:
      plt.xlim([0,self.k*1.01])
      plt.ylim([0,self.cov_q_bar[self.block_size-1]*1.1])

    plt.title(r'C.o.V. of probability of failure $\delta_{pf}$')
    plt.xlabel('Number of simulations')
    plt.ylabel('Coefficient of variation')
    plt.grid(True)
    plt.show()

    # Plot how pf estimate varies with number of simulations plus confidence intervals
    plt.clf()
    plt.plot(x, self.q_bar[idx],'ro')

    plt.title(r'Probability of failure $(p_f)$')
    plt.xlabel('Number of simulations')
    plt.ylabel('Probability of failure')
    plt.grid(True)
    plt.show()

class ImportanceSampling(CrudeMonteCarlo):
  """Importance Sampling

  To decrease the number of simulations and the coefficient of variation,
  other methods can be performed. One commonly applied method is the
  Importance Sampling simulation method (IS).

  :Attributes:
    - analysis_option (AnalysisOption): Option for the structural analysis
    - limit_state (LimitState): Information about the limit state
    - stochastic_model (StochasticModel): Information about the model
  """

  def __init__(self,analysis_options=None,limit_state=None,stochastic_model=None):

    FormAnalysis = Form(analysis_options,limit_state,stochastic_model)
    u = FormAnalysis.getDesignPoint()
    u = np.transpose([u])

    CrudeMonteCarlo.__init__(self,analysis_options,limit_state,stochastic_model,u)


class DistributionAnalysis(MonteCarlo):
  """Distribution Analysis

  To analyze the random variables, used in the limit state function, a
  numerical distribution analysis based on Monte Carlo simulation can be
  performed.

  :Attributes:
    - analysis_option (AnalysisOption): Option for the structural analysis
    - limit_state (LimitState): Information about the limit state
    - stochastic_model (StochasticModel): Information about the model
  """

  def __init__(self,analysis_options=None,limit_state=None,stochastic_model=None):

    MonteCarlo.__init__(self,analysis_options,limit_state,stochastic_model)

    # Set point for crude Monte Carlo / importance sampling
    self.setPoint()

    # Initialize variables # Different
    self.initializeVariables()

    self.k = 0
    while self.k < self.options.getSamples():
      self.block_size = min(self.options.getBlockSize(), self.options.getSamples() - self.k)
      self.k += self.block_size

      # Computation of the random numbers
      self.computeRandomNumbers()

      # Comoute Transformation from u to x space
      self.computeTransformation()

      # Evaluate limit-state function and its gradient
      self.computeLimitState()

      # Collect Data
      self.coumputeDataUpdate()

      # Coumpute percent done
      self.computePercentDone()

    # Compute Distribution Data
    self.computeDistributionData()

    # Show Results # Different
    if self.options.printOutput():
      self.showResults()


  def initializeVariables(self):
    """Initialization of the simulation variables"""
    stdv = self.options.getSimulationStdv()
    samples = self.options.getSamples()
    # Establish covariance matrix, its Cholesky decomposition, and its inverse
    self.covariance = stdv**2 * np.eye(self.nrv)
    self.cholesky_covariance = stdv * np.eye(self.nrv);  # chol_covariance = chol(covariance);

    ng = 1

    self.all_X = np.zeros((self.nrv,samples))
    self.all_G = np.zeros((ng,samples))

    self.done = 0
    self.bins = getBins(self.options.getSamples())


  def coumputeDataUpdate(self):
    """Update data"""
    indx = range((self.k-self.block_size),self.k)
    self.all_X[:,indx] = self.x
    self.all_G[:,indx] = self.G

  def computeDistributionData(self):
    """Compute data for the distributions"""
    x = self.all_G
    x = np.transpose(x)
    self.all_G = x

  def showResults(self):
    """Show results and plots"""
    print ''
    print '=================================================='
    print ''
    print ' RESULTS FROM RUNNING DISTRIBUTION ANALYSIS'
    print ''
    print ' Number of simulations:        ',self.k
    print ' Approximated number of bins   ',self.bins
    print ''
    print '=================================================='
    print ''


    npts = 200
    marg = self.model.getMarginalDistributions()

    for i in range(self.nrv):
      # Plot simulatet distribution
      xi = self.all_X[i]
      minx = np.min(xi)
      maxx = np.max(xi)
      n, bins, patches = plt.hist(xi, self.bins, normed=True, facecolor='green', alpha=0.75)

      # Plot reference distribution
      xr = np.linspace(minx,maxx,npts)
      reference_pdf = pdf(xr,marg[i])

      plt.plot(xr, reference_pdf,'r-')

      name = self.model.getNames()[i]
      string = 'Distribution Analysis for '+name
      plt.title(string)
      plt.xlabel('Random Values')
      plt.ylabel('Probability Density')
      plt.grid(True)

      plt.show()

    xg = self.all_G

    n, bins, patches = plt.hist(xg, self.bins, normed=True, facecolor='green', alpha=0.75)

    plt.title('Distribution Analysis for the Limit State Function')
    plt.xlabel('Random Values')
    plt.ylabel('Probability Density')

    plt.show()
