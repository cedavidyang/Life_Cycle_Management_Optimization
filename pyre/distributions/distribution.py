#!/usr/bin/python -tt
# -*- coding: utf-8 -*-

import numpy as np
import math
import scipy.optimize as opt
import scipy.special as spec
import scipy.stats as stats

def getDistributionType(type):
  listOfDistributions = {1:'Normal',
                         2:'Lognormal',
                         3:'Gamma',
                         4:'ShiftedExponential',
                         5:'ShiftedRayleigh',
                         6:'Uniform',
                         7:'Beta',
                         8:'ChiSquare',
                         11:'TypeIlargestValue',
                         12:'TypeIsmallestValue',
                         13:'TypeIIlargestValue',
                         14:'TypeIIIsmallestValue',
                         15:'Gumbel',
                         16:'Weibull'}
  return listOfDistributions[type]


class Distribution(object):
  """Probability distribution

  Attributes:
    name (str):   Name of the random variable\n
    type (int):   Type of probability distribution\n
                     1:'Normal'\n
                     2:'Lognormal'\n
                     3:'Gamma'\n
                     4:'ShiftedExponential'\n
                     5:'ShiftedRayleigh'\n
                     6:'Uniform'\n
                     7:'Beta'\n
                     8:'ChiSquare'\n
                    11:'TypeIlargestValue'\n
                    12:'TypeIsmallestValue'\n
                    13:'TypeIIlargestValue'\n
                    14:'TypeIIIsmallestValue'\n
                    15:'Gumbel'\n
                    16:'Weibull'\n
    mean (float): Mean or other variable\n
    stdv (float): Standard deviation or other variable\n
    startpoint (float): Start point for seach\n
    p1 (float):   Parameter for the distribution\n
    p2 (float):   Parameter for the distribution\n
    p3 (float):   Parameter for the distribution\n
    p4 (float):   Parameter for the distribution\n
    input_type (any): Change meaning of mean and stdv\n

    Default: all values

  """

  def __init__(self, name=None,type=None,mean=None,stdv=None,startpoint=None,p1=None,p2=None,p3=None,p4=None,input_type=None):
    self.name = name
    self.type = type
    self.mean = mean
    self.stdv = stdv

    self.startpoint = None
    self.setStartPoint(startpoint)

    self.p1 = p1
    self.p2 = p2
    self.p3 = p3
    self.p4 = p4
    self.input_type = input_type
    self.rv = stats.rv_continuous()

  def __repr__(self):
    type = getDistributionType(self.type)
    string = self.name+': '+type+' distribution'
    return string#repr(self.matrix)

  def setName(self, name):
    self.name = name
  
  def getName(self):
    return self.name

  def getType(self):
    return self.type

  def getMean(self):
    return self.mean

  def getStdv(self):
    return self.stdv

  def getStartPoint(self):
    return self.startpoint

  def getP1(self):
    return self.p1

  def getP2(self):
    return self.p2

  def getP3(self):
    return self.p3

  def getP4(self):
    return self.p4

  def getMarginalDistribution(self):
    marg = MarginalDistribution(self.type,self.mean,self.stdv,self.startpoint,self.p1,self.p2,self.p3,self.p4)
    return marg

  def setStartPoint(self,startpoint=None):
    if startpoint==None:
      self.startpoint = self.mean
    else:
      self.startpoint = startpoint
      
  def getRV(self):
    return self.rv
    
  def discretizeRV(self, n_bins=10, retain_percent=0.999, bin_boundaries=None):
    """Discretize retain_percent percentage of continuous distribution into n_bins bins
        
    :Args:
        
        - n_bins (int): number of bins
        - retain_percentage (float, 0<retain_percentage<1): percentage of continuous values considered in the discretization
        
    :Returns:
        
        - bin_boundaries (float list/array): [[lower0, upper0], [lower1, upper1], ..., [lowerN, upperN]], upper(I) = lower(I+1)
        - bin_represent_values (float list/array): [represent_value0, ..., represent_valueN]
        - bin_probabilities (float list/array): [prob0, prob1, ..., probN]
        
    """
    
    if bin_boundaries==None:
        p_bound = 1 - retain_percent
            
        lower_limit = self.rv.ppf(p_bound/2.)
        upper_limit = self.rv.ppf(1. - p_bound/2.)
        bin_edges = np.linspace(lower_limit, upper_limit, num=n_bins+1)
        bin_boundaries = np.array( (bin_edges[:-1], bin_edges[1:]) ).transpose()
        bin_represent_values = ( bin_edges[:-1] + bin_edges[1:] ) / 2.
        bin_probabilities = self.rv.cdf(bin_edges[1:]) - self.rv.cdf(bin_edges[:-1])
        bin_probabilities = bin_probabilities / sum(bin_probabilities)
        
    else:
        bin_represent_values = (bin_boundaries[:,0] + bin_boundaries[:,1]) / 2.
        bin_probabilities = self.rv.cdf(bin_boundaries[:,1]) - self.rv.cdf(bin_boundaries[:,0])
        bin_probabilities = bin_probabilities / sum(bin_probabilities)
    
    return bin_boundaries, bin_represent_values, bin_probabilities


class MarginalDistribution(object):
  """
  """

  def __init__(self, type,mean,stdv,startpoint,p1,p2,p3,p4):
    """

    Arguments:
      - `type`:
      - `mean`:
      - `stdv`:
      - `startpoint`:
      - `p1`:
      - `p2`:
      - `p3`:
      - `p4`:
    """
    self.type = type
    self.mean = mean
    self.stdv = stdv
    self.startpoint = startpoint
    self.p1 = p1
    self.p2 = p2
    self.p3 = p3
    self.p4 = p4
    self.marg = [type,mean,stdv,startpoint,p1,p2,p3,p4]

  def __repr__(self):
        return repr(self.marg)

  def __getitem__(self, key):
    return self.marg[key]

  def __setitem__(self, key, item):
    self.marg[key] = item

  def __len__(self):
    return 1

  def getMarg(self):
    marg = [self.type,self.mean,self.stdv,self.startpoint,self.p1,self.p2,self.p3,self.p4]
    return marg

  def getType(self):
    return self.type

  def getMean(self):
    return self.mean

  def getStdv(self):
    return self.stdv

  def getStartPoint(self):
    return self.startpoint

  def getP1(self):
    return self.p1

  def getP2(self):
    return self.p2

  def getP3(self):
    return self.p3

  def getP4(self):
    return self.p4
