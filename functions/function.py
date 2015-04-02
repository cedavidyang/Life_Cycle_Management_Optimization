#!/usr/bin/python -tt
# -*- coding: utf-8 -*-

import numpy as np
import math
import scipy.special as spec
#from constants import *

def weightedAvgAndStd(values, weights, axis_data=None, axis_weight=None, unbias=None):
    """weighted mean and std
       unbiased version of Std is based on en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    """
    avg = np.average(values, axis=axis_data, weights=weights)
    try:
        tmp_variance = np.average((values-avg)**2, axis=axis_data, weights=weights)
    except ValueError:
        tmp_variance = np.average(((values.T-avg)**2).T, axis=axis_data, weights=weights)
    if unbias==False:
        std = np.sqrt(tmp_variance)
    else:
        v1 = np.sum(weights, axis=axis_weight)
        v2 = np.sum(weights**2, axis=axis_weight)
        weighted_sum = tmp_variance * v1
        variance = weighted_sum / (v1 - (v2/v1))
        std = np.sqrt(variance)
    return avg, std


def weightedSum(values, weights, axis_data, axis_weight=None):
    """Return four sums that are useful to get weighted mean and std
    """
    if axis_weight == None:
        axis_weight = 0
    
    sum1 = np.tensordot(values, weights, axes=(axis_data, axis_weight))
    sum2 = np.tensordot(values**2, weights, axes=(axis_data, axis_weight))
    weight_sum1 = np.sum(weights, axis=axis_weight)
    weight_sum2 = np.sum(weights**2, axis=axis_weight)
    
    return sum1, sum2, weight_sum1, weight_sum2


def weightedAvgAndStdFromSum(sums, unbias=None):
    sum1 = sums[0]
    sum2 = sums[1]
    weight_sum1 = sums[2]
    weight_sum2 = sums[3]
    avg = sum1 / weight_sum1
    if unbias == False:
        var = (sum2 - 2*avg*sum1 + avg**2*weight_sum1) / (weight_sum1)
    else:
        var = (sum2 - 2*avg*sum1 + avg**2*weight_sum1) / (weight_sum1-weight_sum2/weight_sum1)
    std = np.sqrt(var)
    
    return avg, std
    
