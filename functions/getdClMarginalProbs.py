# getMarginalProbabilities of dCL

from ctypes import CDLL, POINTER, c_double
from numpy import ctypeslib

testLib = CDLL('/Users/amadeus/Desktop/testCppXcodeLib/libtestCppXcodeLib.dylib')

_get_dCl_mean_stdv = testLib._Z17_get_dCl_mean_stdv # use nm -g [-m] to check the real name
_get_dCl_mean_stdv.restype = POINTER(c_double)
res = _get_dCl_mean_stdv()
arr = ctypeslib.as_array(res, (2,))