import numpy as np

class DummyObject(object):
    def __init__(self, value):
        self.value = value
    def rvs(self, size = None):
        return self.value * np.ones(size)
    def ppf(self, cdf_array):
        return self.value * np.ones(cdf_array.shape)

class Deterministic(object):
    def __init__(self, name, mu, sigma=None):
        self.name = name
        self.value = mu
        self.rv = DummyObject(mu)
