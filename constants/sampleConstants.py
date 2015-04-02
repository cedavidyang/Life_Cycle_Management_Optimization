""" MC simulation parameters
# n_ele: number of elements
# dt: intervals of points in time
# cm: load distribution factors, cm = [c1, c2...., cm]', should be
#     consistent with sim_type
# nsub: number of samples in primary adaptive sampling
# nadp: number of adaptations
# nk: number of samples to determine kopt
# n_main: number of samples in main sampling
"""

#n_ele = 1
#
#dt = 5
#nsub = 1000
#nadp = 10
#nk = 10
#n_main = 2000

# CE-based importance sampling
#NUM_COMPONENT = 8
#NUM_ADAPTATION = 20
#NUM_PRE_SMP = 1000
#NUM_KOPT_SMP = 200
#NUM_MAIN_SMP = 400

NUM_COMPONENT = 4
NUM_ADAPTATION = 20
NUM_PRE_SMP = 800
NUM_KOPT_SMP = 100
NUM_MAIN_SMP = 200

#NUM_COMPONENT = 1
#NUM_ADAPTATION = 20
#NUM_PRE_SMP = 200
#NUM_KOPT_SMP = 100
#NUM_MAIN_SMP = 200

STD_RATIO1 = 3.0
STD_RATIO2 = 1.5
STD_RATIO3 = 1.2 

INIT_K = 2.0
K_STEPS = 0.05
KRATIO_CRITERIA = 0.5

MAX_NUM_UPDATE = 10

# simulated annealing constants
MIN_ANNEAL_TEMP = 1e-4
MIN_ANNEAL_FUNC = 0.01
ANNEAL_RATE = 0.95
MAX_ITER = 5e6
NUM_DWELL = 10

# time-variant cov
COVT0_COV = 0.1