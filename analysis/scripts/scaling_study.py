# basic imports
import matplotlib.pyplot as plt
import numpy as np
from unyt import MW, GW, km

# osier imports
from osier import CapacityExpansion
from osier.tech_library import nuclear_adv, wind, battery, natural_gas
from osier import total_cost, annual_emission

# pymoo imports
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.pcp import PCP

from multiprocessing.pool import ThreadPool
from pymoo.core.problem import StarmapParallelization


phase_shift = 0  # horizontal shift [radians]
base_shift = 2  # vertical shift [units of demand]
n_hours = 24  # hours per day
total_demand = 185  # [MWh], sets the total demand [units of energy]
rng = np.random.default_rng(seed=1234)

n_days = 2  # days to model
N = n_hours * n_days  # total number of time steps

hours = np.linspace(0, N, N)

demand = (np.sin((hours * np.pi / n_hours * 2 + phase_shift))
          * -1 + np.ones(N) * (base_shift + 1))

noise = np.random.random(N)
demand += noise

demand = demand / demand.sum() * total_demand
wind_speed = np.random.weibull(a=2.5, size=N)
