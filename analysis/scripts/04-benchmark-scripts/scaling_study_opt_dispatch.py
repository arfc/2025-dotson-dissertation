# basic imports
import matplotlib.pyplot as plt
import numpy as np
from unyt import MW, GW, km, hour
import time
import pandas as pd


# osier imports
from osier import CapacityExpansion
import osier.tech_library as lib # nuclear_adv, wind, battery, natural_gas
from osier import total_cost, annual_emission
from osier.utils import synchronize_units

# pymoo imports
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.pcp import PCP

from multiprocessing.pool import ThreadPool
from pymoo.core.problem import StarmapParallelization
import dill as pickle


if __name__ == "__main__":
    phase_shift = 0  # horizontal shift [radians]
    base_shift = 2  # vertical shift [units of demand]
    day_hours = 24  # hours per day
    year_hours = 8760  # hours per year
    total_demand = 185  # [MWh], sets the total demand [units of energy]
    rng = np.random.default_rng(seed=1234)
    squish = 0.5

    tech_list = [lib.wind, lib.natural_gas, lib.battery, lib.nuclear_adv]

    # synchronize the units
    synched_techs = synchronize_units(tech_list, MW, hour)


    n_days_list = [1, 5, 10, 25, 50, 100, 200, 300, 365]
    # n_days_list = [100, 200, 300, 365]
    # n_days_list = [1, 5, 10, 25]
    # n_days_list = [10]
    optimizers = ['optimal','logical']
    # optimizers = ['logical']
    # optimizers = ['optimal']
    time_data = {'n_days':n_days_list,
                 'solve_time_logical':[],
                 'solve_time_optimal':[]
                 }
    population_size = 20
    n_generations = 10
    # populate data
    try:
        for i, n_days in enumerate(n_days_list):  
            for optimizer in optimizers:  
                # n_days = 1  # days to model
                N = day_hours * n_days  # total number of time steps

                hours = np.linspace(0, N, N)
                y = np.sin((hours * np.pi / year_hours + phase_shift)) + np.ones(len(hours)) * (base_shift + 1) * squish

                demand = (np.sin((hours * np.pi / day_hours * 2 + phase_shift))
                        * -1 + np.ones(N) * (base_shift + 1)) * 0.25 + y

                noise = rng.normal(size=N)*0.05
                demand += noise

                demand = demand / demand.max()
                wind_speed = rng.weibull(a=2, size=N)
                wind_speed /= wind_speed.max()

                # run the model
                print(f"======== Starting {optimizer.capitalize()} Dispatch ========")
                print(f"===================== n_days: {n_days} =====================")
                problem = CapacityExpansion(technology_list = synched_techs,
                                        demand=demand*MW,
                                        wind=wind_speed,
                                        upper_bound= 1 / lib.wind.capacity_credit,
                                        objectives = [total_cost, annual_emission],
                                        # model_engine='logical',
                                        model_engine=optimizer,
                                        solver='appsi_highs',
                                        )

                algorithm = NSGA2(pop_size=population_size)
                res = minimize(problem, 
                                algorithm, 
                                termination=("n_gen", n_generations), 
                                seed=1, 
                                verbose=True, 
                                save_history=False)
                time_data[f'solve_time_{optimizer}'].append(res.exec_time)

        df = pd.DataFrame(time_data)
        df.to_csv(snakemake.output.algorithm_times)
        with open("time_data_dict.pkl", "wb") as f:
                pickle.dump(time_data, f)
    except KeyboardInterrupt:
        print("Interrupting simulation. Saving time data.")
        with open("time_data_dict.pkl", "wb") as f:
                pickle.dump(time_data, f)

