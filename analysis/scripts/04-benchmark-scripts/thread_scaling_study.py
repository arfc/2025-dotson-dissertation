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


    population_list = [25, 50, 100, 200]
    n_threads_list = [1, 2, 4, 8, 12]
    time_data = {'population':[],
                 'n_threads':[],
                 'solve_time':[]
                 }
    n_days = 5
    n_generations = 10
    # populate data
    try:
        for i, population_size in enumerate(population_list):  
            for j, n_threads in enumerate(n_threads_list):  
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
                print(f"===================== n_threads: {n_threads} =====================")
                print(f"===================== pop_size: {population_size} =====================")
                pool = ThreadPool(n_threads)
                runner = StarmapParallelization(pool.starmap)
                algorithm = NSGA2(pop_size=population_size)
                problem = CapacityExpansion(technology_list = [lib.wind, lib.natural_gas, lib.battery],
                                            demand=demand*MW,
                                            wind=wind_speed,
                                            upper_bound= 1 / lib.wind.capacity_credit,
                                            objectives = [total_cost, annual_emission],
                                            model_engine='logical',
                                            solver='appsi_highs',
                                            elementwise_runner=runner)
                
                res = minimize(problem, 
                                algorithm, 
                                termination=("n_gen", n_generations), 
                                seed=1, 
                                verbose=True, 
                                save_history=False)
                time_data[f'solve_time'].append(res.exec_time)
                time_data['n_threads'].append(n_threads)
                time_data['population'].append(population_size)
                pool.close()

        df = pd.DataFrame(time_data)
        df.to_csv(snakemake.output.thread_times)
        with open("thread_data_dict.pkl", "wb") as f:
                pickle.dump(time_data, f)
    except KeyboardInterrupt:
        print("Interrupting simulation. Saving time data.")
        with open("thread_data_dict.pkl", "wb") as f:
                pickle.dump(time_data, f)

