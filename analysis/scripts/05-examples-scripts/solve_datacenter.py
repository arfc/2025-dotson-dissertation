# basic imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import dill as pickle
import os
from glob import glob
import time

# pymoo imports
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.termination.max_gen import MaximumGenerationTermination


solver = "appsi_highs"

print(f"Solver set: {solver}")

if __name__ == "__main__":
    
    with open(snakemake.input.dc_problem, "rb") as file:
        problem = pickle.load(file)
    

    checkpoint_list = glob("checkpoint_*.pkl")
    checkpoint_list.sort()
    if len(checkpoint_list) > 0:
        with open(checkpoint_list[-1], 'rb') as f:
            algorithm = pickle.load(f)
            algorithm.termination = MaximumGenerationTermination(200)
            print(f"Loaded {checkpoint_list[-1]}:", algorithm)
    else:
        print("No checkpoints found. Starting new run.")
        algorithm = NSGA2(pop_size=100)
        algorithm.termination = MaximumGenerationTermination(200)
        
    try:
        res = minimize(problem,
                        algorithm,
                        seed=1,
                        copy_algorithm=False,
                        save_history=True,
                        verbose=True)
        # save results
        print('Saving simulation results')
        with open(snakemake.output.dc_results, "wb") as file:
            pickle.dump(res, file)
        # save final checkpoint
        timestr = time.strftime("%Y%m%d-%H%M%S")
        checkpoint_name = f"checkpoint_{timestr}.pkl"
        print(f"Saving final simulation checkpoint to {checkpoint_name}")
        with open(checkpoint_name, "wb") as f:
            pickle.dump(algorithm, f)
    except KeyboardInterrupt:
        # save checkpoint on early termination
        timestr = time.strftime("%Y%m%d-%H%M%S")
        checkpoint_name = f"checkpoint_{timestr}.pkl"
        current_results_name = f"results_{timestr}.pkl"
        print(f"Simulation stopped. Saving checkpoint to {checkpoint_name}")
        with open(checkpoint_name, "wb") as f:
            pickle.dump(algorithm, f)

        
    