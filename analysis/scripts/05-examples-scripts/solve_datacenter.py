# basic imports
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import dill as pickle
from glob import glob
import time

# pymoo imports
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.indicators.hv import HV
from pymoo.termination.max_gen import MaximumGenerationTermination



mpl.use("pgf") 
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True 
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.shadow'] = False

solver = "appsi_highs"

print(f"Solver set: {solver}")

if __name__ == "__main__":
    
    with open(snakemake.input.dc_problem, "rb") as file:
        problem = pickle.load(file)
    
    use_checkpoints = False

    checkpoint_list = glob("checkpoint_*.pkl")
    checkpoint_list.sort()
    if (len(checkpoint_list) > 0) and use_checkpoints:
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
        with open(snakemake.output.dc_results_F, "wb") as file:
            pickle.dump(res.F, file)

        with open(snakemake.output.dc_results_X, "wb") as file:
            pickle.dump(res.X, file)

    except KeyboardInterrupt:
        # save checkpoint on early termination
        timestr = time.strftime("%Y%m%d-%H%M%S")
        checkpoint_name = f"checkpoint_{timestr}.pkl"
        current_results_name = f"results_{timestr}.pkl"
        print(f"Simulation stopped. Saving checkpoint to {checkpoint_name}")
        with open(checkpoint_name, "wb") as f:
            pickle.dump(algorithm, f)

    ref_point = np.max(np.array([opt.F for opt in res.history[0].opt]), axis=0)

    ind = HV(ref_point=ref_point)

    performance_list = []
    for pop in res.history:
        pop_pf = np.array([e.F for e in pop.opt])
        performance_list.append(ind(pop_pf))

    n_evals = np.array([e.evaluator.n_eval for e in res.history])

    plt.title("Convergence", fontsize=18)
    plt.plot(n_evals, performance_list, "--")
    plt.ylabel("Hypervolume", fontsize=18)
    plt.xlabel("Evaluations", fontsize=18)
    plt.ylim(min(performance_list), max(performance_list)*1.005)
    plt.xlim(0, max(n_evals))
    plt.grid()
    plt.tight_layout()
    plt.savefig(snakemake.output.convergence)
        
    