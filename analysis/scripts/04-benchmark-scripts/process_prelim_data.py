import pandas as pd
import dill
import time
from osier import n_mga, get_objective_names, get_tech_names

if __name__ == "__main__":

    print(f"Opening {snakemake.input.four_obj_results}")

    start = time.perf_counter()
    with open(snakemake.input.four_obj_results, 'rb') as file:
        results = dill.load(file)

    end = time.perf_counter()

    print(f"Dataset took {(end-start):.2f} seconds to open")

    # save tech names
    techs = get_tech_names(results.problem.technology_list)
    
    # save objective names
    obj_cols = get_objective_names(results)

    # Pareto front
    F = results.F
    slack = 0.1

    peak_demand = results.problem.max_demand.to_value()

    # N_MGA dataframe
    mdf = n_mga(results, wide_form=True, how='all')
    mdf['optimal'] = 'Sub-opt'
    mdf.iloc[:,4:-1] = mdf.iloc[:,4:-1]*peak_demand

    # Pareto dataframe
    pf_obj = pd.DataFrame(dict(zip(obj_cols, F.T)))
    ddf = pd.DataFrame(dict(zip(techs,results.X.T)))*peak_demand
    pf_obj = pd.concat([pf_obj, ddf], axis=1)
    pf_obj['optimal'] = 'Optimal'

    combined = pd.concat([pf_obj, mdf], axis=0).reset_index(drop=True)
    combined.to_csv(snakemake.output.combined_data, index=False)
