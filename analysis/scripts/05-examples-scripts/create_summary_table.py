import pandas as pd
import numpy as np
import camelot
import yaml

def is_pareto_efficient_dumb(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        is_efficient[i] = np.all(np.any(costs[:i]>c, axis=1)) and np.all(np.any(costs[i+1:]>c, axis=1))
    return is_efficient


if __name__ == "__main__":
    mapping = {'once-through':[f'EG0{i}' for i in range(1,9)],
                'limited-recycle':[f'EG0{i}' if i < 10 else f'EG{i}' for i in range(9,19)],
                'continuous-recycle':[f'EG{i}' for i in range(19,41)]
                }
    summary_tables = camelot.read_pdf(snakemake.input[0], pages="19")
    full_df = pd.read_csv(snakemake.input.metric_data, index_col=0)
    full_df_normed = full_df.div(full_df.max(axis=0), axis=1)
    pareto_optimal_df = full_df.iloc[is_pareto_efficient_dumb(full_df_normed.values),:]

    reactors_df = summary_tables[0].df

    reactors_df.columns = reactors_df.iloc[0,:].values
    reactors_df = reactors_df.iloc[1:,:]


    summary_df = pd.DataFrame(index=[f'EG0{i}' 
                                     if i < 10 
                                     else f'EG{i}' 
                                     for i in range(1,41)]).assign(fuel_cycle_type='', 
                                                                   reactor_type='',
                                                                   set_conclusion='Not promising',
                                                                   pareto_optimal='False')
    
    for fc_type, eg_list in mapping.items():
        for eg in eg_list:
            summary_df.at[eg, 'fuel_cycle_type'] = fc_type

    for reactor in reactors_df['Operations'].values:
        eg_series = reactors_df.loc[reactors_df['Operations']==reactor, 'Used by the Analysis Example for  the Evaluation Group']
        for eg in eg_series:
            eg_list = [group.replace('\n', '').replace(' ','') for group in eg.split(',')]
            for group in eg_list:
                summary_df.at[group, 'reactor_type'] = reactor

    for status, eg_list in snakemake.config['promising_fuelcycles'].items():
        for eg in eg_list:
            summary_df.at[eg, 'set_conclusion'] = " ".join(status.split('_')).capitalize()

    for eg in pareto_optimal_df.index:
        summary_df.at[eg,'pareto_optimal'] = 'True'


    summary_df.to_csv(snakemake.output.summary_data)
    

    summary_df.columns = ["Fuel Cycle Type",
                          "Reactor Type",
                          "EST Conclusion",
                          "Pareto Optimal"]
    
    summary_df.index.name = "EG"
    
    summary_df.style.to_latex(snakemake.output.summary_table, hrules=True)

    