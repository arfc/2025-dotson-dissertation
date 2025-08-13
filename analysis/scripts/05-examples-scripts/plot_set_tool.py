import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob
from pymoo.visualization.pcp import PCP
from pymoo.mcdm.high_tradeoff import HighTradeoffPoints

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

    color_map = {'once-through':'tab:red',
             'limited-recycle':'tab:green',
             'continuous-recycle':'tab:blue'}


    full_df = pd.read_csv(snakemake.input.metric_data, index_col=0)
    full_df.columns = [" ".join(c.split(" ")[:-1]) for c in list(full_df.columns)]


    full_df_normed = full_df.div(full_df.max(axis=0), axis=1)
    pareto_optimal_df = full_df.iloc[is_pareto_efficient_dumb(full_df_normed.values),:]

    # full set plot
    plot = PCP(title=("Objective Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=list(pareto_optimal_df.columns),
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True


    for k, v in mapping.items():
        pareto_optimal_set = [eg for eg in mapping[k] if eg in pareto_optimal_df.index]
        subset = pareto_optimal_df.loc[pareto_optimal_set,:].values
        for i, row in enumerate(subset):
            if i == 0:
                plot.add(row, 
                        alpha=1,
                        color=color_map[k],
                        label=k)
            else:
                plot.add(row, 
                        alpha=1,
                        color=color_map[k]
                        )
    plot.do()

    for label in plot.ax.get_xticklabels():
        label.set_rotation(90)

    plot.save(snakemake.output.full_set_plot)


    # once through set plot
    plot = PCP(title=("Objective Space\n Once-Through Fuel Cycles", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=list(pareto_optimal_df.columns),
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True


    for k, v in mapping.items():
        pareto_optimal_set = [eg for eg in mapping[k] if eg in pareto_optimal_df.index]
        subset = pareto_optimal_df.loc[pareto_optimal_set,:].values
        for i, row in enumerate(subset):
            if (i == 0) and k == 'once-through':
                plot.add(row, 
                        alpha=1,
                        color=color_map[k],
                        label=k,
                        zorder=50)
            elif k == 'once-through':
                plot.add(row, 
                        alpha=1,
                        color=color_map[k],
                        #  label=k
                        zorder=50
                        )
            else:
                plot.add(row, 
                        alpha=0.25,
                        color='grey'
                        )
    plot.do()

    for label in plot.ax.get_xticklabels():
        label.set_rotation(90)

    plt.legend().set_zorder(60)
    plot.save(snakemake.output.once_through_set_plot)

    # "knee" solution (high tradeoff point)
    dm = HighTradeoffPoints()
    I = dm(pareto_optimal_df.values)

    eg_highlight = pareto_optimal_df.index[I[0]]  # EG04
    plot = PCP(title=(eg_highlight, {'pad': 30, 'fontsize':20}),
                n_ticks=10,
                legend=(True, {'loc': "upper left"}),
                labels=list(pareto_optimal_df.columns),
                figsize=(13,6),
                )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True


    for k, v in mapping.items():
        pareto_optimal_set = [eg for eg in mapping[k] if eg in pareto_optimal_df.index]
        subset = pareto_optimal_df.loc[pareto_optimal_set,:].values
        for i, row in enumerate(subset):
            if v[i] == eg_highlight:
                plot.add(row, 
                        alpha=1,
                        color=color_map[k],
                        label=v[i],
                        zorder=50)
            else:
                plot.add(row, 
                        alpha=0.25,
                        color='grey'
                        # label=v[i]
                        )
    plot.do()

    for label in plot.ax.get_xticklabels():
        label.set_rotation(90)

    plot.save(snakemake.output.single_eg_set_plot)


    # non pareto optimal plot
    plot = PCP(title=("Objective Space\n Non-Pareto Optimal", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=list(full_df.columns),
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True


    for k, v in mapping.items():
        pareto_optimal_set = [eg for eg in mapping[k] if eg not in pareto_optimal_df.index]
        subset = full_df.loc[pareto_optimal_set,:].values
        for i, row in enumerate(subset):
            plot.add(row, 
                    alpha=1,
                    color=color_map[k],
                    label=pareto_optimal_set[i])
    plot.do()

    for label in plot.ax.get_xticklabels():
        label.set_rotation(90)

    plot.save(snakemake.output.non_optimal_set_plot)