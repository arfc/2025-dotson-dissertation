import numpy as np
import matplotlib.pyplot as plt
import dill as pickle
from pymoo.visualization.pcp import PCP
import matplotlib as mpl
from mycolorpy import colorlist as mcp
from osier import get_tech_names    

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

if __name__ == "__main__":
    
    with open(snakemake.input.dc_results_X, "rb") as f:
        res_X = pickle.load(f)

    with open(snakemake.input.dc_results_F, "rb") as f:
        res_F = pickle.load(f)
    
    with open(snakemake.input.tech_list, "rb") as f:
        tech_list = pickle.load(f)

    # objective space

    obj_labels=['Total Cost', 'CO2eq', 'EROI']
    obj_markers=['s','o','x']
    obj_colors=['tab:red','tab:green','tab:purple']
    plot = PCP(title=("Objective Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=obj_labels,
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True
    plot.add(res_F, color="grey", alpha=0.3)

    for i, label in enumerate(obj_labels):
        min_idx = np.where(res_F[:,i]==min(res_F[:,i]))[0][0]
        plot.add(res_F[min_idx], 
                 linewidth=2, 
                 label=f"min({label})", 
                 marker=obj_markers[i], 
                 markersize=10,
                 color=obj_colors[i])

    plot.save(snakemake.output.dc_objective_space)


    # design space
    # obj_labels=['Total Cost', 'CO2eq', 'EROI']
    tech_labels = get_tech_names(tech_list)
    tech_labels = [label.replace("_","") for label in tech_labels]
    plot = PCP(title=("Design Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=tech_labels,
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True
    plot.add(res_X, color="grey", alpha=0.3)

    for i, label in enumerate(obj_labels):
        min_idx = np.where(res_F[:,i]==min(res_F[:,i]))[0][0]
        # if i == 2:
        #     continue
        plot.add(res_X[min_idx], linewidth=2, label=f"min({label})", marker=obj_markers[i], markersize=10, color=obj_colors[i])
        # plot.add(res_F[3], linewidth=5, color="tab:green", label=r"Least CO$_2$")
        # plot.add(res_F[6], linewidth=5, color="tab:blue", label="Least Cost")
    plot.save(snakemake.output.dc_design_space)


    # top 10 solutions for each objective - objective
    obj_labels=['Total Cost', 'CO2eq', 'EROI']

    obj_colors = {
                'Total Cost':'Reds_r',
                'CO2eq':'Greens_r',
                'EROI':'Purples_r',
                }
    markers = ['s','o','x']
    n=10

    plot = PCP(title=(f"Objective Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(False, {'loc': "upper left"}),
            labels=obj_labels,
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True
    plot.add(res_F, color="grey", alpha=0.3)
    for i, (obj, colors) in enumerate(obj_colors.items()):
        # get the n_smallest values
        F_i = res_F[:,i]
        idx_smallest = np.argpartition(F_i, n)[:n]

        # generate color map
        cmap = colors
        color=mcp.gen_color(cmap=cmap,n=n+1)

        F_smallest = res_F[idx_smallest]
        for j, f in enumerate(F_smallest):
            plot.add(f, linewidth=2, color=color[j], marker=markers[i])

    plot.save(snakemake.output.dc_objective_space_10)


    # top 10 solutions for each objective - design
    obj_labels=['Total Cost', 'CO2eq', 'EROI']
    plot = PCP(title=(f"Objective Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(False, {'loc': "upper left"}),
            labels=tech_labels,
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True
    plot.add(res_X, color="grey", alpha=0.3)
    for i, (obj, colors) in enumerate(obj_colors.items()):
        # get the n_smallest values
        F_i = res_F[:,i]
        idx_smallest = np.argpartition(F_i, n)[:n]

        # generate color map
        cmap = colors
        color=mcp.gen_color(cmap=cmap,n=n+1)

        X_i = res_X[idx_smallest]
        for j, x in enumerate(X_i):
            plot.add(x, linewidth=2, color=color[j], marker=markers[i])

    plot.save(snakemake.output.dc_design_space_10)