import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import dill
import time
from pymoo.visualization.pcp import PCP
import matplotlib as mpl
from mycolorpy import colorlist as mcp
from tqdm import tqdm

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

from osier import n_mga, get_objective_names, get_tech_names


if __name__ == "__main__":
    print(f"Opening {snakemake.input.combined_data}")
    combined = pd.read_csv(snakemake.input.combined_data)

    obj_labels=['Total Cost', 'CO2eq', 'Land Use', 'Percent \nNonrenewable']
    obj_nice_labels = {'lifecycle_co2_rate':r'Least CO$_2$',
                   'total_cost':'Least Cost',
                   'land_use':'Least Land Use',
                   'percent_nonrenewable':'Highest Renewable'}
    obj_tab_colors = {'lifecycle_co2_rate':'tab:green',
                'total_cost':'tab:red',
                'land_use':'tab:orange',
                'percent_nonrenewable':'tab:blue'}
    opt_df = combined.loc[combined.optimal == "Optimal"]

    ### Objective Space Plot ###
    print("Plotting objective space")
    plot = PCP(title=("Objective Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=obj_labels,
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True

    plot.add(combined[combined.optimal=='Optimal'].iloc[:,:4].values, 
             color="grey", 
             alpha=0.3)
    
    for obj, label in obj_nice_labels.items():
        plot.add(
            opt_df.loc[opt_df[obj] == opt_df[obj].min(), 
                    :'percent_nonrenewable'].values,
            linewidth=5,
            color=obj_tab_colors[obj],
            label=label
        )
    plot.save(snakemake.output.four_obj_objective_space)


    ### Objective Space MGA Plot ###
    print("Plotting objective space with MGA")
    plot = PCP(title=("Objective Space (with MGA Solutions)", 
                      {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=obj_labels,
            figsize=(13,6),
            )

    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True

    plot.add(combined.iloc[:,:4].values, color="grey", alpha=0.3)

    for obj, label in obj_nice_labels.items():
        plot.add(
            opt_df.loc[opt_df[obj] == opt_df[obj].min(), 
                    :'percent_nonrenewable'].values,
            linewidth=5,
            color=obj_tab_colors[obj],
            label=label
        )
    plot.save(snakemake.output.four_obj_objective_mga)


    ### Objective Space Top 5 for Each Objective ###
    n = 5
    print(f"Plotting top {n} performances for each objective...")
    obj_colors = {'lifecycle_co2_rate':'Greens',
                'total_cost':'Reds_r',
                'land_use':'Oranges_r',
                'percent_nonrenewable':'Blues_r'}
    obj_name = {'lifecycle_co2_rate':'CO2eq',
                'total_cost':'Total Cost',
                'land_use':'Land Use',
                'percent_nonrenewable':'Percent Nonrenewable'}
    for objective in tqdm(list(obj_name.keys())):
        plot = PCP(title=(f"Objective Space \n Highlighted: {obj_name[objective]}", 
                          {'pad': 30, 'fontsize':20}),
                n_ticks=10,
                legend=(False, {'loc': "upper left"}),
                labels=obj_labels,
                figsize=(13,6),
                )

        plot.set_axis_style(color="grey", alpha=0.5)
        plot.tight_layout = True

        plot.add(combined.iloc[:,:4].values, 
                 color="grey", 
                 alpha=0.3)
        cmap = obj_colors[objective]
        color=mcp.gen_color(cmap=cmap,n=n+1)

        subset_df = combined.nsmallest(n=n,columns=[objective]).reset_index(drop=True)

        for i in range(n):
            plot.add(subset_df.iloc[i,:4].values, linewidth=3, color=color[i])
        plot.save(f"../docs/figures/04_benchmark_chapter/4_obj_objective_space_{objective}.pgf")


    ### Design Space Plot ###
    print('Plotting design space')
    tech_labels=['Battery',
    'Biomass',
    'Coal\nConventional',
    'Coal\nAdvanced',
    'Natural Gas \nConventional',
    'Natural Gas \nAdvanced',
    'Nuclear',
    'Nuclear\nAdvanced',
    'SolarPanel',
    'WindTurbine']

    plot = PCP(title=("Design Space", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=tech_labels,
            figsize=(13,6),
            )
    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True

    plot.add(opt_df.iloc[:,4:-1].values, 
                color="grey", 
                alpha=0.3)

    for obj, label in obj_nice_labels.items():
            plot.add(
                opt_df.loc[opt_df[obj] == opt_df[obj].min()].iloc[:,4:-1].values,
                linewidth=5,
                color=obj_tab_colors[obj],
                label=label
            )

    plot.save(snakemake.output.four_obj_design_space)

    ### Design Space MGA ###
    print('Plotting design space with MGA')
    plot = PCP(title=("Design Space (with MGA Solutions)", {'pad': 30, 'fontsize':20}),
            n_ticks=10,
            legend=(True, {'loc': "upper left"}),
            labels=tech_labels,
            figsize=(13,6),
            )
    plot.set_axis_style(color="grey", alpha=0.5)
    plot.tight_layout = True

    plot.add(combined.iloc[:,4:-1].values, 
                color="grey", 
                alpha=0.3)

    for obj, label in obj_nice_labels.items():
            plot.add(
                opt_df.loc[opt_df[obj] == opt_df[obj].min()].iloc[:,4:-1].values,
                linewidth=5,
                color=obj_tab_colors[obj],
                label=label
            )

    plot.save(snakemake.output.four_obj_design_mga)


    ### Design Space Top 5 per Objective ###
    print(f"Plotting top {n} designs for each objective... ")
    for objective in tqdm(list(obj_name.keys())):
        plot = PCP(title=(f"Objective Space \n Highlighted: {obj_name[objective]}", 
                          {'pad': 30, 'fontsize':20}),
                n_ticks=10,
                legend=(False, {'loc': "upper left"}),
                labels=tech_labels,
                figsize=(13,6),
                )

        plot.set_axis_style(color="grey", alpha=0.5)
        plot.tight_layout = True

        plot.add(np.hstack(combined.iloc[:,4:-1].values).reshape(26,10), color="grey", alpha=0.3)
        cmap = obj_colors[objective]
        color=mcp.gen_color(cmap=cmap,n=n+1)

        subset_df = combined.nsmallest(n=n,columns=[objective]).reset_index(drop=True)

        for i, d in enumerate(subset_df.iloc[4:-1].values):
            plot.add(d, linewidth=3, color=color[i])
        plot.save(f"../docs/figures/04_benchmark_chapter/4_obj_design_space_{objective}.pgf")
    

    ### Design Space Box Plot ###
    fig, axes = plt.subplots(1,1,figsize=(13,6), facecolor='w', sharex=True, sharey=True)
    
    sb.boxenplot(ax=axes, data=combined.iloc[:,4:-1])

    axes.set_xticklabels(tech_labels, rotation=12.5, size=14)
    axes.set_xticks(range(len(tech_labels)))
    axes.set_xlabel("", size=14)
    # peak wind is at nearly 71 GW
    axes.set_ylim(0,30)
    axes.grid(alpha=0.2, which='major')
    axes.grid(alpha=0.05, which='minor')
    axes.set_title("Design Space (MGA)", fontsize=16)
    plt.tight_layout()
    plt.savefig(snakemake.output.four_obj_design_mga_boxplot)