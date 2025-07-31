import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib as mpl
import dill as pickle
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
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['legend.shadow'] = False

if __name__ == "__main__":

    color_map = {
             'Nuclear':'tab:blue', 
             'Battery':'Purple',
             'NaturalGas Conv':'tab:brown',
             'Battery charge':'orange',
             'Curtailment':'gray',
             'WindTurbine':'tab:cyan',
            }

    dispatch_files = snakemake.input.dispatch_data
    plot_files = snakemake.output.dispatch_plots

    with open(snakemake.input.demand_data, "rb") as f:
        demand_data = pickle.load(f)



    for i, (dispatch, plot) in enumerate(zip(dispatch_files, plot_files)):
        print(f'plotting {dispatch}')
        df = pd.read_csv(dispatch, index_col=0)
        df.columns = [column.replace("_", " ") for column in list(df.columns)]    


        fig, ax = plt.subplots(figsize=(8,6))
        tech_order = [tech for tech in list(color_map.keys()) if tech in df.columns]
        df[tech_order].plot.area(ax=ax,zorder=2, color=color_map)
        ax.plot(demand_data.to_value(), zorder=3, label='Demand', color='k', linestyle='--')
        ax.grid()
        ax.legend(ncols=3, loc='lower left')
        ax.set_xlim(0,167)
        ax.set_ylabel("Energy [MWh]", fontsize=18)
        ax.set_xlabel("Time [hours]", fontsize=18)
        # plt.tight_layout()
        plt.savefig(plot)


    # ================= Combined Plot ===================

    opt_df = pd.read_csv('../data/all_dispatch_results_optimal.csv', index_col=0)
    opt_df.columns = [column.replace("_", " ") for column in list(opt_df.columns)]   
    log_df = pd.read_csv('../data/all_dispatch_results_logical.csv', index_col=0)
    log_df.columns = [column.replace("_", " ") for column in list(log_df.columns)]   
    tech_order = [tech for tech in list(color_map.keys()) if tech in log_df.columns]

    fig, ax = plt.subplots(2,1, figsize=(10,8), sharex=True)
    log_df[tech_order].plot.area(ax=ax[0],zorder=2, color=color_map, legend=False)
    ax[0].plot(demand_data.to_value(), zorder=3,color='k', linestyle='--')
    ax[0].grid()
    # ax[0].set_xlim(0,168)
    ax[0].set_ylabel("Energy [MWh]", fontsize=18)

    opt_df[tech_order].plot.area(ax=ax[1],zorder=2, color=color_map)
    ax[1].plot(demand_data.to_value(), zorder=3, label='Demand', color='k', linestyle='--')
    ax[1].grid()
    ax[1].legend(ncols=3, loc='lower left')
    ax[1].set_xlim(0,167)
    ax[1].set_ylabel("Energy [MWh]", fontsize=18)
    ax[1].set_xlabel("Time [hours]", fontsize=18)

    x_loc, y_loc = 2.5, 1450
    ax[0].text(x_loc,y_loc, "a)", fontsize=14, bbox=dict({'facecolor':'w'}))
    ax[1].text(x_loc,y_loc, "b)", fontsize=14, bbox=dict({'facecolor':'w'}))

    plt.tight_layout()

    plt.savefig("../docs/figures/04_benchmark_chapter/dispatch_comparison_plot_x2.pgf")