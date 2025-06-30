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
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.shadow'] = False

if __name__ == "__main__":

    color_map = {
             'Nuclear':'tab:blue', 
             'Battery':'Purple',
             'NaturalGas_Conv':'tab:brown',
             'Battery_charge':'orange',
             'Curtailment':'gray',
             'WindTurbine':'tab:cyan',
            }

    dispatch_files = snakemake.input.dispatch_data
    plot_files = snakemake.output.dispatch_plots

    with open(snakemake.input.demand_data, "rb") as f:
        demand_data = pickle.load(f)



    for i, (dispatch, plot) in enumerate(zip(dispatch_files, plot_files)):
        df = pd.read_csv(dispatch, index_col=0)
        
        fig, ax = plt.subplots(figsize=(8,6))
        tech_order = [tech for tech in list(color_map.keys()) if tech in df.columns]
        df[tech_order].plot.area(ax=ax,zorder=2, color=color_map)
        ax.plot(demand_data.to_value(), zorder=3, label='Demand', color='k', linestyle='--')
        ax.grid()
        ax.legend(ncols=3, loc='lower left')
        ax.set_xlim(0,168)
        ax.set_ylabel("Energy [MWh]", fontsize=18)
        ax.set_xlabel("Time [hours]", fontsize=18)
        plt.savefig(plot)