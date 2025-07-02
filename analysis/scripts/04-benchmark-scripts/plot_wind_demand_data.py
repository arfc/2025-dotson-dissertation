# basic imports
import pandas as pd
import numpy as np
from unyt import kW, minute, hour, day, MW, kg, lb, kWh, MWh
import dill as pickle

# osier imports
from osier import CapacityExpansion
import osier.tech_library as lib
from osier.equations import total_cost, annual_emission
from osier import technology_dataframe

# import megatonnes from unyt -- must be done after importing osier
from unyt import megatonnes

# pymoo imports
from functools import partial

# plotting imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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
    day_hours = 24
    n_days = 7
    N = day_hours * n_days
    hours = np.linspace(0, N, N)

    with open(snakemake.input.demand_data, "rb") as f:
        demand = pickle.load(f)
    demand = demand/demand.max()


    with open(snakemake.input.windspeed_data, "rb") as f:
        wind_speed = pickle.load(f)

    with open(snakemake.input.windpower_data, "rb") as f:
        wind_power = pickle.load(f)


    # demand plot
    print("Plotting the demand profile")
    fig, ax = plt.subplots()
    ax.plot(hours, demand, color='k', linestyle='--')
    ax.grid()
    ax.set_ylabel("Normalized Demand [-]", fontsize=18)
    ax.set_xlabel("Time [hours]", fontsize=18)
    ax.set_xlim(0,168)
    plt.savefig(snakemake.output.demand_plot)

    # plot wind
    print("Plotting the wind profile")
    fig, ax = plt.subplots()
    ax.plot(hours, wind_power, color='tab:blue', label="Turbine Power")
    ax.plot(hours, wind_speed, color='green', linestyle='-', label="Wind Speed")
    ax.grid()
    ax.set_ylabel("Normalized Wind Data [-]", fontsize=18)
    ax.set_xlabel("Time [hours]", fontsize=18)
    ax.set_xlim(0,47)
    ax.legend(fontsize=12)
    plt.savefig(snakemake.output.wind_plot)