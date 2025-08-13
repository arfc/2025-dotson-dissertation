# basic imports
import matplotlib.pyplot as plt
import numpy as np
from unyt import MW, GW, km, hour
from unyt import m, s, MW, W, kg
from unyt import unyt_array
import time
import pandas as pd
import matplotlib as mpl
mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.shadow'] = False



# osier imports
from osier import DispatchModel, LogicDispatchModel
import osier.tech_library as lib #wind, battery, natural_gas
from osier.utils import synchronize_units
import dill as pickle


cut_in = 3.0*m/s
cut_out = 25.0*m/s
rated = 13.0*m/s
diameter = 103*m
rated_power = 2.75*MW
air_density = 1.225*kg/m**3

power = lambda v: np.min((((0.5*np.pi/4)*(air_density*diameter**2)*v**3).to(MW), rated_power))*MW

def turbine_power(v):
    """
    Calculates the power output of a wind turbine.
    
    Parameters
    ----------
    v : float
        Wind speed in meters per second (m/s).
    Returns
    -------
    power : float
        Power in megawatts (MW).
    """
    
    wind_speed = v*m/s
    
    if wind_speed < cut_in:
        return 0*MW
    elif (wind_speed > rated) and (wind_speed < cut_out):
        return rated_power
    elif wind_speed >= cut_out:
        return 0*MW
    elif wind_speed >= cut_in:
        return power(wind_speed)


if __name__ == "__main__":
    phase_shift = 0  # horizontal shift [radians]
    base_shift = 2  # vertical shift [units of demand]
    day_hours = 24  # hours per day
    year_hours = 8760  # hours per year
    total_demand = 185  # [MWh], sets the total demand [units of energy]
    rng = np.random.default_rng(seed=1234)
    squish = 0.5
    n_days = 7

    natural_gas = lib.natural_gas
    natural_gas.capacity = 0.5e3*MW
    wind = lib.wind
    wind.capacity = 1e3*MW
    battery = lib.battery
    battery.capacity = 0.25e3*MW
    nuclear = lib.nuclear
    nuclear.capacity = 0.5e3*MW

    N = day_hours * n_days  # total number of time steps

    hours = np.linspace(0, N, N)
    y = np.sin((hours * np.pi / year_hours + phase_shift)) + np.ones(len(hours)) * (base_shift + 1) * squish

    demand = (np.sin((hours * np.pi / day_hours * 2 + phase_shift))
            * -1 + np.ones(N) * (base_shift + 1)) * 0.25 + y

    noise = rng.normal(size=N)*0.05
    demand += noise

    demand = demand / demand.max()
    demand *= 1e3 * MW

    with open(snakemake.output.demand_data, "wb") as f:
        pickle.dump(demand, f)


    rng = np.random.default_rng(seed=1234)
    wind_speed = rng.weibull(a=2, size=N)*9
    wind_power = np.array([turbine_power(u) for u in wind_speed]).astype(float)

    wind_speed /= wind_speed.max()
    wind_power /= wind_power.max()

    with open(snakemake.output.windspeed_data, "wb") as f:
        pickle.dump(wind_speed, f)

    with open(snakemake.output.windpower_data, "wb") as f:
        pickle.dump(wind_power, f)

    wind_power *= wind.capacity

    net_demand = demand - wind_power


    file_dict = snakemake.output.dispatch_data
    print(file_dict)
    
    optimizers = ['logical','optimal']

    objective_values_dict = {'scenario':['gas','gas','wind','wind','all','all'],
                             'optimizer':[],
                             'value':[]
                             }
    
    dispatch_models = {
                   'logical': LogicDispatchModel,
                   'optimal': DispatchModel,
                   }

    for i, outfile in enumerate(file_dict):
        if i % 2 == 0:
            opt = 'logical'
            allow_blackout = True
            curtailment = True
        else:
            opt = 'optimal'
            allow_blackout = False
            curtailment = True
        objective_values_dict['optimizer'].append(opt)
        
        print(f"Solving {objective_values_dict['scenario'][i]} with {opt} dispatch")
        if i in [0,1]:
            dem = demand
            tech_list = [natural_gas, nuclear]
        elif i in [2,3]:
            dem = net_demand
            tech_list = [natural_gas, nuclear]
        else:
            dem = net_demand
            tech_list = [natural_gas, nuclear, battery]
        model = dispatch_models[opt](technology_list = tech_list,
                                        net_demand = dem,
                                        allow_blackout = allow_blackout,
                                        curtailment=curtailment,
                                        penalty=1e-10,
                                        solver='appsi_highs')
        print(model)
        model.solve()

        objective_values_dict[f"value"].append(model.objective)
        # breakpoint()
        if i in [0,1]:
            dispatch_df = model.results
        else:
            dispatch_df = pd.concat([model.results, pd.DataFrame({'WindTurbine':wind_power})], axis=1)
        print(f"saving dispatch results to {outfile}")
        dispatch_df.to_csv(outfile)


    objective_df = pd.DataFrame(objective_values_dict)
    objective_df.to_csv(snakemake.output.objective_results)

    objective_df.style.format(precision=5).to_latex(snakemake.output.objective_table, 
                                                    hrules=True)








