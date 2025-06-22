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

solver = "appsi_highs"

def eroi_objective(technology_list, solved_dispatch_model):
    """
    Calculate the objective to ``maximize'' the Energy Return on Investment (EROI)
    by minimizing its inverse.
    """
    weighted_eroi = np.array([t.capacity.to_value() * t.eroi for t in technology_list]).sum()
    obj_value = 1/weighted_eroi
    return obj_value

if __name__ == "__main__":
    n_hours = 24  # hours per day
    n_days = 7  # days to model
    N = n_hours*n_days  # total number of time steps
    phase_shift = np.pi/2  # horizontal shift [radians]
    base_shift = 2  # vertical shift [units of demand]
    hours = np.linspace(0,N,N)

    solar_cp = (-np.sin((hours*np.pi/n_hours*2+phase_shift)))
    rng = np.random.default_rng(1234)

    solar_cp += rng.normal(size=N)*15e-2
    solar_cp[solar_cp<0] = 0

    solar_cp = solar_cp/solar_cp.max()  # rescale

    # plot solar output
    plt.plot(hours, solar_cp, color='gold')
    plt.ylabel('Solar Availability [-]')
    plt.xlabel('Time [hr]')
    plt.grid(alpha=0.2)
    plt.savefig(snakemake.output.solar_plot)

    peak_demand = 1e3*MW 
    demand = np.ones(N)*peak_demand

    # Get technologies
    natural_gas = lib.natural_gas
    natural_gas_adv = lib.natural_gas_adv
    nuclear_ap1000 = lib.nuclear
    nuclear_smr = lib.nuclear_adv
    solar = lib.solar
    battery = lib.battery

    # Get emissions data
    emission_df = pd.read_html("https://www.eia.gov/tools/faqs/faq.php?id=74&t=11")[0].droplevel(1, axis=1).set_index('Unnamed: 0_level_0').iloc[:-1,:-2]
    emission_df.columns = ['Total kWh', 'Total CO2 million mt', 'million short tons','lbs per kWh']
    emission_df.index.name = ''

    # Get EROI data 
    eroi_df = pd.read_csv(snakemake.input.eroi_data, 
                          index_col=0)

    # Get cost data
    cost_df = pd.read_csv(snakemake.input.tech_costs, 
                          index_col=['core_metric_parameter',
                                     'display_name'])


    # for a NG 1-on-1 H Frame design
    natural_gas.co2_rate = (float(emission_df.at['Natural gas', 'lbs per kWh'])*lb/(kW*hour)).to(megatonnes/(MW*hour))
    natural_gas.capital_cost = cost_df.at[('CAPEX','Natural Gas'), 'value'] / kW  / 1e6
    natural_gas.om_cost_variable = cost_df.at[('Variable O&M','Natural Gas'), 'value'] / MWh / 1e6
    natural_gas.om_cost_fixed = cost_df.at[('Fixed O&M','Natural Gas'), 'value'] / kW / 1e6
    natural_gas.eroi = eroi_df.at['Natural Gas (CCGT)', 'EROIstd'] # from Walmsley et al.
    # same design, with 95% CCS
    natural_gas_adv.co2_rate = natural_gas.co2_rate * 0.05
    natural_gas_adv.capital_cost = cost_df.at[('CAPEX','Natural Gas CCS'), 'value'] / kW / 1e6
    natural_gas_adv.om_cost_variable = cost_df.at[('Variable O&M','Natural Gas CCS'), 'value'] / MWh / 1e6
    natural_gas_adv.om_cost_fixed = cost_df.at[('Fixed O&M','Natural Gas CCS'), 'value'] / kW / 1e6
    natural_gas_adv.eroi = eroi_df.at['Natural Gas (CCGT & CCS)', 'EROIstd'] # from Walmsley et al.

    # utility scale solar, with good insolation, middling estimate
    solar.capital_cost = cost_df.at[('CAPEX','UtilityPV'), 'value']/kW / 1e6
    solar.eroi = eroi_df.at['Solar PV (Mono-Si, SE-med)', 'EROIstd'] # from Walmsley et al. 2018
    solar.om_cost_fixed = cost_df.at[('Fixed O&M','UtilityPV'), 'value']/kW / 1e6
    solar.co2_rate = 0*(megatonnes/(MW*hour))
    solar.fuel_cost = 0/(MW*hour)

    battery.capital_cost = cost_df.at[('CAPEX','Battery'), 'value'] / kW / 1e6
    battery.om_cost_fixed = cost_df.at[('Fixed O&M','Natural Gas'), 'value']/ kW / 1e6
    battery.eroi = 10  # actually an 'ESOI,' from Barnhart and Benson 2013
    battery.co2_rate = 0*(megatonnes/(MW*hour))
    battery.fuel_cost = 0/(MW*hour)


    nuclear_ap1000.capital_cost = cost_df.at[('CAPEX','Nuclear'), 'value'] / kW / 1e6
    nuclear_ap1000.om_cost_fixed = cost_df.at[('Fixed O&M','Nuclear'), 'value'] / kW / 1e6
    nuclear_ap1000.om_cost_variable = cost_df.at[('Variable O&M','Nuclear'), 'value']/ MWh / 1e6
    nuclear_ap1000.eroi = eroi_df.at['Nuclear (100% centrifuge)', 'EROIstd']
    nuclear_ap1000.co2_rate = 0*(megatonnes/(MW*hour))

    nuclear_smr.capital_cost = cost_df.at[('CAPEX','Advanced Nuclear'), 'value'] / kW / 1e6
    nuclear_smr.om_cost_fixed = cost_df.at[('Fixed O&M','Advanced Nuclear'), 'value']/ kW / 1e6
    nuclear_smr.om_cost_variable = cost_df.at[('Variable O&M','Advanced Nuclear'), 'value']/ MWh / 1e6
    nuclear_smr.eroi = eroi_df.at['Nuclear (100% centrifuge)', 'EROIstd']
    nuclear_smr.co2_rate = 0*(megatonnes/(MW*hour))

    tech_list = [natural_gas, natural_gas_adv, solar, battery, nuclear_ap1000, nuclear_smr]

    tech_data_df = technology_dataframe(tech_list, cast_to_string=True).T.iloc[[13,14,15,16, 17, 21, ],:]
    tech_data_df.index = tech_data_df.index.str.replace("_", " ")
    tech_data_df.columns = [col.replace("_", " ") for col in tech_data_df.columns]
    tech_data_df.style.format(precision=2).to_latex(snakemake.output.dc_tech_table, hrules=True)

    problem = CapacityExpansion(technology_list=tech_list,
                                demand=demand,
                                solar=solar_cp,
                                upper_bound = 1 / solar.capacity_credit,
                                objectives = [total_cost, 
                                            partial(annual_emission, emission='co2_rate'),
                                            eroi_objective],
                                solver=solver
                                )
    

    with open(snakemake.output.dc_problem, "wb") as file:
        pickle.dump(problem, file)

    