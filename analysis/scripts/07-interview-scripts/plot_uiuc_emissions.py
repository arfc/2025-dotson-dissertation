import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import matplotlib as mpl
from mycolorpy import colorlist as mcp

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
    tables = pd.read_html("https://icap.sustainability.illinois.edu/metric/total-campus-ghg-emissions")

    df = tables[0].iloc[:,:2]
    df = df[::-1].set_index(pd.to_datetime(df['Applicable Date'][::-1]).dt.year)\
                 .drop(columns=['Applicable Date'])
    df.columns = ['UIUC Actual']

    future = tables[-1]
    future = future.set_index(pd.to_datetime(future['Year'], format="%Y"))\
                    .drop(columns=['Year'])\
                    .resample('YE').mean().interpolate(method='linear')
    future.index = future.index.year
    future.columns = ['UIUC Target']
    
    # ==== Illinois emissions ====
    il_df = pd.DataFrame({'year':[2005,2021,2030,2050],
                      'emissions':[283.63e6, 228.01e6,135e6,0]}
                     )
    il_df['year'] = pd.to_datetime(il_df['year'], format="%Y")
    il_df.set_index('year',inplace=True,drop=True)
    il_df = il_df.resample('YE').mean().interpolate('linear')
    il_df.index = il_df.index.year
    il_df.columns = ['IL State Target']


    fig, ax = plt.subplots(figsize=(8,6))
    df.plot(ax=ax, marker='o')
    future.plot(ax=ax, linestyle='--', color='tab:green')
    # plot illinois emissions
    il_df.plot(ax=ax, linestyle='-.', color='tab:purple')
    ax.set_xlim(2008,2050)
    # ax.set_ylim(0,600e3)
    ax.set_yscale('log')
    ax.grid()
    ax.grid(which='minor', alpha=0.4, linestyle='--')
    ax.legend()
    ax.set_xlabel("")
    ax.set_ylabel("Emissions [MTCO2$_{eq}$]", fontsize=18)
    plt.savefig(snakemake.output.uiuc_emissions_plot)
