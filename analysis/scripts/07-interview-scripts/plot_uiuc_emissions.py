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
    df.columns = ['Actual Emissions']

    future = tables[-1]
    future = future.set_index(pd.to_datetime(future['Year'], format="%Y"))\
                    .drop(columns=['Year'])\
                    .resample('YE').mean().interpolate(method='linear')
    future.index = future.index.year
    


    fig, ax = plt.subplots(figsize=(8,6))
    df.plot(ax=ax, marker='o')
    future.plot(ax=ax, linestyle='--', color='tab:green')
    ax.set_xlim(2008,2050)
    # ax.set_ylim(0,600e3)
    ax.set_yscale('log')
    ax.grid()
    ax.legend()
    ax.set_xlabel("")
    ax.set_ylabel("Emissions [MTCO2$_{eq}$]", fontsize=18)
    plt.savefig(snakemake.output.uiuc_emissions_plot)