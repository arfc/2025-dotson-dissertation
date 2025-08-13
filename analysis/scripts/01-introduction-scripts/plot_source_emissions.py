import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import matplotlib as mpl
# # mpl.use("pgf")
# # plt.rcParams['pgf.texsystem'] = 'pdflatex'
# # plt.rcParams['text.usetex'] = True
# # plt.rcParams['pgf.rcfonts'] = False
# plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"


if __name__ == "__main__":
    df = pd.read_html('https://en.wikipedia.org/wiki/Life-cycle_greenhouse_gas_emissions_of_energy_sources')[-2]
    df.columns = ['Technology', 'Class', 'gCO2eq/kWh']
    df = df.loc[~df['Class'].str.contains("with CCS")]


    fig, ax = plt.subplots(figsize=(8, 6), facecolor='w')

    # plot a bar chart
    sb.barplot(ax=ax, x="Technology", y="gCO2eq/kWh", data=df, estimator=np.mean, 
                    errorbar=('ci',65), capsize=.2, color='lightblue',
                    zorder=2)
    # ax.set_ylabel(r'Lifecycle Emissions $\left[\frac{{gCO}_{2}eq}{kWh}\right]$',
    #             fontsize=16)
    ax.set_ylabel(r'Lifecycle Emissions [gCO$_2$eq/kWh]',
                fontsize=16)
    ax.set_xlabel('')
    ax.grid(zorder=0)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14, rotation=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig(snakemake.output.source_emissions_plot)
