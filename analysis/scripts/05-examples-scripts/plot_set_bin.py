import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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

    full_df = pd.read_csv(snakemake.input.metric_data, index_col=0)
    full_df.columns = [" ".join(c.split(" ")[:-1]) for c in list(full_df.columns)]

    fig, ax = plt.subplots(figsize=(9,6))
    (full_df['activity at 100 yrs'].sort_values(ascending=False)/1e6).plot(ax=ax, kind='bar', zorder=3, color='tab:red')
    # ax.grid(zorder=0)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylabel(r"SNF+HLW Activity $\left[\frac{MCi}{GWe-yr}\right]$", fontsize=16)
    ax.set_xlabel("")
    ax.set_ylim(0, 2)

    ax.axvline(x=1.5, color='k')
    ax.text(0.1, 1.85, "D", fontsize=16, color='blue')

    ax.axvline(x=30.5, color='k')
    ax.text(16, 1.85, "C", fontsize=16, color='blue')
    ax.text(35, 1.85, "B", fontsize=16, color='blue')

    # ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig(snakemake.output.set_bin_plot)
    # plt.show()