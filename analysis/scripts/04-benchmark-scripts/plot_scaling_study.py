import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
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
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.shadow'] = False

if __name__ == "__main__":

    algorithm_df = pd.read_csv(snakemake.input.algorithm_times, index_col=0)
    algorithm_df.columns = ['n_days','logical dispatch', 'optimal dispatch']
    thread_df = pd.read_csv(snakemake.input.thread_times, index_col=0)


    #================= Algorithm comparison plot ================
    fig, ax = plt.subplots(figsize=(8,6))
    algorithm_df.plot(ax=ax, 
                      x='n_days', 
                      zorder=5, 
                      style={"logical dispatch":"b-",
                             "optimal dispatch":"r--"}
                    )
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(fontsize=16)
    ax.minorticks_on()
    ax.grid(which='major')
    ax.grid(which='minor', linestyle='--', alpha=0.2)
    ax.set_xlim(1,365)
    ax.set_ylabel("Time [seconds]", fontsize=18)
    ax.set_xlabel("Number of Modeled Days", fontsize=18)
    plt.savefig(snakemake.output.algorithm_scaling_plot)


    #================= Thread scaling plot ================
    fig, ax = plt.subplots(figsize=(8,6))
    sb.lineplot(thread_df.loc[thread_df['n_threads']<=12], 
                x='population',
                y='solve_time', 
                hue='n_threads', 
                ax=ax,
                palette="bright",
                style='n_threads')
    # ax.set_yscale('log')
    ax.legend(title="Number of Threads", fontsize=14, title_fontsize=14)
    ax.minorticks_on()
    ax.grid(which='major')
    ax.grid(which='minor', linestyle='--', alpha=0.2)
    ax.set_xlim(25,200)
    ax.set_ylabel("Time [seconds]", fontsize=18)
    ax.set_xlabel("Population per Generation", fontsize=18)
