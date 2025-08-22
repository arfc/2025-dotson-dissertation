import matplotlib.pyplot as plt
import pandas as pd
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
    df = pd.read_csv(snakemake.input.ptc_data, 
                     usecols=['Year','Ameren','ComEd'], 
                     index_col='Year')

    fig, ax = plt.subplots(figsize=(8,6))
    df.plot(ax=ax, marker='o')
    ax.set_ylim(-3,3)
    ax.axhline(y=0, linestyle='--', color='tab:red')
    ax.set_ylabel(r"$\Delta$ (ARES - PTC) [cent/kWh]", fontsize=18)
    ax.set_xlabel('')
    ax.grid()

    plt.savefig(snakemake.output.ptc_plot)