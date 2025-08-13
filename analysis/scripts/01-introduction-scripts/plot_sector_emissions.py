import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib as mpl
import pandas as pd
import seaborn as sb
import camelot

mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"

if __name__ == "__main__":

    tables = camelot.read_pdf(snakemake.input.epa_inventory,
                          flavor="stream",
                          pages="21")
    
    df = tables[0].df.iloc[5:12,0:8].reset_index(drop=True)
    df.columns = df.iloc[0]
    df = df.iloc[1:,:]
    df.set_index('Economic Sectors', inplace=True)
    for col in df.columns:
        df[col] = df[col].str.replace(",","")

    df = df.astype(float)
    # convert to annual percentages
    df = (df.div(df.sum(axis=0))*100).round(2)


    fig, ax = plt.subplots(figsize=(8,6), facecolor='w')
    df['2022'].plot.bar(ax=ax, legend=False, color='orange', zorder=2)
    ax.grid(which='major', alpha=0.5, zorder=0)
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.grid(which='minor', alpha=0, linestyle='--', zorder=0)
    ax.set_yticklabels(labels=ax.yaxis.get_ticklabels(), fontsize=16)
    xlabels=['Transportation', 'Electric Power', 'Industry', 'Agriculture',
    'Commercial', 'Residential']
    # ax.set_xticklabels(labels=ax.xaxis.get_ticklabels(), fontsize=12, rotation=20)
    ax.set_xticklabels(labels=xlabels, fontsize=11, rotation=0)
    # ax.set_xticklabels(labels=df.index, fontsize=16, rotation=20)
    ax.set_xlabel('')
    for container in ax.containers:
        ax.bar_label(container, fmt='%.1f%%', fontsize=14)

    ax.set_ylabel("Share of Annual CO$_2$ Emissions", fontsize=16)
    ax.set_ylim(0,31)

    plt.tight_layout()
    plt.savefig(snakemake.output.sector_emissions_plot)