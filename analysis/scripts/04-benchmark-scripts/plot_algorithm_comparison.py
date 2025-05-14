from osier import *
import dill
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde

mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'lightgray'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.shadow'] = True

plt.rcParams['text.latex.preamble'] = r'\newcommand{\mathdefault}[1][]{}'


if __name__ == "__main__":

    # open pickle files
    print('Opening datasets... ')
    with open(snakemake.input.unsga3_results, "rb") as file:
        old_results = dill.load(file)

    pf_usnga3 = old_results.F
    with open(snakemake.input.pf_nsga2, "rb") as file:
        pf_nsga2 = dill.load(file)
    with open(snakemake.input.pf_nsga3, "rb") as file:
        pf_nsga3 = dill.load(file)


    

    fig, axScatter = plt.subplots(figsize=(8,6), facecolor='white')
    divider = make_axes_locatable(axScatter)
    axHistx = divider.append_axes("top", 1.2, pad=0., sharex=axScatter)
    axHisty = divider.append_axes("right", 1.2, pad=0., sharey=axScatter)

    # algs_list = [pf_usnga3, pf_nsga3, pf_nsga2, pf_temoa_remap]
    algs_list = [pf_usnga3, pf_nsga3, pf_nsga2]
    markers = ['s', '*', 'o', '^']
    colors = ['tab:blue', 'tab:orange', 'tab:purple', 'tab:red']
    labels = ['Pymoo/UNSGA3', 'DEAP/NSGA3', 'DEAP/NSGA2', 'Temoa+MGA']
    for i,data in enumerate(algs_list):
        x = data[:,0]
        y = data[:,1]
        
        density_x = gaussian_kde(x)
        density_y = gaussian_kde(y)

        axScatter.scatter(x, y, s=30, facecolors=colors[i], label=labels[i], marker=markers[i])
        linex = density_x(np.sort(x))
        liney = density_y(np.sort(y))
        axHistx.plot(np.sort(x), linex, color=colors[i])
        axHistx.fill_between(np.sort(x), 0, linex,alpha=0.3, color=colors[i])
        
        axHisty.plot(liney,np.sort(y), color=colors[i])
        axHisty.fill_between(liney, np.sort(y), np.min(pf_nsga3[:,1]), alpha=0.3, color=colors[i])
            
    plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
            visible=False)

    plt.setp(axHistx.get_yticklabels() + axHisty.get_xticklabels(),
            visible=False)


    axScatter.set_ylim(0,22)
    axScatter.set_xlim(5e3,17.5e3)
    axHisty.tick_params(axis='x', labelrotation=-90)
    axScatter.grid(alpha=0.2, which='major')
    axScatter.grid(alpha=0.05, which='minor')
    axScatter.legend(shadow=True, fontsize=16)

    axScatter.set_xlabel("Total Cost [M\$]", fontsize=18)
    axScatter.set_ylabel("Emissions [MT CO2eq]", fontsize=18)
    axHistx.grid(alpha=0.2)
    axHisty.grid(alpha=0.2)
    # axScatter.set_xscale('log')
    # axScatter.set_yscale('log')
    plt.tight_layout()
    plt.savefig(snakemake.output.alg_comp_plot)
    # plt.show()