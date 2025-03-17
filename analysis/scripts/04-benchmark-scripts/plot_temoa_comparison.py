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
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 16
import seaborn as sb
import pandas as pd
import dill

if __name__ == "__main__":

    print("Opening datasets")

    with open(snakemake.input.unsga3_objective, "rb") as file:
        osier_F = dill.load(file)
    with open(snakemake.input.unsga3_design, "rb") as file:
        osier_X = dill.load(file)
    with open(snakemake.input.temoa_F2, "rb") as file:
        temoa_F2 = dill.load(file)
    with open(snakemake.input.temoa_mga_design, "rb") as file:
        temoa_X = dill.load(file)
    with open(snakemake.input.techs, "rb") as file:
        techs = dill.load(file)


    ### Plot first comparison ###

    df = pd.DataFrame({"Cost":osier_F[:,0], "Carbon":osier_F[:,1]})
    slack = 1.1
    f1 = osier_F[:,0]*slack
    f2 = osier_F[:,1]*slack
    df2 = pd.DataFrame({"Cost":f1, "Carbon":f2})
    color_1 = 'blue'
    color_2 = 'tab:blue'
    color_3 = 'tab:red'
    alpha = 0.5

    fig, ax = plt.subplots(facecolor='w')

    df.sort_values(by="Cost").plot.line(ax=ax,
                                        x='Cost', 
                                        y='Carbon',
                                        legend=False, 
                                        color=color_1, 
                                        marker='o', 
                                        label='osier',
                                        markersize=1)
    df2.sort_values(by="Cost").plot.line(ax=ax,
                                         x='Cost', 
                                         y='Carbon',
                                         legend=False, 
                                         color=color_2, 
                                         alpha=alpha,
                                        label='near-optimal space (osier)')
    ax.fill(np.append(df.sort_values(by="Cost").Cost.values, 
                      df2.sort_values(by="Cost").Cost.values[::-1]),
        np.append(df.sort_values(by="Cost").Carbon.values, 
                  df2.sort_values(by="Cost").Carbon.values[::-1]),
        alpha=alpha, color=color_2)
    ax.scatter(temoa_F2[:, 0], 
               temoa_F2[:, 1], 
               s=20, 
               facecolors='none', 
               edgecolors='red', 
               label='temoa+mga')


    ax.set_xlim(min(temoa_F2[:,0]),max(osier_F[:,0])*slack)
    ax.set_ylim(min(osier_F[:,1]),max(temoa_F2[:,1])*slack)
    ax.grid()

    temoa_min = temoa_F2[:,0].min()

    ax.axvspan(xmin=temoa_min, 
               xmax=temoa_min*1.1, 
               alpha=0.2, 
               hatch='/', 
               label='near-optimal space (Temoa)', 
               color=color_3)

    # inset plot
    axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])

    df.sort_values(by="Cost").plot.line(ax=axins,
                                        x='Cost', 
                                        y='Carbon',
                                        legend=False, 
                                        color=color_1, 
                                        marker='o', 
                                        label='osier',
                                        markersize=6)
    df2.sort_values(by="Cost").plot.line(ax=axins,
                                         x='Cost', 
                                         y='Carbon',
                                         legend=False, 
                                         color=color_2, 
                                         alpha=alpha,
                                         label='Near-optimal space')
    axins.fill(np.append(df.sort_values(by="Cost").Cost.values, 
                         df2.sort_values(by="Cost").Cost.values[::-1]),
        np.append(df.sort_values(by="Cost").Carbon.values, 
                  df2.sort_values(by="Cost").Carbon.values[::-1]),
        alpha=alpha, color=color_2)
    axins.scatter(temoa_F2[:, 0], 
                  temoa_F2[:, 1], 
                  s=20, 
                  facecolors='none', 
                  edgecolors='red', 
                  label='temoa+mga')
    axins.axvspan(xmin=temoa_min, 
                  xmax=temoa_min*1.1, 
                  alpha=0.2, hatch='/', 
                  label='near-optimal space (Temoa)', 
                  color=color_3)

    axins.set_xlim(min(temoa_F2[:,0])-20,5350)
    axins.set_ylim(9,min(temoa_F2[:,1])+3)
    axins.grid()
    axins.set_xlabel('')

    ax.indicate_inset_zoom(axins, edgecolor="k", label='')

    plt.xlabel("Total Cost (M\$)", fontsize=18)
    plt.ylabel("CO2 emissions (MT CO2)", fontsize=18)
    plt.legend(loc='lower right')

    plt.savefig(snakemake.output.temoa_comparison_01)
    # plt.show()