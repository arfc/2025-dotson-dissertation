import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('ggplot')
import pandas as pd

import matplotlib as mpl
# mpl.use("pgf")
# plt.rcParams['pgf.texsystem'] = 'pdflatex'
# plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"


if __name__ == "__main__":
    # set up basic data

    cf_s = 0.2
    cf_w = 0.35
    cf_n = 0.93

    cost_s = 0.811
    cost_w = 1.411

    demand_2019 = 138 # TWh
    nuc_gen_2019 = 98.7 # TWh
    TWh_to_GW = 8.76
    cap_needed = (demand_2019-nuc_gen_2019)/TWh_to_GW

    max_wind = cap_needed/cf_w
    max_solar = cap_needed/cf_s
    
    # plot +
    fig, ax = plt.subplots(figsize=(8,6))

    max_wind = 1.0
    max_solar = 1.0
    cost_s = 1.0
    cost_w = 0.8
    px1 = np.array([max_solar,0])
    px2 = np.array([0, max_wind])
    cx1 = np.array([cost_s*max_solar, 0])
    cx2 = np.array([0, cost_w*max_wind])
    slack = 0.1
    slacked_x1 = cx1*(1+slack)
    slacked_x2 = cx2*(1+slack)
    opt1 = (1,0)
    opt2 = (0.6, 0.4)
    ax.plot(px1, px2, label='x$_1$ + x$_2$ = 1', linestyle='-', color='k')
    ax.plot(cx1, cx2, label='min(c$_1$x$_1$ + c$_2$x$_2$)', color='gray')
    ax.plot(slacked_x1, slacked_x2, label='c$_1$x$_1$ + c$_2$x$_2$ $\leq c_1\cdot$slack', color='gray', linestyle='--')
    ax.set_ylim(0, max_wind*(1+slack))
    ax.set_xlim(0, max_solar*(1+slack))
    # ax.fill_between(px1, px2, max_wind, alpha=0.2, label='Feasible Space')
    # ax.scatter(x=17, y=6.3, c='tab:red', label='CEJA Goal', s=90)
    ax.scatter(x=opt2[0], y=opt2[1], c='k', label='', facecolor='k', linewidth=1, s=90)
    ax.scatter(x=opt1[0],y=opt1[1], color='k', marker='o', s=90, facecolor='w')
    # ax.minorticks_on()
    # ax.grid(which='minor', linestyle='--', color='gray', alpha=0.4)
    ax.grid(which='major', linestyle='-', color='gray', alpha=0.3)
    ax.legend(loc='upper right', fancybox=True, shadow=True, facecolor='w',fontsize=16)
    ax.set_xlabel('x$_1$',fontsize=20)
    ax.set_ylabel('x$_2$',fontsize=20)
    # ax.set_title('What\'s the Optimal Mix of Renewable Energy? \nDemonstration of Modelling to Generate Alternatives',
    #              fontsize=24)
    ax.tick_params(axis='both', which='major', labelsize=14)
    # ax.arrow(x=cx1[0]-1,y=cx1[1]+1, dx=slacked_x1[0]*0.97-cx1[0],dy=slacked_x1[1]-cx1[1], width=0.1, color='k')

    ax.annotate(f'Optimum: {opt1}', opt1, xycoords='data', xytext=(opt1[0]-0.5, opt1[1]+0.1), arrowprops=dict(facecolor='k',
                                                                                                        arrowstyle='->, head_width=0.35'),
                fontsize=14)
    ax.annotate(f'MGA Solution: {opt2}', opt2, xycoords='data', xytext=(opt2[0]+0.1, opt2[1]+.1), arrowprops=dict(facecolor='k',
                                                                                                        arrowstyle='->, head_width=0.35'),
                fontsize=14)
    ax.set_facecolor('w')
    fig.suptitle("Modeling-to-Generate-Alternatives", fontsize=24)
    ax.set_title('Design Space', fontsize=20)
    plt.tight_layout()
    plt.savefig(snakemake.output.standard_mga)
