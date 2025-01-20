import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
from scipy.optimize import nnls
from scipy.optimize import curve_fit

mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"
import seaborn as sb
import pandas as pd

from pymoo.problems import get_problem
from pymoo.util.plotting import plot
from pymoo.visualization.scatter import Scatter

if __name__ == "__main__":
    # F = get_problem("truss2d").pareto_front()
    F = get_problem("bnh").pareto_front()
    a = min(F[:,0])
    b = max(F[:,0])
    f1 = F[:,0]
    f2 = F[:,1]
    shift = 0.75


    """
    Truss 2D example plot
    """
    fig, ax = plt.subplots(figsize=(8,6), facecolor='w', edgecolor='k')
    ax.scatter(F[:,0], F[:,1], label="Solutions", facecolor='none', edgecolor='r')
    ax.plot(F[:,0], F[:,1], label="Pareto-front", color='k', lw=2.5)
    ax.set_xlabel('f1', fontsize=14)
    ax.set_ylabel('f2', fontsize=14)
    ax.fill_between(f1, f2, alpha=0.2, label='Infeasible Space', color='gray')
    ax.fill_between(f1, f2, max(F[:,1]), alpha=0.2, label='Feasible Space', color='tab:blue')

    ax.set_xlim(a,b)
    ax.set_ylim(min(F[:,1]),max(F[:,1]))
    ax.scatter((2-shift)*min(F[:,0]),(2-shift)*min(F[:,1]), label='Ideal', marker="*", s=300, color='tab:blue')
    ax.scatter(0.95*max(F[:,0]),0.95*max(F[:,1]), label='Nadir', marker="s", s=200, color='tab:orange')
    ax.legend(fontsize=14)
    plt.savefig("../docs/figures/truss2d_pareto.pgf")



    """
    Pareto Near-Optimal Space
    """
    slack = 0.2
    alpha = 0.5
    F1 = f1 * (1+slack)
    F2 = f2 * (1+slack)
    fig, ax = plt.subplots(figsize=(8,6), facecolor='w', edgecolor='k')
    subopt = ax.plot(F1, F2, label="", color='lightgray',alpha=alpha)
    opt=ax.plot(f1, f2, label="Pareto-front", color='k', lw=2.5)
    ax.set_xlabel('f1', fontsize=14)
    ax.set_ylabel('f2', fontsize=14)
    plt.fill(np.append(f1, F1[::-1]), np.append(f2, F2[::-1]), 'lightgrey', alpha=alpha, label="Near-optimal space")

    ax.set_xlim(min(F1),max(F1))
    ax.set_ylim(min(F2),max(F2))
    ax.legend(fontsize=14)
    plt.savefig("../docs/figures/near-optimal-pareto.pgf")


    """
    Visualize interior points
    """
    F3 = F*(1+slack)

    rng = np.random.default_rng(seed=1234)
    R = rng.uniform(size=(1000,2))
    # for Truss2d
    # R[:,1] = R[:,1]*15e4
    # R[:,0] = R[:,0]*8e-2
    # R_sub = R[(R[:,1] < 12e4) & (R[:,0] < 0.06)]
    # for BNH
    R[:,1] = R[:,1]*65
    R[:,0] = R[:,0]*165    
    R_sub = R[(R[:,1] < 65) & (R[:,0] < 165)]

    interior_pts = []
    for p in R_sub:
        cond_1 = np.any((p < F3).sum(axis=1) == 2)
        cond_2 = np.any((p > F).sum(axis=1) == 2)
        if cond_1 and cond_2:
            interior_pts.append(p)

    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(*zip(*F), color='black', lw=3, label='Pareto Front')
    ax.plot(*zip(*F3), color='black', lw=1, alpha=0.2)
    ax.scatter(*zip(*R_sub), c='tab:blue', s=3, label='Tested points')
    ax.scatter(*zip(*np.array(interior_pts)), c='tab:red', s=20, label="Alternative solutions")
    ax.fill(np.append(F[:,0], F3[:,0][::-1]), np.append(F[:,1],F3[:,1][::-1]), alpha=0.2, color='gray')
    ax.fill_between(f1*0.98, f2*0.98, alpha=1, color='w')
    ax.set_xlim(min(F1),max(F1))
    ax.set_ylim(min(F2),max(F2))
    ax.set_xlabel('f1', fontsize=14)
    ax.set_ylabel('f2', fontsize=14)
    ax.legend(fontsize=14, shadow=True, loc='upper right')
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.savefig("../docs/figures/nd-mga-paretofront.pgf")