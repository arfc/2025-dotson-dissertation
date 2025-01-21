from pymoo.util.ref_dirs import get_reference_directions
from pymoo.problems import get_problem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"


if __name__ == "__main__":

    np.random.seed(1234)
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=12)
    pf = get_problem("mw8").pareto_front(ref_dirs)
    slack = 0.1
    N = 1000
    pf_slack = pf * (1 + slack)
    R = np.random.rand(N, 3)
    int_pts = []
    for p in R:
        cond1 = np.any((p < pf_slack).sum(axis=1) == 3)
        cond2 = np.any((p > pf).sum(axis=1) == 3)
        if cond1 and cond2:
            int_pts.append(p)

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(132, projection="3d")
    ax2 = fig.add_subplot(133, projection="3d")
    ax3 = fig.add_subplot(131, projection="3d")
    pt_size = 30
    pt_color = "tab:red"
    pf_color = "tab:blue"

    ax.plot_trisurf(*zip(*pf), alpha=1, color=pf_color)
    ax.scatter3D(*zip(*np.array(int_pts)), facecolor=pt_color, s=pt_size)
    ax.plot_trisurf(*zip(*pf_slack), color='tab:green', alpha=0.2)
    ax.view_init(azim=45, elev=45)
    ax.set_xlabel('f1', fontsize=14, labelpad=2)
    ax.set_ylabel('f2', fontsize=14, labelpad=2)
    ax.set_zlabel('f3', fontsize=14, labelpad=2)

    ax3.plot_trisurf(*zip(*pf), alpha=1, color=pf_color)
    ax3.scatter3D(*zip(*np.array(int_pts)), facecolor=pt_color, s=pt_size)
    ax3.plot_trisurf(*zip(*pf_slack), color='tab:green', alpha=1)
    ax3.view_init(azim=45, elev=45)
    ax3.set_xlabel('f1', fontsize=14, labelpad=2)
    ax3.set_ylabel('f2', fontsize=14, labelpad=2)
    ax3.set_zlabel('f3', fontsize=14, labelpad=2)

    ax2.plot_trisurf(*zip(*pf), alpha=1, color=pf_color)
    ax2.scatter3D(*zip(*np.array(int_pts)), facecolor=pt_color, s=pt_size)
    ax2.plot_trisurf(*zip(*pf_slack), color='tab:green', alpha=0.2)

    ax2.view_init(azim=225, elev=-45)

    ax2.set_xlabel('f1', fontsize=14, labelpad=2)
    ax2.set_ylabel('f2', fontsize=14, labelpad=2)
    ax2.set_zlabel('f3', fontsize=14, labelpad=2)

    # turn off tick labels
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.axes.zaxis.set_ticklabels([])

    ax2.axes.xaxis.set_ticklabels([])
    ax2.axes.yaxis.set_ticklabels([])
    ax2.axes.zaxis.set_ticklabels([])

    ax3.axes.xaxis.set_ticklabels([])
    ax3.axes.yaxis.set_ticklabels([])
    ax3.axes.zaxis.set_ticklabels([])

    plt.tight_layout()
    plt.savefig("../docs/figures/3d-mga-paretofront.pgf")
