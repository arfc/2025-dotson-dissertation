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
import seaborn as sb
import pandas as pd

from pymoo.util.ref_dirs import get_reference_directions
from pymoo.visualization.scatter import Scatter

if __name__ == "__main__":

    ref_dirs = get_reference_directions("energy", 3, 90, seed=1)

    x = ref_dirs[:,0]
    y = ref_dirs[:,1]
    z = ref_dirs[:,2]

    from mpl_toolkits import mplot3d
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection='3d')

    ax.scatter3D(x, y, z, alpha=1, s=50);
    ax.set_zlabel("$f_3$", fontsize=14)
    ax.set_ylabel("$f_2$", fontsize=14)
    ax.set_xlabel("$f_1$", fontsize=14)
    ax.view_init(45,45)
    plt.savefig(snakemake.output.reference_directions)