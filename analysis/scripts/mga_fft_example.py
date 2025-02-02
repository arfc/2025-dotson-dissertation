import numpy as np
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
import matplotlib as mpl
import matplotlib.patches as patches
import pandas as pd
import itertools as it
# pymoo imports
from pymoo.problems import get_problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

from osier import distance_matrix, farthest_first, check_if_interior, apply_slack


mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"

# ===============PROBLEM SET UP==================
problem = get_problem("bnh")

seed = 1234
pop_size = 100
n_gen = 200
algorithm = NSGA2(pop_size=pop_size)

res = minimize(problem,
               algorithm,
               ('n_gen', n_gen),
               seed=seed,
               verbose=False,
               save_history=True
               )

PF = problem.pareto_front()
F = res.F
a = min(F[:, 0])
b = max(F[:, 0])
f1 = PF[:, 0]
f2 = PF[:, 1]
slack = 0.2
alpha = 0.5
F1 = f1 * (1 + slack)
F2 = f2 * (1 + slack)

X_hist = np.array([history.pop.get("X")
                  for history in res.history]).reshape(n_gen * pop_size, 2)
F_hist = np.array([history.pop.get("F")
                  for history in res.history]).reshape(n_gen * pop_size, 2)

slack_front = apply_slack(F, slack=slack)

X_hist = np.array([history.pop.get("X")
                  for history in res.history]).reshape(n_gen * pop_size, 2)
F_hist = np.array([history.pop.get("F")
                  for history in res.history]).reshape(n_gen * pop_size, 2)

int_pts = check_if_interior(
    points=F_hist,
    par_front=F,
    slack_front=slack_front)
X_int = X_hist[int_pts]
F_int = F_hist[int_pts]

D = distance_matrix(X=X_int)

n_hist = 1000

n_pts = 6
idxs = farthest_first(X=X_int, D=D, n_points=n_pts, seed=seed)

F_select = F_int[idxs]
X_select = X_int[idxs]

F_df = pd.DataFrame(dict(zip(['f0', 'f1'], F_select.T)))
X_df = pd.DataFrame(dict(zip(['x0', 'x1'], X_select.T)))
mga_df = pd.concat([F_df, X_df], axis=1)


cmap_name = 'tab10'
color1 = mcp.gen_color(cmap=cmap_name, n=n_pts)
cmap = plt.get_cmap(cmap_name, n_pts)

# with plt.style.context('dark_background'):
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

ax[0].set_title("Objective Space", fontsize=16)
ax[0].scatter(F_hist[:n_hist, 0], F_hist[:n_hist, 1], facecolor="none",
              edgecolor="lightgray", alpha=1, label='tested points')
ax[0].scatter(res.F[:, 0], res.F[:, 1], facecolor="none",
              edgecolor="red", label='optimal points')
# ax[0].plot(PF[:,0], PF[:,1], color="g", alpha=1, label='Pareto front', lw=2)
ax[0].scatter(F_select[:, 0], F_select[:, 1],
              c=color1, s=80, label='selected points')

ax[0].fill(np.append(f1, F1[::-1]), np.append(f2, F2[::-1]),
           'lightgrey', alpha=0.3, label="Near-optimal space")
ax[0].legend(fontsize=14)


ax[1].set_title("Design Space", fontsize=16)
ax[1].scatter(X_int[:n_hist, 0], X_int[:n_hist, 1], facecolor="none",
              edgecolor="lightgray", alpha=1, label='tested points')
ax[1].scatter(res.X[:, 0], res.X[:, 1], facecolor="none",
              edgecolor="red", label='optimal points')
ax[0].tick_params(labelsize=14)
ax[1].tick_params(labelsize=14)
ax[1].set_xlim(0, 5)
ax[1].set_ylim(0, 3)
ax[0].set_xlim(0, 160)
ax[0].set_ylim(min(res.F[:, 1]), 60)
# ax[1].legend(loc='upper left')

style = "Simple, tail_width=0.5, head_width=10, head_length=15"
arrows = []
prev = X_select[0]
for i, (c, (x, y)) in enumerate(zip(color1, X_select)):
    ax[1].scatter(x, y, color=c, s=100)
    if i == 0:
        pass
    else:
        kw = dict(arrowstyle=style, color=c, linewidth=3)
        curr = (x, y)
        arrows.append(patches.FancyArrowPatch(prev, curr, **kw))
        prev = curr

for a in arrows:
    ax[1].add_patch(a)

plt.tight_layout()
plt.savefig("../docs/figures/mga-fft-example.pgf")
