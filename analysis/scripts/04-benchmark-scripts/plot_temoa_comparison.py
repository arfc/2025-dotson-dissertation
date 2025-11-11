import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
from osier import n_mga, check_if_interior, distance_matrix, farthest_first


def is_pareto_efficient(costs, return_mask = True):
    """
    Code borrowed from this Stackoverflow post:
    https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python

    Find the pareto-efficient points :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask :return: An array of indices of
    pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        nondominated_point_mask = np.any(costs<costs[next_point_index],
                                         axis=1)
        nondominated_point_mask[next_point_index] = True
        # Remove dominated points
        is_efficient = is_efficient[nondominated_point_mask]  
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient

if __name__ == "__main__":

    print("Opening datasets")
    
    with open(snakemake.input.unsga3_results, "rb") as file:
        data = dill.load(file)
    with open(snakemake.input.unsga3_design, "rb") as file:
        osier_X = dill.load(file)
    with open(snakemake.input.temoa_F2, "rb") as file:
        temoa_F2 = dill.load(file)
    with open(snakemake.input.temoa_mga_design, "rb") as file:
        temoa_X = dill.load(file)
    with open(snakemake.input.techs, "rb") as file:
        techs = dill.load(file)


    # Parameters
    n_gen = 95
    pop_size = 100
    n_vars = 10
    n_objs = 2
    threshold = 14.3
    seed = 42


    # Get data history
    X_hist = np.array([hist.pop.get("X") for hist in data.history]).reshape(
        n_gen * pop_size, n_vars)
    F_hist = np.array([hist.pop.get("F") for hist in data.history]).reshape(
        n_gen * pop_size, n_objs)
    
    # get the Pareto Front
    pareto_idxs = is_pareto_efficient(F_hist)
    F = F_hist[pareto_idxs]
    X = X_hist[pareto_idxs]

    df = pd.DataFrame({"Cost":F[:,0], "Carbon":F[:,1]})
    slack = 1.1
    f2 = F[:,1]*slack
    f1 = F[:,0]*slack
    df2 = pd.DataFrame({"Cost":f1, "Carbon":f2})
    F_slack = np.c_[f1,f2]


    # Get interior points and mga points
    int_points = check_if_interior(F_hist, par_front=F, slack_front=F_slack)
    X_interior = X_hist[int_points]
    F_interior = F_hist[int_points]
    X_hist_list = [list(row) for row in X_interior]
    X_df = pd.DataFrame(dict(Designs = X_hist_list))
    F_df = pd.DataFrame(dict(zip(['Cost','Carbon'], F_interior.T)))
    interior_df = pd.concat([F_df, X_df], axis=1)
    interior_df = interior_df.loc[interior_df['Carbon'] < threshold]

    F_outside = np.array([row for row 
                          in F_hist 
                          if row not in interior_df[['Cost','Carbon']].values])
    F_hist_masked = F_outside[(F_outside[:,0] < 8300) & (F_outside[:,1]<1e40)]


    mga_df = n_mga(results_obj=data, 
                   slack=0.1, 
                   n_points=32, 
                   wide_form=False, 
                   seed=90)
    mga_df = mga_df.loc[mga_df['co2'] < threshold]


    ### Plot first comparison ###
    print(f"Plotting {snakemake.output.temoa_comparison_01}")
    color_1 = 'k'
    color_2 = 'grey'
    # color_1 = 'blue' color_2 = 'tab:blue'
    color_3 = 'tab:red'
    alpha = 0.5
    size=35

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
               s=size, 
               facecolors='none', 
               edgecolors='red', 
               label='temoa+mga')

    rng = np.random.default_rng(seed=seed)
    random_indices = rng.choice(list(range(len(F_hist_masked))), size=500)
    ax.scatter(x=F_hist_masked[random_indices,0], 
                y=F_hist_masked[random_indices,1], 
                fc='None', 
                ec=color_2, 
                alpha=1, 
                s=size, 
                label='osier (tested points)'
                )


    ax.set_xlim(min(temoa_F2[:,0]),max(F[:,0])*slack)
    ax.set_ylim(min(F[:,1]),max(temoa_F2[:,1])*slack)
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

    ax.set_xlabel("Total Cost [M\$]", fontsize=18)
    ax.set_ylabel(r"Emissions [MT CO$_2$eq]", fontsize=18)
    # plt.legend(loc='lower right')
    ax.legend(fontsize=16, loc='lower right', framealpha=1)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(snakemake.output.temoa_comparison_01)
    

    ### Second Comparison Plot ### 
    print(f"Plotting {snakemake.output.temoa_comparison_02}")
    
    F_hist_masked = F_outside[(F_outside[:,0] < 8300) & (F_outside[:,1]<16)]

    fig, ax = plt.subplots(figsize=(8,6), facecolor='w')
    
    rng = np.random.default_rng(seed=69)
    random_indices = rng.choice(list(range(len(F_hist_masked))), size=500)
    ax.scatter(x=F_hist_masked[random_indices,0], 
                y=F_hist_masked[random_indices,1], 
                fc='None', 
                ec=color_1, 
                alpha=1, 
                s=size, 
                label='tested points'
                )
    interior_df.plot.scatter(ax=ax, 
                            x='Cost',
                            y='Carbon', 
                            color='tab:green',
                            s=size,
                            label='sub-optimal'
                            )
    mga_df.plot.scatter(ax=ax, 
                        x='total_cost', 
                        y='co2', 
                        legend=True, 
                        color='tab:red', 
                        marker='o', 
                        label='selected points', 
                        s=size
                        )

    df.sort_values(by="Cost").plot.line(ax=ax,
                                        x='Cost', 
                                        y='Carbon',
                                        legend=False, 
                                        color=color_1, 
                                        marker='.', 
                                        label='Pareto-front', 
                                        lw=2.5)
    df2.sort_values(by="Cost").plot.line(ax=ax,
                                        x='Cost', 
                                        y='Carbon',
                                        legend=False, 
                                        color=color_2, 
                                        alpha=alpha,
                                        label='Near-optimal space')

    patch = ax.fill(np.append(df.sort_values(by="Cost").Cost.values, 
                            df2.sort_values(by="Cost").Cost.values[::-1]),
        np.append(df.sort_values(by="Cost").Carbon.values, 
                    df2.sort_values(by="Cost").Carbon.values[::-1]),
                    alpha=alpha, 
                    color=color_2)

    ax.set_xlim(min(F[:,0]),max(F[:,0])*slack)
    ax.set_ylim(min(F[:,1]),max(F[:,1])*slack)
    ax.grid()

    ax.set_xlabel("Total Cost [M\$]", fontsize=18)
    ax.set_ylabel(r"Emissions [MT CO$_2$eq]", fontsize=18)
    ax.legend(fontsize=16, loc='upper right', framealpha=1)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(snakemake.output.temoa_comparison_02)


    ### Plot design comparison ###
    print(f"Plotting {snakemake.output.temoa_comparison_03}")

    labels = [t.technology_name for t in techs]
    labels = [l.replace('_', '\_') for l in labels]
    peak_demand = 29.251353638620195

    hex_colors = []
    for c in plt.cm.tab10.colors:
        hex_colors.append(mpl.colors.to_hex(c))

    X_select = np.array([np.array(row) for row in mga_df.designs.values])

    fig, axes = plt.subplots(3,1,figsize=(12,14), facecolor='w', sharex=True, sharey=True)
    sb.boxenplot(ax=axes[0], data=(osier_X)*peak_demand, k_depth='trustworthy')
    sb.boxenplot(ax=axes[1], data=(X_select)*peak_demand, k_depth='trustworthy')
    sb.boxenplot(ax=axes[2], data=(temoa_X), k_depth='trustworthy')

    axes[0].minorticks_on()
    axes[1].minorticks_on()
    axes[2].minorticks_on()
    axes[0].set_ylabel('Osier Capacity (GW)', size=18)
    axes[1].set_ylabel('Osier MGA Capacity (GW)', size=18)
    axes[2].set_ylabel('Temoa MGA Capacity (GW)', size=18)

    axes[0].grid(alpha=0.2, which='major',zorder=5)
    axes[0].grid(alpha=0.05, which='minor',zorder=5)
    axes[1].grid(alpha=0.2, which='major',zorder=5)
    axes[2].grid(alpha=0.2, which='major',zorder=5)
    axes[1].grid(alpha=0.05, which='minor',zorder=5)
    axes[2].grid(alpha=0.05, which='minor',zorder=5)
    axes[0].set_ylim(0,27)
    axes[1].set_ylim(0,27)
    axes[2].set_ylim(0,27)

    axes[2].set_xticklabels(labels, rotation=12.5, size=12)
    axes[2].set_xticks(range(len(labels)))
    axes[2].set_xlabel("", size=14)

    x_loc, y_loc = -.31, 24.95
    axes[0].text(x_loc,y_loc, "a)", fontsize=14, bbox=dict({'facecolor':'w'}))
    axes[1].text(x_loc,y_loc, "b)", fontsize=14, bbox=dict({'facecolor':'w'}))
    axes[2].text(x_loc,y_loc, "c)", fontsize=14, bbox=dict({'facecolor':'w'}))

    patch = mpatches.Patch(color='grey', label='manual patch')
    patches = []
    for color, label in zip(hex_colors, labels):
        patches.append(mpatches.Patch(color=color, label=label))

    handles, llabels = axes[0].get_legend_handles_labels()
    handles.extend(patches)

    axes[2].legend(handles=handles, ncol=3, title='Technologies', fontsize=14, title_fontsize=14)

    plt.tight_layout()
    plt.savefig(snakemake.output.temoa_comparison_03)