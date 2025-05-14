import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib as mpl
mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"

if __name__ == "__main__":
    demand_df = pd.read_csv(snakemake.input.load_data, 
                        index_col='time', 
                        parse_dates=True,
                        usecols=['time', 'kw'])
    
    demand_df[demand_df['kw'] < 32e3] = np.nan
    demand_df = demand_df.interpolate('linear')
    demand_df['kw'] = demand_df['kw'] / demand_df['kw'].max()

    grouped = demand_df.groupby(demand_df.index.year)
    frames = []
    for g in grouped.groups:
        df = grouped.get_group(g)
        df.reset_index(inplace=True, drop=True)
        frames.append(df)
    df2 = pd.concat(frames, axis=1)
    df2['demand_mean'] = df2.mean(axis=1)

    df3 = demand_df.copy()
    df3['interval'] = 1
    df3 = df3.sort_values(by=['kw'], ascending=False)
    df3['duration'] = df3['interval'].cumsum()
    df3['percentage'] = df3['duration']*100/len(df3)

    df2_average = df2.copy()
    df2_average['interval'] = 1
    df2_sorted = df2_average.sort_values(by=['demand_mean'], ascending=False)
    df2_sorted['duration'] = df2_sorted['interval'].cumsum()
    df2_sorted['percentage'] = df2_sorted['duration']*100/len(df2_sorted)
    df2_sorted['demand_mean'] = df2_sorted['demand_mean']/df2_sorted['demand_mean'].max()


    fig, axes = plt.subplots(1,2,figsize=(12,6), gridspec_kw={'width_ratios':[2.5,1]}, sharey=True, 
                            edgecolor='k', facecolor='w')
    ax1 = axes[0]
    ax2 = axes[1]
    (df2[['demand_mean']]/df2['demand_mean'].max()).plot(ax=ax1, legend=False)
    sb.lineplot(ax=ax2,x='percentage', y='demand_mean', data=df2_sorted, color='tab:red')

    ax1.minorticks_on()
    ax2.minorticks_on()
    ax1.grid()
    ax2.grid()
    ax1.set_ylim(min(df2_sorted['demand_mean']), max(df2_sorted['demand_mean']))
    ax1.set_xlim(0,8760)
    ax1.set_title('Normalized Demand Curve', fontsize=20)
    ax2.set_title('Load Duration Curve', fontsize=20)
    ax1.tick_params(axis='both', labelsize=14, rotation=0)
    ax2.tick_params(axis='both', labelsize=14, rotation=0)
    ax1.set_ylabel('Demand (--)', fontsize=18)
    ax1.set_xlabel('Time [hours]', fontsize=18)
    ax2.set_xlabel(r'Time [\%]', fontsize=18)
    ax2.set_xlim(-1,101)
    plt.tight_layout()
    plt.savefig(snakemake.output.load_curve)