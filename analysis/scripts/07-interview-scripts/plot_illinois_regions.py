import pandas as pd
import seaborn as sb
import matplotlib as mpl
from mycolorpy import colorlist as mcp
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import us


mpl.use("pgf") 
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True 
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.shadow'] = False

if __name__ == "__main__":
    # gdf = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2024/shp/cb_2024_us_state_20m.zip")

    counties = gpd.read_file(snakemake.input.illinois_data)


    ax = counties.loc[counties['STATE_NAME'] == 'Illinois'].plot(ec='k',alpha=0.2)
    counties.loc[((counties['NAME'] == 'Champaign') & (counties['STATE_NAME']=='Illinois'))].plot(ax=ax, color='tab:blue', ec='k')
    counties.loc[((counties['NAME'].isin(['Will','DuPage'])) & (counties['STATE_NAME']=='Illinois'))].plot(ax=ax, color='tab:blue', ec='k')
    ax.scatter(y=40.15, x=-88.1437, color='tab:red', label="Champaign-Urbana")
    ax.scatter(y=41.75, x=-88.16, color='gold', label="Naperville")
    ax.legend(loc=(0, 1.01), framealpha=0)
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(snakemake.output.illinois_plot)