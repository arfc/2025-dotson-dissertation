import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from pandas.tseries.frequencies import to_offset
import matplotlib.dates as mdates

mpl.use("pgf")
plt.rcParams['pgf.texsystem'] = 'pdflatex'
plt.rcParams['text.usetex'] = True
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['figure.edgecolor'] = 'k'
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['savefig.dpi'] = 400
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = "serif"

if __name__ == "__main__":
    time_columns = ['year', 'month','day', 'hour', 'minute', 'second']
    data_columns = ['value']
    df = pd.read_csv(snakemake.input.mauna_loa_data, 
                    skiprows=158, 
    #                  sep='/s', 
                    delimiter=' ',
                    na_values=[-999.99, -99.99],
                    na_filter=True,
                    usecols=time_columns+data_columns)
    df['timestamp'] = pd.to_datetime(df[time_columns])
    df.set_index('timestamp', inplace=True)
    df.drop(columns=time_columns, inplace=True)

    df.dropna(inplace=True)

    df_monC = df.resample('ME', offset='14D').mean().interpolate(method='cubic')
    df_smooth = df_monC.rolling(11).mean()
    

    fig, axes = plt.subplots(1,1,figsize=(8,6))

    start_year = 1985
    end_year = 2024

    # df.value.plot(ax=axes, label='Raw Data', color='b', marker='.', alpha=0.5)
    # df_monC.value.plot(ax=axes, label='Monthly Data', marker='.', color='r')
    # df_smooth.value.plot(ax=axes, label='Smoothed', marker='.', color='k')

    # df[df.index.year>=start_year].value.plot(ax=axes, label='Raw Data', color='b', 
    #                                    marker='.', alpha=0.5)
    df_monC[df_monC.index.year>=start_year].value.plot(ax=axes, label='Monthly Data', 
                                                marker='.', color='r', lw=1)
    df_smooth[df_smooth.index.year>=start_year].value.plot(ax=axes, label='Smoothed', 
                                                    marker='.', color='k', lw=1)


    axes.legend(fontsize=14)



    axes.minorticks_on()

    axes.tick_params(axis='both',which='both',direction='in',rotation=0, 
                        labelsize=14, width=1.25)
    axes.tick_params(axis='both',which='minor',length=4)
    axes.tick_params(axis='both',which='major',length=8)

    axes.set_ylabel('parts per million (ppm)',fontsize=14)

    axes.set_xlabel('Year', fontsize=14)
    axes.grid()

    axes.set_title('Atmospheric CO$_2$ at Mauna Loa Observatory', fontsize=18)

    co2_min = df[df.index.year==start_year].value.min()
    co2_max = df[df.index.year==end_year].value.max()
    axes.set_ylim(co2_min-1, co2_max+1)
    # axes.set_xlim('01-01-2010', '010')
    plt.tight_layout()
    plt.savefig(snakemake.output.mauna_loa_plot)