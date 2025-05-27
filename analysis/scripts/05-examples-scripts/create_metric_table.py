from osier import ThermalTechnology
from dill import load, dump
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

if __name__ == "__main__":
    files = snakemake.input

    frames = []
    for f in files:
        if "demand_data" in f:
            continue
        elif "4_obj" in f:
            continue
        else: 
            metric_name = f.split('/')[-1].replace('.csv', '')
            print(metric_name)
            df = pd.read_csv(f, index_col=0)
            df.index.name = 'EG'
            df = df.iloc[:,:1]
            # get units from column name
            units = df.columns[0].split(" ")[-1].strip('()')
            print(units)
            df.columns = [f"{metric_name.replace("_", " ")} [{units}]"]
            frames.append(df)


    full_df = pd.concat(frames, axis=1)

    full_df.style.format(precision=2).to_latex(snakemake.output.metric_table, hrules=True)

    # full_df.columns = ["_".join(c.split(" ")[:-1]) for c in list(full_df.columns)]
    full_df.to_csv(snakemake.output.metric_data)

