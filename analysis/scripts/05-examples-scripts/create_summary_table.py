import pandas as pd
import numpy as np
import camelot

if __name__ == "__main__":
    mapping = {'once-through':[f'EG0{i}' for i in range(1,9)],
                'limited-recycle':[f'EG0{i}' if i < 10 else f'EG{i}' for i in range(9,19)],
                'continuous-recycle':[f'EG{i}' for i in range(19,41)]
                }
    summary_tables = camelot.read_pdf(snakemake.input[0], pages="19")
    