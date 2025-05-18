import pandas as pd
import camelot
import numpy as np



if __name__ == "__main__":
    table_pages = snakemake.config['table_pages']
    for name, pages in table_pages.items():
        print(f"processing {name}")
        tables = camelot.read_pdf(snakemake.input[0], pages=pages)
        print(tables)
        frames = []
        col_names = []
        for i, t in enumerate(tables):
            frame = t.df.copy()
            if frame.shape[1] == 4:
                # make sure no empty values
                if frame.replace("", np.nan).isnull().all().all():
                    continue
                else:
                    # check values of column names
                    if "Metric \nData" in frame.iloc[0,:].values:
                        # convert first row into column names
                        col_names = frame.iloc[0,:].values
                        frame.columns = col_names
                        # drop first row
                        frame = frame.iloc[1:,]
                        frames.append(frame)
                    else: 
                        # in theory, this should always follow a table with column names
                        frame.columns = col_names
                        frames.append(frame)
            elif frame.shape[1] == 5:
                # convert first row into column names
                frame.columns = frame.iloc[0,:].values
                # drop first row
                frame = frame.iloc[1:,]
                frames.append(frame)
            else:
                continue
        full_df = pd.concat(frames)
        full_df.columns = [c.replace("\n", "") for c in full_df.columns]
        try:
            full_df = full_df.replace("", np.nan).dropna(axis=0, 
                                                    how='any', 
                                                    subset='EG').reset_index(drop=True)
        except KeyError:
            full_df = full_df.replace("", np.nan).dropna(axis=0,how='all').reset_index(drop=True)   
            
        full_df.iloc[:,1] = full_df.iloc[:,1].astype(float)
        full_df.to_csv(f"../data/{name}.csv", index=False)