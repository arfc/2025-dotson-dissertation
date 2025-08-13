if __name__ == "__main__":
    import pandas as pd

    cols = [
        'core_metric_parameter',
        'core_metric_case',
        'tax_credit_case',
        'crpyears',
        'technology',
        'technology_alias',
        'techdetail',
        'techdetail2',
        'resourcedetail',
        'display_name',
        'default',
        'scale',
        'maturity',
        'scenario',
        'core_metric_variable',
        'units',
        'value'
        ]
    tech_list = ['Nuclear',
                 'UtilityPV',
                 'Utility-Scale Battery Storage',
                 'NaturalGas_FE']
    param_list = ['Fixed O&M',
                  'CAPEX',
                  'Variable O&M',]
    final_techs = {'Utility PV - Class 1':'UtilityPV', 
       'NG 1-on-1 Combined Cycle (H-Frame)':"Natural Gas",
       'NG 1-on-1 Combined Cycle (H-Frame) 95% CCS':"Natural Gas CCS",
       'Nuclear - Large':"Nuclear",
       'Nuclear - Small':"Advanced Nuclear",        
       'Utility-Scale Battery Storage - 4Hr':"Battery",
       }
    
    atbe = pd.read_csv(snakemake.input.atbe, usecols=cols, low_memory=False)

    final_df = atbe.loc[((atbe.display_name.isin(final_techs.keys()))
                        &(atbe.core_metric_parameter.isin(param_list))
                        &(atbe.scenario=='Moderate')
                        &(atbe.core_metric_variable==2032)
                        &(~atbe.tax_credit_case.isin(['ITC',
                                                        'PTC']))
                        &(atbe.crpyears == 30)),['core_metric_parameter',
                                                'display_name',
                                                'value']]\
                    .drop_duplicates(subset=['core_metric_parameter',
                                             'value'])\
                    .reset_index(drop=True)
    
    final_df = final_df.replace(final_techs)

    final_df.to_csv(snakemake.output.tech_costs, index=False)