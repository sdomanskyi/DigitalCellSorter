import pandas as pd

df_expr = pd.read_hdf('PBMC.h5', key='PBMC', mode='r')

print(df_expr)