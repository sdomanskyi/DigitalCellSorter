import os
import pandas as pd
import numpy as np

from DigitalCellSorter import DigitalCellSorter as DigitalCellSorterSubmodule

if __name__ == '__main__':

    #columns = pd.read_hdf('PBMC_processed_CIBERSORT.h5', key='df_markers_expr', mode='r').columns
    #index = pd.read_excel('predictions.xlsx', sheet_name='predictions', header=0, index_col=0).index
    #series = pd.Series(index=columns.get_level_values(1), data=columns.get_level_values(3).str.split(' #', expand=True).get_level_values(0)).loc[index]
    #series.to_excel('se.xlsx')

    #exit()

    source = 'output_TD_coarse_6'

    dataNames = ['BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7', 'BM8', 'PBMC']

    listName = 'CIBERSORT_TD_coarse_6.xlsx'

    # Extract new top-500 markers
    if False:
        for dataName in dataNames:
            DCS = DigitalCellSorterSubmodule.DigitalCellSorter(dataName=dataName, saveDir=os.path.join(source, dataName, ''), geneListFileName=os.path.join('geneLists', listName))
            DCS.extractNewMarkerGenes(top=500)

        for dataName in dataNames:
            saveDir=os.path.join(source, dataName, '')
            df_markers = pd.read_excel(os.path.join(saveDir, 'new_markers.xlsx'), sheet_name='MarkerCellType', index_col='Marker', header=0)
            df_markers.to_hdf('markers.h5', key=dataName, mode='a', complevel=4, complib='zlib')

    # Combine all datasets new top-500 markers
    if False:
        df = pd.DataFrame() 

        for dataName in dataNames:
            df_markers = pd.read_hdf('markers.h5', key=dataName, mode='r')
            df_markers.columns = pd.MultiIndex.from_arrays([[dataName]*len(df_markers.columns), df_markers.columns], names=['batch', 'celltype'])
            print(df_markers.columns, '\n')

            df = pd.concat([df, df_markers], sort=False, axis=1)

        df = df.fillna(0.)
        df = df.groupby(level=['celltype'], axis=1).sum()
        
        print(df)

        df.to_hdf('markers.h5', key='All', mode='a', complevel=4, complib='zlib')

    # Keep those markers ">=3"
    if True:
        df = pd.read_hdf('markers.h5', key='All', mode='r')
        df = df.loc[df.max(axis=1)>=3]

        def set_zeros(series):

            where_max = series.idxmax()
            value = series.loc[where_max]

            series.loc[:] = 0.
            series.loc[where_max] = value

            return series

        df = df.apply(set_zeros, axis=1)

        print((df>=3).sum(axis=0)._values)

        print(df.columns)
        
        df_cut = pd.DataFrame()
        for col in df.columns:
            df_cut = pd.concat([df_cut, df[col][df[col]>0].sort_values(ascending=False)[:50]], sort=False, axis=1)

        df_cut = df_cut.fillna(0.)

        print((df_cut>0).sum(axis=0).values)
        print(df_cut.columns)

        df_cut[df_cut>0] = 1

        df_cut = df_cut.sort_index()

        print(df_cut)

        df_cut.to_excel(source + '.xlsx')
