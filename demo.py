import os
import pandas as pd
import numpy as np

import DigitalCellSorter
from DigitalCellSorter.ReadPrepareDataHCApreviewDataset import PrepareDataOnePatient
    
if __name__ == '__main__':

    DCS = DigitalCellSorter.DigitalCellSorter(nSamplesDistribution=2000, nClusters=20)
    dataName = 'BM1'
    DCS.dataName = dataName
    DCS.saveDir = os.path.join('output_TD_coarse_6', dataName, '')
    DCS.geneListFileName = os.path.join('geneLists', 'CIBERSORT_TD_coarse_6.xlsx')
    print(DCS.geneListFileName)

    PrepareDataOnePatient(os.path.join('data', 'ica_bone_marrow_h5.h5'), dataName, os.path.join('data', ''), useAllData=False, cellsLimitToUse=1000)
    df_expr = pd.read_hdf(os.path.join('data', dataName + '.h5'), key=dataName, mode='r')
    df_expr.columns.names = ['batch', 'cell']

    DCS.process(DCS.prepare(df_expr))