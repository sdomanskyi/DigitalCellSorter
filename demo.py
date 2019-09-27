import os
import pandas as pd
import numpy as np
import DigitalCellSorter.DigitalCellSorter as DigitalCellSorter
from DigitalCellSorter.ReadPrepareDataHCApreviewDataset import PrepareDataOnePatient as HCAtool
    
if __name__ == '__main__':

    for dataName in ['BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7', 'BM8']:

        HCAtool(os.path.join('data', 'ica_bone_marrow_h5.h5'), dataName, os.path.join('data', ''), useAllData=False, cellsLimitToUse=1000)

        df_expr = pd.read_hdf(os.path.join('data', dataName + '.h5'), key=dataName, mode='r')

        DCS = DigitalCellSorter.DigitalCellSorter(nSamplesDistribution=10000, nClusters=20)
        DCS.dataName = dataName
        DCS.saveDir = os.path.join('output', dataName, '')
        DCS.geneListFileName = os.path.join('geneLists', 'CIBERSORT.xlsx')
        DCS.toggleRecordAllExpression = True
        DCS.process(df_expr)

        DCSsub = DigitalCellSorter.DigitalCellSorter(dataName=DCS.dataName, nClusters=10, doQualityControl=False, nSamplesDistribution=10000)

        DCSsub.subclusteringName = 'T cell'
        DCSsub.saveDir = os.path.join('output', dataName, 'subclustering T cell', '')
        DCSsub.geneListFileName = os.path.join('geneLists', 'CIBERSORT_T_SUB.xlsx')
        DCSsub.process(df_expr[DCS.getCells(celltype='T cell')])

        DCSsub.subclusteringName = 'B cell'
        DCSsub.saveDir = os.path.join('output', dataName, 'subclustering B cell', '')
        DCSsub.geneListFileName = os.path.join('geneLists', 'CIBERSORT_B_SUB.xlsx')
        DCSsub.process(df_expr[DCS.getCells(celltype='B cell')])