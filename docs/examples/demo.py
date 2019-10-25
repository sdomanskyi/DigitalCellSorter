import os
import pandas as pd
import numpy as np
import DigitalCellSorter
from DigitalCellSorter import DigitalCellSorter as DigitalCellSorterSubmodule
from DigitalCellSorter import ReadPrepareDataHCApreviewDataset as HCAtools
    
if __name__ == '__main__':

    DCS = DigitalCellSorterSubmodule.DigitalCellSorter()

    DCS.dataName = 'BM1'
    DCS.saveDir = os.path.join('output', DCS.dataName, '')
    DCS.geneListFileName = os.path.join('geneLists', 'CIBERSORT.xlsx')
    DCS.nClusters = 20

    HCAtools.PrepareDataOnePatient(os.path.join('data', 'ica_bone_marrow_h5.h5'), 'BM1', os.path.join('data', ''))

    df_expr = pd.read_hdf(os.path.join('data', DCS.dataName + '.h5'), key=DCS.dataName, mode='r')

    df_expr = DCS.prepare(df_expr)

    DCS.process(df_expr)

    DCSsub = DigitalCellSorterSubmodule.DigitalCellSorter(dataName=DCS.dataName, nClusters=10, doQualityControl=False)

    DCSsub.subclusteringName = 'T cell'
    DCSsub.saveDir = os.path.join('output', DCS.dataName, 'subclustering T cell', '')
    DCSsub.geneListFileName = os.path.join('geneLists', 'CIBERSORT_T_SUB.xlsx')

    DCSsub.process(df_expr[DCS.getCells(celltype='T cell')])

    DCSsub.subclusteringName = 'B cell'
    DCSsub.saveDir = os.path.join('output', DCS.dataName, 'subclustering B cell', '')
    DCSsub.geneListFileName = os.path.join('geneLists', 'CIBERSORT_B_SUB.xlsx')

    DCSsub.process(df_expr[DCS.getCells(celltype='B cell')])