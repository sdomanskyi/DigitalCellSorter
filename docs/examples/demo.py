import os
import pandas as pd
import numpy as np

import DigitalCellSorter
from DigitalCellSorter.ReadPrepareDataHCApreviewDataset import PrepareDataOnePatient
    
if __name__ == '__main__':

    # Create an instance of class DigitalCellSorter. Here we use Default parameter values
    DCS = DigitalCellSorter.DigitalCellSorter()

    # Modify some of the DCS attributes
    DCS.dataName = 'BM1'
    DCS.saveDir = os.path.join(os.path.dirname(__file__), 'output', DCS.dataName, '')
    DCS.nClusters = 20

    # Call function PrepareDataOnePatient to create file BM1.h5 (HDF file of input type 3) in folder 'data'
    PrepareDataOnePatient(os.path.join('data', 'ica_bone_marrow_h5.h5'), 'BM1', os.path.join('data', ''))

    # Load the expression data from h5 file
    df_expr = pd.read_hdf(os.path.join('data', DCS.dataName + '.h5'), key=DCS.dataName, mode='r')

    # Validate the expression data, so that it has correct form
    df_expr = DCS.prepare(df_expr)

    # Process the expression data
    DCS.process(df_expr)


    # Further analysis can be done on cell types of interest, e.g. here 'T cell' and 'B cell'.
    # Let's create a new instance of DigitalCellSorter to run "sub-analysis" with it.
    # It is important to disable Quality control, because the low quality cells have already been identified and filtered with DCS.
    # Parameter dataNam points to the location processed with DCS. 
    DCSsub = DigitalCellSorterSubmodule.DigitalCellSorter(dataName=DCS.dataName, nClusters=10, doQualityControl=False)

    # Modify a few other attributes
    DCSsub.subclusteringName = 'T cell'
    DCSsub.saveDir = os.path.join(os.path.dirname(__file__), 'output', DCS.dataName, 'subclustering T cell', '')
    DCSsub.geneListFileName = os.path.join(os.path.dirname(__file__), 'docs', 'examples', 'CIBERSORT_T_SUB.xlsx')

    # Process subtype 'T cell'
    DCSsub.process(df_expr[DCS.getCells(celltype='T cell')])

    
    # We can reuse DCSsub to analyze cell type 'B cell'. Just modify the attributes
    DCSsub.subclusteringName = 'B cell'
    DCSsub.saveDir = os.path.join(os.path.dirname(__file__), 'output', DCS.dataName, 'subclustering B cell', '')
    DCSsub.geneListFileName = os.path.join(os.path.dirname(__file__), 'docs', 'examples', 'CIBERSORT_B_SUB.xlsx')

    # Process subtype 'B cell'
    DCSsub.process(df_expr[DCS.getCells(celltype='B cell')])