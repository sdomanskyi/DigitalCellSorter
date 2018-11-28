import numpy as np
import pandas as pd
import scripts.DigitalCellSorter as DigitalCellSorter
import scripts.ReadPrepareDataHCApreviewDataset as HCA
import shutil
import os

if __name__ == '__main__':

    patient = 1
    n_clusters = 8
    AllData=True
    AvailableCPUsCount = 27
    data_Folder = 'data'
    N_samples_for_distribution = 10000 #10000
    attemptLoadingSavedTransformedData = True
    cellTypeNameForSubclustering = None #None  'B cells'  'T cells'
    SaveTransformedData = True if cellTypeNameForSubclustering == None else False
    SaveXpcaDataCSV = True if cellTypeNameForSubclustering == None else False

    if cellTypeNameForSubclustering == 'B cells':
        clusterIndex = [1]
        geneListToUse = 'geneLists/G_9_BM_B_cells_subtypes.xlsx'
    elif cellTypeNameForSubclustering == 'T cells':
        clusterIndex = [3,5]
        geneListToUse = 'geneLists/G_10_BM_T_cells_subtypes.xlsx'
    else:
        clusterIndex = None
        geneListToUse = 'geneLists/G_2_Human_cell_markers_BM.xlsx'
    
    HCA.PrepareData('data/ica_bone_marrow_h5.h5', data_Folder, patient, useAllData=AllData, cellsLimitToUse=1000)
                    
    dataName = 'HCA_BM%s_data' % (patient)

    print('\n======================\nDone preparing raw data: %s!\n======================'%(dataName))

    print("\nLoading data of " + dataName)
    dir_expr = '/'.join(['data',dataName,'matrix.mtx'])
    dir_geneNames = '/'.join(['data',dataName,'genes.tsv'])

    with open(dir_expr) as myfile:
        ijv = myfile.readlines()

    header = ijv[:3]

    ijv = np.vstack([np.int_(i.strip('\n').split(' ')) for i in ijv[3:]])
    ijv[:,:2] -= 1

    imax,jmax = np.int_(header[-1].split(' ')[:2])
    df_expr = np.zeros([imax,jmax])

    df_expr[ijv[:,0],ijv[:,1]] = ijv[:,2]
    df_expr = pd.DataFrame(df_expr,index=pd.read_csv(dir_geneNames,delimiter='\t',header=None).values[:,1])
    print('\n======================\nDone loading raw data!\n======================')

    pathToRemove = data_Folder + '/' + dataName
    if os.path.exists(pathToRemove):
        shutil.rmtree(pathToRemove)

    DigitalCellSorter.DigitalCellSorter().Process(df_expr, 
                                                    dataName, 
                                                    saveDir = 'demo_output/' + dataName + '/', 
                                                    geneListFileName = geneListToUse,
                                                    N_samples_for_distribution = N_samples_for_distribution,
                                                    SaveTransformedData = SaveTransformedData,
                                                    attemptLoadingSavedTransformedData = attemptLoadingSavedTransformedData,
                                                    SaveXpcaDataCSV = SaveXpcaDataCSV,
                                                    AvailableCPUsCount = AvailableCPUsCount,
                                                    clusterIndex=clusterIndex,
                                                    clusterName=cellTypeNameForSubclustering,
                                                    n_clusters=n_clusters)