import numpy as np
import pandas as pd
import scripts.DigitalCellSorter as DigitalCellSorter
import scripts.ReadPrepareDataHCApreviewDataset as HCA
import shutil
import os


if __name__ == '__main__':

    typeToProcess = 'BM' #'PBMC' 'BM'

    if typeToProcess=='BM':
        patient = 1
        dataName = 'HCA_BM%s_data' % (patient)
        prepareData = True
        AllData=True
    elif typeToProcess=='PBMC':
        dataName = 'filtered_matrices_mex'
        prepareData = False

    n_clusters = 8

    AvailableCPUsCount = 11
    data_Folder = 'data'
    N_samples_for_distribution = 10000
    cellTypeNameForSubclustering = None #None  'B cells'  'T cells'

    if cellTypeNameForSubclustering == 'B cells':
        clusterIndex = [0]
        geneListToUse = 'geneLists/CIBERSORT_B_SUB.xlsx'
    elif cellTypeNameForSubclustering == 'T cells':
        clusterIndex = [2,6,7]
        geneListToUse = 'geneLists/CIBERSORT_T_SUB.xlsx'
    else:
        clusterIndex = None
        geneListToUse = 'geneLists/CIBERSORT.xlsx'
    
    if prepareData:
        HCA.PrepareData('data/ica_bone_marrow_h5.h5', data_Folder, patient, useAllData=AllData, cellsLimitToUse=1000)
        
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

    if prepareData:
        pathToRemove = data_Folder + '/' + dataName
        if os.path.exists(pathToRemove):
            shutil.rmtree(pathToRemove)

    DigitalCellSorter.DigitalCellSorter().Process(df_expr, 
                                                    dataName, 
                                                    saveDir = 'demo_output/' + dataName + '/', 
                                                    geneListFileName = geneListToUse,
                                                    N_samples_for_distribution = N_samples_for_distribution,
                                                    AvailableCPUsCount = AvailableCPUsCount,
                                                    clusterIndex=clusterIndex,
                                                    clusterName=cellTypeNameForSubclustering,
                                                    n_clusters=n_clusters)