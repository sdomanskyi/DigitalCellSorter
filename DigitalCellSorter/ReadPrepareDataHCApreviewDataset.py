'''The user can download the file ica_bone_marrow_h5.h5 from https://preview.data.humancellatlas.org/ (Raw Counts Matrix - Bone Marrow) and place in folder data. 
The file is ~485Mb and contains all 378000 cells from 8 bone marrow donors (BM1-BM8).
'''

import os
import copy
import h5py
import numpy as np
import pandas as pd

def PrepareDataOnePatient(filename, patient, saveFolderName,  
                          useAllData=True, cellsLimitToUse=1000, randomlySample=True, randomSeed=0):

    '''Prepare data from Human Cell Atlas (HCA) h5 data file.
        
    Parameters:
        filename: str
            Path and name of the file to store binary data in

        patient: str
            Identifier of the patient: 'BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7' or 'BM8'

        saveFolderName: str
            Path where to save prepared data file

        useAllData: boolean, Default True
            Whether to use all data or a subset

        cellsLimitToUse: int, Default 1000
            Number of cells to use if useAllData=False

        randomlySample: boolean, Default True
            Whether to sample cell randomly of pick top number

        randomSeed: int, Default 0
            Random seed

    Returns:
        None
        
    Usage:
        PrepareDataOnePatient(os.path.join('data', 'ica_bone_marrow_h5.h5'), 'BM1', os.path.join('data', ''), useAllData=False, cellsLimitToUse=5000)
    '''

    print('------ Reading values of %s --------------------------------------------'%(patient))

    #Open and Prepare data file
    h5File = h5py.File(filename, 'r')
    barcodes = np.array(list(map(lambda s: s.decode('UTF-8'), copy.deepcopy(h5File['GRCh38/barcodes'][()]))))
    gene_names = np.array(list(map(lambda s: s.decode('UTF-8'), copy.deepcopy(h5File['GRCh38/gene_names'][()]))))
    indptr = copy.deepcopy(h5File['GRCh38/indptr'][()])
    patients = np.array([cell[6:9] for cell in barcodes])
    patient_cells = np.where(patients==patient)[0]
        
    if not useAllData:
        np.random.seed(randomSeed)
        patient_cells = np.random.choice(patient_cells, replace=False, size=cellsLimitToUse) if randomlySample else patient_cells[:cellsLimitToUse]

    columns = pd.MultiIndex.from_arrays([patients[patient_cells], barcodes[patient_cells]], names=['batch', 'cell'])
    df = pd.DataFrame(data=np.zeros((len(gene_names), len(patient_cells))), index=gene_names, columns=columns)

    print('Combining cells into DataFrame. Data size: %s genes, %s cells' % df.shape)
    for cell, index in enumerate(patient_cells):
        df.iloc[copy.deepcopy(h5File['GRCh38/indices'][indptr[index]:indptr[index + 1]]), cell] = copy.deepcopy(h5File['GRCh38/data'][indptr[index]:indptr[index + 1]])
    
    print('Compressing %s to Pandas HDF'%(patient))
    df.to_hdf(os.path.join(saveFolderName, '%s.h5'%(patient)), key=patient, mode='a', complevel=4, complib='zlib')

    print('-------------------------------------------------------------------------\n')

    return None

