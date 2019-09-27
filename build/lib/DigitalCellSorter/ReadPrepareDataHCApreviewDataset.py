import copy
import h5py
import os
import numpy as np
import pandas as pd

def PrepareDataOnePatient(filename, patient, saveFolderName, doQC=True,  useAllData=True, cellsLimitToUse=1000, randomlySample=True, randomSeed=0):

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

    columns = pd.MultiIndex.from_arrays([patients[patient_cells], barcodes[patient_cells]], names=['patient', 'cell'])
    df = pd.DataFrame(data=np.zeros((len(gene_names), len(patient_cells))), index=gene_names, columns=columns)

    print('Combining cells into DataFrame. Data size: %s genes, %s cells' % df.shape)
    for cell, index in enumerate(patient_cells):
        df.iloc[copy.deepcopy(h5File['GRCh38/indices'][indptr[index]:indptr[index + 1]]), cell] = copy.deepcopy(h5File['GRCh38/data'][indptr[index]:indptr[index + 1]])
    
    print('Compressing %s to Pandas HDF'%(patient))
    df.to_hdf(os.path.join(saveFolderName, 'HCA_%s.h5'%(patient)), key=patient, mode='a', complevel=4, complib='zlib')

    print('-------------------------------------------------------------------------\n')

