import os
import copy
import h5py
import numpy as np
import pandas as pd

from scipy import io
import urllib.request

from .GenericFunctions import write, read, extractFromZipOfGz

def getHCAdataByURL(url, extractPath, extractData = True):

    '''Download and extract data from Human Cell Atlas Portal

    Parameters:
        url: str
            URL of the data of interest

        extractPath: str
            Path where to save and extract data to

        extractData: boolean, Default True
            Whether to extract downloaded data

    Returns:
        None
        
    Usage:
        getHCAdataByURL(url, extractPath)
    '''

    # Create directories for extract path
    if not os.path.exists(extractPath):
        print('Creating directories:\n', extractPath)

        os.makedirs(extractPath)
    else:
        print('Extraction directory already exists. To re-extract files remove the directory:\n', extractPath)

        return

    # Download the *.gz file
    if not os.path.isfile(os.path.join(os.path.dirname(extractPath), os.path.basename(url))):
        print('Downloading file:', url)
        urllib.request.urlretrieve(url, os.path.join(os.path.dirname(extractPath), os.path.basename(url)))

    if extractData:
        # Extract file to extract path
        extractFromZipOfGz(extractPath + '.zip', removeDownloadedZipFile=True)
    
        for tempName in ['cells', 'genes', 'barcodes']:
            write(pd.read_csv(os.path.join(extractPath, tempName + '.tsv'), delimiter='\t', index_col=0, header=0), 
                    os.path.join(extractPath, tempName))

    return

def recordFilesOfIndividualDonors(extractPath, organName=None, donorIDcolumn='donor_organism.provenance.document_id', 
                                          organColumn='derived_organ_parts_label', useHogoGeneNames=True):

    '''Record h5 files of HCA individual donors in a dataset

    Parameters:
        extractPath: str
            Path of directories where HCA matrix files were downloaded and extracted.
            See function getHCAdataByURL() for detail.

        organName: str, Default None
            Name of the organ name. E.g. 'pancreas', 'bone marrow', etc.

        donorIDcolumn: str, Default donor_organism.provenance.document_id
            Column with unique IDs of donors in the file. Another option is
            'specimen_from_organism.provenance.document_id' IDs at samples level
            is needed.

        organColumn: str, Default 'derived_organ_parts_label'
            'derived_organ_label'

            'derived_organ_parts_label'
            This option is ignored when organName parameter is None.

        useHogoGeneNames: boolean, Default True
            Whether to use HUGO gene names.

    Returns:
        list
            List of donor IDs

    Usage:
        recordFilesOfIndividualDonors(extractPath, organName='retina')
    '''

    print('extractPath:', extractPath)

    if not os.path.isfile(os.path.join(extractPath, 'cells' + '.pklz')):

        print('Data files not found. Download and prepare them by calling function getHCAdataByURL()')

        return

    cells = read(os.path.join(extractPath, 'cells'))

    if organName is None:
        subsetOfCells = cells
    else:
        subsetOfCells = cells.iloc[np.where(cells[organColumn]==organName)[0]]
    
    donorIDs = np.unique(subsetOfCells[donorIDcolumn].values).tolist()

    dataSparseMatrix = None
    
    for donorID in donorIDs:
        print('Donor:', donorID)

        if not os.path.isfile(os.path.join(extractPath, 'dfDonorID %s'%(donorID) + '.h5')):

            if dataSparseMatrix is None:
                print('Reading data from *.mtx file')
                dataSparseMatrix = io.mmread(os.path.join(extractPath, 'matrix.mtx')).tocsr()

            # Identifiers of cells of single Bone Marrow donor 
            if organName is None:
                idxDonorID = np.where(cells[donorIDcolumn]==donorID)[0]
            else:
                idxDonorID = np.where((cells[organColumn]==organName) & (cells[donorIDcolumn]==donorID))[0]

            # Ensembl to HUGO conversion
            genes = read(os.path.join(extractPath, 'genes'))['featurename']

            # Read all data and keep only selected cells
            dfDonorID = pd.DataFrame(data=dataSparseMatrix[:,idxDonorID].todense(),
                                index=genes.index if not useHogoGeneNames else genes.values,
                                columns=pd.MultiIndex.from_arrays([[donorID]*len(idxDonorID), cells.index[idxDonorID].values]))

            print(dfDonorID.shape)

            # Merge duplicates that might have appeared after gene name conversion
            dfDonorID = dfDonorID.groupby(level=0, axis=0).sum()

            print(dfDonorID.shape)

            # Drop all-zero genes
            dfDonorID = dfDonorID.loc[dfDonorID.sum(axis=1)>0.]

            # Write data of single Bone Marrow donor 
            print('Recording data of %s'%(donorID), dfDonorID.shape)
            dfDonorID.to_hdf(os.path.join(extractPath, 'dfDonorID %s'%(donorID) + '.h5'), key=donorID, mode='a', complevel=4, complib='zlib')

    return donorIDs

def getDataframeByDonorID(extractPath, donorID):

    '''Get pandas.DataFrame by Donor ID

    Parameters:
        extractPath: str
            Path of directories where HCA matrix files were downloaded and extracted.
            See function getHCAdataByURL() for detail.

        donorID: str
            Donor ID.

    Returns:
        pandas.DataFrame
            Matrix corresponding to the Donor ID

    Usage:
        getDataframeByDonorID(extractPath, donorID)
    '''

    fileName = 'dfDonorID %s'%(donorID) + '.h5'
    fullFileName = os.path.join(extractPath, fileName)

    if os.path.isfile(fullFileName):
        print('Reading compressed DataFrame at:', fullFileName)

        dfDonorID = pd.read_hdf(fullFileName, key=donorID, mode='r')
    else:
        print('File %s not found'%(fullFileName))

        return

    return dfDonorID

def PrepareDataOnePatient_PREVIEW_DATASET(filename, patient, saveFolderName,  
                          useAllData=True, cellsLimitToUse=1000, randomlySample=True, randomSeed=0):

    '''Prepare data from Human Cell Atlas (HCA) preview dataset h5 data file. 
    The user can download the file ica_bone_marrow_h5.h5 from https://preview.data.humancellatlas.org/ 
    (Raw Counts Matrix - Bone Marrow) and place in folder data. 
    The file is ~485Mb and contains all 378000 cells from 8 bone marrow donors (BM1-BM8).
    Note: this data file is no longer available at HCA data server, however, some users may have a 
    copy of it and need to extract data from it.
        
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
