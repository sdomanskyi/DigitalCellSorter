import sys
sys.path.append("..")

import os
import urllib.request
import DigitalCellSorter
import DigitalCellSorter.ReadPrepareDataHCA as prep

if __name__ == '__main__':

    print('This is a large dataset demo.\nFor the "5k PBMC demo" run "python -m DigitalCellSorter"\n')

    here = os.path.dirname(__file__)

    url = "https://data.humancellatlas.org/project-assets/project-matrices/cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.mtx.zip"
    extractPath = os.path.join(here, 'data', os.path.splitext(os.path.basename(url))[0])

    ## Download and unpack data
    #prep.getHCAdataByURL(url, extractPath)

    # Call function recordFilesOfIndividualDonors to load the data from HCA Data Portal
    #id = prep.recordFilesOfIndividualDonors(extractPath, organName='bone marrow')[0]
    id = '085e737d-adb5-4597-bd54-5ebeda170038'

    # Get the data. The file will be downloaded from github if not found locally
    try:
        if not os.path.exists(extractPath):
            os.makedirs(extractPath)

        if not os.path.isfile(os.path.join(extractPath, 'dfDonorID %s.h5' % id)):
            print('Downloading 110 Mb data file (50000 cells)')
            urllib.request.urlretrieve('https://github.com/sdomanskyi/DigitalCellSorter/raw/master/data/dfDonorID %s.h5' % id, extractPath)
    except Exception as exception:
        print('Could not download the file\n', exception)
        exit()

    # Load gene expression data from h5 file
    df_expr = prep.getDataframeByDonorID(extractPath, id)
    df_expr.columns.names = ['batch', 'cell']

    # Create an instance of class DigitalCellSorter. 
    # Here we use Default parameter values for most of the parameters
    DCS = DigitalCellSorter.DigitalCellSorter(dataName='BM1', geneNamesType = 'ensembl',
                                              saveDir=os.path.join(here, 'output', 'BM1', ''),
                                              geneListFileName='CIBERSORT_LM22_7')

    # Validate the expression data, so that it has correct form
    DCS.prepare(df_expr)

    # Delete df_expr as now DCS contains the master copy of it
    del df_expr

    # Process the expression data, i.e. quality control, dimensionality reduction, clustering
    DCS.process()

    # Load marker genes and annotate cells
    DCS.annotate()

    # Make plots of annotated data
    DCS.visualize()
    
    # Make CD19 gene expression plot
    for name in DCS.getHugoName('CD19'):
        DCS.makeIndividualGeneExpressionPlot(name)
            
    # Make CD33 gene expression plot
    for name in DCS.getHugoName('CD33'):
        DCS.makeIndividualGeneExpressionPlot(name)
    
    # Further analysis can be done on cell types of interest, e.g. here 'T cell' and 'B cell'.
    # Let's create a new instance of DigitalCellSorter to run "sub-analysis" with it.
    # It is important to disable Quality control, because the low quality cells have 
    # already been identified and filtered with DCS.
    # Parameter dataName points to the location processed with DCS. 
    DCSsub = DigitalCellSorter.DigitalCellSorter(dataName='BM1', 
                                                 nClusters=10, 
                                                 doQualityControl=False,
                                                 layout='PHATE',
                                                 subclusteringName='T cell')

    # Modify a few other attributes
    DCSsub.saveDir = os.path.join(here, 'output', 'BM1', 'subclustering T cell', '')
    DCSsub.geneListFileName = os.path.join(here, 'docs', 'examples', 'CIBERSORT_T_SUB.xlsx')

    # Get index of T cells
    indexOfTcells = DCS.getCells(celltype='T cell')

    # Get expression of these T cells using their index
    df_expr = DCS.getExprOfCells(indexOfTcells)

    # Insert expression data into DCSsub
    DCSsub.prepare(df_expr)

    # Process subtype 'T cell'
    DCSsub.process(dataIsNormalized=True)

    # Load marker genes and annotate cells
    DCSsub.annotate()

    # Make plots of annotated data
    DCSsub.visualize()
