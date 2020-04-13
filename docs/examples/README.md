# Details of ```demo.py``` script

[![DOI](https://badge.fury.io/gh/sdomanskyi%2FDigitalCellSorter.svg)](https://badge.fury.io/gh/sdomanskyi%2FDigitalCellSorter)
[![DOI](https://badge.fury.io/py/DigitalCellSorter.svg)](https://pypi.org/project/DigitalCellSorter)
[![DOI](https://readthedocs.org/projects/digital-cell-sorter/badge/?version=latest)](https://digital-cell-sorter.readthedocs.io/en/latest/?badge=latest)

The ready-to-read documentation is available at https://digital-cell-sorter.readthedocs.io/.
The documentation of our software is built with [**Sphinx**](https://www.sphinx-doc.org/ "Sphinx") at 
[**ReadTheDocs.org**](https://readthedocs.org/).

[Data preparation](#usage)
[Main cell types](#main-cell-types)
[Cell sub-types](#cell-sub-types)
[Output](#output)

## Data preparation

In the demo, folder ```data``` is intentionally left empty. 
The data file (cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.mtx.zip) is about 2.4Gb in size and
will be downloaded with the ```demo.py``` script.

Note that the HCA BM1 data contains ~50000 sequenced cells, requiring more than 60Gb of RAM (we recommend to use High Performance Computers).
If you want to run our example on a regular PC or a laptop, you can use a randomly chosen small number (~5000) of cells.

In our example  we use the data of BM1 only, howerver all 8 bone marrow samples are downloaded.
Create a local variable that points to location where the ```demo.py``` script was executed:

    here = os.path.dirname(__file__)

The following URL points to data file at Data Portal of Human Cell Atlas:

    url = "https://data.humancellatlas.org/project-assets/project-matrices/cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.mtx.zip"

Create a variable to point to where the data will be downloaded and extracted:

    extractPath = os.path.join(here, 'data', os.path.splitext(os.path.basename(url))[0])

Import a module from DigitalCellSorter and call a function to download and unpack the data:

    import DigitalCellSorter.ReadPrepareDataHCA as prep

    prep.getHCAdataByURL(url, extractPath)
    
Call function recordFilesOfIndividualDonors to load the data from HCA Data Portal. 
Note that this data file and the individual files are large and will take up to 4Gb of disk space:

    id = prep.recordFilesOfIndividualDonors(extractPath, organName='bone marrow')[0]

Load gene expression data from h5 file:

    df_expr = prep.getDataframeByDonorID(extractPath, id)
    df_expr.columns.names = ['batch', 'cell']


## Main cell types

In these instructions we have already created an instance of ```DigitalCellSorter``` class (see section **Loading the package**) .
Let's modify some of the ```DCS``` attributes:

	DCS.dataName = 'BM1'
	DCS.saveDir = os.path.join(here, 'output', 'BM1', '')
    DCS.geneListFileName = 'CIBERSORT_LM22_7'

Now we are ready to load the data, ```prepare```(validate) it and ```process```.

Validate the expression data, so that it has correct form:

    DCS.prepare(df_expr)

Delete df_expr as now DCS contains the master copy of it:

    del df_expr

Process the expression data, i.e. quality control, dimensionality reduction, clustering:

	DCS.process(df_expr)

Load markers and annotate the processe data:

	DCS.annotate()

Then generate all the default plots by:

	DCS.visualize()

Make CD19 and CD33 gene expression plots:

    for name in DCS.getHugoName('CD19'):
        DCS.makeIndividualGeneExpressionPlot(name)
            
    for name in DCS.getHugoName('CD33'):
        DCS.makeIndividualGeneExpressionPlot(name)


## Cell sub-types

Further analysis can be done on cell types of interest, e.g. here 'T cell'.
Let's create a new instance of DigitalCellSorter to run "sub-analysis" with it:

    DCSsub = DigitalCellSorter.DigitalCellSorter(dataName='BM1', 
                                                nClusters=10, 
                                                doQualityControl=False,
                                                layout='PHATE',
                                                subclusteringName='T cell')

It is important to disable Quality control, because the low quality cells have already been identified and filtered with ```DCS```.
Also ```dataName``` parameter points to the location processed with ```DCS```. 
Next modify a few other attributes and process cell type 'T cell':

    DCSsub.saveDir = os.path.join(here, 'output', 'BM1', 'subclustering T cell', '')
    DCSsub.geneListFileName = os.path.join('here', 'docs', 'examples', 'CIBERSORT_T_SUB.xlsx')

Get index of T cells:

    indexOfTcells = DCS.getCells(celltype='T cell')

Get expression of these T cells using their index:

    df_expr = DCS.getExprOfCells(indexOfTcells)

Insert expression data into DCSsub:

    DCSsub.prepare(df_expr)

Process subtype 'T cell':

    DCSsub.process(dataIsNormalized=True)

Load marker genes and annotate cells:

    DCSsub.annotate()

Make plots of annotated data:

    DCSsub.visualize()


This way the 2D layout with annotated clusters (left) of T cell sub-types and the corresponding voting matrix (right) 
are generated by the function ```process()```:

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/subclustering T cell/BM1_clusters_by_clusters_annotated.png?raw=true" width="400"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/subclustering T cell/BM1_scores_matrix.png?raw=true" height="400"/>
</p>


## Output

All the output files are saved in ```output``` directory inside the directory where the ```demo.py``` script is. 
If you specify any other directory, the results will be generetaed in it.
If you do not provide any directory the results will appear in the root where the script was executed.



