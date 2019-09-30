# Digital Cell Sorter
Identification of hematological cell types from heterogeneous single cell RNA-seq data.

[Polled Digital Cell Sorter (p-DCS): Automatic identification of hematological cell types from single cell RNA-sequencing clusters](
https://doi.org/10.1186/s12859-019-2951-x 
"Polled Digital Cell Sorter (p-DCS): Automatic identification of hematological cell types from single cell RNA-sequencing clusters")
Sergii Domanskyi, Anthony Szedlak, Nathaniel T Hawkins, Jiayin Wang, Giovanni Paternostro & Carlo Piermarocchi,
 *BMC Bioinformatics* volume 20, Article number: 369 (**2019**)


## Getting Started

These instructions will get you a copy of the project up and running on your machine for data analysis, development or testing purposes.

### Prerequisites

The code runs in Python >= 3.7 environment. 

It is highly recommended to install Anaconda.
Installers are available at https://www.anaconda.com/distribution/

It uses packages ```numpy```, ```pandas```, ```matplotlib```, ```scikit-learn```, ```scipy```, 
```mygene```, ```fftw```, ```pynndescent```, ```networkx```, ```python-louvain```, ```fitsne```
and a few other standard Python packages. Most of these packages are installed with installation of the 
latest release of ```DigitalCellSorter```:

	pip install DigitalCellSorter

Alternatively, you can install this module directly from GitHub using:

	pip install git+https://github.com/sdomanskyi/DigitalCellSorter

Installing ```fftw``` from the ```conda-forge``` channel can be achieved by adding ```conda-forge``` to your channels with:

	conda config --add channels conda-forge

Once the conda-forge channel has been enabled, ```fftw``` can be installed with:

	conda install fftw

### Loading the package

In your script import the package:

	import DigitalCellSorter.DigitalCellSorter as DigitalCellSorter

Create and instance of a ```DigitalCellSorter```, here, for simplicity, we use Default parameter values:

	DCS = DigitalCellSorter.DigitalCellSorter()


<details><summary>During the initialization the following parameters can be specified (click me)</summary><p>

```dataName```: name used in output files, Default ''

```geneListFileName```: marker cell type list name, Default None

```mitochondrialGenes```: list of mitochondrial genes for quality conrol routine, Default None

```sigmaOverMeanSigma```: threshold to consider a gene constant, Default 0.3

```nClusters```: number of clusters, Default 5

```nComponentsPCA```: number of pca components, Default 100

```nSamplesDistribution```: number of random samples to generate, Default 10000

```saveDir```: directory for output files, Default is current directory

```makeMarkerSubplots```:  whether to make subplots on markers, Default True

```makePlots```: whether to make all major plots, Default True

```votingScheme```: voting shceme to use instead of the built-in, Default None

```availableCPUsCount```: number of CPUs available, Default os.cpu_count()

```zScoreCutoff```: zscore cutoff when calculating Z_mc, Default 0.3

```clusterName```: parameter used in subclustering, Default None

```doQualityControl```: whether to remove low quality cells, Default True

```doBatchCorrection```: whether to correct data for batches, Default False

</p></details>

These and other parameters can be modified after initialization using, e.g.:

	DCS.toggleMakeStackedBarplot = False



### Gene Expression Data Format

<details open><summary>The input gene expression data is expected in one of the following formats:</summary><p>

1. Condensed matrix in a form ```('cell', 'gene', 'expression count')```, while cell and gene names may be supplied in the 
additional lists, in ```tsv```, ```csv```, ```xlsx```, etc. format. 
If there are batches in the data the matrix has to be of the form 
```('batch', 'cell', 'gene', 'expression count')``` or the batches information can be supplied as a separate list.

2. Spreadsheet in ```tsv```, ```csv```, ```xlsx```, etc. format where rows are genes, columns are cells with gene expression counts.
Cell and gene names may be supplied in the additional lists, in ```tsv```, ```csv```, ```xlsx```, etc. format.
If there are batches in the data the spreadsheet the header should have two rows: ```'batch'``` and ```'cell'``` or the batches information can be supplied as a separate matrix 
in a form ```('batch', 'cell')```.

3. ```Pandas DataFrame``` where ```axis 0``` is genes and ```axis 1``` are cells.
If the are batched in the data then the index of ```axis 1``` should have two levels, e.g. ```('batch', 'cell')```, 
with the first level indicating patient, batch or expreriment where that cell was sequenced, and the
second level containing cell barcodes for identification.

4. ```Pandas Series ``` where index should have two levels, e.g. ```('cell', 'gene')```. If there are batched in the data
the first level should be indicating patient, batch or expreriment where that cell was sequenced, the second level cell barcodes for 
identification and the third level gene names.

</p></details>

Any of the data types outlined above need to be prepared/validated with a function ```prepare()```. 
Let us demonstrate this on the input of type 1:

	df_expr = DCS.prepare(data='data/testData/data.tsv', 
				genes='data/testData/genes.tsv', 
				cells='data/testData/barcodes.tsv',
				batches=None)

### Other Data

```geneLists/markersDCS.xlsx```: An excel book with marker data. Rows are markers and columns are cell types. 
'1' means that the gene is a marker for that cell type, and '0' otherwise.
This gene marker file included in the package is used by Default. 
If you use your own file it has to be prepared in the same format.

```geneLists/Human.MitoCarta2.0.xls```: An excel book with human mitochondrial genes, created within work 
[MitoCarta2.0: an updated inventory of mammalian mitochondrial proteins](https://doi.org/10.1093/nar/gkv1003 "MitoCarta2.0")
Sarah E. Calvo, Karl R. Clauser, Vamsi K. Mootha, *Nucleic Acids Research*, Volume 44, Issue D1, 4 January 2016, Pages D1251–D1257.


## Functionality

The main class for cell sorting functions and producing output images is DigitalCellSorter

<details open><summary>The class includes tools for:</summary><p>

  1. Pre-preprocessing single cell mRNA sequencing data (gene expression data)
     1. Cleaning: filling in missing values, zemoving all-zero genes and cells, converting gene index to a desired convention, etc.
     2. Normalizing: rescaling all cells expression, log-transforming, etc.

  2. Quality control
  3. Batch effects correction
  4. Dimensionality reduction
  5. Clustering (Hierarchical, K-Means, knn-graph-based, etc.)
  6. **Annotating cell types**
  7. **Vizualization**
       1. t-SNE layout plot
       2. Quality Control histogram plot
       3. Marker expression t-SNE subplot
       4. Marker-centroids expression plot
       5. Voting results matrix plot
       6. Cell types stacked barplot
       7. Histogram null distribution plot
       8. New markers plot
       9. Sankey diagram (a.k.a. river plot)
  
  8. Post-processing functions, e.g. extract cells of interest, find significantly expressed genes, 
plot marker expression of the cells of interest, etc.

</p></details>


The ```process()``` function will produce all necessary files for post-analysis of the data. 

<details open><summary>The visualization tools include:</summary><p>
 
- ```makeMarkerExpressionPlot()```: a heatmap that shows all markers and their expression levels in the clusters, 
in addition this figure contains relative (%) and absolute (cell counts) cluster sizes

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_voting.png?raw=true" width="1000"/>
</p>

- ```getIndividualGeneExpressionPlot()```:  t-SNE layout colored by individual gene's expression

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/marker_subplots/BM1_CD19_(B4_CVID3_CD19).png?raw=true" width="400"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/marker_subplots/BM1_CD4_(CD4_CD4mut).png?raw=true" width="400"/>
</p>

- ```makeVotingResultsMatrixPlot()```: z-scores of the voting results for each input cell type and each cluster, 
in addition this figure contains relative (%) and absolute (cell counts) cluster sizes

<p align="middle">
 <img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_matrix_voting.png?raw=true" height="700"/>
</p>

- ```makeHistogramNullDistributionPlot()```: null distribution for each cluster and each cell type illustrating 
the "machinery" of the Digital Cell Sorter

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_null_distributions.png?raw=true" width="800"/>
</p>

- ```makeQualityControlHistogramPlot()```: Quality control histogram plots

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/QC_plots/BM1_number_of_genes_histogram.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/QC_plots/BM1_count_depth_histogram.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/QC_plots/BM1_fraction_of_mitochondrialGenes_histogram.png?raw=true" width="250"/>
</p>

- ```makeTSNEplot()```: t-SNE layouts colored by number of unique genes expressed, 
number of counts measured, and a faraction of mitochondrial genes..

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_number_of_genes.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_count_depth.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_fraction_of_mitochondrialGenes.png?raw=true" width="250"/>
</p>

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_is_quality_cell.png?raw=true" width="500"/>
</p>

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_clusters.png?raw=true" width="375"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_patients.png?raw=true" height="375"/>
</p>


- ```makeStackedBarplot()```: plot with fractions of various cell types

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_clusters_annotated.png?raw=true" width="500"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_subclustering_stacked_barplot_.png?raw=true" height="500"/>
</p>


- ```makeSankeyDiagram()```: river plot to compare various results

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/Sankey_example.png?raw=true" width="800"/>
</p>

- ```getAnomalyScoresPlot()```: plot with anomaly scores per cell

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_anomaly_score All.png?raw=true" width="750"/>
</p>

Calculate and plot anomaly scores for an arbitrary cell type or cluster:

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_anomaly_score B cell.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_anomaly_score T cell.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_clusters_by_anomaly_score Cluster2.png?raw=true" width="250"/>
</p>


- ```makePlotOfNewMarkers()```: genes significantly expressed in the annotated cell types

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/output/BM1/BM1_new_markers.png?raw=true" width="1000"/>
</p>

</p></details>


## Demo

### Usage


The function ```process()``` takes takes as an input parameter a pandas DataFrame validated by function ```process()```:

	DCS.process(df_expr) 

We have made an example execution file ```demo.py``` that shows how to use ```DigitalCellSorter```.

In the demo, folder ```data``` is intentionally left empty. The reader can download the file ```ica_bone_marrow_h5.h5``` 
from https://preview.data.humancellatlas.org/ (Raw Counts Matrix - Bone Marrow) and place in folder ```data```. 
The file is ~485Mb and contains all 378000 cells from 8 bone marrow donors (BM1-BM8). 
In our example, the data of BM1 is prepared by 
function ```PrepareDataOnePatient()``` in module ```ReadPrepareDataHCApreviewDataset```.
To load this function one can do:

	from DigitalCellSorter.ReadPrepareDataHCApreviewDataset import PrepareDataOnePatient as HCAtool

Now call this function to create a ```BM1.h5``` file (HDF file of input type 3) in the ```data``` folder:

	HCAtool(os.path.join('data', 'ica_bone_marrow_h5.h5'), 'BM1', os.path.join('data', '') 

In these instructions we have already created an instance of ```DigitalCellSorter```, 
let's modify some of the ```DCS``` attributes:

	DCS.dataName = 'BM1'
	DCS.saveDir = os.path.join('output', dataName, '')
	DCS.geneListFileName = os.path.join('geneLists', 'CIBERSORT.xlsx')
	DCS.nClusters = 20

Now we are ready to ```load``` the data, ```validate``` it and ```process```:

	df_expr = pd.read_hdf(os.path.join('data', 'BM1.h5'), key='BM1', mode='r')

	df_expr = DCS.prepare(df_expr)
	
	DCS.process(df_expr)

Further analysis can be done on cell types of interest, e.g. here 'T cell' and 'B cell'.
Let's create a new instance of DigitalCellSorter to run "sub-analysis" with it:

    DCSsub = DigitalCellSorter.DigitalCellSorter(dataName=DCS.dataName, 
                                                nClusters=10, 
                                                doQualityControl=False)

Here it was important to disable Quality control, because the low quality cells have already been identified with ```DCS```.
Also ```dataName``` parameter points to the location processed with ```DCS```. 
Now modify a few other attributes and process cell type 'T cell':

    DCSsub.subclusteringName = 'T cell'
    DCSsub.saveDir = os.path.join('output', dataName, 'subclustering T cell', '')
    DCSsub.geneListFileName = os.path.join('geneLists', 'CIBERSORT_T_SUB.xlsx')
	
    DCSsub.process(df_expr[DCS.getCells(celltype='T cell')])

We can reuse the ```DCSsub``` to analyze cell type 'B cell':

    DCSsub.subclusteringName = 'B cell'
    DCSsub.saveDir = os.path.join('output', dataName, 'subclustering B cell', '')
    DCSsub.geneListFileName = os.path.join('geneLists', 'CIBERSORT_B_SUB.xlsx')
	
    DCSsub.process(df_expr[DCS.getCells(celltype='B cell')])

For a complete script see:

	python demo.py

### Output

All the output files are saved in ```output``` directory. If you specify any other directory, the results will be generetaed in it.
If you do not provide any directory the results will appear in the root where the script was executed.
