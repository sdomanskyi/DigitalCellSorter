# Digital Cell Sorter

[![DOI](https://badge.fury.io/gh/sdomanskyi%2FDigitalCellSorter.svg)](https://badge.fury.io/gh/sdomanskyi%2FDigitalCellSorter)
[![DOI](https://badge.fury.io/py/DigitalCellSorter.svg)](https://pypi.org/project/DigitalCellSorter)
[![DOI](https://readthedocs.org/projects/digital-cell-sorter/badge/?version=latest)](https://digital-cell-sorter.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3749878.svg)](https://doi.org/10.5281/zenodo.3749878) 

Digital Cell Sorter (DCS): a single cell RNA-seq analysis toolkit for clustering, cell type identification, and anomaly detection.

> **Note:** We are currently preparing a manuscript describing the toolkit located this repository.
> If you want to access the package detailed in our latest publication of Polled Digital Cell Sorter
> go to https://zenodo.org/record/2603265 and download the package (v1.1).


> **The latest publication describing the methodology of cell types identification:**
> [Polled Digital Cell Sorter (p-DCS): Automatic identification of hematological cell types from single cell RNA-sequencing clusters](
> https://doi.org/10.1186/s12859-019-2951-x 
> "Polled Digital Cell Sorter (p-DCS): Automatic identification of hematological cell types from single cell RNA-sequencing clusters")
> Sergii Domanskyi, Anthony Szedlak, Nathaniel T Hawkins, Jiayin Wang, Giovanni Paternostro & Carlo Piermarocchi, 
> *BMC Bioinformatics* volume 20, Article number: 369 (**2019**)


The documentation is available at https://digital-cell-sorter.readthedocs.io/.

- [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Loading the package](#loading-the-package)
  * [Gene Expression Data Format](#gene-expression-data-format)
  * [Other Data](#other-data)
- [Functionality](#functionality)
  * [Overall](#overall)
  * [Visualization](#visualization)
- [Demo](#demo)
  * [Usage](#usage)
    + [Main cell types](#main-cell-types)
    + [Cell sub-types](#cell-sub-types)
  * [Output](#output)

## Getting Started

These instructions will get you a copy of the project up and running on your machine for data analysis, development or testing purposes.

### Prerequisites

#### Environment setup
The software runs in Python >= 3.7

It is highly recommended to install Anaconda.
Installers are available at https://www.anaconda.com/distribution/
Whether you already had Anaconda installed or just installed it we recommend to
update all packages by running:

	conda update conda

With conda, you can create, export, list, remove, and update environments that 
have different versions of Python and/or packages installed in them. 
Switching or moving between environments is called activating the environment.

	conda create --name DCS
	conda activate DCS

Now, in your new environment, the packages can be installed or updated without affecting
your other environments. Note, environments use is not necessary, and the 
default ```(base)``` is used if you dont set up any other. For more information see
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

#### Installation of the DigitalCellSorter package

Install ```DigitalCellSorter``` with ```pip```. Most of the dependencies packages 
are automatically installed with installation of the  latest release 
of ```DigitalCellSorter```:

	pip install DigitalCellSorter

Alternatively, you can clone and install this module directly from GitHub using:

	pip install git+https://github.com/sdomanskyi/DigitalCellSorter

Similarly, one can create a local copy of this project for development purposes, and 
install the package from the cloned directory:

	git clone https://github.com/sdomanskyi/DigitalCellSorter
	python setup.py install

Our software uses packages ```numpy```, ```pandas```, ```matplotlib```, 
```scikit-learn```, ```scipy```, ```mygene```, ```fftw```, 
```fitsne```, ```adjustText``` and a few other standard Python packages. 
Some of the packages used in ```DigitalCellSorter``` are not installed by default, 
and should by installed by separately if using certain functionality with 
Digital Cell Sorter. For example, for network-based clustering
install packages ```pynndescent```, ```networkx```, ```python-louvain```. 
Other packages that have to be installed separately are ```fitsne```, ```umap```, 
```phate``` and ```orca```. The detailed instructions are below.

#### t-SNE
With datasets containing less than 2000 cells ```sklearn.manifold.TSNE``` is used.
For large datasets Fast Fourier Transform-accelerated Interpolation-based t-SNE (FIt-SNE)
implemented by **KlugerLab** is used (https://github.com/KlugerLab/FIt-SNE).
To use FIt-SNE the following need to be installed. First update ```cython``` by

	pip install --upgrade cython

Then install ```fftw``` from the ```conda-forge``` channel 
add ```conda-forge``` to your channels, and install ```fftw```:

	conda config --add channels conda-forge
	conda install fftw

The next installation step is platform specific. To install FI-tSNE for Linux:

	pip install fitsne

On macOS Mojave C++ compiler has to be specified explicitly:

	env CC=clang CXX=clang++ pip install fitsne

On Windows the FI-tSNE wrapper and executable are already 
included with ```DigitalCellSorter```. 

#### Other layouts

To use UMAP layout

	pip install umap-learn

To use PHATE 

	pip install phate

> Note, if neither ```fitsne```, ```umap``` nor ```phate``` are installed 
> ```DigitalCellSorter``` defaults to PCA two largest principal components for 
> visualization layout.

#### Interactive HTML figures
To use Sankey diagrams that are part of Digital Cell Sorter 
install ```plotly``` and ```orca```:

    conda install -c plotly plotly-orca
    conda install -c anaconda psutil

```orca	``` is necessary to convert Sankey diagrams to static images.
If for any reason ```orca``` is unavailable the Sankey diagrams will be saved as 
ineractive HTML figure, that can be opened in a browser (Chrome, Firefox etc.) and 
saved as static image. The visualization of ```DigitalCellSorter``` are implemented 
with ```matplotlib```, allowing all the figures to be saved in either raster or 
vactor format. Since ```plotly``` can convert simple ```matplotlib``` figures 
(scatter, line, bar plots, but not heatmaps, splines or other complex patch objects) to
ineractive HTML format ```DigitalCellSorter``` can attempt to save any of its figures
as HTML. This is particulatly useful with ```Projection``` plots, even though the color
bars are not rendered in HTML figures.

### Loading the package

In your script import the package:

	import DigitalCellSorter

Create an instance of class ```DigitalCellSorter```. Here, for simplicity, we use Default parameter values:

	DCS = DigitalCellSorter.DigitalCellSorter()

During the initialization a number of parameters can be specified. For detailed list see documentation.
Many of these parameters are transfered to DCS attributes thus can be modified after initialization using, e.g.:

	DCS.toggleMakeStackedBarplot = False



### Gene Expression Data Format

The input gene expression data is expected in one of the following formats:

1. Spreadsheet of comma-separated values ```csv``` containing condensed matrix in a form ```('cell', 'gene', 'expr')```. 
If there are batches in the data the matrix has to be of the form ```('batch', 'cell', 'gene', 'expr')```. Columns order can be arbitrary.

<details closed><summary>Examples:</summary><p>

| cell | gene | expr |
|------|------|------|
| C1   | G1   | 3    |
| C1   | G2   | 2    |
| C1   | G3   | 1    |
| C2   | G1   | 1    |
| C2   | G4   | 5    |
| ...  | ...  | ...  |

or:

| batch  | cell | gene | expr |
|--------|------|------|------|
| batch0 | C1   | G1   | 3    |
| batch0 | C1   | G2   | 2    |
| batch0 | C1   | G3   | 1    |
| batch1 | C2   | G1   | 1    |
| batch1 | C2   | G4   | 5    |
| ...    | ...  | ...  | ...  |

</p></details>


2. Spreadsheet of comma-separated values ```csv``` where rows are genes, columns are cells with gene expression counts.
If there are batches in the data the spreadsheet the first row should be ```'batch'``` and the second ```'cell'```.

<details closed><summary>Examples:</summary><p>

| cell  | C1     | C2     | C3     | C4     |
|-------|--------|--------|--------|--------|
| G1    |        | 3      | 1      | 7      |
| G2    | 2      | 2      |        | 2      |
| G3    | 3      | 1      |        | 5      |
| G4    | 10     |        | 5      | 4      |
| ...   | ...    | ...    | ...    | ...    |

or:

| batch | batch0 | batch0 | batch1 | batch1 |
|-------|--------|--------|--------|--------|
| cell  | C1     | C2     | C3     | C4     |
| G1    |        | 3      | 1      | 7      |
| G2    | 2      | 2      |        | 2      |
| G3    | 3      | 1      |        | 5      |
| G4    | 10     |        | 5      | 4      |
| ...   | ...    | ...    | ...    | ...    |

</p></details>

3. ```Pandas DataFrame``` where ```axis 0``` is genes and ```axis 1``` are cells.
If the are batched in the data then the index of ```axis 1``` should have two levels, e.g. ```('batch', 'cell')```, 
with the first level indicating patient, batch or expreriment where that cell was sequenced, and the
second level containing cell barcodes for identification.

<details closed><summary>Examples:</summary><p>

    df = pd.DataFrame(data=[[2,np.nan],[3,8],[3,5],[np.nan,1]], 
                      index=['G1','G2','G3','G4'], 
                      columns=pd.MultiIndex.from_arrays([['batch0','batch1'],['C1','C2']], names=['batch', 'cell']))    


</p></details>

4. ```Pandas Series ``` where index should have two levels, e.g. ```('cell', 'gene')```. If there are batched in the data
the first level should be indicating patient, batch or expreriment where that cell was sequenced, the second level cell barcodes for 
identification and the third level gene names.

<details closed><summary>Examples:</summary><p>

    se = pd.Series(data=[1,8,3,5,5], 
                   index=pd.MultiIndex.from_arrays([['batch0','batch0','batch1','batch1','batch1'],
                                                    ['C1','C1','C1','C2','C2'],
                                                    ['G1','G2','G3','G1','G4']], names=['batch', 'cell', 'gene']))


</p></details>

Any of the data types outlined above need to be prepared/validated with a function ```prepare()```. 
Let us demonstrate this on the input of type 1:

	df_expr = DCS.prepare('data/testData/dataFileCondensedWithBatches.tsv')

### Other Data

```markersDCS.xlsx```: An excel book with marker data. Rows are markers and columns are cell types. 
'1' means that the gene is a marker for that cell type, '-1' means that this gene is not expressed in this cell type, and '0' otherwise.
This gene marker file included in the package is used by Default. 
If you use your own file it has to be prepared in the same format (including the two-line header). Note that only the first worksheet will be read,
and its name can be arbitrary. The first column should contain gene names. The second row should contain cell types, and the first row how 
those cell types are grouped. If any of the cell types need to be skipped, have "NA" in the corresponding cell of the first row of that cell type.

<details closed><summary>Example:</summary><p>

|A       |B            |C             |D           |E          |F                |G                         |H                           |I                        |J                         |K                  |L               |M                 |...      |
|--------|-------------|--------------|------------|-----------|-----------------|--------------------------|----------------------------|-------------------------|--------------------------|-------------------|----------------|------------------|---------|
|        |B cells      |B cells       |B cells     |T cells    |T cells          |T cells                   |T cells                     |T cells                  |T cells                   |T cells            |NK cells        |NK cells          |...      |
|Marker  |B cells naive|B cells memory|Plasma cells|T cells CD8|T cells CD4 naive|T cells CD4 memory resting|T cells CD4 memory activated|T cells follicular helper|T cells regulatory (Tregs)|T cells gamma delta|NK cells resting|NK cells activated|...      |
|ABCB4   |1            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ABCB9   |0            |0             |1           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ACAP1   |0            |0             |0           |0          |1                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ACHE    |0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ACP5    |0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ADAM28  |1            |1             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ADAMDEC1|0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ADAMTS3 |0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ADRB2   |0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|AIF1    |0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|AIM2    |0            |1             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ALOX15  |0            |0             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ALOX5   |0            |1             |0           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|AMPD1   |0            |0             |1           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|ANGPT4  |0            |0             |1           |0          |0                |0                         |0                           |0                        |0                         |0                  |0               |0                 |...      |
|...     |...          |...           |...         |...        |...              |...                       |...                         |...                      |...                       |...                |...             |...               |...      |

</p></details>

```Human.MitoCarta2.0.csv```: An ```csv``` spreadsheet with human mitochondrial genes, created within work 
[MitoCarta2.0: an updated inventory of mammalian mitochondrial proteins](https://doi.org/10.1093/nar/gkv1003 "MitoCarta2.0")
Sarah E. Calvo, Karl R. Clauser, Vamsi K. Mootha, *Nucleic Acids Research*, Volume 44, Issue D1, 4 January 2016.


## Functionality

### Overall

The main class, DigitalCellSorter, includes tools for:

  1. **Pre-preprocessing**
  2. **Quality control**
  3. **Batch effects correction**
  4. **Cells anomaly score evaluation**
  4. **Dimensionality reduction**
  5. **Clustering**
  6. **Annotating cell types**
  7. **Vizualization**  
  8. **Post-processing**.


### Visualization

Function ```visualize()``` will produce most of the necessary files for post-analysis of the data. 

See examples of the visualization tools below.


<details closed><summary>The visualization tools include:</summary><p>
 
- ```makeMarkerExpressionPlot()```: a heatmap that shows all markers and their expression levels in the clusters, 
in addition this figure contains relative (%) and absolute (cell counts) cluster sizes

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_marker_expression.png?raw=true" width="1000"/>
</p>

- ```getIndividualGeneExpressionPlot()```:  2D layout colored by individual gene's expression

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/marker_subplots/BM1_CD19_(B4_CVID3_CD19).png?raw=true" width="400"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/marker_subplots/BM1_CD33_(SIGLEC-3_CD33_p67_SIGLEC3).png?raw=true" width="400"/>
</p>

- ```makeVotingResultsMatrixPlot()```: z-scores of the voting results for each input cell type and each cluster, 
in addition this figure contains relative (%) and absolute (cell counts) cluster sizes

<p align="middle">
 <img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_scores_matrix.png?raw=true" height="700"/>
</p>

- ```makeHistogramNullDistributionPlot()```: null distribution for each cluster and each cell type illustrating 
the "machinery" of the Digital Cell Sorter

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_null_distributions.png?raw=true" width="800"/>
</p>

- ```makeQualityControlHistogramPlot()```: Quality control histogram plots

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/QC_plots/BM1_number_of_genes_histogram.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/QC_plots/BM1_count_depth_histogram.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/QC_plots/BM1_fraction_of_mitochondrialGenes_histogram.png?raw=true" width="250"/>
</p>

- ```makeProjectionPlot()```: 2D layout colored by number of unique genes expressed, 
number of counts measured, and a faraction of mitochondrial genes..

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_number_of_genes.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_count_depth.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_fraction_of_mitochondrialGenes.png?raw=true" width="250"/>
</p>

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_is_quality_cell.png?raw=true" width="500"/>
</p>

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_clusters.png?raw=true" width="375"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_patients.png?raw=true" width="375"/>
</p>

Effect of batch correction demostrated on combining BM1, BM2, BM3 and processing the data jointly without (left) and with (right) batch correction option:

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/BM123_no_corr_clusters_by_patients.png?raw=true" width="375"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/BM123_with_corr_clusters_by_patients.png?raw=true" width="375"/>
</p>

- ```makeStackedBarplot()```: plot with fractions of various cell types

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_clusters_annotated.png?raw=true" width="500"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_subclustering_stacked_barplot_BM1.png?raw=true" height="500"/>
</p>


- ```makeSankeyDiagram()```: river plot to compare various results 

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/Sankey_example.png?raw=true" width="800"/>
</p>

- ```getAnomalyScoresPlot()```: plot with anomaly scores per cell

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score_All.png?raw=true" width="750"/>
</p>

Calculate and plot anomaly scores for an arbitrary cell type or cluster:

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score_B_cells.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score_T_cells.png?raw=true" width="250"/>
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score_cluster_7.0.0.png?raw=true" width="250"/>
</p>


- ```getIndividualGeneTtestPlot()```: Produce heatmap plot of t-test p-Values calculated gene-pair-wise
        from the annotated clusters

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_ttest_CD4_(CD4_CD4mut).png?raw=true" width="500"/>
</p>


- ```makePlotOfNewMarkers()```: genes significantly expressed in the annotated cell types

<p align="middle">
	<img src="https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_new_markers.png?raw=true" width="1000"/>
</p>

</p></details>


## Demo

### Usage

We have made an example execution file ```demo.py``` that shows how to use ```DigitalCellSorter```.

In the demo, folder ```data``` is intentionally left empty. 
The data file (cc95ff89-2e68-4a08-a234-480eca21ce79.homo_sapiens.mtx.zip) is about 2.4Gb in size and
will be downloaded with the ```demo.py``` script.

> Previously the HCA preview data was consolidated in file ```ica_bone_marrow_h5.h5``` and downloadable  
> from https://preview.data.humancellatlas.org/ (Raw Counts Matrix - Bone Marrow). 
> That file was ~485Mb and containing 378000 cells from 8 bone marrow donors (BM1-BM8). 

See details of the script ```demo.py``` at:

> [Example walkthrough of demo.py script](https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/ "Examples")


To execute the complete script ```demo.py``` run:

	python demo.py

*Note that the HCA BM1 data contains ~50000 sequenced cells, requiring more than 60Gb of RAM (we recommend to use High Performance Computers).
If you want to run our example on a regular PC or a laptop, you can use a randomly chosen number of cells:

    df_expr.sample(n=5000, axis=1)


### Output

All the output files are saved in ```output``` directory inside the directory where the ```demo.py``` script is. 
If you specify any other directory, the results will be generetaed in it.
If you do not provide any directory the results will appear in the root where the script was executed.
