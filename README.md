# DigitalCellSorter
Identification of hematological cell types from heterogeneous single cell RNA-seq data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The code runs in Python>=3.3 environment. It uses packages numpy, pandas, matplotlib, scikit-learn, scipy, mygene, biopython, progressbar3, pickle, xlrd, openpyxl, h5py, asdf, shutil, gzip, multiprocessing. One convenient way to install these packages is through pip. After pip is installed, you can install most packages by

```
pip3 install 'package name'
```

### Input Data Format

The input data is expected in MatrixMarket IJV format, e.g. there should be three files ```matrix.mtx```, ```genes.tsv``` and ```barcodes.tsv```.

In the demo folder "data" is intentionally left empty. The reader can download the file ```ica_bone_marrow_h5.h5``` from https://preview.data.humancellatlas.org/ (Raw Counts Matrix - Bone Marrow) and place in folder ```data```. The file is ~485Mb and contains all 378000 cells from 8 bone marrow donors. In our example, the MatrixMarket files are prepared automatically on the fly for each donor separately using ```ReadPrepareDataHCApreviewDataset.PrepareData()```. For detail see ```demo.py```.

The data is converted to a pandas dataframe as input to the main function. The dataframe has genes as rows and cells as columns, and each grid is the expression of that gene in the cell. 

### Other Data

geneLists/G_2_Human_cell_markers_BM.xlsx: An excel sheet with marker data. Rows are markers and columns are cell types. A '1' digit means that the gene is a marker for that cell type, and '0' otherwise.

## Method

The main class for cell sorting functions and producing output images is DigitalCellSorter located in scripts/DigitalCellSorter.py. It takes the following parameters that you can customize.

- df_expr: expression data in panda dataframe format

- dataName: name used in output files

- sigma_over_mean_sigma: threshold when keeping only genes with large enough standard deviation, default is 0.3.

- n_clusters: number of clusters, default is 11

- n_components_pca: number of principal components to keep for clustering, default is 100
 
- zscore_cutoff: zscore cutoff when calculating Zmc, default is 0.3

- saveDir: directory to save all outputs, default is None

- marker_expression_plot: whether to produce marker expression plot, default is True

- tSNE_cluster_plot: whether to produce cluster plot, default is True

- save_processed_data: whether to save processed data, default is True

- marker_subplot: whether to make subplots on markers, default is True

- votingScheme: voting function to predict celltypes for each cluster, default is None. If it is None, then our cell sorting algorithm from the paper would be used. If you want to your own voting function, make sure it returns a dictionary with cluster labels as keys and predicted celltypes as values. You can modify code to change inputs to this function if needed. See the code and comments for more clarifications.

The Process function will produce the following outputs. Below explains what they are and shows some sample outputs produced using our sample data.
 
- dataName_clusters.png: an image that shows the clustering of cells and identified cell type of each cluster. 

 <img src="https://github.com/wangjiayin1010/DigitalCellSorter/blob/master/demo_output/HCA_BM1_data/HCA_BM1_data_clusters2D.png" align="center" height="500" width="500" >
 
- dataName_voting.png: a heatmap that shows all markers and their expression levels in the clusters

<img src="https://github.com/wangjiayin1010/DigitalCellSorter/blob/master/demo_output/aml035pre_voting.png" align="center">

- dataName_voting.xlsx: an excel sheet that shows the voting results

- dataName_expression_labeled.tar.gz: a zip file that contains processed expression data

- marker_subplots: a directory that contains subplots of each marker and their expression levels in the clusters. For example below is the subplot of CD33, a myeloid marker.

<img src="https://github.com/wangjiayin1010/DigitalCellSorter/blob/master/demo_output/HCA_BM1_data/marker_subplots/HCA_BM1_data_CD33_CD33.png" align="center" height="500" width="500" >

## Demo

### Usage

We've made an example execution file demo.py that shows how to use DigitalCellSorter. The data set we are analyzing, the healthy and AML data from 10X, are in MatrixMarket IJV format (genes.tsv and matrix.mtx under data/aml350pre/), so it first converts them to a pandas dataframe as input for Process function. Then it creates a DigitalCellSorter class and calls its Process function, which does all the projection, clustering, identification and image producing.

You can run it using 

```
python DCS_demo_on_HCA_BM.py
```
Note that to use your own data, you would also need to convert them to a pandas dataframe. You can also customize the parameters in Process function as listed above. For example, you can change it to

```
dcs.Process(n_clusters=10, n_components_pca=50, saveDir='demo_output/', marker_subplot = False)
```

To see our example work, you just need to download everything, go to the directory then run 

```
python DCS_demo_on_HCA_BM.py
```

### Output

All the outputs are saved in ```demo_output/``` directory. 

