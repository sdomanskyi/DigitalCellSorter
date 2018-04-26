# DigitalCellSorter
Cell type clustering and identification in heterogeneous single cell RNA-seq data

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The code runs in Python>=3.3 environment. It uses packages numpy, pandas, matplotlib, scikit-learn, scipy. One convenient way to install these packages is through pip. After pip is installed, you can install most packages by

```
pip3 install 'package name'
```

### Input Data Format

Input data are located in data/ folder. Each set of data should be an individual subfolder with a name you specify. In our example, we use data from data/aml035pre/ as input.


Your data should be in the same directory hierarchy in order for the code to work. 

### Other Data

-geneLists/cd_marker_handbook.xlsx: an excel with marker information. Rows are markers and columns are cell types. A '+' sign means that the gene is a marker for that cell type.

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

The Process function will produce the following outputs. Below are outputs produced using our sample data.
 
- dataName_clusters.png: shows the clustering of cells and identified cell type of each cluster. 

 <img src="https://github.com/wangjiayin1010/DigitalCellSorter/blob/master/demo_output/aml035pre_clusters.png" align="left" height="600" width="600" >
 

- dataName_voting.png: a heatmap that shows all markers and their expression levels in the clusters

<img src="https://github.com/wangjiayin1010/DigitalCellSorter/blob/master/demo_output/aml035pre_voting.png" align="left" height="600" width="1200" >

- dataName_voting.xlsx: an excel that shows the voting results

- dataName_expression_labeled.tar.gz: a zip file that contains processed expression data

- marker_subplots: a directory that contains subplots of each marker and their expression levels in the clusters


## Demo

### Usage

We've made an example execution file demo.py that shows how to use the DigitalCellSorter. The data set we are analyzing, the healthy and AML data from 10X, are in MatrixMarket IJV format (genes.tsv and matrix.mtx), so it first converts them to a pandas dataframe as input for Process function. Then it creates a DigitalCellSorter class and calls its Process function, which does all the projection, clustering, identification and image producing.

You can run it using 

```
python demo.py
```
Note that you would also need to convert your data to a pandas dataframe. You can also customize the parameters in Process function. For example, you can change it to

```
dcs.Process(n_clusters=10, n_components_pca=50, saveDir='demo_output/', marker_subplot = False)
```

To see our example work, you just need to download everything, go to the directory then run 

```
python demo.py
```

For other purposes, you can add your own data data/ folder, change the function call in run.py then run the same command.

### Output

All the outputs are saved in demo/output directory. 

