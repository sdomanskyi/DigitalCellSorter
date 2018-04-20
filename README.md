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

Input data are located in data/ folder. Each set of data should be an individual subfolder with a name you specify. In our example, we use data from data/aml035pre/ as input. This subfolder contains two required data files:

- genes.tsv: a tsv file where the second column of each row is the gene name. 

- matrix.mtx: #TODO

Your data should be in the same directory hierarchy in order for the code to work. 

### Other Data

-geneLists/cd_marker_handbook.xlsx: an excel with marker information. Rows are markers and columns are cell types. A '+' sign means that the gene is a marker for that cell type.

## Running the tests

The main class for cell sorting functions and producing output images is DigitalCellSorter located in scripts/DigitalCellSorter.py. The __init__ method takes in dataName, which is the subfolder name in data/ containing your input data. It processes the raw data. The Process function does all the clustering, identification and image producing.

The example execution file is run.py. It creates a DigitalCellSorter class and calls its Process function. You can run it using 

```
python run.py
```
Note that you can customize the parameters in Process function. 

- sigma_over_mean_sigma: threshold when keeping only genes with large enough standard deviation, default is 0.3.

- n_clusters: number of clusters, default is 11

- n_components_pca: number of principal components to keep for clustering, default is 100
 
- zscore_cutoff: zscore cutoff when calculating Zmc, default is 0.3

- saveDir: directory to save all outputs, default is None

You can change them by changing the parameters in the call to Process in run.py. For example, you can change it to

```
dcs.Process(n_clusters=10, n_components_pca=50, saveDir='dcs_output/')
```

To see our example work, you just need to download everything, go to the directory then run 

```
python run.py
```

For other purposes, you can add your own data data/ folder, change the function call in run.py then run the same command.
