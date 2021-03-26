- 1.3.7.6
   * Minor bug fix of fail-safe fitsne package import

- 1.3.7.5
   * Upgraded function makeViolinPlot to enable combinations of columns

- 1.3.7.3-4
   * Added option of the fixed random seed (option *random_state*) for the whole prcessing pipeline

- 1.3.7.2
   * Added documentation of the function makeViolinPlot

- 1.3.7.1
   * Added makeViolinPlot visualization function enabling exploratory analysis of the numeric values distributions for comparison between categorical variables

   * Updated visualization functions returning

- 1.3.7
   * Added a function to import data from kallisto-bustools and cellranger
   * Updated documentation 

- 1.3.6
   * Added quick-demo materials
   
- 1.3.5
   * Miscellaneous code improvements and bug fixes

- 1.3.4.0-1.3.4.11
   * Integrated plotly offline figure saving (when orca is unavailable)
   * Added Quality Control pre-cut

- 1.3.2
   * Added Hopfield landscape visuzlization capability
   * Added network of underlying biological gene-gene interaction to the Hopfield annotaiton scheme

- 1.3.1
   * Minor API modifications

- 1.3.0
   * Modified pDCS algorithm for cell type identification to account for markers that should not be expressed in a given cell type (negative markers)
   * Modified pDCS celltype/marker matrix normalization
   * Modified pDCS algorithm account for low quality scores
   * Added Hopfield classifier for cell type annotation
   * Added ratio method for cell type annotation
   * Added options for consensus cell type annotation
   * Added cell markers pie summary function and plot
   * Added t-test for individual gene plot
   * Added several new user functions, for efficient and flexible extraction of cells, genes, clusters, etc.
   * Added anomaly score calculation and visualization
   * Refactored function for extraction of new markers based on cell type annotations to separate it from function process() of class DigitalCellSorter
   * Optimized implementation (for higher performance) of various function of this package
   * Detailed visualization functions API
   * Incorporated different clustering methods in addition to the widely-utilized hierarchical clustering 
   * Incorporated several types of high-dimensional data projection methods, such as efficient t-SNE, UMAP and simple PCA components.
   * Extended options for input data format
   * Included a set of functions to load data from Human Cell Atlas (HCA) and prepare it for processing

- 1.2.3
   * API updates, documentation updates

- 1.2.1
   * Minor updates, reshaped DigitalCellSorter into a stand-alone package

- 1.2.0
   * More features, better runtime efficiency

- 1.1
   * Updated method, signature matrices

- 1.0
   * First Release
