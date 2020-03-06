
- 1.3.0
   * Modified core algorithm for cell type identification to account for markers that should not be expressed in a given cell type (negative markers)
   * Modified core algorithm to better normalize celltype/marker matrix
   * Modified core algorithm account for low quality votes and abundunt markers
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
