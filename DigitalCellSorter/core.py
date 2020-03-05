'''The main class for cell sorting functions and producing output images is DigitalCellSorter. The class includes tools for:

  1. **Pre-preprocessing** of single cell mRNA sequencing data (gene expression data)
       1. Cleaning: filling in missing values, zemoving all-zero genes and cells, converting gene index to a desired convention, etc.
       2. Normalizing: rescaling all cells expression, log-transforming, etc.

  2. **Quality control**

  3. **Batch effects correction**

  4. **Cells anomaly score evaluation**

  5. **Dimensionality reduction**

  6. **Clustering** (Hierarchical, K-Means, knn-graph-based, etc.)

  7. **Annotating cell types**

  8. **Vizualization**
       1. 2D layout (projection) plot
       2. Quality Control histogram plot
       3. Marker expression projection subplot
       4. Marker-centroids expression plot
       5. Voting results matrix plot
       6. Cell types stacked barplot
       7. Anomaly scores plot
       8. Histogram null distribution plot
       9. New markers plot
       10. Sankey diagram (a.k.a. river plot)
  
  9. **Post-processing** functions, e.g. extract cells of interest, find significantly expressed genes, 
  plot marker expression of the cells of interest, etc.
'''

import os
import platform
import copy
import multiprocessing

import numpy as np
import pandas as pd

import scipy.stats
import scipy.signal
import scipy.cluster.hierarchy
import scipy.spatial.distance
from scipy.interpolate import UnivariateSpline

import sklearn.metrics
import sklearn.cluster
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralCoclustering
from sklearn.ensemble import IsolationForest

from . import GeneNameConverter, geneLists
from .Combat import combat
from .VisualizationFunctions import VisualizationFunctions
from .GenericFunctions import read, write
from .VisualizationFunctions import cm

class DigitalCellSorter(VisualizationFunctions):

    '''Class of Digital Cell Sorter with methods for processing single cell mRNA-seq data.
    Includes analyses and visualization tools.

    Parameters:
        dataName: str, Default 'dataName'
            Name used in output files

        geneNamesType: str, Default 'alias'
            Input gene name convention

        geneListFileName: str, Default None
            Name of the marker genes file

        mitochondrialGenes: list, Default None
            List of mitochondrial genes to use in quality control

        sigmaOverMeanSigma: float, Default 0.3
            Threshold to consider a gene constant

        nClusters: int, Default 5
            Number of clusters

        nFineClusters: int, Default 3
            Number of fine clusters to determine with Spectral Co-clustering routine.
            This option is ignored is doFineClustering is False.

        doFineClustering: boolean, Default True
            Whether to do fine clustering or not

        minSizeForFineClustering: int, Default 50
            Minimum number of cells required to do fine clustering of a cluster.
            This option is ignored is doFineClustering is False.

        clusteringFunction: function, Default AgglomerativeClustering
            Clustering function to use.
            Other options: KMeans, {k_neighbors:40}, etc.

        nComponentsPCA: int, Default 100
            Number of pca components

        nSamplesDistribution: int, Default 10000
            Number of random samples in distribution

        saveDir: str, Default os.path.join('')
            Directory for output files

        makeMarkerSubplots: boolean, Default True
            Whether to make subplots on markers

        makePlots: boolean, Default True
            Whether to make all major plots

        votingScheme: function, Default None
            Voting shceme to use instead of the built-in

        availableCPUsCount: int, Default os.cpu_count()
            Number of CPUs available

        zScoreCutoff: float, Default 0.3
            Z-Score cutoff when calculating Z_mc

        minimumScoreForUnknown: float, Default 0.3
            Z-Score cutoff when assigning label "Unknown"

        clusterName: str, Default None
            Parameter used in subclustering

        doQualityControl: boolean, Default True
            Whether to remove low quality cells

        doBatchCorrection: boolean, Default False
            Whether to correct data for batches
    
    Usage:
        DCS = DigitalCellSorter.DigitalCellSorter()

        df_data = DCS.Clean(df_data)
    '''

    gnc = GeneNameConverter.GeneNameConverter(dictDir=os.path.join('DigitalCellSorter', 'pickledGeneConverterDict', 'ensembl_hugo_entrez_alias_dict.pythdat'))

    def __init__(self, df_expr=None, dataName='dataName', geneNamesType='alias', geneListFileName=None, dataIsNormalized=False, mitochondrialGenes=None,
                sigmaOverMeanSigma=0.3, nClusters=5, nFineClusters=3, doFineClustering=True, minSizeForFineClustering=50, 
                clusteringFunction=AgglomerativeClustering, nComponentsPCA=200, nSamplesDistribution=10000, 
                saveDir=os.path.join(''), makeMarkerSubplots=False, availableCPUsCount=min(12, os.cpu_count()), zScoreCutoff=0.3,
                subclusteringName=None, doQualityControl=True, doBatchCorrection=True, makePlots=True,
                minimumNumberOfMarkersPerCelltype=10, minimumScoreForUnknown=0.3, layout='TSNE'):

        '''Initialization function. Automatically called when an instance on Digital Cell Sorter is created'''

        defaultGeneList = 'CIBERSORT' # 'CIBERSORT' 'markersDCS'

        self.saveDir = saveDir

        if self.saveDir!=os.path.join('') and not os.path.exists(self.saveDir):
            os.makedirs(self.saveDir)

        self.dataName = dataName

        self.df_expr = df_expr

        self.geneNamesType = geneNamesType
        
        self.defualtGeneListsDir = os.path.join(os.path.dirname(__file__), 'geneLists')

        self.geneListFileName = os.path.join(self.defualtGeneListsDir, defaultGeneList + '.xlsx') if geneListFileName is None else geneListFileName

        self.dataIsNormalized = dataIsNormalized

        self.mitochondrialGenes = mitochondrialGenes

        self.sigmaOverMeanSigma = sigmaOverMeanSigma
        self.nClusters = nClusters
        self.doFineClustering = doFineClustering
        self.nFineClusters = nFineClusters
        self.minSizeForFineClustering = minSizeForFineClustering
        self.nComponentsPCA = nComponentsPCA
        self.zScoreCutoff = zScoreCutoff
        self.minimumNumberOfMarkersPerCelltype = minimumNumberOfMarkersPerCelltype

        self.layout = layout

        self.minimumScoreForUnknown = minimumScoreForUnknown

        self.nSamplesDistribution = nSamplesDistribution
        self.availableCPUsCount = availableCPUsCount

        self.clusteringFunction = clusteringFunction

        self.subclusteringName = subclusteringName

        self.toggleRecordAllExpression = True

        self.toggleRemoveLowQualityCells  =  doQualityControl
        self.toggleDoBatchCorrection = doBatchCorrection

        self.toggleMakeHistogramNullDistributionPlot = makePlots
        self.toggleMakeVotingResultsMatrixPlot = makePlots
        self.toggleMakeMarkerExpressionPlot = makePlots
        self.toggleMakeProjectionPlotsQC = makePlots
        self.toggleMakeProjectionPlotClusters = makePlots
        self.toggleMakeProjectionPlotAnnotatedClusters = makePlots
        self.toggleMakeProjectionPlotBatches = makePlots
        self.toggleMakeStackedBarplot = makePlots
        self.toggleAnomalyScoresProjectionPlot = False

        self.toggleGetMarkerSubplots = makeMarkerSubplots

        super().__init__(saveDir=saveDir, dataName=dataName)

        return None

    # Main functions of class #################################################################################################################################
    def prepare(self, obj):

        '''Prepare pandas.DataFrame for input to function process()
        If input is pd.DataFrame validate the input whether it has correct structure.

        Parameters:
            obj: str, pandas.DataFrame, pandas.Series
                Expression data in a form of pandas.DataFrame, pandas.Series, or name and path to a csv file with data

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            dDCS.preapre('data.csv')
        '''

        if type(obj) is pd.Series:
            print('Received data in a form of pandas.Series', flush=True)
            print('Validating it and converting it to pandas.DataFrame', flush=True)

            if not 'cell' in obj.index.names:
                print('Column "cell" not found. Returning None', flush=True)
                return None

            if not 'gene' in obj.index.names:
                print('Column "gene" not found. Returning None', flush=True)
                return None

            if not 'batch' in obj.index.names:
                print('Column "batch" not found. Assuming one batch in the data.', flush=True)
                batch = np.array(['batch0']*len(obj.index.get_level_values('cell')))
            else:
                batch = obj.index.get_level_values('batch')

            obj.index = pd.MultiIndex.from_arrays([batch, obj.index.get_level_values('cell'), obj.index.get_level_values('gene')], names=['batch', 'cell', 'gene'])

            obj = obj.unstack(level='gene').T

            self.df_expr = obj

            return None

        elif type(obj) is pd.DataFrame:
            print('Received data in a form of pandas.DataFrame', flush=True)
            print('Validating pandas.DataFrame', flush=True)

            try:
                obj.index.name='gene'
            except:
                print('DataFrame index format is not understood. Returning None', flush=True)
                return None

            if not 'cell' in obj.columns.names:
                print('Columns level "cell" not found. Returning None', flush=True)
                return None

            if not 'batch' in obj.columns.names:
                print('Columns level "batch" not found. Assuming one batch in the data.', flush=True)
                batch = np.array(['batch0']*len(obj.columns.get_level_values('cell')))
            else:
                batch = obj.columns.get_level_values('batch')

            obj.columns = pd.MultiIndex.from_arrays([batch, obj.columns.get_level_values('cell')], names=['batch', 'cell'])

            self.df_expr = obj

            return None

        elif type(obj) is str:
            columns = pd.read_csv(obj, header=0, index_col=None, nrows=5).columns.values.tolist()

            if ('cell' in columns) and ('gene' in columns) and ('expr' in columns):
                print('Received data in a form of condensed matrix. Reading data', flush=True)
                df_expr = pd.read_csv(obj, header=0, index_col=None)

                print('Converting it to pandas.DataFrame')

                if not 'cell' in df_expr.columns:
                    print('The column with "cell" identifiers is not found. Returning None', flush=True)
                    return None

                if not 'gene' in df_expr.columns:
                    print('The column with "gene" identifiers is not found. Returning None', flush=True)
                    return None

                if not 'expr' in df_expr.columns:
                    print('The column with expression values is not found. Returning None', flush=True)
                    return None

                if not 'batch' in df_expr.columns:
                    print('The column with "batch" identifiers is not found. Assuming one batch in the data', flush=True)
                    df_expr['batch'] = np.array(['batch0']*len(df_expr))

                df_expr = df_expr.set_index(['batch', 'cell', 'gene'])['expr'].unstack(level='gene').T

                return df_expr
            else:
                print('Received data in a form of matrix. Reading data', flush=True)
                df_expr = pd.read_csv(obj, header=None, index_col=0)

                print('Converting it to pandas.DataFrame', flush=True)

                if not 'cell' in df_expr.index:
                    print('The row with "cell" identifiers is not found. Returning None', flush=True)
                    return None

                if not 'batch' in df_expr.index:
                    print('The row with "batch" identifiers is not found. Assuming one batch in the data', flush=True)
                    df_expr.loc['batch'] = np.array(['batch0']*df_expr.shape[1])

                df_expr = df_expr.T.set_index(['batch', 'cell']).T

                df_expr.index.name = 'gene'

                self.df_expr = df_expr

                return None
        else:
            print('Unknown input data format. Returning None', flush=True)

        return None

    def convert(self, nameFrom='alias', nameTo='hugo'):

        '''Convert index to hugo names, if any names in the index are
        duplicated, remove duplicates

        Parameters:
            nameFrom: str, Default 'alias'
                Gene name type to convert from

            nameTo: str, Default 'hugo'
                Gene name type to convert to

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.convertIndex()
        '''

        if nameTo=='hugo' and nameFrom=='alias':
            reversed = self.createReverseDictionary(self.gnc.conversionDict[nameTo][nameFrom])
            self.df_expr.index = [reversed[gene][0] if (gene in reversed.keys()) else gene for gene in self.df_expr.index]
        else:
            self.df_expr.index = self.gnc.Convert(list(self.df_expr.index), nameFrom, nameTo, returnUnknownString=False)

        len_total, len_unique = len(self.df_expr.index), len(np.unique(self.df_expr.index))
        if len_total != len_unique:
            unique_items = np.unique(self.df_expr.index, return_counts=True)
            #self.df_expr = self.df_expr.loc[~self.df_expr.index.duplicated(keep='first')]
            self.df_expr = self.df_expr.groupby(level=0, axis=0).sum()
            print('Merged %s duplicated items in the index of size %s' % (len_total - len_unique, len_total), flush=True)
            print(unique_items[0][unique_items[1] > 1], flush=True)

        return None

    def clean(self):

        '''Clean pandas.DataFrame: validate index,
        replace missing with zeros, remove all-zero rows and columns

        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.clean()
        '''

        # Replace any NaN(s) with zeros.
        self.df_expr.fillna(0.0, inplace=True)
        print('Replaced missing values with zeros. Data size: %s genes, %s cells' % self.df_expr.shape, flush=True)

        # Check is any names in the index are duplicated, remove duplicates
        len_total, len_unique = len(self.df_expr.index), len(np.unique(self.df_expr.index))
        if len_total != len_unique:
            unique_items = np.unique(self.df_expr.index, return_counts=True)
            #self.df_expr = self.df_expr.loc[~self.df_expr.index.duplicated(keep='first')]
            self.df_expr = self.df_expr.groupby(level=0, axis=0).sum()
            print('Merged %s duplicated items in the index of size %s' % (len_total - len_unique, len_total), flush=True)
            print(unique_items[0][unique_items[1] > 1], flush=True)

        # Keep only cells with at least one expressed gene.
        self.df_expr = self.df_expr.T[self.df_expr.sum(axis=0) > 0].T
        print('Removed all-zero cells. Data size: %s genes, %s cells' % self.df_expr.shape, flush=True)

        # Keep only genes expressed in at least one cell.
        self.df_expr = self.df_expr[self.df_expr.sum(axis=1) > 0]
        print('Removed all-zero genes. Data size: %s genes, %s cells' % self.df_expr.shape, flush=True)

        return None

    def normalize(self, median=None):

        '''Normalize pandas.DataFrame: rescale all cells, log-transform data,
        remove constant genes, sort index

        Parameters:
            median: float, Default None
                Scale factor, if not provided will be computed as median across all cells in data

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.normalize()
        '''

        # Scale all cells.
        median = np.median(np.sum(self.df_expr,axis=0)).astype(float) if median is None else median
        print('Rescaling all cells by "sum of values = %s".' % (median), flush=True)
        self.df_expr = self.df_expr.apply(lambda q: q * median / np.sum(q),axis=0)

        print('Log-transforming data.', flush=True)
        # Replace zeros with minimum value.
        MIN = np.min(self.df_expr.values[self.df_expr.values > 0.])
        if MIN <= 0.:
            raise ValueError
        self.df_expr = self.df_expr.replace(0., MIN)

        # Take log2 of expression.
        self.df_expr = np.log2(self.df_expr)
        self.df_expr -= np.min(self.df_expr.values)

        # Keep only genes expressed in at least one cell.
        self.df_expr = self.df_expr[self.df_expr.sum(axis=1) > 0]
        print('Removed all-zero genes. Data size: %s genes, %s cells' % self.df_expr.shape, flush=True)  

        if not self.sigmaOverMeanSigma is None:
            # Keep only those genes with large enough standard deviation.
            self.df_expr = self.df_expr[np.std(self.df_expr, axis=1) / np.mean(np.std(self.df_expr.values)) > self.sigmaOverMeanSigma]
            print('Removed constant genes. Data size: %s genes, %s cells' % self.df_expr.shape, flush=True)

        # Sort rows by gene name
        self.df_expr = self.df_expr.sort_index()

        return None
    
    def project(self, PCAonly=False, do_fast_tsne=True):

        '''Project pandas.DataFrame to lower dimensions

        Parameters:
            PCAonly: boolean, Default False
                Perform Principal component analysis only

            do_fast_tsne: boolean, Default True
                Do FI-tSNE instead of "exact" tSNE

        Returns:
            pandas.DataFrame
                Processed data
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            xPCA, PCs, tSNE = DCS.project()
        '''

        print('Performing PC projection from %s to %s features...' % (self.df_expr.shape[0], self.nComponentsPCA), flush=True)
        _PCA = PCA(n_components=self.nComponentsPCA)

        idx = np.argsort(np.var(self.df_expr.values.T, axis=0)/np.mean(self.df_expr.values.T, axis=0))[-2000:]
        X_pca = _PCA.fit_transform(self.df_expr.values.T[:, idx]).T

        print('Explained variance:', np.round(np.sum(_PCA.explained_variance_ratio_) * 100., 2), "%", flush=True)
        
        if not PCAonly:
            if self.layout=='TSNE':
                if X_pca.T.shape[0] < 2000:
                    do_fast_tsne = False

                print('Performing tSNE projection from %s to %s features...' % (self.nComponentsPCA,2), flush=True)
                if do_fast_tsne:
                    if platform.system() == "Windows":
                        from .FastTSNE import fast_tsne
                        print('Windows detected. Using FastTSNE submodule wrapper', flush=True)
                        X_projection2 = fast_tsne(X_pca.T, perplexity = 30, seed=42).T
                    else:
                        import fitsne
                        X_projection2 = fitsne.FItSNE(np.array(X_pca.T, order='C')).T
                else:
                    X_projection2 = TSNE(n_components=2).fit_transform(X_pca.T).T
            elif self.layout=='UMAP':
                print('Performing UMAP projection from %s to %s features...' % (self.nComponentsPCA,2), flush=True)
                import umap
                X_projection2 = umap.UMAP(random_state=42).fit_transform(X_pca.T).T
            elif self.layout=='PCA':
                print('Using PC1 and PC2 for layout', flush=True)
                X_projection2 = X_pca[[0,1],:]

            return X_pca, _PCA.components_, X_projection2

        return X_pca, _PCA.components_
    
    def cluster(self, X_pca, clusteringFunction=AgglomerativeClustering):

        '''Cluster PCA-reduced data into a desired number of clusters

        Parameters:
            X_pca: 2d numpy.array
                PCA-reduced expression data

            clusteringFunction: function, Default AgglomerativeClustering
                Note: the function should have .fit method and same input and output.
                For Network-based clustering pass a dictionary 
                {'k_neighbors':40, metric:'euclidean', 'clusterExpression':True},
                thsi way the best number of clusters will be determined automatically

        Returns:
            1d numpy.array
                Cell cluster index labels
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            index = DCS.cluster(xPCA, 10)
        '''

        if type(clusteringFunction) is dict:

            import pynndescent
            import networkx as nx
            import community
            
            try:
                k_neighbors = clusteringFunction[k_neighbors]
            except:
                k_neighbors = 40

            try:
                metric = clusteringFunction[metric]
            except:
                metric = 'euclidean'
            
            try:
                clusterExpression = clusteringFunction[clusterExpression]
            except:
                clusterExpression = False

            data = self.df_expr.values.T if clusterExpression else X_pca.T

            print('Searching for %s nearest neighbors'%(k_neighbors), flush=True)
            knn = pynndescent.NNDescent(data, metric=metric, n_neighbors=k_neighbors).query(data, k=k_neighbors)

            print('k(=%s) nearest neighbors found. Constructing a NetworkX graph'%(k_neighbors), flush=True)
            A = np.zeros((len(knn[0]),len(knn[0])))
            for i in range(len(knn[0])):
                A[i, knn[0][i]] = knn[1][i]

            G = nx.from_numpy_array(A)

            print('Clustering the graph', flush=True)
            cellClusterIndex = pd.Series(community.best_partition(G)).sort_index().values
        else:
            data = X_pca

            cellClusterIndex = clusteringFunction(n_clusters=self.nClusters).fit(data.T).labels_.astype(float).astype(str)
            print(np.unique(cellClusterIndex, return_counts=True))

            if self.doFineClustering:
                fineCellClusterIndex = cellClusterIndex.copy()

                clusters = np.unique(cellClusterIndex)

                for cluster in clusters:
                    cellsOfCluster = np.where(cellClusterIndex==cluster)[0]
                    subData = data[:, cellsOfCluster]
                    subData = subData[subData.sum(axis=1)>0.]
            
                    if len(cellsOfCluster)>=self.minSizeForFineClustering:
                        try:
                            model = SpectralCoclustering(n_clusters=self.nFineClusters, random_state=0)
                            model.fit(subData.T)
                            tempCellClusterIndex = np.zeros((subData.T.shape[0],))
                            tempCellClusterIndex[:] = np.nan
                            for i, subCluster in enumerate(model.biclusters_[0]):
                                tempCellClusterIndex[np.where(subCluster)[0]] = i
                        except:
                            print('Exception', subData.shape, cellsOfCluster.shape)
                            tempCellClusterIndex = np.zeros((subData.T.shape[0],))
                    else:
                        print('Small cluster', len(cellsOfCluster))
                        tempCellClusterIndex = np.zeros((subData.T.shape[0],))
            
                    fineCellClusterIndex[cellsOfCluster] = np.array([str(cluster) + '.' + label for label in tempCellClusterIndex.astype(int).astype(str)])

                cellClusterIndex = fineCellClusterIndex

        return cellClusterIndex
    
    def vote(self, df_markers_expr, df_marker_cell_type, votingScheme):
        
        '''Produce cluster voting results

        Parameters:
            df_markers_expr: pandas.DataFrame
                Data with marker genes by cells expression

            df_marker_cell_type: pandas.DataFrame
                Data with marker genes by cell types

            votingScheme: function
                Voting function

        Returns:
            dictionary
                Voting results, a dictionary in form of:
                {cluster label: assigned cell type}
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            voting = DCS.vote(df_markers_expr, df_marker_cell_type)
        '''
        
        return votingScheme(df_markers_expr, df_marker_cell_type)
    
    def process(self, visualize=True, cellTypeIdentificationOnly=False):

        '''Main function

        Parameters:
            visualize: boolean, Default True
                Whether to make the plots

            cellTypeIdentificationOnly: boolean, Default False
                Whether to perform cell type identification without recalculating clustering and 2D layout.
                This option is useful to re-run cell types analysis with a new cell type/marker list.
                Processed data will be loaded from the previous calculation. If the data is not found this
                option is ignored.

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()
        '''

        np.random.seed(0)
        
        tempFile = os.path.join(self.saveDir, self.dataName + '_processed.h5')

        if (not cellTypeIdentificationOnly) or (not self.KeyInFile('df_expr', tempFile)) or (not self.KeyInFile('df_clusters', tempFile)):

            if self.saveDir!=os.path.join('') and not os.path.exists(self.saveDir):
                os.makedirs(self.saveDir)
        
            ##############################################################################################
            # Convert index to hugo names, clean data
            ##############################################################################################
            #self.convert(nameFrom=self.geneNamesType)
            if not self.dataIsNormalized:
                self.convert(nameFrom=self.geneNamesType)
                self.clean()
        
            ##############################################################################################
            # Calculate QC measures
            ##############################################################################################
            if self.mitochondrialGenes is None:
                self.mitochondrialGenes = pd.read_csv(os.path.join(self.defualtGeneListsDir, 'Human.MitoCarta2.0.csv'), index_col=None, header=0)['Symbol'].values.squeeze().tolist()
            print('Calculating quality control measures (count depth, number of genes, fraction of mitochondrial genes) for each cell', flush=True)
            df_QC = pd.concat([(self.df_expr).sum(axis=0), (self.df_expr > 0).sum(axis=0), (self.df_expr.loc[self.df_expr.index.intersection(self.mitochondrialGenes)] > 0).sum(axis=0) / (self.df_expr > 0).sum(axis=0)], axis=1, sort=False)
            df_QC.columns = ['count_depth', 'number_of_genes', 'fraction_of_mitochondrialGenes']
            df_QC.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='QC', mode='a', complevel=4, complib='zlib')

            ##############################################################################################
            # Normalize
            ##############################################################################################
            if not self.dataIsNormalized:
                self.normalize()
          
                ##############################################################################################
                # Correct for batch effects
                ##############################################################################################
                if self.toggleDoBatchCorrection:
                    print('ComBat transformation', flush=True)
                    cells = self.df_expr.columns.get_level_values('cell').values.copy()
                    patients = self.df_expr.columns.get_level_values('batch').values.copy()

                    if len(np.unique(patients))==1:
                        print('Only one batch provided. Batch correction is not necessary', flush=True)
                    else:
                        print('Reading positions of zeros', flush=True)
                        where_zeros = self.df_expr == 0.
                        values = combat(pd.DataFrame(data=self.df_expr.values.copy(), index=self.df_expr.index.values.copy(), columns=cells), pd.Series(data=patients, index=cells)).values
                        self.df_expr = pd.DataFrame(data=values, index=self.df_expr.index, columns=self.df_expr.columns)
                        print('Setting original zeros back to zeros', flush=True)
                        self.df_expr[where_zeros] = 0.
                        print('Setting negative values to zeros', flush=True)
                        self.df_expr[self.df_expr < 0.0] = 0.

            ##############################################################################################
            # Compress and record DataFrame to dense matrix
            ##############################################################################################
            if self.toggleRecordAllExpression:
                print('Recording compressed DataFrame', flush=True)
                self.df_expr.replace(0, np.nan).T.stack().to_frame().to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='a', complevel=4, complib='zlib')

            self.df_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r').unstack().T.fillna(0.)
            self.df_expr.index = self.df_expr.index.get_level_values(1)

            ##############################################################################################
            # Calculate PCA and 2D projection
            ##############################################################################################
            print('Preparing xpca, pcs, 2D projection of df_expr', flush=True)
            X_pca, PCs, X_projection = self.project()
            print('Recording xpca, pcs, 2D projection of df_expr', flush=True)
            pd.DataFrame(data=X_pca, columns=self.df_expr.columns).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_xpca', mode='a', complevel=4, complib='zlib')
            pd.DataFrame(data=X_projection, columns=self.df_expr.columns).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection_pre_QC' if self.toggleRemoveLowQualityCells else 'df_projection', mode='a', complevel=4, complib='zlib')

            ##############################################################################################
            # Remove low quality cells
            ##############################################################################################
            index = self.getIndexOfGoodQualityCells(self.saveDir, self.dataName)
            if self.toggleRemoveLowQualityCells:
                self.df_expr = self.df_expr[self.df_expr.columns.intersection(index).sort_values()]
                print('Removed low quality cells. Data size: %s genes, %s cells' % self.df_expr.shape, flush=True)
                X_pca, PCs = self.project(PCAonly=True)
                pd.DataFrame(data=X_pca, columns=self.df_expr.columns).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_xpca', mode='a', complevel=4, complib='zlib')
                #pd.DataFrame(data=PCs, columns=self.df_expr.index).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_pcs', mode='a', complevel=4, complib='zlib')
                df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection_pre_QC', mode='r')
                df_projection[self.df_expr.columns].to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='a', complevel=4, complib='zlib')

            ##############################################################################################
            # Cluster data, append cluster index to expression DataFrame
            ##############################################################################################
            print('Calculating clustering of PCA data', flush=True)
            df_xpca = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_xpca', mode='r')
            cellClusterIndexLabel = self.cluster(df_xpca.values, 
                                                 clusteringFunction=self.clusteringFunction)
            df_clusters = pd.DataFrame(data=cellClusterIndexLabel, index=self.df_expr.columns)
            df_clusters.columns = ['cluster']
            df_clusters.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='a', complevel=4, complib='zlib')

        print('Loading processed data', flush=True)
        self.df_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r').unstack(level=2, fill_value=0).T
        self.df_expr.index = self.df_expr.index.get_level_values(-1)

        df_clusters = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='r')
        self.df_expr = pd.concat([df_clusters, self.df_expr.T], sort=False, axis=1).reset_index().set_index(['batch', 'cell', 'cluster']).T

        ##############################################################################################
        # Get dictionary to map from markers to cell types. Select genes from the marker list only
        ##############################################################################################
        df_marker_cell_type = self.readMarkerFile()
        df_marker_cell_type = df_marker_cell_type.loc[df_marker_cell_type.index.intersection(self.df_expr.index)]
        print('Markers/celltypes:', df_marker_cell_type.shape, flush=True)
        # Remove underepresented cell types
        df_marker_cell_type = df_marker_cell_type[df_marker_cell_type.columns[(df_marker_cell_type>0.).sum(axis=0) >= self.minimumNumberOfMarkersPerCelltype]]
        df_marker_cell_type = df_marker_cell_type[(df_marker_cell_type!=0).sum(axis=1) >= 1]
        df_marker_cell_type.index.name = 'Marker'
        df_marker_cell_type.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_marker_cell_type', mode='a', complevel=4, complib='zlib')
        df_markers_expr = self.df_expr.loc[self.df_expr.index.intersection(df_marker_cell_type.index)]
        print('Selected genes from the marker list. DataFrame shape:', df_markers_expr.shape, flush=True)
        
        ##############################################################################################
        # Vote on cell type.  Use the voting scheme described in the paper if
        # no other scheme is given. Update marker expression with "True" labels and record it
        ##############################################################################################
        votingResults = self.vote(df_markers_expr, df_marker_cell_type.T, self.dcsVotingScheme)
        print(votingResults, flush=True)
        df_markers_expr = df_markers_expr.T.reset_index().T
        df_markers_expr.loc['label'] = np.array([votingResults[i] for i in df_markers_expr.loc['cluster']])
        df_markers_expr = df_markers_expr.T.set_index(['batch', 'cell', 'cluster', 'label']).T.apply(pd.to_numeric)
        df_markers_expr.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='a', complevel=4, complib='zlib')
        
        ##############################################################################################
        # Create a colormap for cell types (based on a marker-celltype matrix)
        ##############################################################################################
        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'w') as temp_file:
            df_marker_cell_type = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_marker_cell_type', mode='r')
            cellTypes = np.sort(df_marker_cell_type.columns.values.tolist())
            for cellTypeIndex in range(len(cellTypes)):
                temp_file.write(cellTypes[cellTypeIndex] + '\t' + str(cm.jet(cellTypeIndex / len(cellTypes))) + '\n')
            temp_file.write('Unknown' + '\t' + str(cm.jet(1.0)) + '\n')

        ##############################################################################################
        # Make all plots following the analysis
        ##############################################################################################
        if visualize:
            self.visualize()

        return None
    
    def visualize(self):

        '''A convenient aggregate of visualization tools of this class.

        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.visualize()
        '''

        ##############################################################################################
        # Plot voting results matrix
        ##############################################################################################
        if self.toggleMakeVotingResultsMatrixPlot:
            print('Making voting results matrix plot')
            self.makeVotingResultsMatrixPlot()

        ##############################################################################################
        # Plot null distributions
        ##############################################################################################
        if self.toggleMakeHistogramNullDistributionPlot:
            print('Making null distributions plot')
            self.makeHistogramNullDistributionPlot()
        
        ##############################################################################################
        # Plot mean marker expression
        ##############################################################################################
        if self.toggleMakeMarkerExpressionPlot:
            print('Making marker expression plot')
            self.makeMarkerExpressionPlot()

        ##############################################################################################
        # Make 2D projection plots
        ##############################################################################################
        if self.toggleMakeProjectionPlotsQC:
            self.getProjectionPlotsQC()

        if self.toggleMakeProjectionPlotClusters:
            self.getProjectionPlotByClusters()

        if self.toggleMakeProjectionPlotBatches:
            self.getProjectionPlotByBatches()

        if self.toggleMakeProjectionPlotAnnotatedClusters:
            self.getProjectionPlotAnnotated()

        if self.toggleAnomalyScoresProjectionPlot:
            self.getAnomalyScoresPlot()
           
        ##############################################################################################
        # Make stacked barplot of cell type fractions
        ##############################################################################################
        if self.toggleMakeStackedBarplot:        
            self.makeStackedBarplot(self.subclusteringName)
        
        ##############################################################################################
        # Make projection plots showing relative expression of different markers (one for each marker)
        ##############################################################################################
        if self.toggleGetMarkerSubplots:
            self.getMarkerSubplots()

        return None


    # User functions of class ###########################################################################################################################
    def getProjectionPlotAnnotated(self):

        '''Produce projection plot colored by cell types
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getProjectionPlotAnnotated()
        '''

        print('Making projection plot by clusters "True" labels')

        df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r')

        df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')
        df_projection = df_projection[df_markers_expr.groupby(level=['batch', 'cell'], sort=False, axis=1).count().columns]

        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}

        self.makeProjectionPlot(df_projection.values, df_markers_expr.columns.get_level_values('label'), 'by_clusters_annotated',
                            colormap=colormap, legend=False)

        return None
    
    def getProjectionPlotByBatches(self):

        '''Produce projection plot colored by batches
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getProjectionPlotByBatches()
        '''

        print('Making projection plot by batches')

        df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')

        self.makeProjectionPlot(df_projection.values, df_projection.columns.get_level_values('batch').values, 'by_patients')

        return None
    
    def getProjectionPlotByClusters(self):

        '''Produce projection plot colored by clusters
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getProjectionPlotByClusters()
        '''

        print('Making projection plot by clusters')

        df_clusters = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='r')
        df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')

        self.makeProjectionPlot(df_projection.values, np.array(['Cluster #%s' % (label[0]) for label in df_clusters.values]), 'by_clusters')

        return None

    def getProjectionPlotsQC(self):

        '''Produce Quality Control projection plots
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getProjectionPlotsQC()
        '''

        print('Making projection plots of QC', flush=True)

        df_QC = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='QC', mode='r')
        if self.toggleRemoveLowQualityCells:
            df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection_pre_QC', mode='r')
        else:
            df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')

        goodQUalityCells = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r').columns

        self.makeProjectionPlot(df_projection.values, np.array([cell in goodQUalityCells for cell in df_projection.columns]), 'by_is_quality_cell', legend=False, labels=True)
        self.makeProjectionPlot(df_projection.values, df_QC['number_of_genes'].values, 'by_number_of_genes', legend=False, labels=False, colorbar=True)
        self.makeProjectionPlot(df_projection.values, df_QC['count_depth'].values, 'by_count_depth', legend=False, labels=False, colorbar=True)
        self.makeProjectionPlot(df_projection.values, df_QC['fraction_of_mitochondrialGenes'].values, 'by_fraction_of_mitochondrialGenes', legend=False, labels=False, colorbar=True)

        return None

    def getMarkerSubplots(self):

        '''Produce subplots on each marker and its expression on all clusters
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getMarkerSubplots()
        '''

        df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')
        df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r')
        df_projection = df_projection[pd.MultiIndex.from_arrays([df_markers_expr.columns.get_level_values('batch'), df_markers_expr.columns.get_level_values('cell')])]
        hugo_cd_dict = dict(zip(df_markers_expr.index.values.tolist(), self.gnc.Convert(list(df_markers_expr.index), 'hugo', 'alias', returnUnknownString=False)))
        self.internalMarkerSubplots(df_markers_expr, df_projection.values, hugo_cd_dict)

        return

    def getAnomalyScoresPlot(self, cells='All'):

        '''Make anomaly scores plot

        Parameters:
            cells: pandas.MultiIndex, Default 'All'
                Index of cells of interest

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            cells = DCS.getCells(celltype='T cell')

            DCS.makeAnomalyScoresPlot(cells)
        '''

        df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')

        if type(cells) is str:
            if cells=='All':
                cells = df_projection.columns
            else:
                print('"cells" value "%s" is unknown'%(cells))
                raise ValueError

        df_sel = self.getExprOfCells(cells)

        scores = self.getAnomalyScores(df_sel, df_sel)
        print('Loaded expression data')

        scores = pd.DataFrame(index=cells, data=scores).reindex(df_projection.columns).values.T[0]

        self.makeProjectionPlot(df_projection.values, scores, suffix='by_anomaly_score', legend=False, labels=False, colorbar=True)

        return None
    
    def getIndividualGeneTtestPlot(self, gene, analyzeBy='label'):

        '''Produce individual gene t-test plot of the two-tailed p-value.

        Parameters:
            gene: str
                Name of gene of interest

            analyzeBy: str, Default 'cluster'
                What level of lablels to include.
                Other possible options are 'label' and 'celltype'

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.getIndividualGeneTtestPlot('SDC1')
        '''

        df_markers_expr = self.getExprOfGene(gene, analyzeBy=analyzeBy)

        if df_markers_expr is None:

            return None

        groups = np.unique(df_markers_expr.columns.get_level_values(analyzeBy).values)

        ttestStatistic = pd.DataFrame(index=groups, columns=groups)
        ttestpValue = pd.DataFrame(index=groups, columns=groups)
        for groupA in groups:
            for groupB in groups:
                A = df_markers_expr.xs(key=groupA, level=analyzeBy, axis=1).values.squeeze()
                B = df_markers_expr.xs(key=groupB, level=analyzeBy, axis=1).values.squeeze()
                ttestStatistic.loc[groupA, groupB], ttestpValue.loc[groupA, groupB] = scipy.stats.ttest_ind(A[np.where(A!=0)], B[np.where(B!=0)])

        alt = self.gnc.Convert([gene], 'hugo', 'alias', returnUnknownString=False)[0]
        alt = [alt] if type(alt) is str else alt
        hugo_cd_dict = '%s\n(%s)'%(gene, ('\n').join(list(alt)))

        self.makeTtestPlot(ttestStatistic, ttestpValue, label=hugo_cd_dict)

        return None    

    def getIndividualGeneExpressionPlot(self, gene, hideClusterLabels=False, outlineClusters=True):

        '''Produce individual gene expression plot on a 2D layout

        Parameters:
            gene: str
                Name of gene of interest

            hideClusterLabels: boolean, Default False
                Whether to hide the clusters labels

            outlineClusters: boolean, Default True
                Whether to outline the clusters with circles

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.getIndividualGeneExpressionPlot('CD4')
        '''

        df_markers_expr = self.getExprOfGene(gene, analyzeBy='cluster')

        if df_markers_expr is None:

            return None

        df_projection = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_projection', mode='r')
        
        hugo_cd_dict = {gene: self.gnc.Convert([gene], 'hugo', 'alias', returnUnknownString=False)[0]}

        self.internalMarkerSubplots(df_markers_expr, 
                                df_projection.values, 
                                hugo_cd_dict, 
                                NoFrameOnFigures=True, 
                                HideClusterLabels=hideClusterLabels, 
                                outlineClusters=outlineClusters)

        return None

    def getAnomalyScores(self, trainingSet, testingSet, printResults=False):

        '''Function to get anomaly score of cells based on some reference set

        Parameters:
            trainingSet: pandas.DataFrame
                With cells to trail isolation forest on

            testingSet: pandas.DataFrame
                With cells to score

            printResults: boolean, Default False
                Whether to print results

        Returns:
            1d numpy.array
                Anomaly score(s) of tested cell(s)
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            cutoff = DCS.checkCellsByIsolationForest(df_expr.iloc[:, 5:], df_expr.iloc[:, :5])
        '''

        instanceIsolationForest = IsolationForest(behaviour='new',
                                                  max_samples=np.min([trainingSet.shape[1]-1, 100]), 
                                                  random_state=np.random.RandomState(None), 
                                                  contamination='auto')

        instanceIsolationForest.fit(trainingSet.values.T)

        if type(testingSet) is pd.Series:
            testingSet = testingSet.to_frame()

        scores = instanceIsolationForest.score_samples(testingSet.values.T)

        scores *= -1.

        results = list(zip(testingSet.columns, scores))

        if printResults:
            for cell in results:
                print('Cell:', cell[0], 'Score:', cell[1])

        return scores

    def getHugoName(self, gene, printAliases=False):
        
        '''Get gene hugo name(s).

        Parameters:
            gene: str
                'hugo' or 'alias' name of a gene

            printAliases: boolean, Default False
                Print aliases of the found hugos

        Returns:
            list
                List containing the input name with any found hugo names

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.getHugoName('CD138')
        '''
        
        conversionDict = self.gnc.conversionDict['hugo']['alias']

        se = pd.Series(conversionDict).apply(','.join).str.split(',', expand=True)

        names = se.index[np.apply_along_axis(lambda arr: gene in arr, axis=1, arr=se.values)].values

        if printAliases:
            try:
                print(gene, ':', conversionDict[gene])
            except:
                pass

            for name in names:
                print(name, ':', conversionDict[name])

        return [gene] + names.tolist()

    def getExprOfGene(self, gene, analyzeBy='cluster'):

        '''Get expression of a gene.
        Run this function only after function process()

        Parameters:
            cells: pandas.MultiIndex
                Index of cells of interest

            analyzeBy: str, Default 'cluster'
                What level of lablels to include.
                Other possible options are 'label' and 'celltype'

        Returns:
            pandas.DataFrame
                With expression of the cells of interest

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getExprOfGene('SDC1')
        '''

        fileName = os.path.join(self.saveDir, self.dataName + '_processed.h5')

        if not os.path.isfile(fileName):
            print('Processed data file not found')
            return
        
        try:
            df_markers_expr = pd.read_hdf(fileName, key='df_expr', mode='r').xs(key=gene, axis=0, level=2).T
        except:
            print('Gene %s not in index'%(gene))
            return

        df_markers_expr.index = [gene]

        labelled = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r').columns
        
        columns = pd.MultiIndex.from_arrays([labelled.get_level_values('batch'), labelled.get_level_values('cell')])
        df_markers_expr = df_markers_expr.reindex(columns, axis=1).fillna(0.)

        if analyzeBy=='celltype':
            analyzeBy = 'label'
            columns = labelled.to_series().reset_index().set_index(['batch', 'cell'])[analyzeBy].loc[df_markers_expr.columns].reset_index().set_index(['batch', 'cell', analyzeBy]).index

            df_markers_expr.columns = pd.MultiIndex.from_arrays([columns.get_level_values('batch'), 
                                                                 columns.get_level_values('cell'), 
                                                                 columns.get_level_values(analyzeBy).str.split(' #', expand=True).get_level_values(0)], names=['batch', 'cell','celltype'])
        else:
            df_markers_expr.columns = labelled.to_series().reset_index().set_index(['batch', 'cell'])[analyzeBy].loc[df_markers_expr.columns].reset_index().set_index(['batch', 'cell', analyzeBy]).index

        return df_markers_expr
    
    def getExprOfCells(self, cells):

        '''Get expression of a set of cells.
        Run this function only after function process()

        Parameters:
            cells: pandas.MultiIndex
                Index of cells of interest

        Returns:
            pandas.DataFrame
                With expression of the cells of interest

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getExprOfCells(cells)
        '''

        df_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r')
        df_expr.index.names = ['batch', 'cell', 'gene']
        df_expr = df_expr.unstack(level='gene', fill_value=0).T
        df_expr.index = df_expr.index.get_level_values('gene')

        return df_expr[cells]

    def getCells(self, celltype=None, clusterIndex=None, clusterName=None):

        '''Get cell annotations in a form of pandas.Series

        Parameters:
            celltype: str, Default None
                Cell type to extract

            clusterIndex: int, Default None
                Cell type to extract

            clusterName: str, Default None
                Cell type to extract

        Returns:
            pandas.MultiIndex
                Index of labelled cells
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            labels = DCS.getCells()
        '''

        df = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r').columns.to_frame()

        if not clusterIndex is None:

            se = df.set_index(['batch', 'cell'])['cluster']

            condition = se==clusterIndex

            if not condition.any():
                print('No cells found')
                return None

            selectedCells = se[condition].index

            print('Selected %s cells from cluster: %s'%(len(selectedCells), clusterIndex))

            return selectedCells

        if not clusterName is None:

            se = df.set_index(['batch', 'cell'])['cluster']

            condition = se==eval(clusterName.split('#')[1])

            if not condition.any():
                print('No cells found')
                return None

            selectedCells = se[condition].index

            print('Selected %s cells from cluster: %s'%(len(selectedCells), clusterName))

            return selectedCells

        se = df.set_index(['batch', 'cell'])['label']

        if celltype is None and clusterIndex is None:

            print('Available cell types:', np.unique(se.str.split(' #', expand=True)[0].values))
            print('Note: To get a specific cell type call this function with a specified cell type')

            return se

        condition = se.str.find(celltype)!=-1

        if not condition.any():
            print('No cells found')
            return None

        selectedCells = se[condition].index

        print('Selected %s cells of type: %s'%(len(selectedCells), celltype))

        return selectedCells

    def getIndexOfGoodQualityCells(self, saveDir, dataName, count_depth_cutoff=0.5, number_of_genes_cutoff=0.5, mitochondrial_genes_cutoff_upper_bound=3.0):

        '''Get index of sells that satisfy the QC criteria

        Parameters:
            saveDir: str
                Directory for output files

            dataName: str 
                Name used in output files

            count_depth_cutoff: float, Default 0.5
                Fraction of median to take as count depth cutoff

            number_of_genes_cutoff: float, Default 0.5
                Fraction of median to take as number of genes cutoff

            mitochondrial_genes_cutoff: float, Default 3.0
                The cutoff is median + standard_deviation * this_parameter

        Returns:
            pandas.Index
                Index of cells
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            index = DCS.get_index_of_good_quality_cells(saveDir, dataName)
        '''

        df_QC = pd.read_hdf(os.path.join(saveDir, dataName + '_processed.h5'), key='QC', mode='r')

        plotsDir = os.path.join(saveDir, 'QC_plots', '')

        if not os.path.exists(plotsDir):
            os.makedirs(plotsDir)

        # Calculate cutoffs
        cutoff_count_depth = self.getQualityControlCutoff(df_QC['count_depth'], count_depth_cutoff, plotPathAndName=os.path.join(plotsDir, '%s_count_depth' % (dataName)), MakeHistogramPlot=True)
        cutoff_number_of_genes = self.getQualityControlCutoff(df_QC['number_of_genes'], number_of_genes_cutoff, plotPathAndName=os.path.join(plotsDir, '%s_number_of_genes' % (dataName)), MakeHistogramPlot=True)
        cutoff_fraction_of_mitochondrialGenes = self.getQualityControlCutoff(df_QC['fraction_of_mitochondrialGenes'], mitochondrial_genes_cutoff_upper_bound, plotPathAndName=os.path.join(plotsDir, '%s_fraction_of_mitochondrialGenes' % (dataName)), mito=True, MakeHistogramPlot=True)

        df_QC['isGoodQuality'] = np.zeros(len(df_QC)).astype(bool)

        mask = (df_QC['count_depth'] > cutoff_count_depth) & \
            (df_QC['number_of_genes'] > cutoff_number_of_genes) & \
            (df_QC['fraction_of_mitochondrialGenes'] < cutoff_fraction_of_mitochondrialGenes)

        return df_QC['isGoodQuality'].index[mask].sort_values()

    def getQualityControlCutoff(self, se, cutoff, mito=False, plotPathAndName=None, MakeHistogramPlot=True):

        '''Function to calculate QC quality cutoff

        Parameters:
            se: pandas.Series 
                With data to analyze

            cutoff: float
                Parameter for calculating the quality control cutoff

            mito: boolean, Default False
                Whether the analysis of mitochondrial genes fraction

            plotPathAndName: str, Default None
                Text to include in the figure title and file name

            MakeHistogramPlot: boolean, Default True
                Whether to make a histogram plot

        Returns:
            float
                Cutoff value
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            cutoff = DCS.getCutoff(se)
        '''

        subset = se[~np.isnan(se.values)].values

        if mito:
            median = np.round(np.median(subset),4)
            std = np.round(np.std(subset),4)

            hist_of_subset = scipy.stats.rv_histogram(np.histogram(subset, bins=100, range=(0, 0.3)))
            hist_data = hist_of_subset._hpdf / 100
            hist_bins = hist_of_subset._hbins

            xs = np.linspace(hist_bins[0], hist_bins[-1], 1000)
            spline_data = np.vstack((xs, UnivariateSpline(hist_bins, hist_data[:-1], k=5, s=0)(xs))).T
            smoothed = scipy.signal.savgol_filter(spline_data.T[1], 101, 3)

            x1 = spline_data.T[0][np.where(spline_data.T[0]>median)[0][0]]
            y1 = smoothed[np.where(spline_data.T[0]>median)[0][0]]

            x2 = spline_data.T[0][np.where(spline_data.T[0]>(median+1.5*std))[0][0]]
            y2 = smoothed[np.where(spline_data.T[0]>(median+1.5*std))[0][0]]

            cutoff = min(median + cutoff * std, x1 - y1 * (x2-x1)/(y2-y1))
        else:
            cutoff = int(cutoff * np.median(subset))

        if MakeHistogramPlot:
            self.makeQualityControlHistogramPlot(subset, cutoff, plotPathAndName=plotPathAndName, mito=mito)
            
        return cutoff

    def alignSeries(self, se1, se2, tagForMissing):

        '''Align two pandas.Series

        Parameters:
            se1: pandas.Series
                Series with the first set of items

            se2: pandas.Series 
                Series with the second set of items

            tagForMissing: str, Default 'Missing'
                Label to assign to non-overlapping items

        Returns:
            pandas.DataFrame
                Contains two aligned pandas.Series
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df = alignSeries(pd.Index(['A', 'B', 'C', 'D']).to_series(), pd.Index(['B', 'C', 'D', 'E', 'F']).to_series())
        '''
        
        se1.index.name = 'index'
        se2.index.name = 'index'

        append = lambda se1, se2: pd.concat([se1, pd.Series(index=se2.index.difference(se1.index), data=[tagForMissing] * len(se2.index.difference(se1.index)))], axis=0, sort=False)

        se1 = append(se1, se2)
        se2 = append(se2, se1)

        se1.name = 'se1'
        se2.name = 'se2'

        return pd.concat((se1, se2.loc[se1.index]), axis=1, sort=True)

    def getCountsDataframe(self, se1, se2, tagForMissing='Missing'):

        '''Get a pandas.DataFrame with cross-counts (overlaps) between two pandas.Series

        Parameters:
            se1: pandas.Series
                Series with the first set of items

            se2: pandas.Series 
                Series with the second set of items

            tagForMissing: str, Default 'Missing'
                Label to assign to non-overlapping items

        Returns:
            pandas.DataFrame
                Contains counts
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df = DCS.getCountsDataframe(se1, se2)
        '''

        df = self.alignSeries(se1, se2, tagForMissing)

        counts = {group[0]:{k:len(v) for k, v in group[1].groupby(by='se1').groups.items()} for group in df.reset_index().drop(columns=['index']).set_index('se2').groupby('se2')}

        df = pd.DataFrame.from_dict(counts).fillna(0.0).astype(int)

        moveTag = lambda df: pd.concat([df.iloc[np.where(df.index != tagForMissing)[0]], df.iloc[np.where(df.index == tagForMissing)[0]]], axis=0, sort=False) if tagForMissing in df.index else df

        return moveTag(moveTag(df.T).T)
    
    def getNewMarkerGenes(self, cluster=None, top=100, zScoreCutoff=None, removeUnknown=False):

        '''Extract new marker genes based on the cluster annotations

        Parameters:
            cluster: int, Default None
                Cluster #, if provided genes of only this culster will be returned

            top: int, Default 100
                Upper bound for number of new markers per cell type

            zScoreCutoff: float, Default 0.3
                Lower bound for a marker z-score to be significant

            removeUnknown: boolean, Default False
                Whether to remove type "Unknown"

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.extractNewMarkerGenes()
        '''

        if zScoreCutoff is None:
            zScoreCutoff = self.zScoreCutoff

        tempFile = os.path.join(self.saveDir, self.dataName + '_processed.h5')

        if not self.KeyInFile('df_clusters', tempFile):

            print('Clusters key not found in HDF. Execute function process() first')

            return None

        print('Loading processed data', flush=True)
        self.df_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r').unstack(level=2, fill_value=0).T
        self.df_expr.index = self.df_expr.index.get_level_values(-1)

        df_clusters = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='r')
        self.df_expr = pd.concat([df_clusters, self.df_expr.T], sort=False, axis=1).reset_index().set_index(['batch', 'cell', 'cluster']).T

        df_gene_cluster_centroids = self.df_expr.groupby(level=['cluster'], sort=True, axis=1).mean()

        df_gene_cluster_centroids_merged = pd.DataFrame()
        df_votingResults = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_voting.xlsx'), sheet_name='z-scores', index_col='cluster')[['# cells in cluster', 'Predicted cell type']]

        groups = pd.concat([df_votingResults['# cells in cluster'], df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0]], sort=False, axis=1).reset_index().set_index(0).groupby(by=0)

        for group in groups:
            se_cluster_size = group[1].reset_index().set_index('cluster')['# cells in cluster']
            se_type = (df_gene_cluster_centroids[se_cluster_size.index] * se_cluster_size.values).sum(axis=1) / se_cluster_size.values.sum()
            se_type.name = group[0]
            df_gene_cluster_centroids_merged = pd.concat([df_gene_cluster_centroids_merged, se_type], sort=False, axis=1)

        df_gene_cluster_centroids_merged = df_gene_cluster_centroids_merged.apply(self.zScoreOfSeries, axis=1)

        if not cluster is None:

            clusterGenes = df_gene_cluster_centroids[cluster].sort_values(ascending=False)
            #clusterGenes = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_gene_cluster_centroids', mode='r')[cluster].sort_values(ascending=False)
            clusterGenes = clusterGenes[clusterGenes >= zScoreCutoff].iloc[:top]

            return {'genes':clusterGenes.index.values.tolist(), 'zscore':clusterGenes.values.tolist()}

        if removeUnknown and 'Unknown' in df_gene_cluster_centroids_merged.columns:
            df_gene_cluster_centroids_merged.drop(columns=['Unknown'], inplace=True)

        df_new_marker_genes = df_gene_cluster_centroids_merged.T

        df_new_marker_genes[df_new_marker_genes<0.] = 0.
        print('New marker cell type shape:', df_new_marker_genes.shape)

        df_new_marker_genes = df_new_marker_genes.T.loc[df_new_marker_genes.max(axis=0)>=zScoreCutoff].T
        print('New marker cell type shape:', df_new_marker_genes.shape)

        df_new_marker_genes = df_new_marker_genes.loc[df_new_marker_genes.max(axis=1)>=zScoreCutoff]
        print('New marker cell type shape:', df_new_marker_genes.shape)

        df_new_marker_list = pd.DataFrame()
        for i, celltype in enumerate(df_new_marker_genes.index.values):
            all = df_new_marker_genes.loc[celltype]
            sel = all[all>1.0].sort_values()[:top]
            print('Candidates of %s:'%(celltype), len(sel))
            df_new_marker_list = pd.concat([df_new_marker_list, sel], sort=False, axis=1)
        df_new_marker_list = df_new_marker_list.fillna(0.).sort_index()
        df_new_marker_list[df_new_marker_list>0.] = 1.0

        if 'Unknown' in df_new_marker_list.columns:
            df_new_marker_list.drop(columns=['Unknown'], inplace=True)

        df_new_marker_list.index.name = 'Marker'
        fileName = os.path.join(self.saveDir, 'new_markers.xlsx')
        print('Recording voting results to:', fileName)
        writer = pd.ExcelWriter(fileName)
        df_new_marker_list.to_excel(writer, 'MarkerCellType')
        df_temp = pd.Series(data=df_new_marker_list.columns, index=df_new_marker_list.columns)
        df_temp.name = 'CellTypeGrouped'
        df_temp.index.name = 'CellType'
        df_temp.to_excel(writer, 'CellTypesGrouped')
        writer.save()

        df_marker_cell_type = pd.read_excel(os.path.join(self.saveDir, os.path.basename(self.geneListFileName)), index_col=0, header=0)
        df_new_marker_genes = df_new_marker_genes[df_new_marker_list.index]

        self.makePlotOfNewMarkers(df_marker_cell_type, df_new_marker_genes)

        return None

    
    # Supporting functions of class ###########################################################################################################################
    @classmethod
    def convertColormap(cls, colormap):

        '''Convert colormap from the form (1.,1.,1.,1.) to 'rgba(255,255,255,1.)'

        Parameters:
            colormap: dictionary
                Colormap to convert

        Returns:
            dictionary
                Converted colomap
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            colormap = DCS.convertColormap(colormap)
        '''

        return {k:'rgba' + str(tuple(list((np.array(v[:3]) * 255).astype(int)) + [v[3]])) for k, v in colormap.items()}

    @classmethod
    def zScoreOfSeries(cls, se):
    
        '''Calculate z-score of pandas.Series and modify the Series in place
    
        Parameters:
            se: pandas.Series
                Series to process

        Returns:
            pandas.Series
                Processed series

        Usage:
            se = zScoreOfSeries(se)
        '''

        se.iloc[:] = scipy.stats.zscore(se.values)

        return se
       
    @classmethod
    def KeyInFile(cls, key, file):

        '''Check is a key exists in a HDF file.
    
        Parameters:
            key: str
                Key name to check

            file: str
                HDF file name to check

        Returns:
            boolean
                True if the key is found
                False otherwise

        Usage:
            KeyInFile('df_expr', 'data/file.h5')
        '''

        if not os.path.isfile(file):
            return False 

        with pd.HDFStore(file) as file:
            return True if "/" + key.strip("/") in file.keys() else False

        return

    @classmethod
    def calculateV(cls, args):

        '''Calculate the vting scores (celltypes by clusters)
    
        Parameters:
            args: tuple
                Tuple of sub-arguments

                df_M: pandas.DataFrame
                    Marker cell type DataFrame

                df_markers_expr: pandas.DataFrame
                    Markers expression DataFrame

                cluster_index: 1d numpy.array
                    Clustering index

                zscore_cutoff: float
                    Z-Score cutoff for a given marker to be significant

        Returns:
            pandas.DataFrame 
                Contains voting scores per celltype per cluster 

        Usage:
            df = getV((df_marker_celltype, df_markers_expr, cluster_index, 0.3))
        '''

        df_M, df_markers_expr, cluster_index, zscore_cutoff = args

        df_expr_cluster_centroids = pd.DataFrame(data=df_markers_expr.values, index=df_markers_expr.index, columns=cluster_index).groupby(level=0, sort=True, axis=1).mean()

        df_Z = pd.DataFrame(data=scipy.stats.zscore(df_expr_cluster_centroids.values, axis=1), index=df_expr_cluster_centroids.index, columns=df_expr_cluster_centroids.columns)

        r = np.mean(df_expr_cluster_centroids.values, axis=1) / np.std(df_expr_cluster_centroids.values, axis=1)

        where_significant_high = df_Z > (zscore_cutoff)*r[:,None]
        df_Z[:] = 0.
        df_Z[where_significant_high] = 1.

        df_V = df_M.dot(df_Z)

        G = ((df_M>0.).dot(df_Z))
        R = ((df_M<0.).dot(df_Z))
        A = (df_Z>0.).sum(axis=0)

        df_V *= (G + R) / A

        df_V = df_V.stack()
        df_V[np.isnan(df_V)] = 0.

        return pd.Series(index=df_V.index, data=df_V.values)

    def dcsVotingScheme(self, df_markers_expr, df_marker_cell_type):
        
        '''Produce cluster voting results

        Parameters:
            df_markers_expr: pandas.DataFrame 
                Data with marker genes by cells expression

            df_marker_cell_type: pandas.DataFrame 
                Data with marker genes by cell types

        Returns:
            dictionary
                Voting results, a dictionary in form of:
                {cluster label: assigned cell type}
        
        Usage:
            Function should be called internally only
        '''

        #def norm(s):
            
        #    pos_sum = np.abs(s.iloc[np.where(s>0.)[0]].sum())
        #    neg_sum = np.abs(s.iloc[np.where(s<0.)[0]].sum())

        #    s.iloc[np.where(s>0.)[0]] /= pos_sum if pos_sum>0. else 1.
        #    s.iloc[np.where(s<0.)[0]] /= neg_sum if neg_sum>0. else 1.

        #    return s

        #df_marker_cell_type[df_marker_cell_type<0.] = 0.

        ### Normalize markers by the number of cell types expressing that marker (i.e. the "specificity")
        ##df_marker_cell_type = df_marker_cell_type.apply(lambda q: q / (np.sum(q[q==1]) if np.sum(q[q==1])>0. else 1.) , axis=0)
        ## Normalize cell types so that the absolute number of known markers in a given cell type is irrelevant
        #df_marker_cell_type = df_marker_cell_type.apply(norm, axis=1)

        df_positive = df_marker_cell_type.copy()
        df_negative = df_marker_cell_type.copy()
        df_positive[df_marker_cell_type<.0] = 0.
        df_negative[df_marker_cell_type>.0] = 0.
        df_positive /= (df_positive>0.).sum(axis=0).fillna(1.).replace(0., 1.)
        df_negative /= (df_negative<0.).sum(axis=0).fillna(1.).replace(0., 1.)

        df_marker_cell_type = df_positive + df_negative

        print(df_marker_cell_type)

        # Align df_markers_expr and df_marker_cell_type
        df_markers_expr.sort_index(inplace=True, axis=0)
        df_marker_cell_type.sort_index(inplace=True, axis=1)

        # Calculate score (Vkc) for the best clustering
        df_V = self.calculateV((df_marker_cell_type, df_markers_expr, df_markers_expr.columns.get_level_values('cluster').values, self.zScoreCutoff)).unstack()
        df_V.index.name = None

        df_markers_expr_copy = df_markers_expr.copy()

        uclusters = np.unique(df_markers_expr_copy.columns.get_level_values(2))
        uindex = pd.Series(df_markers_expr_copy.columns.get_level_values(2)).replace(dict(zip(uclusters, range(len(uclusters))))).values
        df_markers_expr_copy.columns = pd.MultiIndex.from_arrays([df_markers_expr_copy.columns.get_level_values(0),
                                                                  df_markers_expr_copy.columns.get_level_values(1),
                                                                  uindex], names=df_markers_expr_copy.columns.names)

        # Generate random cluster index
        randomClusterIndex = np.vstack([np.random.choice(df_markers_expr_copy.columns.get_level_values('cluster'), size=df_markers_expr_copy.shape[1], replace=False) for i in range(self.nSamplesDistribution)]).astype(int)

        # Generate random cluster configurations and calculate scores (Pkc) of those
        print('Generating null distribution')

        def process_batch(batch_range):

            print('\t', batch_range)

            tuples =  [(df_marker_cell_type, df_markers_expr_copy, randomClusterIndex[i], self.zScoreCutoff) for i in batch_range]

            pool = multiprocessing.Pool(processes = self.availableCPUsCount)
            temp_random_df_V = pd.concat(pool.map(self.calculateV, tuples), sort=False, axis=1)
            pool.close()
            pool.join()

            return temp_random_df_V

        batch_size = self.availableCPUsCount * 80
        N_batches = int(self.nSamplesDistribution / batch_size)

        random_df_V = pd.concat([process_batch(range(i*batch_size,(i+1)*batch_size)) for i in range(N_batches)], sort=False, axis=1)

        random_df_V.index = pd.MultiIndex.from_arrays([random_df_V.index.get_level_values(0),
                                                pd.Series(random_df_V.index.get_level_values(1)).replace(dict(zip(range(len(uclusters)), uclusters))).values], 
                                                names=[random_df_V.index.names[0], 'cluster'])

        min_value = np.nanmin(random_df_V.values)
        max_value = np.nanmax(random_df_V.values)

        print('Min:', min_value, '\t', 'Max:', max_value)

        Nbins = 300

        # Calculate null distribution histograms data for plots
        df_null_distributions = random_df_V.apply(lambda s: scipy.stats.rv_histogram(np.histogram(s, bins=Nbins, range=(min_value,max_value)))._hpdf / Nbins, axis=1).apply(pd.Series).T
        df_null_distributions.columns.names = ['CellType', 'Cluster']
        df_null_distributions.index = [min_value] + (scipy.stats.rv_histogram(np.histogram(random_df_V.iloc[0], bins=Nbins, range=(min_value,max_value)))._hbins).tolist()

        # Calculate z-score (Lkc) for the best clustering
        print('Processing voting results')
        sigma = pd.Series(data=np.std(random_df_V.values, axis=1), index=random_df_V.index).replace(0., np.inf).unstack()
        mean = pd.Series(data=np.mean(random_df_V.values, axis=1), index=random_df_V.index).unstack()
        df_L = (df_V - mean) / sigma

        df_L.index.name = None

        # Determine winning score
        winning_score = np.array([df_V.iloc[ind[0], ind[1]] for ind in np.flip(np.array(list(enumerate(np.argmax(df_L.values, axis=0)))), axis=1)])

        # Determine winning cell type
        wtypes = df_L.index[np.argmax(df_L.values, axis=0)].values.copy()

        # Properly rename winning cell type
        T = df_L.index[np.argmax(df_L.values, axis=0)].values
        T[(df_L.values < self.minimumScoreForUnknown).all(axis=0)] = 'Unknown'
        T = pd.Index(T) + ' #0'
        for i in range(len(T)):
            T = T.where(~T.duplicated(), T.str.replace(' #%s' % (i), ' #%s' % (i + 1)))
        predicted_celltype = T.values

        # Determine cluster sizes
        cluster_sizes = df_markers_expr.groupby(level='cluster', sort=True, axis=1).count().iloc[0].values

        # Save score columns names
        columns_scores = df_V.index.values.copy().tolist()

        # Update the DataFrames
        df_V.loc['Winning score'] = df_L.loc['Winning score'] = winning_score
        df_V.loc['Predicted cell type'] = df_L.loc['Predicted cell type'] = predicted_celltype
        df_V.loc['# cells in cluster'] = df_L.loc['# cells in cluster'] = cluster_sizes

        # Determine all markers hits for each cluster
        df_custer_centroids = df_markers_expr.groupby(level='cluster', sort=True, axis=1).mean()
        
        r = np.mean(df_custer_centroids.values, axis=1) / np.std(df_custer_centroids.values, axis=1)
        #df_custer_centroids_sig = df_custer_centroids.apply(self.zScoreOfSeries, axis=1) > (self.zScoreCutoff)*0.3/r[:,None]
        df_custer_centroids_sig = df_custer_centroids.apply(self.zScoreOfSeries, axis=1) > (self.zScoreCutoff)*r[:,None]
        #df_custer_centroids_sig = df_custer_centroids.apply(self.zScoreOfSeries, axis=1) > (self.zScoreCutoff)*(3. - r[:,None])/3.

        dict_supporting_markers = {cluster:df_custer_centroids_sig.index[df_custer_centroids_sig[cluster] > 0].values for cluster in df_custer_centroids_sig}
        all_markers_in_cluster = [(' // ').join(dict_supporting_markers[cluster].tolist()) for cluster in df_V.columns]
        df_V.loc['All markers'] = df_L.loc['All markers'] = all_markers_in_cluster

        all_markers = {celltype:df_marker_cell_type.T.index[df_marker_cell_type.T[celltype] > 0].values for celltype in columns_scores}
        all_markers['Unknown'] = np.array([])

        supporting_markers = [(' // ').join((np.intersect1d(all_markers[df_V.loc['Predicted cell type'].loc[cluster].split(' #')[0]],
                                                            df_V.loc['All markers'].loc[cluster].split(' // '))).tolist()) for cluster in df_V.columns]
        df_V.loc['Supporting markers'] = df_L.loc['Supporting markers'] = supporting_markers
        
        # Record the DataFrames
        fileName = os.path.join(self.saveDir, self.dataName + '_voting.xlsx')
        print('Recording voting results to:', fileName)
        writer = pd.ExcelWriter(fileName)
        df_L.columns.name = 'cluster'
        df_L.T.to_excel(writer, 'z-scores', columns=['Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers'] + columns_scores)
        df_V.columns.name = 'cluster'
        df_V.T.to_excel(writer, 'Voting scores', columns=['Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers'] + columns_scores)
        df_custer_centroids.T.to_excel(writer, 'Cluster centroids')
        df_null_distributions.to_excel(writer, 'Null distributions')
        df_marker_cell_type.to_excel(writer, 'Marker cell type weight matrix')
        writer.save()

        return df_V.loc['Predicted cell type'].to_dict()
    
    def dcsVotingSchemeEEH(self, df_markers_expr, df_marker_cell_type):
        
        '''Produce cluster voting results

        Parameters:
            df_markers_expr: pandas.DataFrame 
                Data with marker genes by cells expression

            df_marker_cell_type: pandas.DataFrame 
                Data with marker genes by cell types

        Returns:
            dictionary
                Voting results, a dictionary in form of:
                {cluster label: assigned cell type}
        
        Usage:
            Function should be called internally only
        '''

        df_marker_cell_type[~(df_marker_cell_type==1.)] = 0.

        means = []
        temp_df = pd.DataFrame(data=df_markers_expr.values, index=df_markers_expr.index, columns=df_markers_expr.columns.get_level_values('cluster'))
        for i in range(100):
            print(i, end=' ', flush=True)
            np.random.shuffle(temp_df.columns.values)
            means.append(temp_df.groupby(level=0, sort=True, axis=1).mean().stack().values)
        means = np.array(means)
        print()

        df_V = df_markers_expr.groupby(level=2, sort=True, axis=1).mean().stack()

        df_R = (pd.DataFrame(data=means, columns=df_V.index).apply(self.zScoreOfSeries, axis=0) > self.zScoreCutoff)*1.
        df_R = df_R.xs(key='0.0', level=1, axis=1)
        print(df_R)

        def norm(s):
            
            pos_sum = np.abs(s.iloc[np.where(s>0.)[0]].sum())
            neg_sum = np.abs(s.iloc[np.where(s<0.)[0]].sum())

            s.iloc[np.where(s>0.)[0]] /= pos_sum if pos_sum>0. else 1.
            s.iloc[np.where(s<0.)[0]] /= neg_sum if neg_sum>0. else 1.

            return s

        # Normalize markers by the number of cell types expressing that marker (i.e. the "specificity")
        df_marker_cell_type = df_marker_cell_type.apply(lambda q: q / (np.sum(q[q==1]) if np.sum(q[q==1])>0. else 1.) , axis=0)
        # Normalize cell types so that the absolute number of known markers in a given cell type is irrelevant
        df_marker_cell_type = df_marker_cell_type.apply(norm, axis=1)
        df_marker_cell_type = df_marker_cell_type[df_R.columns]
        print(df_marker_cell_type)


        df_Pr = df_marker_cell_type.dot(df_R.T)
        print(df_Pr)

        v0 = df_marker_cell_type.dot(df_V.to_frame().xs(key='0.0', level=1, axis=0))
        (v0.T - np.mean(df_Pr, axis=1)) / np.std(df_Pr, axis=1)



        df_L = (df_V - np.mean(means, axis=0)) / np.std(means, axis=0)
        df_L = df_L.unstack().fillna(0.).replace(np.inf, 0.)

        #df_L[df_L<0.] = 0.
        #df_L = df_L.replace(0., np.nan)
        print('\n', df_L)


        genes1 = df_marker_cell_type.columns[df_marker_cell_type.loc['T Cell']==1].values
        genes2 = df_marker_cell_type.columns[df_marker_cell_type.loc['Erythrocyte']==1].values
        genesIn0 = (df_L['0.0']>0.).index[(df_L['0.0']>0.)].values
        w1 = len(np.intersect1d(genesIn0, genes1)) # 244 / 67 = 3.6
        w2 = len(np.intersect1d(genesIn0, genes2)) # 17 / 5 = 3.4

        df_Score = df_marker_cell_type[df_L.index].dot(df_L)/df_marker_cell_type.sum(axis=1).values[:, None]
        print(df_Score)

        res = list(zip(df_Score.index[np.argmax(df_Score.values, axis=0)].values, df_markers_expr.groupby(level=2, sort=True, axis=1).count().iloc[0].values))
        print(res)

    def createReverseDictionary(self, inputDictionary):

        """Efficient way to create a reverse dictionary from a dictionary.
        Utilizes Pandas.Dataframe.groupby and Numpy arrays indexing.
    
        Parameters: 
            inputDictionary: dictionary
                Dictionary to reverse

        Returns:
            dictionary
                Reversed dictionary

        Usage:
            revDict = createReverseDictionary(Dict)
        """

        keys, values = np.array(list(inputDictionary.keys())), np.array(list(inputDictionary.values()))
        df = pd.DataFrame(np.array([[keys[i], value] for i in range(len(keys)) for value in values[i]]))
        dfGrouped = df.groupby(df.columns[1])
        keys, values = list(dfGrouped.indices.keys()), list(dfGrouped.indices.values())
        GOs = df.values.T[0]

        return dict(zip(keys, [GOs[value].tolist() for value in values]))

    def readMarkerFile(self):

        '''Read markers file, prepare markers

        Parameters: 
            None

        Returns:
            pandas.DataFrame
                Celltype/markers matrix

        Usage:
            df_marker_cell_type = readMarkerFile()
        '''

        df_marker_cell_type = pd.read_excel(self.geneListFileName, index_col=0, header=[0,1]).replace(0., np.nan)
        df_marker_cell_type.columns.names = ['CellTypeGrouped', 'CellType']

        reversed = self.createReverseDictionary(self.gnc.conversionDict['hugo']['alias'])

        df_marker_cell_type.index = [reversed[gene][0] if (gene in reversed.keys()) else gene for gene in df_marker_cell_type.index]

        df_marker_cell_type = df_marker_cell_type.drop(columns=[col for col in df_marker_cell_type.columns if col[0]=='NA'])
        df_marker_cell_type = df_marker_cell_type.groupby(level='CellTypeGrouped', axis=1).max(skipna=False)

        df_marker_cell_type = df_marker_cell_type.fillna(0.)
        print('Markers/celltypes:', df_marker_cell_type.shape, flush=True)

        # Merge duplicates that might have appeared after gene name conversion
        df_marker_cell_type = df_marker_cell_type.groupby(level=0, axis=0).sum()
        df_marker_cell_type[df_marker_cell_type>0.] = 1.
        df_marker_cell_type[df_marker_cell_type<0.] = -1.
        print('Markers/celltypes:', df_marker_cell_type.shape, flush=True)

        def norm(s):
            
            pos_sum = np.abs(s.iloc[np.where(s>0.)[0]].sum())
            neg_sum = np.abs(s.iloc[np.where(s<0.)[0]].sum())

            s.iloc[np.where(s>0.)[0]] /= pos_sum if pos_sum>0. else 1.
            s.iloc[np.where(s<0.)[0]] /= neg_sum if neg_sum>0. else 1.

            return s

        df_marker_cell_type = df_marker_cell_type.apply(norm, axis=0)

        return df_marker_cell_type
