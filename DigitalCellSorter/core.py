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
       1. t-SNE layout plot
       2. Quality Control histogram plot
       3. Marker expression t-SNE subplot
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
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.ensemble import IsolationForest

from . import GeneNameConverter, geneLists
from .Combat import combat
from .VisualizationFunctions import VisualizationFunctions
from .GenericFunctions import read, write, timeMark, getStartTime, getElapsedTime
from .VisualizationFunctions import cm

class DigitalCellSorter(VisualizationFunctions):

    '''Class of Digital Cell Sorter with methods for processing single cell mRNA-seq data.
    Includes analyses and visualization tools.

    Parameters:
        dataName: str, Default 'dataName'
            Name used in output files

        geneListFileName: str, Default None
            Name of the marker genes file

        mitochondrialGenes: list, Default None
            List of mitochondrial genes to use in quality control

        sigmaOverMeanSigma: float, Default 0.3
            Threshold to consider a gene constant

        nClusters: 
            Number of clusters, Default 5

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

    def __init__(self, dataName='dataName', geneListFileName=None, mitochondrialGenes=None,
                sigmaOverMeanSigma=0.3, nClusters=5, clusteringFunction=AgglomerativeClustering, nComponentsPCA=200, nSamplesDistribution=10000, 
                saveDir=os.path.join(''), makeMarkerSubplots=False, availableCPUsCount=os.cpu_count(), zScoreCutoff=0.3,
                subclusteringName=None, doQualityControl=True, doBatchCorrection=True, makePlots=True,
                minimumNumberOfMarkersPerCelltype=10, minimumScoreForUnknown=0.3):

        '''Initialization function. Automatically called when an instance on Digital Cell Sorter is created'''

        defaultGeneList = 'CIBERSORT' # 'CIBERSORT' 'markersDCS'

        self.saveDir = saveDir
        self.dataName = dataName
        
        self.defualtGeneListsDir = os.path.join(os.path.dirname(__file__), 'geneLists')

        self.geneListFileName = os.path.join(self.defualtGeneListsDir, defaultGeneList + '.xlsx') if geneListFileName is None else geneListFileName

        self.mitochondrialGenes = mitochondrialGenes

        self.sigmaOverMeanSigma = sigmaOverMeanSigma
        self.nClusters = nClusters
        self.nComponentsPCA = nComponentsPCA
        self.zScoreCutoff = zScoreCutoff
        self.minimumNumberOfMarkersPerCelltype = minimumNumberOfMarkersPerCelltype

        self.minimumScoreForUnknown = minimumScoreForUnknown

        self.nSamplesDistribution = nSamplesDistribution
        self.availableCPUsCount = availableCPUsCount

        self.clusteringFunction = clusteringFunction

        self.subclusteringName = subclusteringName

        self.toggleRecordAllExpression = True

        self.toggleDoQualityControl  =  doQualityControl
        self.toggleDoBatchCorrection = doBatchCorrection

        self.toggleMakeHistogramNullDistributionPlot = makePlots
        self.toggleMakeVotingResultsMatrixPlot = makePlots
        self.toggleMakeMarkerExpressionPlot = makePlots
        self.toggleMakeTSNEplotQC = makePlots
        self.toggleMakeTSNEplotClusters = makePlots
        self.toggleMakeTSNEplotAnnotatedClusters = makePlots
        self.toggleMakeTSNEplotBatches = makePlots
        self.toggleMakeStackedBarplot = makePlots
        self.toggleAnomalyScoresTSNEplot = False

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
            pandas.DataFrame
                Pandas DataFrame
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_expr = DCS.preapre('data.csv')
        '''

        if type(obj) is pd.Series:
            print('Received data in a form of pandas.Series')
            print('Validating it and converting it to pandas.DataFrame')

            if not 'cell' in obj.index.names:
                print('Column "cell" not found. Returning None')
                return None

            if not 'gene' in obj.index.names:
                print('Column "gene" not found. Returning None')
                return None

            if not 'batch' in obj.index.names:
                print('Column "batch" not found. Assuming one batch in the data.')
                batch = np.array(['batch0']*len(obj.index.get_level_values('cell')))
            else:
                batch = obj.index.get_level_values('batch')

            obj.index = pd.MultiIndex.from_arrays([batch, obj.index.get_level_values('cell'), obj.index.get_level_values('gene')], names=['batch', 'cell', 'gene'])

            obj = obj.unstack(level='gene').T

            return obj

        elif type(obj) is pd.DataFrame:
            print('Received data in a form of pandas.DataFrame')
            print('Validating pandas.DataFrame')

            try:
                obj.index.name='gene'
            except:
                print('DataFrame index format is not understood. Returning None')
                return None

            if not 'cell' in obj.columns.names:
                print('Columns level "cell" not found. Returning None')
                return None

            if not 'batch' in obj.columns.names:
                print('Columns level "batch" not found. Assuming one batch in the data.')
                batch = np.array(['batch0']*len(obj.columns.get_level_values('cell')))
            else:
                batch = obj.columns.get_level_values('batch')

            obj.columns = pd.MultiIndex.from_arrays([batch, obj.columns.get_level_values('cell')], names=['batch', 'cell'])

            return obj

        elif type(obj) is str:
            columns = pd.read_csv(obj, header=0, index_col=None, nrows=5).columns.values.tolist()

            if ('cell' in columns) and ('gene' in columns) and ('expr' in columns):
                print('Received data in a form of condensed matrix. Reading data')
                df_expr = pd.read_csv(obj, header=0, index_col=None)

                print('Converting it to pandas.DataFrame')

                if not 'cell' in df_expr.columns:
                    print('The column with "cell" identifiers is not found. Returning None')
                    return None

                if not 'gene' in df_expr.columns:
                    print('The column with "gene" identifiers is not found. Returning None')
                    return None

                if not 'expr' in df_expr.columns:
                    print('The column with expression values is not found. Returning None')
                    return None

                if not 'batch' in df_expr.columns:
                    print('The column with "batch" identifiers is not found. Assuming one batch in the data')
                    df_expr['batch'] = np.array(['batch0']*len(df_expr))

                df_expr = df_expr.set_index(['batch', 'cell', 'gene'])['expr'].unstack(level='gene').T

                return df_expr
            else:
                print('Received data in a form of matrix. Reading data')
                df_expr = pd.read_csv(obj, header=None, index_col=0)

                print('Converting it to pandas.DataFrame')

                if not 'cell' in df_expr.index:
                    print('The row with "cell" identifiers is not found. Returning None')
                    return None

                if not 'batch' in df_expr.index:
                    print('The row with "batch" identifiers is not found. Assuming one batch in the data')
                    df_expr.loc['batch'] = np.array(['batch0']*df_expr.shape[1])

                df_expr = df_expr.T.set_index(['batch', 'cell']).T

                df_expr.index.name = 'gene'

                return df_expr
        else:
            print('Unknown input data format. Returning None')

        return None

    def convert(self, df_expr, nameFrom='alias', nameTo='hugo'):

        '''Convert index to hugo names, if any names in the index are
        duplicated, remove duplicates

        Parameters:
            df_expr: pandas.DataFrame
                Data to process

            nameFrom: str, Default 'alias'
                Gene name type to convert from

            nameTo: str, Default 'hugo'
                Gene name type to convert to

        Returns:
            pandas.DataFrame
                Processed data

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_expr = DCS.convertIndex(df_expr)
        '''

        df_expr.index = self.gnc.Convert(list(df_expr.index), 'alias', 'hugo', returnUnknownString=False)
        len_total, len_unique = len(df_expr.index), len(np.unique(df_expr.index))
        if len_total != len_unique:
            unique_items = np.unique(df_expr.index, return_counts=True)
            df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
            print('Removed %s duplicated items in the index of size %s' % (len_total - len_unique, len_total))
            print(unique_items[0][unique_items[1] > 1])

        return df_expr

    def clean(self, df_expr):

        '''Clean pandas.DataFrame: validate index,
        replace missing with zeros, remove all-zero rows and columns

        Parameters:
            df_expr: pandas.DataFrame
                Data to clean

        Returns:
            pandas.DataFrame
                Processed data
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_expr = DCS.clean(df_expr)
        '''

        # Check is any names in the index are duplicated, remove duplicates
        len_total, len_unique = len(df_expr.index), len(np.unique(df_expr.index))
        if len_total != len_unique:
            unique_items = np.unique(df_expr.index, return_counts=True)
            df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]
            print('Removed %s duplicated items in the index of size %s' % (len_total - len_unique, len_total))
            print(unique_items[0][unique_items[1] > 1])

        # Replace any NaN(s) with zeros.
        df_expr.fillna(0.0, inplace=True)
        print('Replaced missing values with zeros. Data size: %s genes, %s cells' % df_expr.shape)

        # Keep only cells with at least one expressed gene.
        df_expr = df_expr.T[df_expr.sum(axis=0) > 0].T
        print('Removed all-zero cells. Data size: %s genes, %s cells' % df_expr.shape)

        # Keep only genes expressed in at least one cell.
        df_expr = df_expr[df_expr.sum(axis=1) > 0]
        print('Removed all-zero genes. Data size: %s genes, %s cells' % df_expr.shape)

        return df_expr

    def normalize(self, df_expr, median=None):

        '''Normalize pandas.DataFrame: rescale all cells, log-transform data,
        remove constant genes, sort index

        Parameters:
            df_expr: pandas.DataFrame
                Data to normalize

            median: float, Default None
                Scale factor, if not provided will be computed as median across all cells in data

        Returns:
            pandas.DataFrame
                Processed data
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_expr = DCS.normalize(df_expr, 0.3)
        '''

        # Scale all cells.
        median = np.median(np.sum(df_expr,axis=0)).astype(float) if median is None else median
        print('Rescaling all cells by "sum of values = %s".' % (median))
        df_expr = df_expr.apply(lambda q: q * median / np.sum(q),axis=0)

        print('Log-transforming data.')
        # Replace zeros with minimum value.
        MIN = np.min(df_expr.values[df_expr.values > 0.])
        if MIN <= 0.:
            raise ValueError
        df_expr = df_expr.replace(0., MIN)

        # Take log2 of expression.
        df_expr = np.log2(df_expr)
        df_expr -= np.min(df_expr.values)

        # Keep only genes expressed in at least one cell.
        df_expr = df_expr[df_expr.sum(axis=1) > 0]
        print('Removed all-zero genes. Data size: %s genes, %s cells' % df_expr.shape)  

        if not self.sigmaOverMeanSigma is None:
            # Keep only those genes with large enough standard deviation.
            df_expr = df_expr[np.std(df_expr, axis=1) / np.mean(np.std(df_expr.values)) > self.sigmaOverMeanSigma]
            print('Removed constant genes. Data size: %s genes, %s cells' % df_expr.shape)

        # Sort rows by gene name
        df_expr = df_expr.sort_index()

        return df_expr
    
    def project(self, df_expr, PCAonly=False, do_fast_tsne=True):

        '''Project pandas.DataFrame to lower dimensions

        Parameters:
            df_expr: pandas.DataFrame
                Data to process

            PCAonly: boolean, Default False
                Perform Principal component analysis only

            do_fast_tsne: boolean, Default True
                Do FI-tSNE instead of "exact" tSNE

        Returns:
            pandas.DataFrame
                Processed data
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            xPCA, PCs, tSNE = DCS.project(df_expr, 100)
        '''

        print('Performing PC projection from %s to %s features...' % (df_expr.shape[0], self.nComponentsPCA))
        _PCA = PCA(n_components=self.nComponentsPCA)
        X_pca = _PCA.fit_transform(df_expr.values.T).T

        print('Explained variance:', np.round(np.sum(_PCA.explained_variance_ratio_) * 100., 2), "%")
        
        if not PCAonly:
            print('Performing tSNE projection from %s to %s features...' % (self.nComponentsPCA,2))
            if do_fast_tsne:
                if platform.system() == "Windows":
                    from .FastTSNE import fast_tsne
                    print('Windows detected. Using FastTSNE submodule wrapper')
                    X_tsne2 = fast_tsne(X_pca.T, perplexity = 30, seed=42).T
                else:
                    import fitsne
                    X_tsne2 = fitsne.FItSNE(np.array(X_pca.T, order='C')).T
            else:
                X_tsne2 = TSNE(n_components=2).fit_transform(X_pca.T).T

            return X_pca, _PCA.components_, X_tsne2

        return X_pca, _PCA.components_
    
    def cluster(self, X_pca, df_expr=None, clusteringFunction=AgglomerativeClustering):

        '''Cluster PCA-reduced data into a desired number of clusters

        Parameters:
            X_pca: 2d numpy.array
                PCA-reduced expression data

            df_expr: 2d numpy.array, Default None
                Gene expression data 

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

            data = df_expr.values.T if clusterExpression else X_pca.T

            print('Searching for %s nearest neighbors'%(k_neighbors))
            knn = pynndescent.NNDescent(data, metric=metric, n_neighbors=k_neighbors).query(data, k=k_neighbors)

            print('k(=%s) nearest neighbors found. Constructing a NetworkX graph'%(k_neighbors))
            A = np.zeros((len(knn[0]),len(knn[0])))
            for i in range(len(knn[0])):
                A[i, knn[0][i]] = knn[1][i]

            G = nx.from_numpy_array(A)

            print('Clustering the graph')
            cellClusterIndex = pd.Series(community.best_partition(G)).sort_index().values
        else:
            cellClusterIndex = clusteringFunction(n_clusters=self.nClusters).fit(X_pca.T).labels_

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
    
    def process(self, df_expr, visualize=True):

        '''Main function

        Parameters:
            df_expr: pandas.DataFrame
                Gene expression data

            visualize: boolean, Default True
                Whether to make the plots

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process(df_expr)
        '''

        np.random.seed(0)
        
        if self.saveDir!=os.path.join('') and not os.path.exists(self.saveDir):
            os.makedirs(self.saveDir)

        ##############################################################################################
        # Convert index to hugo names, clean data
        ##############################################################################################
        df_expr = self.convert(df_expr)
        df_expr = self.clean(df_expr)
        timeMark()
        
        ##############################################################################################
        # Calculate QC measures
        ##############################################################################################
        if self.toggleDoQualityControl:
            if self.mitochondrialGenes is None:
                self.mitochondrialGenes = pd.read_csv(os.path.join(self.defualtGeneListsDir, 'Human.MitoCarta2.0.csv'), index_col=None, header=0)['Symbol'].values.squeeze().tolist()
            print('Calculating quality control measures (count depth, number of genes, fraction of mitochondrial genes) for each cell')
            df_QC = pd.concat([(df_expr).sum(axis=0), (df_expr > 0).sum(axis=0), (df_expr.loc[self.mitochondrialGenes] > 0).sum(axis=0) / (df_expr > 0).sum(axis=0)], axis=1, sort=False)
            df_QC.columns = ['count_depth', 'number_of_genes', 'fraction_of_mitochondrialGenes']
            df_QC.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='QC', mode='a', complevel=4, complib='zlib')
            timeMark()

        ##############################################################################################
        # Normalize
        ##############################################################################################
        df_expr = self.normalize(df_expr)
        timeMark()
          
        ##############################################################################################
        # Correct for batch effects
        ##############################################################################################
        if self.toggleDoBatchCorrection:
            print('ComBat transformation')
            cells = df_expr.columns.get_level_values('cell').values.copy()
            patients = df_expr.columns.get_level_values('batch').values.copy()

            if len(np.unique(patients))==1:
                print('Only one batch provided. Batch correction is unnecessary')
            else:
                print('Reading positions of zeros')
                where_zeros = df_expr == 0.
                values = combat(pd.DataFrame(data=df_expr.values.copy(), index=df_expr.index.values.copy(), columns=cells), pd.Series(data=patients, index=cells)).values
                df_expr = pd.DataFrame(data=values, index=df_expr.index, columns=df_expr.columns)
                print('Setting original zeros back to zeros')
                df_expr[where_zeros] = 0.
                print('Setting negative values to zeros')
                df_expr[df_expr < 0.0] = 0.
                timeMark()

        ##############################################################################################
        # Compress and record DataFrame to dense matrix
        ##############################################################################################
        if self.toggleRecordAllExpression:
            print('Recording compressed DataFrame')
            df_expr.replace(0, np.nan).T.stack().to_frame().to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='a', complevel=4, complib='zlib')
            timeMark()

        ##############################################################################################
        # Calculate PCA and tSNE
        ##############################################################################################
        print('Preparing xpca, pcs, tsne of df_expr')
        X_pca, PCs, X_tsne2 = self.project(df_expr)
        print('Recording xpca, pcs, tsne of df_expr')
        pd.DataFrame(data=X_pca, columns=df_expr.columns).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_xpca', mode='a', complevel=4, complib='zlib')
        pd.DataFrame(data=PCs, columns=df_expr.index).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_pcs', mode='a', complevel=4, complib='zlib')
        pd.DataFrame(data=X_tsne2, columns=df_expr.columns).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne_pre_QC' if self.toggleDoQualityControl else 'df_tsne', mode='a', complevel=4, complib='zlib')
        timeMark()

        ##############################################################################################
        # Remove low quality cells
        ##############################################################################################
        if self.toggleDoQualityControl:
            df_expr = df_expr[df_expr.columns.intersection(self.getIndexOfGoodQualityCells(self.saveDir, self.dataName)).sort_values()]
            print('Removed low quality cells. Data size: %s genes, %s cells' % df_expr.shape)
            X_pca, PCs = self.project(df_expr, PCAonly=True)
            pd.DataFrame(data=X_pca, columns=df_expr.columns).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_xpca', mode='a', complevel=4, complib='zlib')
            pd.DataFrame(data=PCs, columns=df_expr.index).to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_pcs', mode='a', complevel=4, complib='zlib')
            df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne_pre_QC', mode='r')
            df_tsne[df_expr.columns].to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='a', complevel=4, complib='zlib')
            timeMark()

        ##############################################################################################
        # Cluster data, append cluster index to expression DataFrame
        ##############################################################################################
        print('Calculating clustering of PCA data')
        df_xpca = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_xpca', mode='r')
        sendDfExpr = False
        if type(self.clusteringFunction) is dict:
            try:
                if clusteringFunction[clusterExpression]:
                    sendDfExpr = True
            except:
                pass
        cellClusterIndexLabel = self.cluster(df_xpca.values, 
                                             df_expr=df_expr if sendDfExpr else None, 
                                             clusteringFunction=self.clusteringFunction)
        df_clusters = pd.DataFrame(data=cellClusterIndexLabel, index=df_expr.columns)
        df_clusters.columns = ['cluster']
        df_clusters.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='a', complevel=4, complib='zlib')
        df_expr = pd.concat([df_clusters, df_expr.T], sort=False, axis=1).reset_index().set_index(['batch', 'cell', 'cluster']).T
        timeMark()

        ##############################################################################################
        # Calculate average of each gene in each cluster
        ##############################################################################################
        df_gene_cluster_centroids = df_expr.groupby(level=['cluster'], sort=True, axis=1).mean()
        df_gene_cluster_centroids.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_gene_cluster_centroids', mode='a', complevel=4, complib='zlib')

        ##############################################################################################
        # Get dictionary to map from markers to cell types. Select genes from the marker list only
        ##############################################################################################
        df_marker_cell_type = pd.read_excel(self.geneListFileName, index_col=0, header=[0,1])
        df_marker_cell_type.columns.names = ['CellTypeGrouped', 'CellType']
        df_marker_cell_type.index = self.gnc.Convert(list(df_marker_cell_type.index), 'alias', 'hugo', returnUnknownString=False)
        df_marker_cell_type = df_marker_cell_type.drop(columns=[col for col in df_marker_cell_type.columns if col[0]=='NA'])
        df_marker_cell_type = df_marker_cell_type.groupby(level='CellTypeGrouped', axis=1).max().loc[df_marker_cell_type.index.intersection(df_expr.index)]
        df_marker_cell_type = df_marker_cell_type[df_marker_cell_type.columns[df_marker_cell_type.sum(axis=0) > self.minimumNumberOfMarkersPerCelltype]]
        df_marker_cell_type = df_marker_cell_type[df_marker_cell_type.sum(axis=1) > 0]
        df_marker_cell_type.index.name = 'Marker'
        df_marker_cell_type.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_marker_cell_type', mode='a', complevel=4, complib='zlib')
        df_marker_cell_type.to_excel(os.path.join(self.saveDir, os.path.basename(self.geneListFileName)))
        df_markers_expr = df_expr.loc[df_expr.index.intersection(df_marker_cell_type.index)]
        print('Selected genes from the marker list. DataFrame shape:', df_markers_expr.shape)
        timeMark()
            
        ##############################################################################################
        # Vote on cell type.  Use the voting scheme described in the paper if
        # no other scheme is given. Update marker expression with "True" labels and record it
        ##############################################################################################
        votingResults = self.vote(df_markers_expr, df_marker_cell_type.T, self.dcsVotingScheme)
        print(votingResults)
        df_markers_expr = df_markers_expr.T.reset_index().T
        df_markers_expr.loc['label'] = np.array([votingResults[int(i)] for i in df_markers_expr.loc['cluster']])
        df_markers_expr = df_markers_expr.T.set_index(['batch', 'cell', 'cluster', 'label']).T.apply(pd.to_numeric)
        df_markers_expr.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='a', complevel=4, complib='zlib')
        timeMark()

        ##############################################################################################
        # Extract new marker genes from the annotation results
        ##############################################################################################
        df_gene_cluster_centroids = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_gene_cluster_centroids', mode='r')
        df_gene_cluster_centroids_merged = pd.DataFrame()
        df_votingResults = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_voting.xlsx'), sheet_name='z-scores', index_col='cluster')[['# cells in cluster', 'Predicted cell type']]
        groups = pd.concat([df_votingResults['# cells in cluster'], df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0]], sort=False, axis=1).reset_index().set_index(0).groupby(by=0)
        for group in groups:
            se_cluster_size = group[1].reset_index().set_index('cluster')['# cells in cluster']
            se_type = (df_gene_cluster_centroids[se_cluster_size.index] * se_cluster_size.values).sum(axis=1) / se_cluster_size.values.sum()
            se_type.name = group[0]
            df_gene_cluster_centroids_merged = pd.concat([df_gene_cluster_centroids_merged, se_type], sort=False, axis=1)
        df_gene_cluster_centroids_merged = df_gene_cluster_centroids_merged.apply(self.zScoreOfSeries, axis=1)
        df_gene_cluster_centroids_merged.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_new_marker_genes', mode='a', complevel=4, complib='zlib')
        self.getNewMarkerGenes()

        ##############################################################################################
        # Create a colormap for cell types (based on a marker-celltype file)
        ##############################################################################################
        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'w') as temp_file:
            df_marker_cell_type = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_marker_cell_type', mode='r')
            cellTypes = np.sort(df_marker_cell_type.columns.values.tolist())
            for cellTypeIndex in range(len(cellTypes)):
                temp_file.write(cellTypes[cellTypeIndex] + '\t' + str(cm.jet(cellTypeIndex / len(cellTypes))) + '\n')
            temp_file.write('Unknown' + '\t' + str(cm.jet(1.0)) + '\n')

        ##############################################################################################
        # Make all plots to conclude analysis
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

            DCS.process(df_expr)

            DCS.visualize()
        '''

        ##############################################################################################
        # Plot null distributions
        ##############################################################################################
        if self.toggleMakeHistogramNullDistributionPlot:
            print('Making null distributions plot')
            self.makeHistogramNullDistributionPlot()
            timeMark()

        ##############################################################################################
        # Plot voting results matrix
        ##############################################################################################
        if self.toggleMakeVotingResultsMatrixPlot:
            print('Making voting results matrix plot')
            self.makeVotingResultsMatrixPlot()
            timeMark()
        
        ##############################################################################################
        # Plot mean marker expression
        ##############################################################################################
        if self.toggleMakeMarkerExpressionPlot:
            print('Making marker expression plot')
            self.makeMarkerExpressionPlot()
            timeMark()

        ##############################################################################################
        # Make tSNE plots
        ##############################################################################################
        if self.toggleMakeTSNEplotQC and self.toggleDoQualityControl:
            self.getTSNEplotsQC()
            timeMark()    

        if self.toggleMakeTSNEplotClusters:
            self.getTSNEplotByClusters()
            timeMark()

        if self.toggleMakeTSNEplotBatches:
            self.getTSNEplotByBatches()
            timeMark()

        if self.toggleMakeTSNEplotAnnotatedClusters:
            self.getTSNEplotAnnotated()
            timeMark()

        if self.toggleAnomalyScoresTSNEplot:
            self.getAnomalyScoresPlot()
            timeMark()
           
        ##############################################################################################
        # Make stacked barplot of cell type fractions
        ##############################################################################################
        if self.toggleMakeStackedBarplot:        
            self.makeStackedBarplot(self.subclusteringName)
            timeMark()
        
        ##############################################################################################
        # Make tSNE plots showing relative expression of different markers (one
        # for each marker)
        ##############################################################################################
        if self.toggleGetMarkerSubplots:
            self.getMarkerSubplots()
            timeMark()

        return None


    # User functions of class ###########################################################################################################################
    def getTSNEplotAnnotated(self):

        '''Produce t-SNE plot colored by cell types
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process(df_expr)

            DCS.getTSNEplotAnnotated()
        '''

        print('Making tSNE plot by clusters "True" labels')

        df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r')

        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
        df_tsne = df_tsne[df_markers_expr.groupby(level=['batch', 'cell'], sort=False, axis=1).count().columns]

        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}

        self.makeTSNEplot(df_tsne.values, df_markers_expr.columns.get_level_values('label'), 'by_clusters_annotated',
                            colormap=colormap, legend=False)

        return None
    
    def getTSNEplotByBatches(self):

        '''Produce t-SNE plot colored by batches
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process(df_expr)

            DCS.getTSNEplotByBatches()
        '''

        print('Making tSNE plot by batches')

        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')

        self.makeTSNEplot(df_tsne.values, df_tsne.columns.get_level_values('batch').values, 'by_patients')

        return None
    
    def getTSNEplotByClusters(self):

        '''Produce t-SNE plot colored by clusters
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process(df_expr)

            DCS.getTSNEplotByClusters()
        '''

        print('Making tSNE plot by clusters')

        df_clusters = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='r')
        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')

        self.makeTSNEplot(df_tsne.values, np.array(['Cluster #%s' % (label[0]) for label in df_clusters.values]), 'by_clusters')

        return None

    def getTSNEplotsQC(self):

        '''Produce Quality Control t-SNE plots
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process(df_expr)

            DCS.getTSNEplotsQC()
        '''

        print('Making tSNE plots of QC')

        df_QC = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='QC', mode='r')
        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne_pre_QC', mode='r')

        goodQUalityCells = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r').columns

        self.makeTSNEplot(df_tsne.values, df_QC['number_of_genes'].values, 'by_number_of_genes', legend=False, labels=False, colorbar=True)
        self.makeTSNEplot(df_tsne.values, df_QC['count_depth'].values, 'by_count_depth', legend=False, labels=False, colorbar=True)
        self.makeTSNEplot(df_tsne.values, df_QC['fraction_of_mitochondrialGenes'].values, 'by_fraction_of_mitochondrialGenes', legend=False, labels=False, colorbar=True)
        self.makeTSNEplot(df_tsne.values, np.array([cell in goodQUalityCells for cell in df_QC.index]), 'by_is_quality_cell', legend=False, labels=True)

        return None

    def getMarkerSubplots(self):

        '''Produce subplots on each marker and its expression on all clusters
              
        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process(df_expr)

            DCS.getMarkerSubplots()
        '''

        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
        df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r')
        df_tsne = df_tsne[pd.MultiIndex.from_arrays([df_markers_expr.columns.get_level_values('batch'), df_markers_expr.columns.get_level_values('cell')])]
        hugo_cd_dict = dict(zip(df_markers_expr.index.values.tolist(), self.gnc.Convert(list(df_markers_expr.index), 'hugo', 'alias', returnUnknownString=False)))
        self.internalMarkerSubplots(df_markers_expr, df_tsne.values, hugo_cd_dict)

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

        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')

        if type(cells) is str:
            if cells=='All':
                cells = df_tsne.columns
            else:
                print('"cells" value "%s" is unknown'%(cells))
                raise ValueError

        df_sel = self.getExprOfCells(cells)
        timeMark()

        scores = self.getAnomalyScores(df_sel, df_sel)
        print('Loaded expression data')
        timeMark()

        scores = pd.DataFrame(index=cells, data=scores).reindex(df_tsne.columns).values.T[0]

        self.makeTSNEplot(df_tsne.values, scores, suffix='by_anomaly_score', legend=False, labels=False, colorbar=True)

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

        '''Produce individual gene expression plot on a tSNE layout

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

        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
        
        hugo_cd_dict = {gene: self.gnc.Convert([gene], 'hugo', 'alias', returnUnknownString=False)[0]}

        self.internalMarkerSubplots(df_markers_expr, 
                                df_tsne.values, 
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

        try:
            df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r').xs(key=gene, axis=0, level=2).T
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

            DCS.process(df_expr)

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

        df_QC['isGoodQuality'].loc[(df_QC['count_depth'] > cutoff_count_depth) & \
            (df_QC['number_of_genes'] > cutoff_number_of_genes) & \
            (df_QC['fraction_of_mitochondrialGenes'] < cutoff_fraction_of_mitochondrialGenes)] = True

        return df_QC['isGoodQuality'][df_QC['isGoodQuality']].index.sort_values()

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

        def alignSeries(se1, se2, tagForMissing):

            '''Unit test:
                se1, se2 = alignSeries(pd.Index(['A', 'B', 'C', 'D']).to_series(), pd.Index(['B', 'C', 'D', 'E', 'F']).to_series())
            '''

            append = lambda se1, se2: pd.concat([se1, pd.Series(index=se2.index.difference(se1.index), data=[tagForMissing] * len(se2.index.difference(se1.index)))], axis=0, sort=False)

            se1 = append(se1, se2)
            se2 = append(se2, se1)

            se1.name = 'se1'
            se2.name = 'se2'

            return se1, se2.loc[se1.index]

        se1.index.name = 'index'
        se2.index.name = 'index'

        df = pd.concat(alignSeries(se1, se2, tagForMissing), axis=1, sort=True)

        counts = {group[0]:{k:len(v) for k, v in group[1].groupby(by='se1').groups.items()} for group in df.reset_index().drop(columns=['index']).set_index('se2').groupby('se2')}

        df = pd.DataFrame.from_dict(counts).fillna(0.0).astype(int)

        moveTag = lambda df: pd.concat([df.iloc[np.where(df.index != tagForMissing)[0]], df.iloc[np.where(df.index == tagForMissing)[0]]], axis=0, sort=False) if tagForMissing in df.index else df

        return moveTag(moveTag(df.T).T)
    
    def getNewMarkerGenes(self, cluster=None, top=100, zScoreCutoff=0.3, removeUnknown=False):

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

        if not cluster is None:

            clusterGenes = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_gene_cluster_centroids', mode='r')[cluster].sort_values(ascending=False)
            clusterGenes = clusterGenes[clusterGenes >= zScoreCutoff].iloc[:top]

            return {'genes':clusterGenes.index.values.tolist(), 'zscore':clusterGenes.values.tolist()}

        df_gene_cluster_centroids_merged = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_new_marker_genes', mode='r')

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
    def convertColormap(self, colormap):

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

    def zScoreOfSeries(self, se):
    
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
        
    def calculateV(self, args):

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

        where_significant = df_Z > zscore_cutoff
        df_Z[where_significant] = 1.0
        df_Z[~where_significant] = 0.0

        df_V = df_M.dot(df_Z).stack()

        df_V[np.isnan(df_V)] = 0.

        return df_V

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

        # Normalize columns (cell types) so that the absolute number of known
        # markers in a given cell type is irrelevant
        df_marker_cell_type = df_marker_cell_type.apply(lambda q: q / np.sum(q), axis=0)
        # Normalize rows (markers) by the number of cell types expressing that marker (i.e. the "specificity")
        df_marker_cell_type = df_marker_cell_type.apply(lambda q:q / np.sum(q > 0), axis=1)

        # Align df_markers_expr and df_marker_cell_type
        df_markers_expr.sort_index(inplace=True, axis=0)
        df_marker_cell_type.sort_index(inplace=True, axis=1)

        # Generate random cluster index
        randomClusterIndex = np.vstack([np.random.choice(df_markers_expr.columns.get_level_values('cluster'), size=df_markers_expr.shape[1], replace=False) for i in range(self.nSamplesDistribution)])

        # Calculate score (Vkc) for the best clustering
        df_V = self.calculateV((df_marker_cell_type, df_markers_expr, df_markers_expr.columns.get_level_values('cluster'), self.zScoreCutoff)).unstack()
        df_V.index.name = None

        # Generate random cluster configurations and calculate scores (Pkc) of
        # those
        print('Generating null distribution')
        startTime = getStartTime()
        pool = multiprocessing.Pool(processes = self.availableCPUsCount)
        random_df_V = pd.concat(pool.map(self.calculateV, [(df_marker_cell_type, df_markers_expr, randomClusterIndex[i], self.zScoreCutoff) for i in range(self.nSamplesDistribution)]), sort=False, axis=1)
        pool.close()
        pool.join()
        getElapsedTime(startTime)

        # Calculate null distribution histograms data for plots
        df_null_distributions = random_df_V.apply(lambda s: scipy.stats.rv_histogram(np.histogram(s, bins=100, range=(0,1)))._hpdf / 100., axis=1).apply(pd.Series).T
        df_null_distributions.columns.names = ['CellType', 'Cluster']

        # Calculate z-score (Lkc) for the best clustering
        print('Processing voting results')
        df_L = (df_V - pd.Series(data=np.mean(random_df_V.values, axis=1), index=random_df_V.index).unstack()) / \
           pd.Series(data=np.std(random_df_V.values, axis=1), index=random_df_V.index).replace(0., np.inf).unstack()
        #df_L[df_L<0.] = 0.
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
        df_custer_centroids_sig = df_custer_centroids.apply(self.zScoreOfSeries,axis=1) > self.zScoreCutoff
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
        df_L.T.to_excel(writer, 'z-scores', columns=['Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers'] + columns_scores)
        df_V.T.to_excel(writer, 'Voting scores', columns=['Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers'] + columns_scores)
        df_custer_centroids.T.to_excel(writer, 'Cluster centroids')
        df_null_distributions.to_excel(writer, 'Null distributions')
        df_marker_cell_type.to_excel(writer, 'Marker cell type weight matrix')
        writer.save()

        return df_V.loc['Predicted cell type'].to_dict()
