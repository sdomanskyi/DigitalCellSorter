'''The main class of DigitalCellSorter. The class includes tools for:

  1. **Pre-preprocessing** of single cell RNA sequencing data

  2. **Quality control**

  3. **Batch effects correction**

  4. **Cells anomaly score evaluation**

  5. **Dimensionality reduction**

  6. **Clustering**

  7. **Annotation of cell types**

  8. **Vizualization**
  
  9. **Post-processing**

'''

import os
import platform
import copy
import multiprocessing
import warnings

import numpy as np
import pandas as pd

import scipy.stats
import scipy.signal
import scipy.linalg
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

    '''Class of Digital Cell Sorter with methods for processing single cell 
    RNA-seq data. Includes analyses and visualization tools.

    Parameters:
        df_expr: pandas.DataFrame, Defauld None
            Gene expression in a form of a table, where genes are rows, and
            cells/batches are columns

        dataName: str, Default 'dataName'
            Name used in output files

        geneNamesType: str, Default 'alias'
            Input gene name convention

        geneListFileName: str, Default None
            Name of the marker genes file

        mitochondrialGenes: list, Default None
            List of mitochondrial genes to use in quality control

        sigmaOverMeanSigma: float, Default 0.1
            Threshold to consider a gene constant

        nClusters: int, Default 10
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
            Note: the function should have .fit method and same input and output.
            For Network-based clustering pass a dictionary 
            {'k_neighbors':40, metric:'euclidean', 'clusterExpression':True},
            this way the best number of clusters will be determined automatically

        nComponentsPCA: int, Default 200
            Number of pca components

        nSamples_pDCS: int, Default 3000
            Number of random samples in distribution for pDCS annotation method

        nSamples_Hopfield: int, Default 500
            Number of repetitions for Hopfield annotation method

        saveDir: str, Default os.path.join('')
            Directory for output files

        makeMarkerSubplots: boolean, Default False
            Whether to make subplots on markers

        makePlots: boolean, Default True
            Whether to make all major plots

        availableCPUsCount: int, Default min(12, os.cpu_count())
            Number of CPUs used in pDCS method

        zScoreCutoff: float, Default 0.3
            Z-Score cutoff when setting expression of a cluster as significant

        thresholdForUnknown: float, Default 0.3
            Threshold when assigning label "Unknown".
            This option is used only with a combination of 2 or more annotation methods

        thresholdForUnknown_pDCS: float, Default 0.1
            Threshold when assigning label "Unknown" in pDCS method
            
        thresholdForUnknown_ratio: float, Default 0.1
            Threshold when assigning label "Unknown" in ratio method
            
        thresholdForUnknown_Hopfield: float, Default 0.1
            Threshold when assigning label "Unknown" in Hopfield method

        annotationMethod: str, Default 'ratio-pDCS-Hopfield'
            Metod to use for annotation of cell types to clusters. Options are:
                'pDCS': main DCS voting scheme with null testing

                'ratio': simple voting score

                'Hopfield': Hopfield Network classifier

                'pDCS-ratio': 'pDCS' adjusted with 'ratio'

                'pDCS-Hopfield': 'pDCS' adjusted with 'Hopfield'

                'ratio-Hopfield': 'ratio' adjusted with 'Hopfield'

                'pDCS-ratio-Hopfield': 'pDCS' adjusted with 'ratio' and 'Hopfield'

        subclusteringName: str, Default None
            Parameter used in for certain labels on plots

        doQualityControl: boolean, Default True
            Whether to remove low quality cells

        doBatchCorrection: boolean, Default False
            Whether to correct data for batches

        minimumNumberOfMarkersPerCelltype: int, Default 10
            Minimum number of markers per cell type to keep that cell type in annotation options

        nameForUnknown: str, Default 'Unassigned'
            Name to use for clusters where label assignment yielded uncertain results

        nameForLowQC: str, Default 'Failed QC'
            Name to use for cell that do not pass quality control

        layout: str, Default 'TSNE'
            Projection layout used in visualization. Options are:
                'TSNE': t-SNE layout
                L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms. 
                Journal of Machine Learning Research 15(Oct):3221-3245, 2014.

                'PCA': use two largest principal components

                'UMAP': use uniform manifold approximation,
                McInnes, L., Healy, J., UMAP: Uniform Manifold Approximation and Projection for 
                Dimension Reduction, ArXiv e-prints 1802.03426, 2018

                'PHATE': use potential of heat diffusion for affinity-based transition embedding,
                Moon, K.R., van Dijk, D., Wang, Z. et al. Visualizing structure and 
                transitions in high-dimensional biological data. Nat Biotechnol 37, 1482–1492 (2019). 
    
    Usage:
        DCS = DigitalCellSorter.DigitalCellSorter()

        df_data = DCS.Clean(df_data)
    '''

    def __init__(self, df_expr = None, dataName = 'dataName', species = 'Human', geneNamesType = 'alias', geneListFileName = None, 
                 mitochondrialGenes = None, sigmaOverMeanSigma = 0.01, nClusters = 10, nFineClusters = 3, doFineClustering = True, 
                 minSizeForFineClustering = 50, clusteringFunction = AgglomerativeClustering, nComponentsPCA = 200, nSamples_pDCS = 3 * 10 ** 3, 
                nSamples_Hopfield = 200, saveDir = os.path.join(''), makeMarkerSubplots = False, availableCPUsCount = min(6, os.cpu_count()), 
                zScoreCutoff = 0.3, subclusteringName = None, doQualityControl = True, doBatchCorrection = False, makePlots = True, 
                useUnderlyingNetwork = True, minimumNumberOfMarkersPerCelltype = 10, nameForUnknown = 'Unassigned', nameForLowQC = 'Failed QC',
                matplotlibMode = None, countDepthCutoffQC = 0.5, numberOfGenesCutoffQC = 0.5, mitochondrialGenesCutoffQC = 1.5, 
                excludedFromQC = None, countDepthPrecutQC = 500, numberOfGenesPrecutQC = 250, precutQC = False, minSubclusterSize = 25,
                thresholdForUnknown_pDCS = 0., thresholdForUnknown_ratio = 0., thresholdForUnknown_Hopfield = 0., thresholdForUnknown = 0.2, 
                layout = 'TSNE', safePlotting = True, HopfieldTemperature = 0.1, annotationMethod = 'ratio-pDCS-Hopfield',
                useNegativeMarkers = True, removeLowQualityScores = True):

        '''Initialization function that is automatically called when an instance on Digital Cell Sorter is created'''

        self.gnc = GeneNameConverter.GeneNameConverter(dictDir=os.path.join(os.path.dirname(__file__)), species = species)

        self.species = species
        self.dataName = dataName
        self.saveDir = saveDir
        self.safePlotting = safePlotting

        if not df_expr is None:
            self.prepare(df_expr)
        else:
            self._df_expr = None
        
        self.defaultGeneListFileName = 'CIBERSORT_LM22'
        self.defaultGeneListsDir = os.path.join(os.path.dirname(__file__), 'geneLists')
        self.geneListFileName = geneListFileName
        self.mitochondrialGenes = mitochondrialGenes

        self.useNegativeMarkers = useNegativeMarkers
        self.removeLowQualityScores = removeLowQualityScores

        self.geneNamesType = geneNamesType

        self.countDepthCutoffQC = countDepthCutoffQC
        self.numberOfGenesCutoffQC = numberOfGenesCutoffQC
        self.mitochondrialGenesCutoffQC = mitochondrialGenesCutoffQC
        self.countDepthPrecutQC = countDepthPrecutQC
        self.numberOfGenesPrecutQC = numberOfGenesPrecutQC
        self.precutQC = precutQC
        self.excludedFromQC = excludedFromQC

        self.minSubclusterSize = minSubclusterSize
        self.sigmaOverMeanSigma = sigmaOverMeanSigma
        self.nClusters = nClusters
        self.doFineClustering = doFineClustering
        self.nFineClusters = nFineClusters
        self.minSizeForFineClustering = minSizeForFineClustering
        self.nComponentsPCA = nComponentsPCA
        self.zScoreCutoff = zScoreCutoff
        self.minimumNumberOfMarkersPerCelltype = minimumNumberOfMarkersPerCelltype

        self.useUnderlyingNetwork = useUnderlyingNetwork

        self.annotationMethod = annotationMethod

        self.nameForUnknown = nameForUnknown
        self.nameForLowQC = nameForLowQC

        self.layout = layout

        self.thresholdForUnknown_pDCS = thresholdForUnknown_pDCS
        self.thresholdForUnknown_ratio = thresholdForUnknown_ratio
        self.thresholdForUnknown_Hopfield = thresholdForUnknown_Hopfield
        self.thresholdForUnknown = thresholdForUnknown

        self.HopfieldTemperature = HopfieldTemperature

        self.nSamples_pDCS = nSamples_pDCS
        self.nSamples_Hopfield = nSamples_Hopfield

        self.availableCPUsCount = availableCPUsCount

        self.clusteringFunction = clusteringFunction

        self.subclusteringName = subclusteringName

        self.doQualityControl = doQualityControl
        self.toggleRemoveLowQualityCells = doQualityControl
        self.toggleDoBatchCorrection = doBatchCorrection

        self.toggleMakeHistogramNullDistributionPlot = False
        self.toggleMakeVotingResultsMatrixPlot = makePlots
        self.toggleMakeMarkerExpressionPlot = makePlots
        self.toggleMakeProjectionPlotsQC = makePlots
        self.toggleMakeProjectionPlotClusters = makePlots
        self.toggleMakeProjectionPlotAnnotatedClusters = makePlots
        self.toggleMakeProjectionPlotBatches = makePlots
        self.toggleMakeStackedBarplot = makePlots

        self.toggleAnomalyScoresProjectionPlot = False

        self.toggleGetMarkerSubplots = makeMarkerSubplots

        super().__init__(saveDir=saveDir, dataName=dataName, safePlotting=safePlotting, matplotlibMode=matplotlibMode)

        return None
    
    @property
    def saveDir(self):

        return self._saveDir

    @saveDir.setter
    def saveDir(self, value):

        self._saveDir = value

        if self.saveDir != os.path.join('') and not os.path.exists(self.saveDir):
            os.makedirs(self.saveDir)

        self._fileHDFpath = os.path.join(self.saveDir, self.dataName + '_processed.h5')

        return

    @property
    def fileHDFpath(self):

        return self._fileHDFpath
    
    @fileHDFpath.setter
    def fileHDFpath(self, value):

        warnings.warn("Direct change of HDF file path is not allowed. " + 
                      "Modify 'saveDir' in order to change 'fileHDFpath'", UserWarning)

        return

    @property
    def df_expr(self):

        if self._df_expr is None:

            print('Expression data is not assigned')

        return self._df_expr
    
    @df_expr.setter
    def df_expr(self, value):

        if value is None:
            self._df_expr = value
        else:
            warnings.warn("Direct change of 'df_expr' is not allowed. " + 
                          "Use function 'prepare' in order to set 'df_expr' " +
                          "or, alternatively, initialize class instance with 'df_expr'. " +
                          "'df_expr' can be removed by assigning None to it", UserWarning)

        return

    @property
    def geneListFileName(self):

        return self._geneListFileName

    @geneListFileName.setter
    def geneListFileName(self, value):

        if value is None:
            self._geneListFileName = os.path.join(self.defaultGeneListsDir, self.defaultGeneListFileName + '.xlsx')
        else:
            if os.path.isfile(value):
                self._geneListFileName = value
            else:
                if os.path.isfile(os.path.join(self.defaultGeneListsDir, value + '.xlsx')):
                    self._geneListFileName = os.path.join(self.defaultGeneListsDir, value + '.xlsx')
                else:
                    print('Gene list file not found. Uning default list: %s'%(self.defaultGeneListFileName), flush=True)

                    self._geneListFileName = os.path.join(self.defaultGeneListsDir, self.defaultGeneListFileName + '.xlsx')

        return


    # Primary functions of class #######################################################
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
                batch = np.array(['batch0'] * len(obj.index.get_level_values('cell')))
            else:
                batch = obj.index.get_level_values('batch')

            obj.index = pd.MultiIndex.from_arrays([batch, obj.index.get_level_values('cell'), obj.index.get_level_values('gene')], names=['batch', 'cell', 'gene'])

            obj = obj.unstack(level='gene').T

            self._df_expr = obj

            print('Done', flush=True)

            return None

        elif type(obj) is pd.DataFrame:
            print('Received data in a form of pandas.DataFrame', flush=True)
            print('Validating pandas.DataFrame', flush=True)

            try:
                obj.index.name = 'gene'
            except Exception as exception:
                print(exception)
                print('DataFrame index format is not understood. Returning None', flush=True)
                return None

            if len(obj.columns.names) > 1:
                if not 'cell' in obj.columns.names:
                    print('Columns level "cell" not found. Returning None', flush=True)
                    return None
            else:
                obj.columns.names = ['cell']

            if not 'batch' in obj.columns.names:
                print('Columns level "batch" not found. Assuming one batch in the data.', flush=True)
                batch = np.array(['batch0'] * len(obj.columns.get_level_values('cell')))
            else:
                batch = obj.columns.get_level_values('batch')

            obj.columns = pd.MultiIndex.from_arrays([batch, obj.columns.get_level_values('cell')], names=['batch', 'cell'])

            self._df_expr = obj

            print('Done', flush=True)

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
                    df_expr['batch'] = np.array(['batch0'] * len(df_expr))

                self._df_expr = df_expr.set_index(['batch', 'cell', 'gene'])['expr'].unstack(level='gene').T

                print('Done', flush=True)

                return None
            else:
                print('Received data in a form of matrix. Reading data', flush=True)
                df_expr = pd.read_csv(obj, header=None, index_col=0)

                print('Converting it to pandas.DataFrame', flush=True)

                if not 'cell' in df_expr.index:
                    print('The row with "cell" identifiers is not found. Returning None', flush=True)
                    return None

                if not 'batch' in df_expr.index:
                    print('The row with "batch" identifiers is not found. Assuming one batch in the data', flush=True)
                    df_expr.loc['batch'] = np.array(['batch0'] * df_expr.shape[1])

                df_expr = df_expr.T.set_index(['batch', 'cell']).T

                df_expr.index.name = 'gene'

                self._df_expr = df_expr

                print('Done', flush=True)

                return None
        else:
            print('Unknown input data format. Returning None', flush=True)

        return None

    def convert(self, nameFrom = None, nameTo = None, **kwargs):

        '''Convert index to hugo names, if any names in the index are
        duplicated, remove duplicates

        Parameters:
            nameFrom: str, Default 'alias'
                Gene name type to convert from

            nameTo: str, Default 'hugo'
                Gene name type to convert to

            Any parameters that function 'mergeIndexDuplicates' can accept

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.convertIndex()
        '''

        if nameFrom is None:
            nameFrom = self.geneNamesType

        if nameTo is None:
            nameTo = 'hugo'

        if nameTo != nameFrom:
            self._df_expr.index = self.gnc.Convert(list(self._df_expr.index), nameFrom, nameTo, returnUnknownString=False)

            self._df_expr = self.mergeIndexDuplicates(self.df_expr, **kwargs)

        return None

    def clean(self):

        '''Clean pandas.DataFrame: validate index, remove index duplicates,
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
        self._df_expr.fillna(0.0, inplace=True)
        print('Replaced missing values with zeros. Data size: %s genes, %s cells' % self._df_expr.shape, flush=True)

        # Check is any names in the index are duplicated, remove duplicates
        self._df_expr = self.mergeIndexDuplicates(self._df_expr) 

        # Keep only cells with at least one expressed gene.
        self._df_expr = self._df_expr.T[self._df_expr.sum(axis=0) > 0].T
        print('Removed all-zero cells. Data size: %s genes, %s cells' % self._df_expr.shape, flush=True)

        # Keep only genes expressed in at least one cell.
        self._df_expr = self._df_expr[self._df_expr.sum(axis=1) > 0]
        print('Removed all-zero genes. Data size: %s genes, %s cells' % self._df_expr.shape, flush=True)

        return None

    def normalize(self, median = None):

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
        median = np.median(np.sum(self._df_expr,axis=0)).astype(float) if median is None else median
        print('Rescaling all cells by "sum of values = %s".' % (median), flush=True)
        self._df_expr = self._df_expr.apply(lambda q: q * median / np.sum(q),axis=0)

        print('Log-transforming data.', flush=True)
        # Replace zeros with minimum value.
        MIN = np.min(self._df_expr.values[self._df_expr.values > 0.])
        if MIN <= 0.:
            raise ValueError
        self._df_expr = self._df_expr.replace(0., MIN)

        # Take log2 of expression.
        self._df_expr = np.log2(self._df_expr)
        self._df_expr -= np.min(self._df_expr.values)

        # Keep only genes expressed in at least one cell.
        self._df_expr = self._df_expr[self._df_expr.sum(axis=1) > 0]
        print('Removed all-zero genes. Data size: %s genes, %s cells' % self._df_expr.shape, flush=True)  

        if not self.sigmaOverMeanSigma is None:
            # Keep only those genes with large enough standard deviation.
            self._df_expr = self._df_expr[np.std(self._df_expr, axis=1) / np.mean(np.std(self._df_expr.values)) > self.sigmaOverMeanSigma]
            print('Removed constant genes. Data size: %s genes, %s cells' % self._df_expr.shape, flush=True)

        # Sort rows by gene name
        self._df_expr = self._df_expr.sort_index()

        return None
    
    def project(self, PCAonly = False, do_fast_tsne = True):

        '''Project pandas.DataFrame to lower dimensions

        Parameters:
            PCAonly: boolean, Default False
                Perform Principal component analysis only

            do_fast_tsne: boolean, Default True
                Do FI-tSNE instead of "exact" tSNE
                This option is ignored if layout is not 'TSNE'

        Returns:
            tuple
                Processed data
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            xPCA, PCs, tSNE = DCS.project()
        '''

        print('Preparing xpca, PCs, 2D projection of df_expr', flush=True)

        print('Performing PC projection from %s to %s features...' % (self._df_expr.shape[0], self.nComponentsPCA), flush=True)
        _PCA = PCA(n_components=self.nComponentsPCA)

        idx = np.argsort(np.var(self._df_expr.values.T, axis=0) / np.mean(self._df_expr.values.T, axis=0))[-2000:]
        X_pca = _PCA.fit_transform(self._df_expr.values.T[:, idx]).T

        print('Explained variance:', np.round(np.sum(_PCA.explained_variance_ratio_) * 100., 2), "%", flush=True)
        
        print('Recording xpca, PCs, 2D projection of df_expr', flush=True)
        columns=self._df_expr.columns
        columns = pd.MultiIndex.from_arrays([columns.get_level_values('batch'), columns.get_level_values('cell')], names=['batch', 'cell'])
        pd.DataFrame(data=X_pca, columns=columns).to_hdf(self.fileHDFpath, key='df_xpca', mode='a', complevel=4, complib='zlib')
        pd.DataFrame(data=_PCA.components_).to_hdf(self.fileHDFpath, key='df_pcs', mode='a', complevel=4, complib='zlib')

        if not PCAonly:
            if self.layout == 'TSNE':
                if X_pca.T.shape[0] < 2000:
                    do_fast_tsne = False

                print('Performing tSNE projection from %s to %s features...' % (self.nComponentsPCA,2), flush=True)
                if do_fast_tsne:
                    if platform.system() == "Windows":
                        from .FastTSNE import fast_tsne
                        print('Windows detected. Using FastTSNE submodule wrapper', flush=True)
                        X_projection2 = fast_tsne(X_pca.T, perplexity = 30, seed=42).T
                    else:
                        try:
                            import fitsne
                        except Exception as exception:
                            print(exception)
                            print('FI-tSNE package could not be imported. Install the package following \
                            instructions of DigitalCellSorter installation on GitHub. Defaulting to PCA layout')
                            self.layout = 'PCA'

                        if self.layout == 'TSNE':
                            X_projection2 = fitsne.FItSNE(np.array(X_pca.T, order='C')).T
                else:
                    X_projection2 = TSNE(n_components=2).fit_transform(X_pca.T).T

            elif self.layout == 'UMAP':
                print('Performing UMAP projection from %s to %s features...' % (self.nComponentsPCA, 2), flush=True)
                try:
                    import umap
                except Exception as exception:
                    print(exception)
                    print('UMAP package could not be imported. Install the package following \
                    instructions of DigitalCellSorter installation on GitHub. Defaulting to PCA layout')
                    self.layout = 'PCA'

                if self.layout == 'UMAP':
                    X_projection2 = umap.UMAP(random_state=42).fit_transform(X_pca.T).T

            elif self.layout == 'PHATE':
                print('Performing PHATE projection from %s to %s features...' % (self.nComponentsPCA, 2), flush=True)
                try:
                    import phate
                except Exception as exception:
                    print(exception)
                    print('PHATE package could not be imported. Install the package following \
                    instructions of DigitalCellSorter installation on GitHub. Defaulting to PCA layout')
                    self.layout = 'PCA'

                if self.layout == 'PHATE':
                    X_projection2 = phate.PHATE().fit_transform(X_pca.T).T

            if self.layout == 'PCA':
                print('Using PC1 and PC2 for layout', flush=True)
                X_projection2 = X_pca[[0,1],:]

            print('Recording 2D projection of df_expr', flush=True)

            pd.DataFrame(data=X_projection2, columns=columns).to_hdf(self.fileHDFpath, key='df_projection_pre_QC', mode='a', complevel=4, complib='zlib')
            pd.DataFrame(data=X_projection2, columns=columns).to_hdf(self.fileHDFpath, key='df_projection', mode='a', complevel=4, complib='zlib')

            return X_pca, _PCA.components_, X_projection2

        return X_pca, _PCA.components_
    
    def cluster(self):

        '''Cluster PCA-reduced data into a desired number of clusters

        Parameters:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.cluster()
        '''

        print('Calculating clustering of PCA data', flush=True)
        df_xpca = pd.read_hdf(self.fileHDFpath, key='df_xpca')

        if type(self.clusteringFunction) is dict:
            try:
                import pynndescent
                import networkx as nx
                import community
            except Exception as exception:
                print(exception)
                print('Install packages pynndescent, networkx and python-louvain following \
                instructions of DigitalCellSorter installation on GitHub. Defaulting to agglomerative clustering method')
                self.clusteringFunction = AgglomerativeClustering

        if type(self.clusteringFunction) is dict:

            import pynndescent
            import networkx as nx
            import community
            
            try:
                k_neighbors = self.clusteringFunction[k_neighbors]
            except Exception as exception:
                print(exception)
                k_neighbors = 40

            try:
                metric = self.clusteringFunction[metric]
            except Exception as exception:
                print(exception)
                metric = 'euclidean'
            
            try:
                clusterExpression = self.clusteringFunction[clusterExpression]
            except Exception as exception:
                print(exception)
                clusterExpression = False

            data = self._df_expr.values.T if clusterExpression else df_xpca.values.T

            print('Searching for %s nearest neighbors' % (k_neighbors), flush=True)
            knn = pynndescent.NNDescent(data, metric=metric, n_neighbors=k_neighbors).query(data, k=k_neighbors)

            print('k(=%s) nearest neighbors found. Constructing a NetworkX graph' % (k_neighbors), flush=True)
            A = np.zeros((len(knn[0]),len(knn[0])))
            for i in range(len(knn[0])):
                A[i, knn[0][i]] = knn[1][i]

            G = nx.from_numpy_array(A)

            print('Clustering the graph', flush=True)
            cellClusterIndex = pd.Series(community.best_partition(G)).sort_index().values
        else:
            data = df_xpca.values

            cellClusterIndex = self.clusteringFunction(n_clusters=self.nClusters).fit(data.T).labels_.astype(float).astype(str)
            print(np.unique(cellClusterIndex, return_counts=True))

            if self.doFineClustering:
                fineCellClusterIndex = cellClusterIndex.copy()

                clusters = np.unique(cellClusterIndex)

                for cluster in clusters:
                    cellsOfCluster = np.where(cellClusterIndex == cluster)[0]
                    subData = data[:, cellsOfCluster]
                    subData = subData[subData.sum(axis=1) > 0.]
            
                    if len(cellsOfCluster) >= self.minSizeForFineClustering:
                        try:
                            model = SpectralCoclustering(n_clusters=self.nFineClusters, random_state=0)
                            model.fit(subData.T)
                            tempCellClusterIndex = np.zeros((subData.T.shape[0],))
                            tempCellClusterIndex[:] = np.nan
                            for i, subCluster in enumerate(model.biclusters_[0]):
                                tempCellClusterIndex[np.where(subCluster)[0]] = i
                        except Exception as exception:
                            print(exception)
                            print('Exception', subData.shape, cellsOfCluster.shape)
                            tempCellClusterIndex = np.zeros((subData.T.shape[0],))
                    else:
                        print('Small cluster', len(cellsOfCluster))
                        tempCellClusterIndex = np.zeros((subData.T.shape[0],))

                    subSizes = np.unique(tempCellClusterIndex, return_counts=True)[1]

                    if (subSizes < self.minSubclusterSize).any():
                        print('Fine clusters too small (less than %s)' % (self.minSubclusterSize), subSizes)
                        tempCellClusterIndex = np.zeros((subData.T.shape[0],))
            
                    fineCellClusterIndex[cellsOfCluster] = np.array([str(cluster) + '.' + label for label in tempCellClusterIndex.astype(int).astype(str)])

                cellClusterIndex = fineCellClusterIndex

        df_clusters = pd.DataFrame(data=cellClusterIndex, index=self._df_expr.columns)
        df_clusters.columns = ['cluster']

        df_clusters.to_hdf(self.fileHDFpath, key='df_clusters', mode='a', complevel=4, complib='zlib')

        self._df_expr = pd.concat([df_clusters, self._df_expr.T], sort=False, axis=1).reset_index().set_index(['batch', 'cell', 'cluster']).T

        return
    
    def annotate(self, mapNonexpressedCelltypes = True):
        
        '''Produce cluster voting results, annotate cell types, 
        and update marker expression with cell type labels

        Parameters:
            mapNonexpressedCelltypes: boolean, Default True
                If True then cell types coloring will be consistent across all datasets,
                regardless what cell types are annotated in all datasets
                for a given input marker list file.

        Returns:
            dictionary
                Voting results, a dictionary in form of:
                {cluster label: assigned cell type}
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            results = DCS.annotate(df_markers_expr, df_marker_cell_type)
        '''

        self.loadExpressionData()

        if self._df_expr is None:

            return

        self.prepareMarkers(expressedGenes=self._df_expr.index if mapNonexpressedCelltypes else None)

        df_marker_cell_type = pd.read_hdf(self.fileHDFpath, key='df_marker_cell_type')

        df_markers_expr = self._df_expr.loc[self._df_expr.index.intersection(df_marker_cell_type.index)]
        print('Selected genes common with the marker list. Data shape:', df_markers_expr.shape, flush=True)

        methodsToUse = self.annotationMethod.split('-')

        if 'pDCS' in methodsToUse:
            print('\nCalculating results by %s method' % ('pDCS'), flush = True)
            annotationResults_pDCS = list(self.annotateWith_pDCS_Scheme(df_markers_expr.copy(), df_marker_cell_type.T.copy()))

        if 'ratio' in methodsToUse:
            print('\nCalculating results by %s method' % ('ratio'), flush = True)
            annotationResults_ratio = list(self.annotateWith_ratio_Scheme(df_markers_expr.copy(), df_marker_cell_type.T.copy()))

        if 'Hopfield' in methodsToUse:
            print('\nCalculating results by %s method' % ('Hopfield'), flush = True)
            annotationResults_Hopfield = list(self.annotateWith_Hopfield_Scheme(df_markers_expr.copy(), df_marker_cell_type.T.copy()))
            #self.annotateWith_Hopfield_Scheme(df_markers_expr.copy(), df_marker_cell_type.T.copy())
            #return

        if ('pDCS' in methodsToUse) and ('ratio' in methodsToUse) and ('Hopfield' in methodsToUse):
            annotationResults = annotationResults_pDCS
            annotationResults[2] *= annotationResults_ratio[2]
            annotationResults[2] *= annotationResults_Hopfield[2]
            annotationResults[2] = np.power(annotationResults[2], 1./3.)
            annotationResults[2][annotationResults[2] < self.thresholdForUnknown] = 0.

        elif ('pDCS' in methodsToUse) and ('ratio' in methodsToUse) and (not 'Hopfield' in methodsToUse):
            annotationResults = annotationResults_pDCS
            annotationResults[2] *= annotationResults_ratio[2]
            annotationResults[2] = np.power(annotationResults[2], 1./2.)
            annotationResults[2][annotationResults[2] < self.thresholdForUnknown] = 0.

        elif ('pDCS' in methodsToUse) and (not 'ratio' in methodsToUse) and ('Hopfield' in methodsToUse):
            annotationResults = annotationResults_pDCS
            annotationResults[2] *= annotationResults_Hopfield[2]
            annotationResults[2] = np.power(annotationResults[2], 1./2.)
            annotationResults[2][annotationResults[2] < self.thresholdForUnknown] = 0.

        elif (not 'pDCS' in methodsToUse) and ('ratio' in methodsToUse) and ('Hopfield' in methodsToUse):
            annotationResults = annotationResults_ratio
            annotationResults[2] *= annotationResults_Hopfield[2]
            annotationResults[2] = np.power(annotationResults[2], 1./2.)
            annotationResults[2][annotationResults[2] < self.thresholdForUnknown] = 0.

        elif ('pDCS' in methodsToUse) and (not 'ratio' in methodsToUse) and (not 'Hopfield' in methodsToUse):
            annotationResults = annotationResults_pDCS

        elif (not 'pDCS' in methodsToUse) and ('ratio' in methodsToUse) and (not 'Hopfield' in methodsToUse):
            annotationResults = annotationResults_ratio

        elif (not 'pDCS' in methodsToUse) and (not 'ratio' in methodsToUse) and ('Hopfield' in methodsToUse):
            annotationResults = annotationResults_Hopfield

        else:
            print('Annotation method %s not supported' % (self.annotationMethod), flush=True)

            return

        annotationResults = self.recordAnnotationResults(*tuple(annotationResults))

        print(annotationResults, flush=True)
        print('Updating the expression data with annotation results', flush=True)

        df_markers_expr = df_markers_expr.T.reset_index().T

        df_markers_expr.loc['label'] = np.array([annotationResults[i] for i in df_markers_expr.loc['cluster']])

        df_markers_expr = df_markers_expr.T.set_index(['batch', 'cell', 'cluster', 'label']).T.apply(pd.to_numeric)

        df_markers_expr.to_hdf(self.fileHDFpath, key='df_markers_expr', mode='a', complevel=4, complib='zlib')
        
        return annotationResults
    
    def process(self, dataIsNormalized = False, cleanData = True):

        '''Process data before using any annotation of visualization functions

        Parameters:
            dataIsNormalized: boolean, Default False
                Whether DCS.df_expr is normalized or not

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()
        '''

        # Convert index to hugo names, clean data
        if not dataIsNormalized:
            self.convert()
        
        if cleanData:
            self.clean()
        
        # Calculate QC measures
        if self.doQualityControl:
            self.calculateQCmeasures()

        # Normalize and then correct for batch effects
        if not dataIsNormalized:
            self.normalize()

        if self.toggleDoBatchCorrection:
            self.batchEffectCorrection()

        # Calculate PCA and 2D projection
        self.project()

        # Remove low quality cells
        if self.doQualityControl:
            self.qualityControl()

        # Cluster data, append cluster index to expression DataFrame
        self.cluster()

        # Compress and record DataFrame to dense matrix
        self.recordExpressionData()
        
        return

    def visualize(self):

        '''Aggregate of visualization tools of this class.

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
            print('Making voting results matrix plot', flush=True)
            self.makeAnnotationResultsMatrixPlot()

        ##############################################################################################
        # Plot null distributions
        ##############################################################################################
        if self.toggleMakeHistogramNullDistributionPlot:
            print('Making null distributions plot', flush=True)
            self.makeHistogramNullDistributionPlot()
        
        ##############################################################################################
        # Plot mean marker expression
        ##############################################################################################
        if self.toggleMakeMarkerExpressionPlot:
            print('Making marker expression plot', flush=True)
            self.makeMarkerExpressionPlot()

        ##############################################################################################
        # Make 2D projection plots
        ##############################################################################################
        if self.toggleMakeProjectionPlotsQC and self.doQualityControl:
            self.makeProjectionPlotsQualityControl()

        if self.toggleMakeProjectionPlotClusters:
            self.makeProjectionPlotByClusters()

        if self.toggleMakeProjectionPlotBatches:
            self.makeProjectionPlotByBatches()

        if self.toggleMakeProjectionPlotAnnotatedClusters:
            self.makeProjectionPlotAnnotated()

        if self.toggleAnomalyScoresProjectionPlot:
            self.makeAnomalyScoresPlot()
           
        ##############################################################################################
        # Make stacked barplot of cell type fractions
        ##############################################################################################
        if self.toggleMakeStackedBarplot:        
            self.makeStackedBarplot(self.subclusteringName)
        
        ##############################################################################################
        # Make projection plots showing relative expression of different markers (one for each
        # marker)
        ##############################################################################################
        if self.toggleGetMarkerSubplots:
            self.makeMarkerSubplots()

        return None


    # Plotting user functions of class #################################################
    def makeProjectionPlotAnnotated(self, **kwargs):

        '''Produce projection plot colored by cell types
              
        Parameters:
            Any parameters that function 'makeProjectionPlot' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.makeProjectionPlotAnnotated()
        '''

        print('Making projection plot by clusters with annotated labels')

        df_markers_expr = pd.read_hdf(self.fileHDFpath, key='df_markers_expr')

        df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')
        df_projection = df_projection[df_markers_expr.groupby(level=['batch', 'cell'], sort=False, axis=1).count().columns]

        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}

        kwargs.setdefault('suffix', 'by_clusters_annotated')
        kwargs.setdefault('colormap', colormap)
        kwargs.setdefault('legend', False)

        fig = self.makeProjectionPlot(df_projection.values, df_markers_expr.columns.get_level_values('label'), **kwargs)

        return fig
    
    def makeProjectionPlotByBatches(self, **kwargs):

        '''Produce projection plot colored by batches
              
        Parameters:
            Any parameters that function 'makeProjectionPlot' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.makeProjectionPlotByBatches()
        '''

        print('Making projection plot by batches')

        df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')

        batches = df_projection.columns.get_level_values('batch').values

        #wh1 = np.where(batches == np.unique(batches)[0])
        #wh2 = np.where(batches == np.unique(batches)[1])
        #wh3 = np.where(batches == np.unique(batches)[2])
        #batches[wh1] = 'BM1'
        #batches[wh2] = 'BM2'
        #batches[wh3] = 'BM3'

        kwargs.setdefault('legend', True)
        kwargs.setdefault('labels', False)
        kwargs.setdefault('suffix', 'by_patients')

        fig = self.makeProjectionPlot(df_projection.values, batches, **kwargs)

        return fig
    
    def makeProjectionPlotByClusters(self, **kwargs):

        '''Produce projection plot colored by clusters
              
        Parameters:
            Any parameters that function 'makeProjectionPlot' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.makeProjectionPlotByClusters()
        '''

        print('Making projection plot by clusters')

        df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')
        df_clusters = pd.read_hdf(self.fileHDFpath, key='df_clusters').reindex(df_projection.columns)

        kwargs.setdefault('suffix', 'by_clusters')

        labels = np.array(['Cluster #%s' % (label[0]) if label==label else label for label in df_clusters.values])

        fig = self.makeProjectionPlot(df_projection.values, labels, **kwargs)

        return fig

    def makeProjectionPlotsQualityControl(self, **kwargs):

        '''Produce Quality Control projection plots
              
        Parameters:
            Any parameters that function 'makeProjectionPlot' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.makeProjectionPlotsQualityControl()
        '''

        def reduce(v, size = 100):

            bins =  np.linspace(np.min(v), np.max(v), num=size)

            return bins[np.digitize(v, bins) - 1]

        print('Making projection plots of QC', flush=True)

        df_QC = pd.read_hdf(self.fileHDFpath, key='QC')
        if self.toggleRemoveLowQualityCells:
            df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection_pre_QC')
        else:
            df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')

        goodQUalityCells = pd.read_hdf(self.fileHDFpath, key='df_projection').columns

        kwargs.setdefault('legend', False)
        kwargs.setdefault('labels', False)
        kwargs.setdefault('colorbar', True)

        kwargs.update(suffix='by_number_of_genes')
        self.makeProjectionPlot(df_projection.values, reduce(df_QC['number_of_genes'].values), **kwargs)

        kwargs.update(suffix='by_count_depth')
        self.makeProjectionPlot(df_projection.values, reduce(df_QC['count_depth'].values), **kwargs)
        
        kwargs.update(suffix='by_fraction_of_mitochondrialGenes')
        self.makeProjectionPlot(df_projection.values, reduce(df_QC['fraction_of_mitochondrialGenes'].values), **kwargs)

        kwargs.update(suffix='by_is_quality_cell')
        kwargs.update(colorbar=False)

        self.makeProjectionPlot(df_projection.values, df_projection.columns.isin(goodQUalityCells.get_level_values('cell'), level='cell'), **kwargs)

        return

    def makeMarkerSubplots(self, **kwargs):

        '''Produce subplots on each marker and its expression on all clusters
              
        Parameters:
            Any parameters that function 'internalMakeMarkerSubplots' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.makeMarkerSubplots()
        '''

        df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')
        df_markers_expr = pd.read_hdf(self.fileHDFpath, key='df_markers_expr')
        df_projection = df_projection[pd.MultiIndex.from_arrays([df_markers_expr.columns.get_level_values('batch'), df_markers_expr.columns.get_level_values('cell')])]

        hugo_cd_dict = dict(zip(df_markers_expr.index.values.tolist(), self.gnc.Convert(list(df_markers_expr.index), 'hugo', 'alias', returnUnknownString=False)))
        
        kwargs.setdefault('analyzeBy', 'label')
        kwargs.setdefault('NoFrameOnFigures', True)
        kwargs.setdefault('HideClusterLabels', False)
        kwargs.setdefault('outlineClusters', True)

        self.internalMakeMarkerSubplots(df_markers_expr, df_projection.values, hugo_cd_dict, **kwargs)

        return None

    def makeAnomalyScoresPlot(self, cells = 'All', suffix = '', noPlot = False, **kwargs):

        '''Make anomaly scores plot

        Parameters:
            cells: pandas.MultiIndex, Default 'All'
                Index of cells of interest

            Any parameters that function 'makeProjectionPlot' can accept

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            cells = DCS.getCells(celltype='T cell')

            DCS.makeAnomalyScoresPlot(cells)
        '''

        if self._df_expr is None:
            self.loadExpressionData()

        if type(cells) is str:
            if cells == 'All':
                cells = pd.MultiIndex.from_arrays([self._df_expr.columns.get_level_values('batch'), 
                                                   self._df_expr.columns.get_level_values('cell') ], names=['batch', 'cell']) 
            else:
                print('"cells" value "%s" is unknown' % (cells))
                raise ValueError

        df_sel = self.getExprOfCells(cells)

        df_sel.columns = pd.MultiIndex.from_arrays([df_sel.columns.get_level_values('batch'), 
                                                    df_sel.columns.get_level_values('cell') ], names=['batch', 'cell'])

        scores = self.getAnomalyScores(df_sel, df_sel)
        scores = pd.DataFrame(index=df_sel.columns, data=scores)
        scores.columns = ['score']
        scores.sort_values(by='score', ascending=False).to_excel(os.path.join(self.saveDir, 'anomalyScores%s.xlsx' % (suffix)), merge_cells=False)
        
        if noPlot:

            return

        df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')
        scores = scores.reindex(df_projection.columns).values.T[0]

        kwargs.setdefault('suffix', 'by_anomaly_score')
        kwargs.setdefault('legend', False)
        kwargs.setdefault('labels', False)
        kwargs.setdefault('colorbar', True)

        fig = self.makeProjectionPlot(df_projection.values, scores, **kwargs)

        return fig
    
    def makeIndividualGeneTtestPlot(self, gene, analyzeBy = 'label', **kwargs):

        '''Produce individual gene t-test plot of the two-tailed p-value.

        Parameters:
            gene: str
                Name of gene of interest

            analyzeBy: str, Default 'label'
                What level of lablels to include.
                Other possible options are 'label' and 'celltype'

            Any parameters that function 'makeTtestPlot' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeIndividualGeneTtestPlot('SDC1')
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
                ttestStatistic.loc[groupA, groupB], ttestpValue.loc[groupA, groupB] = scipy.stats.ttest_ind(A[np.where(A != 0)], B[np.where(B != 0)])

        alt = self.gnc.Convert([gene], 'hugo', 'alias', returnUnknownString=False)[0]
        alt = [alt] if type(alt) is str else alt

        kwargs.setdefault('label', '%s\n(%s)' % (gene, ('\n').join(list(alt))))

        fig = self.makeTtestPlot(ttestStatistic, ttestpValue, **kwargs)

        return fig    

    def makeIndividualGeneExpressionPlot(self, genes, **kwargs):

        '''Produce individual gene expression plot on a 2D layout

        Parameters:
            gene: str, or list-like
                Name of gene of interest. E.g. 'CD4, CD33', 'PECAM1', ['CD4', 'CD33']

            hideClusterLabels: boolean, Default False
                Whether to hide the clusters labels

            outlineClusters: boolean, Default True
                Whether to outline the clusters with circles

            Any parameters that function 'internalMakeMarkerSubplots' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeIndividualGeneExpressionPlot('CD4')
        '''

        kwargs.setdefault('analyzeBy', 'cluster')
        kwargs.setdefault('NoFrameOnFigures', True)
        kwargs.setdefault('HideClusterLabels', False)
        kwargs.setdefault('outlineClusters', True)

        df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection')

        if type(genes) is str:
            genes = genes.replace(' ', '').split(',')

        for gene in genes:
            df_markers_expr = self.getExprOfGene(gene, analyzeBy=kwargs['analyzeBy'])

            if not df_markers_expr is None:
                hugo_cd_dict = {gene: self.gnc.Convert([gene], 'hugo', 'alias', returnUnknownString=False)[0]}

                cells = pd.MultiIndex.from_arrays([df_markers_expr.columns.get_level_values('batch'),
                                                   df_markers_expr.columns.get_level_values('cell')], 
                                                  names=['batch', 'cell'])

                self.internalMakeMarkerSubplots(df_markers_expr, 
                                                df_projection.reindex(cells, axis=1).values, 
                                                hugo_cd_dict, **kwargs)

        return None

    def makeHopfieldLandscapePlot(self, meshSamplingRate = 1000, plot3D = True, reuseData = False, **kwargs):

        '''Make and plot Hopfield landscape

        Parameters:
            meshSamplingRate: int, Default 1000
                Defines quality of sampling around attractor states
    
            plot3D: boolean, Default False
                Whether to plot 2D or 3D figure

            reuseData: boolean, Default False
                Whether to attempt using precalculated data.

            Any parameters that function 'HopfieldLandscapePlot' or 'HopfieldLandscapePlot3D' can accept

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.makeHopfieldLandscapePlot()
        '''

        dataPath = os.path.join(self.saveDir, 'Hopfiled_mode_1')

        if reuseData and (not os.path.exists(dataPath)):
            reuseData = False

        if not reuseData:
            df_attrs = self.readMarkerFile()

            self.propagateHopfield(xi=df_attrs.values, typesNames=df_attrs.columns.values, path=dataPath, meshSamplingRate=meshSamplingRate, mode=1)

        if plot3D:
            fig = self.HopfieldLandscapePlot3D(trPath=dataPath, **kwargs)
        else:
            fig = self.HopfieldLandscapePlot(trPath=dataPath, colorbarva=0.15, fontsize=12, **kwargs)

        return fig

    
    # Extraction user functions of class ###############################################
    def getAnomalyScores(self, trainingSet, testingSet, printResults = False):

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

            cutoff = DCS.getAnomalyScores(df_expr.iloc[:, 5:], df_expr.iloc[:, :5])
        '''

        print('Calculating anomaly scores of the selected cells', flush=True)

        instanceIsolationForest = IsolationForest(n_jobs=self.availableCPUsCount, max_samples=np.min([trainingSet.shape[1] - 1, 256]))

        print('Training isolation forest')
        instanceIsolationForest.fit(trainingSet.values.T)

        if type(testingSet) is pd.Series:
            testingSet = testingSet.to_frame()

        print('Scoring samples')
        scores = instanceIsolationForest.score_samples(testingSet.values.T)

        scores *= -1.

        results = list(zip(testingSet.columns, scores))

        if printResults:
            for cell in results:
                print('Cell:', cell[0], 'Score:', cell[1])

        return scores

    def getHugoName(self, gene, printAliases = False):
        
        '''Get gene hugo name(s).

        Parameters:
            gene: str
                'hugo' or 'alias' name of a gene

        Returns:
            str
                Hugo name if found, otherwise input name

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.getHugoName('CD138')
        '''

        return self.gnc.Convert(gene, 'alias', 'hugo', returnUnknownString=False)

    def getExprOfGene(self, gene, analyzeBy = 'cluster'):

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

        fileName = self.fileHDFpath

        if not os.path.isfile(fileName):
            print('Processed data file not found')
            return
        
        try:
            df_markers_expr = pd.read_hdf(fileName, key='df_expr').xs(key=gene, axis=0, level=-1).T
        except Exception as exception:
            print(exception)
            print('Gene %s not in index' % (gene))
            return

        df_markers_expr.index = [gene]

        df_markers_expr.columns = pd.MultiIndex.from_arrays([df_markers_expr.columns.get_level_values('batch'), df_markers_expr.columns.get_level_values('cell')])

        labelled = pd.read_hdf(self.fileHDFpath, key='df_markers_expr').columns
        
        columns = pd.MultiIndex.from_arrays([labelled.get_level_values('batch'), labelled.get_level_values('cell')])
        df_markers_expr = df_markers_expr.reindex(columns, axis=1).fillna(0.)

        if analyzeBy == 'celltype':
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
                2-level Index of cells of interest, must include levels 'batch' and 'cell'

        Returns:
            pandas.DataFrame
                With expression of the cells of interest

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.process()

            DCS.getExprOfCells(cells)
        '''

        if self._df_expr is None:
            self.loadExpressionData()

        index = pd.MultiIndex.from_arrays([self._df_expr.columns.get_level_values('batch'), 
                                   self._df_expr.columns.get_level_values('cell') ], names=['batch', 'cell'])

        if cells.reindex(index)[1] is None:
            print('Selected expression of all cells')

            return self.df_expr

        return self._df_expr[self._df_expr.columns[cells.reindex(index)[1] > -1]]

    def getCells(self, celltype = None, clusterIndex = None, clusterName = None):

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

        df = pd.read_hdf(self.fileHDFpath, key='df_markers_expr').columns.to_frame()

        if not clusterIndex is None:

            se = df.set_index(['batch', 'cell'])['cluster']

            condition = se == clusterIndex

            if not condition.any():
                print('No cells found')
                return None

            selectedCells = se[condition].index

            print('Selected %s cells from cluster: %s' % (len(selectedCells), clusterIndex))

            return selectedCells

        if not clusterName is None:

            se = df.set_index(['batch', 'cell'])['cluster']

            condition = se == eval(clusterName.split('#')[1])

            if not condition.any():
                print('No cells found')
                return None

            selectedCells = se[condition].index

            print('Selected %s cells from cluster: %s' % (len(selectedCells), clusterName))

            return selectedCells

        se = df.set_index(['batch', 'cell'])['label']

        if celltype is None and clusterIndex is None:

            print('Available cell types:', np.unique(se.str.split(' #', expand=True)[0].values))
            print('Note: To get a specific cell type call this function with a specified cell type')

            return se

        condition = se.str.find(celltype) != -1

        if not condition.any():
            print('No cells found')
            return None

        selectedCells = se[condition].index

        print('Selected %s cells of type: %s' % (len(selectedCells), celltype))

        return selectedCells

    def getIndexOfGoodQualityCells(self, QCplotsSubDir = 'QC_plots', **kwargs):

        '''Get index of sells that satisfy the QC criteria

        Parameters:
            count_depth_cutoff: float, Default 0.5
                Fraction of median to take as count depth cutoff

            number_of_genes_cutoff: float, Default 0.5
                Fraction of median to take as number of genes cutoff

            mitochondrial_genes_cutoff: float, Default 3.0
                The cutoff is median + standard_deviation * this_parameter

            Any parameters that function 'makeQualityControlHistogramPlot' can accept

        Returns:
            pandas.Index
                Index of cells
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            index = DCS.getIndexOfGoodQualityCells()
        '''

        df_QC = pd.read_hdf(self.fileHDFpath, key='QC')

        if not self.excludedFromQC is None:
            df_QC = df_QC.loc[df_QC.index.difference(self.excludedFromQC)]

        plotsDir = os.path.join(self.saveDir, QCplotsSubDir, '')

        kwargs.setdefault('MakeHistogramPlot', True)

        # Calculate cutoffs
        cutoff_count_depth = self.getQualityControlCutoff(df_QC['count_depth'], self.countDepthCutoffQC, plotPathAndName=os.path.join(plotsDir, '%s_count_depth' % (self.dataName)), precut=self.countDepthPrecutQC, **kwargs)
        cutoff_number_of_genes = self.getQualityControlCutoff(df_QC['number_of_genes'], self.numberOfGenesCutoffQC, plotPathAndName=os.path.join(plotsDir, '%s_number_of_genes' % (self.dataName)), precut=self.numberOfGenesPrecutQC, **kwargs)
        cutoff_fraction_of_mitochondrialGenes = self.getQualityControlCutoff(df_QC['fraction_of_mitochondrialGenes'], self.mitochondrialGenesCutoffQC, plotPathAndName=os.path.join(plotsDir, '%s_fraction_of_mitochondrialGenes' % (self.dataName)), mito=True, **kwargs)

        df_QC['isGoodQuality'] = np.zeros(len(df_QC)).astype(bool)

        mask = (df_QC['count_depth'] > cutoff_count_depth) & \
            (df_QC['number_of_genes'] > cutoff_number_of_genes) & \
            (df_QC['fraction_of_mitochondrialGenes'] < cutoff_fraction_of_mitochondrialGenes)

        index = df_QC['isGoodQuality'].index[mask]

        if not self.excludedFromQC is None:
            index = pd.concat([pd.Series(index=index), pd.Series(index=self.excludedFromQC)]).index

        return index.sort_values()

    def getQualityControlCutoff(self, se, cutoff, precut = 1., mito = False, MakeHistogramPlot = True, **kwargs):

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

            Any parameters that function 'makeQualityControlHistogramPlot' can accept

        Returns:
            float
                Cutoff value
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            cutoff = DCS.getQualityControlCutoff(se)
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

            x1 = spline_data.T[0][np.where(spline_data.T[0] > median)[0][0]]
            y1 = smoothed[np.where(spline_data.T[0] > median)[0][0]]

            wh = np.where(spline_data.T[0] > (median + cutoff * std))[0][0]

            x2 = spline_data.T[0][wh]
            y2 = smoothed[wh]

            cutoff = min(median + 3. * std, x1 - y1 * (x2 - x1) / (y2 - y1))
            cutoff = max(median + 1. * std, cutoff)
        else:
            if self.precutQC:
                cutoff = int(cutoff * np.median(subset[subset > precut]))
                cutoff = min(cutoff, 2 * precut)
            else:
                cutoff = int(cutoff * np.median(subset))

        kwargs.setdefault('mito', mito)

        if MakeHistogramPlot:
            self.makeQualityControlHistogramPlot(subset, cutoff, **kwargs)
            
        return cutoff

    def getCountsDataframe(self, se1, se2, tagForMissing = 'N/A'):

        '''Get a pandas.DataFrame with cross-counts (overlaps) between two pandas.Series

        Parameters:
            se1: pandas.Series
                Series with the first set of items

            se2: pandas.Series 
                Series with the second set of items

            tagForMissing: str, Default 'N/A'
                Label to assign to non-overlapping items

        Returns:
            pandas.DataFrame
                Contains counts
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df = DCS.getCountsDataframe(se1, se2)
        '''

        df = self.alignSeries(se1, se2, tagForMissing)

        counts = {group[0]:{k:len(v) for k, v in group[1].groupby(by='se1').groups.items()} for group in df.reset_index(drop=True).set_index('se2').groupby('se2')}

        df = pd.DataFrame.from_dict(counts).fillna(0.0).astype(int)

        moveTag = lambda df: pd.concat([df.iloc[np.where(df.index != tagForMissing)[0]], df.iloc[np.where(df.index == tagForMissing)[0]]], axis=0, sort=False) if tagForMissing in df.index else df

        return moveTag(moveTag(df.T).T)
    
    def getNewMarkerGenes(self, cluster = None, top = 100, zScoreCutoff = None, removeUnknown = False, **kwargs):

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

            Any parameters that function 'makePlotOfNewMarkers' can accept

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.extractNewMarkerGenes()
        '''

        if zScoreCutoff is None:

            zScoreCutoff = self.zScoreCutoff

        if self._df_expr is None:

            self.loadExpressionData()

        if self._df_expr is None:

            return

        df_gene_cluster_centroids = self._df_expr.groupby(level=['cluster'], sort=True, axis=1).mean()

        df_gene_cluster_centroids_merged = pd.DataFrame()
        df_votingResults = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='z-scores', index_col='cluster')[['# cells in cluster', 'Predicted cell type']]

        groups = pd.concat([df_votingResults['# cells in cluster'], df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0]], sort=False, axis=1).reset_index().set_index(0).groupby(by=0)

        for group in groups:
            se_cluster_size = group[1].reset_index().set_index('cluster')['# cells in cluster']
            se_type = (df_gene_cluster_centroids[se_cluster_size.index] * se_cluster_size.values).sum(axis=1) / se_cluster_size.values.sum()
            se_type.name = group[0]
            df_gene_cluster_centroids_merged = pd.concat([df_gene_cluster_centroids_merged, se_type], sort=False, axis=1)

        df_gene_cluster_centroids_merged = df_gene_cluster_centroids_merged.apply(self.zScoreOfSeries, axis=1)

        if not cluster is None:

            clusterGenes = df_gene_cluster_centroids[cluster].sort_values(ascending=False)
            clusterGenes = clusterGenes[clusterGenes >= zScoreCutoff].iloc[:top]

            return {'genes':clusterGenes.index.values.tolist(), 'zscore':clusterGenes.values.tolist()}

        if removeUnknown and self.nameForUnknown in df_gene_cluster_centroids_merged.columns:
            df_gene_cluster_centroids_merged.drop(columns=[self.nameForUnknown], inplace=True)

        df_new_marker_genes = df_gene_cluster_centroids_merged.T

        df_new_marker_genes[df_new_marker_genes < 0.] = 0.
        print('New marker cell type shape:', df_new_marker_genes.shape)

        df_new_marker_genes = df_new_marker_genes.T.loc[df_new_marker_genes.max(axis=0) >= zScoreCutoff].T
        print('New marker cell type shape:', df_new_marker_genes.shape)

        df_new_marker_genes = df_new_marker_genes.loc[df_new_marker_genes.max(axis=1) >= zScoreCutoff]
        print('New marker cell type shape:', df_new_marker_genes.shape)

        df_new_marker_list = pd.DataFrame()
        for i, celltype in enumerate(df_new_marker_genes.index.values):
            all = df_new_marker_genes.loc[celltype]
            sel = all[all > 1.0].sort_values()[:top]
            print('Candidates of %s:' % (celltype), len(sel))
            df_new_marker_list = pd.concat([df_new_marker_list, sel], sort=False, axis=1)

        df_new_marker_list = df_new_marker_list.fillna(0.).sort_index()
        df_new_marker_list[df_new_marker_list > 0.] = 1.0

        if self.nameForUnknown in df_new_marker_list.columns:
            df_new_marker_list.drop(columns=[self.nameForUnknown], inplace=True)

        df_new_marker_list.index.name = 'Marker'
        fileName = os.path.join(self.saveDir, 'new_markers.xlsx')
        print('Recording results to:', fileName)
        writer = pd.ExcelWriter(fileName)
        df_new_marker_list.to_excel(writer, 'MarkerCellType')
        df_temp = pd.Series(data=df_new_marker_list.columns, index=df_new_marker_list.columns)
        df_temp.name = 'CellTypeGrouped'
        df_temp.index.name = 'CellType'
        df_temp.to_excel(writer, 'CellTypesGrouped')
        writer.save()

        df_marker_cell_type = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), 
                                   sheet_name='Marker cell type weight matrix', index_col=0).T

        df_new_marker_genes = df_new_marker_genes[df_new_marker_list.index]

        self.makePlotOfNewMarkers(df_marker_cell_type, df_new_marker_genes, **kwargs)

        return None

    
    # Annotation functions of class ####################################################
    @classmethod
    def calculateV(cls, args):

        '''Calculate the voting scores (celltypes by clusters)
    
        Parameters:
            args: tuple
                Tuple of sub-arguments

                df_M: pandas.DataFrame
                    Marker cell type DataFrame

                df_X: pandas.DataFrame
                    Markers expression DataFrame

                cluster_index: 1d numpy.array
                    Clustering index

                cutoff: float
                    Significance cutoff, i.e. a threshold for a given marker to be significant

                giveSignificant: boolean
                    Whether to return the significance matrix along with the scores

                removeLowQCscores: boolean
                    Whether to remove low quality scores, 
                    i.e. those with less than 10% of markers that a re supporting

        Returns:
            pandas.DataFrame 
                Contains voting scores per celltype per cluster 

        Usage:
            Function is used internally.

            df = calculateV((df_M, df_X, cluster_index, 0.3, False, True))
        '''

        df_M, df_X, cluster_index, cutoff, giveSignificant, removeLowQCscores = args

        df_Y = pd.DataFrame(data=df_X.values, index=df_X.index, columns=cluster_index).groupby(level=0, sort=True, axis=1).mean()

        if True:
            df_Z = (df_Y.copy() - np.mean(df_Y.values, axis=1)[:, None]) / np.std(df_Y.values, axis=1)[:, None]
            df_Z = 1. * (df_Z > cutoff)
        else:
            df_Z = df_Y > ((1. + cutoff) * np.mean(df_Y.values, axis=1)[:, None])

        df_V = pd.Series(df_M.dot(df_Z).stack().fillna(0.))

        # Remove scores with small number of supporting markers
        if removeLowQCscores:
            df_Q = (1. * (df_M > 0.)).dot(1. * (df_Z > 0.)) < (0.1 * (df_M > 0.).sum(axis=1)[:, None]).astype(int)
            #df_Q = (1. * (df_M > 0.)).dot(1. * (df_Z > 0.)) < (0.1 * (1. * (df_Z > 0.)).sum(axis=0)[None, :]).astype(int)
            df_V[pd.Series(df_Q.stack().fillna(0.))] = 0.

        if giveSignificant:
            return df_V, df_Z

        return df_V

    def annotateWith_pDCS_Scheme(self, df_markers_expr, df_marker_cell_type):
        
        '''Produce cluster annotation results

        Parameters:
            df_markers_expr: pandas.DataFrame 
                Data with marker genes by cells expression

            df_marker_cell_type: pandas.DataFrame 
                Data with marker genes by cell types

        Returns:
            tuple
        
        Usage:
            Function should be called internally only
        '''

        df_positive = df_marker_cell_type.copy()
        df_negative = df_marker_cell_type.copy()
        df_positive[df_marker_cell_type < .0] = 0.
        df_negative[df_marker_cell_type > .0] = 0.
        df_positive /= (df_positive > 0.).sum(axis=0).fillna(1.).replace(0., 1.).replace(0, 1.)
        df_negative /= (df_negative < 0.).sum(axis=0).fillna(1.).replace(0., 1.).replace(0, 1.)
        df_marker_cell_type = df_positive + df_negative

        print(df_marker_cell_type.shape)

        # Align df_markers_expr and df_marker_cell_type
        df_markers_expr.sort_index(inplace=True, axis=0)
        df_marker_cell_type.sort_index(inplace=True, axis=1)

        # Calculate score (Vkc) for the best clustering
        temp = self.calculateV((df_marker_cell_type, 
                                df_markers_expr, 
                                df_markers_expr.columns.get_level_values('cluster').values, 
                                self.zScoreCutoff, True, self.removeLowQualityScores))
        df_V = temp[0].unstack()
        df_V.index.name = None
        df_Z_best = temp[1]

        dict_expressed_markers = {cluster:df_Z_best.index[df_Z_best[cluster] > 0].values for cluster in df_Z_best}

        np.random.seed(0)

        df_markers_expr_copy = df_markers_expr.copy()

        uclusters = np.unique(df_markers_expr_copy.columns.get_level_values(2).astype(str))
        uindex = pd.Series(df_markers_expr_copy.columns.get_level_values(2).astype(str)).replace(dict(zip(uclusters, range(len(uclusters))))).values
        df_markers_expr_copy.columns = pd.MultiIndex.from_arrays([df_markers_expr_copy.columns.get_level_values(0),
                                                                    df_markers_expr_copy.columns.get_level_values(1),
                                                                    uindex], names=df_markers_expr_copy.columns.names)

        batch_size = self.availableCPUsCount * 80

        if self.nSamples_pDCS < batch_size:
            self.nSamples_pDCS = batch_size
            print('Distribution size too small. Increasing it to minimum batch size: %s' % (self.nSamples_pDCS))

        N_batches = int(self.nSamples_pDCS / batch_size)

        if self.nSamples_pDCS != batch_size * N_batches:
            self.nSamples_pDCS = batch_size * N_batches
            print('Adjusting distribution size: %s' % (self.nSamples_pDCS))

        # Generate random cluster index
        randomClusterIndex = np.vstack([np.random.choice(df_markers_expr_copy.columns.get_level_values('cluster'), size=df_markers_expr_copy.shape[1], replace=False) for i in range(self.nSamples_pDCS)]).astype(int)

        # Generate random cluster configurations and calculate scores (Pkc) of those
        print('Generating null distribution')

        def process_batch(batch_range):

            print('\t', batch_range)

            tuples = [(df_marker_cell_type, df_markers_expr_copy, randomClusterIndex[i], self.zScoreCutoff, False, False) for i in batch_range]

            pool = multiprocessing.Pool(processes = self.availableCPUsCount)
            temp_random_df_V = pd.concat(pool.map(self.calculateV, tuples), sort=False, axis=1)
            pool.close()
            pool.join()

            return temp_random_df_V

        random_df_V = pd.concat([process_batch(range(i * batch_size,(i + 1) * batch_size)) for i in range(N_batches)], sort=False, axis=1)

        random_df_V.index = pd.MultiIndex.from_arrays([random_df_V.index.get_level_values(0),
                                                pd.Series(random_df_V.index.get_level_values(1)).replace(dict(zip(range(len(uclusters)), uclusters))).values], 
                                                names=[random_df_V.index.names[0], 'cluster'])

        random_df_V = random_df_V.replace(0., np.nan)

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
        sigma = pd.Series(data=np.nanstd(random_df_V.values, axis=1), index=random_df_V.index).replace(0., np.inf).unstack()
        mean = pd.Series(data=np.nanmean(random_df_V.values, axis=1), index=random_df_V.index).unstack()
        df_L = (df_V - mean) / sigma
        
        df_L[df_L < self.thresholdForUnknown_pDCS] = 0.
        df_L[df_V < 0.] = 0.

        return (df_marker_cell_type, df_markers_expr, df_L, df_V, dict_expressed_markers, df_null_distributions)

    def annotateWith_ratio_Scheme(self, df_markers_expr, df_marker_cell_type):
        
        '''Produce cluster annotation results

        Parameters:
            df_markers_expr: pandas.DataFrame 
                Data with marker genes by cells expression

            df_marker_cell_type: pandas.DataFrame 
                Data with marker genes by cell types

        Returns:
            tuple
        
        Usage:
            Function should be called internally only
        '''

        df_positive = df_marker_cell_type.copy()
        df_negative = df_marker_cell_type.copy()
        df_positive[df_marker_cell_type < .0] = 0.
        df_negative[df_marker_cell_type > .0] = 0.
        df_positive /= (df_positive > 0.).sum(axis=0).fillna(1.).replace(0., 1.).replace(0, 1.)
        df_negative /= (df_negative < 0.).sum(axis=0).fillna(1.).replace(0., 1.).replace(0, 1.)
        df_marker_cell_type = df_positive + df_negative
        print(df_marker_cell_type.shape)

        # Align df_markers_expr and df_marker_cell_type
        df_markers_expr.sort_index(inplace=True, axis=0)
        df_marker_cell_type.sort_index(inplace=True, axis=1)

        # Calculate score (Vkc) for the best clustering
        temp = self.calculateV((df_marker_cell_type, 
                                df_markers_expr, 
                                df_markers_expr.columns.get_level_values('cluster').values, 
                                self.zScoreCutoff, True, self.removeLowQualityScores))
        df_V = temp[0].unstack()
        df_V.index.name = None
        df_Z_best = temp[1]

        dict_expressed_markers = {cluster:df_Z_best.index[df_Z_best[cluster] > 0].values for cluster in df_Z_best}

        df_marker_cell_type_full = self.readMarkerFile().T.loc[df_marker_cell_type.index]
        df_marker_cell_type_full[df_marker_cell_type_full < .0] = 0.
        df_marker_cell_type_full /= (df_marker_cell_type_full > 0.).sum(axis=0).fillna(1.).replace(0., 1.)

        df_L = df_V / df_marker_cell_type_full.sum(axis=1).values[:, None]

        df_L[df_L < self.thresholdForUnknown_ratio] = 0.
        df_L[df_V < 0.] = 0.

        return (df_marker_cell_type, df_markers_expr, df_L, df_V, dict_expressed_markers, None)

    def annotateWith_Hopfield_Scheme(self, df_markers_expr, df_marker_cell_type):

        '''Produce cluster annotation results 
                
        Parameters:
                df_markers_expr: pandas.DataFrame
                    Markers expression DataFrame

                df_marker_cell_type: pandas.DataFrame
                    Marker cell type DataFrame

        Returns:
            tuple

        Usage:
            Function should be called internally only
        '''
        
        df_marker_cell_type = df_marker_cell_type.T

        np.random.seed(0)

        # The scheme does not account for negative markers.
        # Therefore remove them
        df_marker_cell_type[df_marker_cell_type > 0.] = 1.
        df_marker_cell_type[df_marker_cell_type < 0.] = 0.

        # Keep markers that have +1 or -1 for at least one cell type
        df_marker_cell_type = df_marker_cell_type.loc[np.abs(df_marker_cell_type).sum(axis=1) > 0]

        # Check and keep only cell types with at least a minium number of +1 markers
        df_marker_cell_type = df_marker_cell_type[df_marker_cell_type.columns[np.abs(df_marker_cell_type).sum(axis=0) > self.minimumNumberOfMarkersPerCelltype]]

        print(df_marker_cell_type.shape)

        # Remove unused markers
        df_marker_cell_type = df_marker_cell_type.loc[np.abs(df_marker_cell_type).sum(axis=1) > 0]

        # Normalize marker/celltype matrix in both direction by counts.
        # Note that order of the following two lines is irrelevant
        df_marker_cell_type[df_marker_cell_type > 0.] /= (df_marker_cell_type[df_marker_cell_type > 0.] > 0).sum(axis=1).replace(0, 1.).values[:, None]
        df_marker_cell_type[df_marker_cell_type > 0.] /= (df_marker_cell_type[df_marker_cell_type > 0.] > 0).sum(axis=0).replace(0, 1.).values
        df_marker_cell_type[df_marker_cell_type < 0.] /= (df_marker_cell_type[df_marker_cell_type < 0.] < 0).sum(axis=1).replace(0, 1.).values[:, None]
        df_marker_cell_type[df_marker_cell_type < 0.] /= (df_marker_cell_type[df_marker_cell_type < 0.] < 0).sum(axis=0).replace(0, 1.).values

        print(df_marker_cell_type.shape)

        # Order gene expression index in the same order as markers index
        df_markers_expr = df_markers_expr.loc[df_marker_cell_type.index]

        df_expr_cluster_centroids = df_markers_expr.groupby(level='cluster', axis=1, sort=False).mean()
        df_Z = df_expr_cluster_centroids.copy().apply(self.zScoreOfSeries, axis=1)

        df_sig = df_Z > self.zScoreCutoff

        dict_expressed_markers = {cluster:df_sig.index[df_sig[cluster] > 0].values for cluster in df_sig}

        if self.useUnderlyingNetwork:
            try:
                underlyingNetwork = self.getSubnetworkOfPCN(df_marker_cell_type.index.values, min_shared_first_targets=1)

                print('Recording underlying network', flush=True)
                underlyingNetwork.to_hdf(self.fileHDFpath, key='df_underlyingNetwork', mode='a', complevel=4, complib='zlib')

                underlyingNetwork = underlyingNetwork.values
            except Exception as exception:
                print(exception)
                print('Not using underlying network connection')
                underlyingNetwork = None
        else:
            underlyingNetwork = None

        print('Generating Hopfield networks:', flush=True)

        listOfSeries = []

        for i in range(self.nSamples_Hopfield):
            if i % np.int(self.nSamples_Hopfield / 10) == 0:
                print('\n%s%%' % (np.int(100 * i / self.nSamples_Hopfield)), end='\t', flush=True)

            # Run Hopfield network dynamics
            hop = self.propagateHopfield(sigma=df_sig.values * 2. - 1., 
                                         xi=df_marker_cell_type.values, 
                                         typesNames=df_marker_cell_type.columns.values, underlyingNetwork=underlyingNetwork,
                                         clustersNames=df_sig.columns.values, id = i, T = self.HopfieldTemperature, 
                                         path = os.path.join(self.saveDir, 'HopfieldTrajectories'),
                                         recordTrajectories = False)

            # Determine top scores
            res = np.argmax(hop, axis=0)

            # Best votes that below threshold are set to Unknown
            res[np.where(np.max(hop, axis=0) < 10. ** -12)] = len(df_marker_cell_type.columns)

            listOfSeries.append(res)

        print(flush=True)

        r = np.vstack(listOfSeries)

        df_L = pd.DataFrame(index=list(range(len(df_marker_cell_type.columns) + 1)), columns=df_sig.columns).fillna(0.)
        for i in range(len(df_sig.columns)):
            w = np.unique(r[:,i], return_counts=True)
            for j in range(len(w[0])):
                df_L.iloc[w[0][j], i] = w[1][j] / len(r)

        df_L[df_L < self.thresholdForUnknown_Hopfield] = 0.
        df_L.iloc[-1, :] = 0.9 * self.thresholdForUnknown_Hopfield

        # Best results are chosen as highest frequency cell type for each cluster
        se_celltypes = pd.Series(data = np.argmax(df_L.values, axis=0).astype(int).astype(str), index = df_sig.columns)

        # Define cell names dictionary
        celltypeNames = dict(zip(np.array(range(len(df_marker_cell_type.columns))).astype(int).astype(str), df_marker_cell_type.columns))
        celltypeNames.update({str(len(celltypeNames)): self.nameForUnknown})

        # Convert cell names from index to words
        se_celltypes = se_celltypes.replace(to_replace=celltypeNames)

        df_L.index = pd.Series(df_L.index.astype(str)).replace(to_replace=celltypeNames).values

        df_L.drop(self.nameForUnknown, axis=0, inplace=True)

        return (df_marker_cell_type.T, df_markers_expr, df_L, df_L.copy(), dict_expressed_markers)

    def recordAnnotationResults(self, df_marker_cell_type, df_markers_expr, df_L, df_V, dict_expressed_markers, df_null_distributions = None):

        '''Record cell type annotation results to spreadsheets.

        Parameters:
            df_marker_cell_type: pandas.DataFrame
                Markers to cell types table

            df_markers_expr: pandas.DataFrame
                Markers expression in each cluster

            df_L: pandas.DataFrame
                Annotation scores along with other information

            df_V: pandas.DataFrame
                Annotation scores along with other information

            dict_expressed_markers: dictionary
                Dictionary of markers signigicantly expressed in each cluster

            df_null_distributions: pandas.DataFrame, Default None
                Table with null distributions

        Returns:
            None

        Usage:
            This function is intended to be used internally only
        '''

        df_L = df_L.fillna(0.)

        # Determine cluster sizes
        se_cluster_sizes = df_markers_expr.groupby(level='cluster', sort=True, axis=1).count().iloc[0]
        cluster_sizes = se_cluster_sizes.values
        
        df_L = df_L[se_cluster_sizes.index]
        df_V = df_V[se_cluster_sizes.index]

        # Rename winning cell types
        T = df_L.index[np.argmax(df_L.values, axis=0)].values
        T[(df_L.values <= 0.).all(axis=0)] = self.nameForUnknown
        T = pd.Index(T) + ' #0'
        for i in range(len(T)):
            T = T.where(~T.duplicated(), T.str.replace(' #%s' % (i), ' #%s' % (i + 1)))
        predicted_celltype = T.values

        # Save score columns names
        columns_scores = df_V.index.values.copy().tolist()
        
        # Determine winning score
        winning_score = np.array([df_V.iloc[ind[0], ind[1]] for ind in np.flip(np.array(list(enumerate(np.argmax(df_L.values, axis=0)))), axis=1)])

        # Update the DataFrames
        df_V.loc['Winning score'] = df_L.loc['Winning score'] = winning_score
        df_V.loc['Predicted cell type'] = df_L.loc['Predicted cell type'] = predicted_celltype
        df_V.loc['# cells in cluster'] = df_L.loc['# cells in cluster'] = cluster_sizes

        # Determine all markers hits for each cluster
        df_custer_centroids = df_markers_expr.groupby(level='cluster', sort=True, axis=1).mean()

        all_markers_in_cluster = [(' // ').join(dict_expressed_markers[cluster].tolist()) for cluster in df_V.columns]
        df_V.loc['All markers'] = df_L.loc['All markers'] = all_markers_in_cluster

        all_supp_markers = {celltype:df_marker_cell_type.T.index[df_marker_cell_type.T[celltype] > 0].values for celltype in columns_scores}
        all_cont_markers = {celltype:df_marker_cell_type.T.index[df_marker_cell_type.T[celltype] < 0].values for celltype in columns_scores}
        all_supp_markers[self.nameForUnknown] = np.array([])
        all_cont_markers[self.nameForUnknown] = np.array([])

        supporting_markers = [(' // ').join((np.intersect1d(all_supp_markers[df_V.loc['Predicted cell type'].loc[cluster].split(' #')[0]],
                                                            df_V.loc['All markers'].loc[cluster].split(' // '))).tolist()) for cluster in df_V.columns]
        contradicting_markers = [(' // ').join((np.intersect1d(all_cont_markers[df_V.loc['Predicted cell type'].loc[cluster].split(' #')[0]],
                                                            df_V.loc['All markers'].loc[cluster].split(' // '))).tolist()) for cluster in df_V.columns]

        df_V.loc['Supporting markers'] = df_L.loc['Supporting markers'] = supporting_markers
        df_V.loc['Contradicting markers'] = df_L.loc['Contradicting markers'] = contradicting_markers
        
        # Record the results to spreadsheets
        fileName = os.path.join(self.saveDir, self.dataName + '_annotation.xlsx')
        print('Recording voting results to:', fileName)
        writer = pd.ExcelWriter(fileName)
        df_L.columns.name = 'cluster'
        df_L.T.to_excel(writer, 'z-scores', columns=['Predicted cell type', 
                                                     '# cells in cluster', 
                                                     'Winning score', 
                                                     'Supporting markers', 
                                                     'Contradicting markers', 
                                                     'All markers'] + columns_scores)
        df_V.columns.name = 'cluster'
        df_V.T.to_excel(writer, 'Voting scores', columns=['Predicted cell type', 
                                                          '# cells in cluster', 
                                                          'Winning score', 
                                                          'Supporting markers', 
                                                          'Contradicting markers', 
                                                          'All markers'] + columns_scores)

        df_custer_centroids.T.to_excel(writer, 'Cluster centroids')

        if not df_null_distributions is None:
            df_null_distributions.to_excel(writer, 'Null distributions')

        df_marker_cell_type.to_excel(writer, 'Marker cell type weight matrix')

        writer.save()

        return df_V.loc['Predicted cell type'].to_dict()

    def propagateHopfield(cls, sigma = None, xi = None, T = 0.2, tmax = 200, fractionToUpdate = 0.5, mode = 4, meshSamplingRate = 200,
                            underlyingNetwork = None, typesNames = None, clustersNames = None, printInfo = False,
                            recordTrajectories = True, id = None, printSwitchingFraction = False, path = None):

        '''Function is used internally to propagate Hopfield network 
        over a set number of time steps
                    
        Parameters:
            sigma: pandas.DataFrame, Default None
                Markers expression

            xi: pandas.DataFrame, Default None
                Marker cell type DataFrame

            T: float, Default 0.2
                Noise (Temperature) parameter

            tmax: int, Default 200
                Number of step to iterate through

            fractionToUpdate: float, Default 0.5
                Fraction of nodes to randomly update at each iteration

            mode: int, Default 4
                Options are:
                    1: non-onthogonalized, non-weighted attractors
                    2: onthogonalized, non-weighted attractors
                    3: onthogonalized, weighted attractors
                    4: onthogonalized, weighted attractors, asymetric and diluted dynamics

            meshSamplingRate: int, Default 100
                Visualization parameter to control the quality of the color mesh near the attractors

            underlyingNetwork: 2d numpy.array, Default None
                Network of underlying connections between genes

            typesNames: list-like, Default None
                Names of cell types

            clustersNames: list-like, Default None
                Names or identifiers of the clusters

            printInfo: boolean, Default False
                Whether to print detailes

            recordTrajectories: boolean, Default True
                Whether to record trajectories data to files

            id: int, Default None
                Identifier of this function call

            printSwitchingFraction: boolean, Default False
                Whether to print fraction of clusters that switch theie 
                maximum overlapping attractor

            path: str, Default None
                Path for saving trajectories data

        Returns:
            2d numpy.array 
                Overlaps 

        Usage:
            result = propagateHopfield(sigma=sigma, xi=df_attrs)
        '''
          
        if xi is None:
            print('xi is None')

            return

        if not path is None:           
            if not os.path.exists(path):
                os.makedirs(path)

        if mode==1 or mode==2:
            xi[xi > 0.] = 1.
            xi[xi <= 0.] = -1.

        if printInfo:
            print('Unique nodes across Hopfield network:', len(np.unique(xi, axis=0)))

        if mode == 1:
            A = xi.T
            J = (xi).dot(A) / (xi.shape[0]) # m,m
        else:
            Q = (xi.T).dot(xi)              # n,n

            try:
                Q_inv = np.linalg.inv(Q)    # n,n
            except Exception as exception:
                print(exception)
                print('Cannot invert matrix Q')

                return

            A = (Q_inv).dot(xi.T)           # n,m <-- xi_inv
            J = (xi).dot(A)                 # m,m

        if not underlyingNetwork is None:
            if J.shape == underlyingNetwork.shape:
                #print('Udsing network of underlying interactions')
                J *= underlyingNetwork
            else:
                print('Network of underlying interactions shape is inconsistent with J. Ignoring it')

        getHopfieldEnergy = lambda s: -0.5 * s[None, :].dot(J).dot(s)

        if recordTrajectories:
            pca = PCA(n_components=None).fit(A)

            attrs = pca.transform(A)
            variance = pca.explained_variance_ratio_

            if id == 0 or id is None:
                if xi.shape[1] == 2:
                    dim = 2
                    lim = 15
                    g = np.linspace(-lim, lim, num=100)
                    mesh_coords = np.meshgrid(*tuple([g]*dim))
                    mesh_coords = np.vstack([g.flatten() for g in mesh_coords]).T
                    
                    mesh = np.zeros((mesh_coords.shape[0], xi.shape[1]))
                    mesh[:, :len(mesh_coords.T)] = mesh_coords

                    mesh_sigma = pca.inverse_transform(mesh)

                    if False:
                        mesh_sigma[mesh_sigma > 0.] = 1.
                        mesh_sigma[mesh_sigma <= 0.] = -1.
                        digital_mesh = pca.transform(mesh_sigma)
                        mesh_coords = digital_mesh
                    
                    mesh_energy = np.array([getHopfieldEnergy(s) for s in mesh_sigma])
                    write(np.hstack([mesh_coords[:,:2], mesh_energy]), os.path.join(path, 'mesh'))
                else:
                    print('Generating samples in the vicinity of attractors')
                    temp = []
                    for i in range(xi.shape[1]):
                        for j in range(25):
                            st = A[i]
                            temp.append(np.array([st * ((np.random.rand(st.shape[0]) > 0.025*j) * 2. - 1.) for i in range(meshSamplingRate)]))
                    mesh_sigma = np.vstack(temp)

                    print('Calculating energy of the samples')
                    if False:
                        mesh_energy = np.array([getHopfieldEnergy(s) for s in mesh_sigma])
                    else:
                        print('Calculating chunks')

                        def func(mesh_sigma):

                            return -0.5 * ((J).dot(mesh_sigma.T) * mesh_sigma.T).sum(axis=0)[:,None]

                        mesh_energy = np.vstack([func(item) for item in np.split(mesh_sigma, meshSamplingRate)])

                    print('\nProjecting samples on PCs')
                    mesh_coords = pca.transform(mesh_sigma)

                    print('Recording mesh data')
                    write(np.hstack([mesh_coords, mesh_energy]), os.path.join(path, 'mesh'))

        if recordTrajectories:
            write((np.vstack([attrs, variance]), typesNames), os.path.join(path, 'attrs'))

        if sigma is None:
            print('Sigma is None')

            return

        if recordTrajectories:
            trajectories = np.zeros((tmax, sigma.shape[1], xi.shape[1]))
                    
        initial = np.argmax((A).dot(sigma), axis=0)
        initial[np.where(np.max((A).dot(sigma), axis=0) < 10. ** -12)] = -1

        if printInfo:
            sIn = sigma.copy()

        for t in range(tmax): 
            h = (J).dot(sigma)

            #print(np.round((100.*(sigma==0.)*1.).sum().sum()/(sigma.shape[0]*sigma.shape[1]), 2), end=', ', flush=True)
                
            if recordTrajectories:
                if id == 0 or id is None:
                    trajectories[t,:,:] = np.hstack([pca.transform(sigma.T)])

            sigmaNew = ((1. / (1. + np.exp(-2. * h / T))) > np.random.rand(*sigma.shape)) * 2. - 1.

            # Make dynamics asymmetric and diluted: Deactivate nodes that are being turned off
            if mode == 4:
                sigmaNew[(sigmaNew - sigma) < -1.] = 0.

            f = 1. if t >= 100 else max(0.05, float(t) / 100.)

            whereToUpdate = (np.random.rand(*sigma.shape) <= f * fractionToUpdate)
            sigma[whereToUpdate] = sigmaNew[whereToUpdate] * np.abs(sigma[whereToUpdate])

        if printInfo:
            print('In positive:', ((sIn>0)*1.).sum(axis=0))
            print('Out positive:', ((sigma>0)*1.).sum(axis=0))
            print('Positive common in-out:', (((sIn>0.)*(sigma>0.))*1.).sum(axis=0))

        final = np.argmax((A).dot(sigma), axis=0)
        final[np.where(np.max((A).dot(sigma), axis=0) < 10. ** -12)] = -1

        if printSwitchingFraction:
            print(np.round(100. * (final != initial).sum() / (len(initial)), 2).astype(int), end=' ', flush=True)

        if recordTrajectories:
            if id == 0 or id is None:
                write(trajectories, os.path.join(path, 'trajectories%s')%(id))
                write((initial, final, typesNames, clustersNames), os.path.join(path, 'additional'))

        overlap = (A).dot(sigma)
        overlap[overlap <= 0.] = 0.

        return overlap


    # Other functions of class #########################################################
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
            DCS = DigitalCellSorter.DigitalCellSorter()

            se = DCS.zScoreOfSeries(se)
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
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.KeyInFile('df_expr', 'data/file.h5')
        '''

        if not os.path.isfile(file):
            return False 

        with pd.HDFStore(file) as file:
            return True if "/" + key.strip("/") in file.keys() else False

        return

    def getSubnetworkOfPCN(self, subnetworkGenes, min_shared_first_targets = 30):

        '''Extract subnetwork of PCN network

        Parameters:
            subnetworkGenes: list-like
                Set of genes that the subnetwork should contain

            min_shared_first_targets: int, Default 30
                Number of minimum first shared targets to connect two nodes

        Returns:
            pandas.DataFrame
                Adjacency matrix

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_subnetwork = DCS.getSubnetworkOfPCN(genes)
        '''

        print('\nReading PCN (directed, unweighted) network from file')
        PCN = pd.read_csv(os.path.join('data', 'PCN_%s.txt.gz' % (self.species)), compression='gzip', header=0, index_col=[0], delimiter='\t')['Target']

        PCN.index = self.gnc.Convert(PCN.index.tolist(), 'alias', 'hugo', returnUnknownString=False)
        PCN[:] = self.gnc.Convert(PCN.values.tolist(), 'alias', 'hugo', returnUnknownString=False)

        print('Calculating targets of PCN network genes')
        targets = PCN.groupby(level=0).agg('unique').apply(set).to_dict()

        del PCN

        index = pd.Index(set(subnetworkGenes).intersection(set(targets.keys()))).sort_values()
        print('Number of PCN genes: %s'%(len(index)))

        print('Calculating number of common targets for %s genes'%(len(index)))
        data = np.zeros((len(index), len(index))).astype(int)

        for iA in range(len(index)):
            #print(iA, end=', ', flush=True)
            geneA = index[iA]

            for iB in range(iA, len(index)):
                geneB = index[iB]
                direct = 1 if (geneB in targets[geneA]) else 0
                data[iA,iB] = data[iB,iA] = len(targets[geneA].intersection(targets[geneB])) + direct

        targetsOverlaps = pd.DataFrame(index=index, columns=index, data=data)
        #print(targetsOverlaps)

        subnetwork = 1.*(targetsOverlaps >= min_shared_first_targets)
        print(subnetwork.shape)

        subnetwork = subnetwork.reindex(subnetworkGenes, axis=0).reindex(subnetworkGenes, axis=1).fillna(0.)
        print(subnetwork.shape)

        con = subnetwork.sum().sum()
        tot = (len(subnetwork)-1)*len(subnetwork)
        print('\nConnections in subnetwork: %s, out of %s possible (%s%%)\n'%(con, tot, np.round(100.*con/tot, 0)))

        return subnetwork

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

            df = DCS.alignSeries(pd.Index(['A', 'B', 'C', 'D']).to_series(), pd.Index(['B', 'C', 'D', 'E', 'F']).to_series())
        '''
        
        se1.index.name = 'index'
        se2.index.name = 'index'

        append = lambda se1, se2: pd.concat([se1, pd.Series(index=se2.index.difference(se1.index), data=[tagForMissing] * len(se2.index.difference(se1.index)))], axis=0, sort=False)

        se1 = append(se1, se2)
        se2 = append(se2, se1)

        se1.name = 'se1'
        se2.name = 'se2'

        return pd.concat((se1, se2.loc[se1.index]), axis=1, sort=True)

    def createReverseDictionary(self, inputDictionary):

        '''Efficient way to create a reverse dictionary from a dictionary.
        Utilizes Pandas.Dataframe.groupby and Numpy arrays indexing.
    
        Parameters: 
            inputDictionary: dictionary
                Dictionary to reverse

        Returns:
            dictionary
                Reversed dictionary

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            revDict = DCS.createReverseDictionary(Dict)
        '''

        keys, values = np.array(list(inputDictionary.keys())), np.array(list(inputDictionary.values()))
        df = pd.DataFrame(np.array([[keys[i], value] for i in range(len(keys)) for value in values[i]]))
        dfGrouped = df.groupby(df.columns[1])
        keys, values = list(dfGrouped.indices.keys()), list(dfGrouped.indices.values())
        GOs = df.values.T[0]

        return dict(zip(keys, [GOs[value].tolist() for value in values]))

    def readMarkerFile(self, mergeFunction = 'mean', mergeCutoff = 0.25):

        '''Read markers file, prepare markers

        Parameters: 
            mergeCutoff: str, Default 'mean'
                Function used for grouping of the cell sub-types. Options are:
                    'mean': average of the values
                    'max': maxium of the values, effectively a logiacal OR function

            mergeCutoff: float, Default 0.25
                Values below cutoff are set to zero. 
                This option is used if mergeCutoff is 'mean'

        Returns:
            pandas.DataFrame
                Celltype/markers matrix

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_marker_cell_type = DCS.readMarkerFile()
        '''

        df_marker_cell_type = pd.read_excel(self.geneListFileName, index_col=0, header=[0,1]).replace(np.nan, 0.).replace(0, 0.)
        df_marker_cell_type.columns.names = ['CellTypeGrouped', 'CellType']

        df_marker_cell_type.index = self.gnc.Convert(list(df_marker_cell_type.index), 'alias', 'hugo', returnUnknownString=False)

        df_marker_cell_type = df_marker_cell_type.drop(columns=[col for col in df_marker_cell_type.columns if col[0] == 'NA'])
        df_marker_cell_type = df_marker_cell_type.groupby(level='CellTypeGrouped', axis=1).agg(mergeFunction).fillna(0.)

        where_positive = df_marker_cell_type >= mergeCutoff
        where_negative = df_marker_cell_type <= -mergeCutoff
        df_marker_cell_type *= 0.
        df_marker_cell_type[where_positive] = 1.
        df_marker_cell_type[where_negative] = -1.

        df_marker_cell_type = df_marker_cell_type.loc[np.abs(df_marker_cell_type).sum(axis=1) > 0]
        print('Markers/celltypes:', df_marker_cell_type.shape, flush=True)

        # Merge duplicates that might have appeared after gene name conversion
        df_marker_cell_type = df_marker_cell_type.groupby(level=0, axis=0).sum()
        df_marker_cell_type[df_marker_cell_type > 0.] = 1.
        df_marker_cell_type[df_marker_cell_type < 0.] = -1.

        def norm(s_input):
            
            s = s_input.copy()

            pos_sum = np.abs(s.iloc[np.where(s > 0.)].sum())
            neg_sum = np.abs(s.iloc[np.where(s < 0.)].sum())
            
            s.iloc[np.where(s > 0.)] /= pos_sum if pos_sum > 0. else 1.
            s.iloc[np.where(s < 0.)] /= neg_sum if neg_sum > 0. else 1.

            return s

        df_marker_cell_type = df_marker_cell_type.apply(norm, axis=0)

        print('Markers/celltypes:', df_marker_cell_type.shape, flush=True)

        df_marker_cell_type.columns = df_marker_cell_type.columns.str.strip()

        if not self.useNegativeMarkers:
            df_marker_cell_type[df_marker_cell_type < 0.] = 0.
            df_marker_cell_type = df_marker_cell_type.loc[np.abs(df_marker_cell_type).sum(axis=1) > 0]
            print('Removed negative markers. Markers/celltypes:', df_marker_cell_type.shape, flush=True)

        return df_marker_cell_type
    
    def mergeIndexDuplicates(cls, df_expr, method = 'average', printDuplicates = False):

        '''Merge index duplicates

        Parameters: 
            df_expr: pandas.DataFrame
                Gene expression table

            method: str, Default None
                How to deal with index duplicates. Option are:
                    'average': average values of duplicates

                    'first': keep only first of duplicates, discard rest

        Returns:
            pandas.DataFrame
                Gene expression table

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            df_expr = DCS.mergeIndexDuplicates(df_expr)
        '''

        # Check is any names in the index are duplicated, remove duplicates
        len_total, len_unique = len(df_expr.index), len(np.unique(df_expr.index))

        if len_total != len_unique:
            unique_items = np.unique(df_expr.index, return_counts=True)

            if method == 'average':
                df_expr = df_expr.groupby(level=0, axis=0).sum()

            elif method == 'first':
                df_expr = df_expr.loc[~df_expr.index.duplicated(keep='first')]

            else:
                print('Unknown method')

                return df_expr

            print('Merged %s duplicated items in the index of size %s' % (len_total - len_unique, len_total), flush=True)
            
            if printDuplicates:
                print(unique_items[0][unique_items[1] > 1], flush=True)

        return df_expr

    def recordExpressionData(self):

        '''Record expression data from the internal HDF storage.

        Parameters: 
            None

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.recordExpressionData()
        '''

        print('Recording compressed DataFrame', flush=True)

        self._df_expr.replace(0, np.nan).T.stack().to_frame().to_hdf(self.fileHDFpath, key='df_expr', mode='a', complevel=4, complib='zlib')

        return

    def loadAnnotatedLabels(self, detailed = False, includeLowQC = True, infoType = 'label'):

        '''Load cell annotations resulted from function 'annotate'

        Parameters: 
            detailed: boolean, Default False
                Whether to give cluster- or celltype- resolution data

            includeLowQC: boolean, Default False
                Whether to include low quality cells in the output

        Returns:
            pandas.Series

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.loadAnnotatedLabels()
        '''

        if self.KeyInFile('df_markers_expr', self.fileHDFpath):
            try:
                se = pd.read_hdf(self.fileHDFpath, key='df_markers_expr').T.reset_index().set_index(['batch', 'cell'])['label']
            except Exception as exception:
                print(exception)
                print('Error reading annotation results', flush=True)

                return

            if not detailed:
                se = se.str.split(' #', expand=True)[0]

        else:
            print('Annotation results not found', flush=True)

            return

        if includeLowQC:
            if self.KeyInFile('df_projection_pre_QC', self.fileHDFpath):
                allCells = pd.read_hdf(self.fileHDFpath, key='df_projection_pre_QC').columns
                allBathes = allCells.get_level_values('batch')
            elif self.KeyInFile('df_projection', self.fileHDFpath):
                allCells = pd.read_hdf(self.fileHDFpath, key='df_projection').columns
                allBathes = allCells.get_level_values('batch')
            else:
                print('Labelled data not found', flush=True)

                return

            if infoType=='label':
                se = se.reindex(allCells, fill_value=self.nameForLowQC)
            elif infoType=='batch':
                se_low = pd.Series(data=allBathes, index=allCells)
                se.iloc[:] = se.index.get_level_values('batch').values
                se = pd.concat([se, se_low.loc[se_low.index.difference(se.index)]])

        return se
    
    def loadExpressionData(self):

        '''Load processed expression data from the internal HDF storage.

        Parameters: 
            None

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.loadExpressionData()
        '''

        if self.KeyInFile('df_expr', self.fileHDFpath):
            print('Loading processed data', flush=True)

            self._df_expr = pd.read_hdf(self.fileHDFpath, key='df_expr').unstack(level=-1, fill_value=0.).T
            self._df_expr.index = self._df_expr.index.get_level_values(-1)
            self._df_expr.sort_index(inplace=True)
        else:
            print('Processed data not found')

        return

    def prepareMarkers(self, expressedGenes = None, createColormapForCelltypes = True):

        '''Get dictionary of markers for each cell types.
        
        Parameters:
            expressedGenes: pandas.Index, Default None
                If not None then the marker DataFrame will be intersected with this index,
                i.e. all non-expressed genes will be filtered from the marker file

            createColormapForCelltypes: boolean, Default True
                Create (or update) a colormap for cell types
                based on a marker-celltype matrix. This will make coloring of cell clusters
                consistent across all plots.

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.prepareMarkers()
        '''

        df_marker_cell_type = self.readMarkerFile()

        # Remove non-expressed markers
        if not expressedGenes is None:
            df_marker_cell_type = df_marker_cell_type.loc[df_marker_cell_type.index.intersection(expressedGenes)]
            print('Removed non-expressed genes. Markers/celltypes shape:', df_marker_cell_type.shape, flush=True)

        # Remove underepresented cell types
        df_marker_cell_type = df_marker_cell_type[df_marker_cell_type.columns[(df_marker_cell_type > 0.).sum(axis=0) >= self.minimumNumberOfMarkersPerCelltype]]
        df_marker_cell_type = df_marker_cell_type[(df_marker_cell_type != 0).sum(axis=1) >= 1]
        df_marker_cell_type.index.name = 'Marker'

        df_marker_cell_type.to_hdf(self.fileHDFpath, key='df_marker_cell_type', mode='a', complevel=4, complib='zlib')

        # Create (or update) a colormap for cell types
        if createColormapForCelltypes:
            with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'w') as tempMapFile:
                df_marker_cell_type = pd.read_hdf(self.fileHDFpath, key='df_marker_cell_type')

                cellTypes = np.sort(df_marker_cell_type.columns.values.tolist())

                for cellTypeIndex in range(len(cellTypes)):
                    tempMapFile.write(cellTypes[cellTypeIndex] + '\t' + str(cm.jet(cellTypeIndex / len(cellTypes))) + '\n')

                tempMapFile.write(self.nameForUnknown + '\t' + str(cm.jet(1.0)) + '\n')

        return
    
    def calculateQCmeasures(self):

        '''Calculate Quality Control (QC) measures

        Parameters:
            None

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.calculateQCmeasures()
        '''

        if self.mitochondrialGenes is None:
            if self.species == 'Human':
                mitoGenesPath = os.path.join(self.defaultGeneListsDir, 'Human.MitoCarta2.0.csv')
            elif self.species == 'Mouse':
                mitoGenesPath = os.path.join(self.defaultGeneListsDir, 'Mouse.MitoCarta2.0.csv')
            else:
                print('Mitochondrial genes for scpecies %s not found' % (self.species))
                print('Implemented lists are for Human and Mouse')

                raise NotImplementedError

            self.mitochondrialGenes = pd.read_csv(mitoGenesPath, index_col=None, header=0)['Symbol'].values.squeeze().tolist()

        print('Calculating quality control measures (count depth, number of genes, fraction of mitochondrial genes) for each cell', flush=True)

        df_QC = pd.concat([(self._df_expr).sum(axis=0), 
                           (self._df_expr > 0).sum(axis=0), 
                           (self._df_expr.loc[self._df_expr.index.intersection(self.mitochondrialGenes)] > 0).sum(axis=0) / \
                              (self._df_expr > 0).sum(axis=0)], axis=1, sort=False)

        df_QC.columns = ['count_depth', 'number_of_genes', 'fraction_of_mitochondrialGenes']

        df_QC.to_hdf(self.fileHDFpath, key='QC', mode='a', complevel=4, complib='zlib')

        return

    def qualityControl(self, **kwargs):

        '''Remove low quality cells

        Parameters:
            None

        Returns:
            Any parameters that function 'getIndexOfGoodQualityCells' can accept

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.qualityControl()
        '''

        index = self.getIndexOfGoodQualityCells(**kwargs)

        if self.toggleRemoveLowQualityCells:
            self._df_expr = self._df_expr[self._df_expr.columns.intersection(index).sort_values()]
            self._df_expr = self._df_expr[self._df_expr.sum(axis=1) > 0]
            print('Removed low quality cells. Data size: %s genes, %s cells' % self._df_expr.shape, flush=True)

            df_projection = pd.read_hdf(self.fileHDFpath, key='df_projection_pre_QC')
            df_projection[self._df_expr.columns].to_hdf(self.fileHDFpath, key='df_projection', mode='a', complevel=4, complib='zlib')

            self.project(PCAonly=True)

        return

    def batchEffectCorrection(self, method = 'COMBAT'):

        '''Batch effect correction.

        Parameters:
            method: str, Default 'COMBAT'
                Stein, C.K., Qu, P., Epstein, J. et al. Removing batch effects from purified 
                plasma cell gene expression microarrays with modified ComBat. 
                BMC Bioinformatics 16, 63 (2015)

        Returns:
            None

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.batchEffectCorrection()
        '''

        if method == 'COMBAT':
            cells = self._df_expr.columns.get_level_values('cell').values.copy()
            patients = self._df_expr.columns.get_level_values('batch').values.copy()

            if len(np.unique(patients)) == 1:
                print('Only one batch provided. Batch correction is not necessary', flush=True)
            else:
                print('ComBat transformation', flush=True)
                print('Reading positions of zeros', flush=True)
                where_zeros = self._df_expr == 0.
                values = combat(pd.DataFrame(data=self._df_expr.values.copy(), index=self._df_expr.index.values.copy(), columns=range(len(cells))), pd.Series(data=patients, index=range(len(cells)))).values
                self._df_expr = pd.DataFrame(data=values, index=self._df_expr.index, columns=self._df_expr.columns)
                print('Setting original zeros back to zeros', flush=True)
                self._df_expr[where_zeros] = 0.
                print('Setting negative values to zeros', flush=True)
                self._df_expr[self._df_expr < 0.0] = 0.

        else:
            print('Batch effect correction method unknown')

        return
