import os
import platform
import copy
import pickle
import gzip
import multiprocessing
import time

import numpy as np
import pandas as pd

import scipy.stats
import scipy.signal
import scipy.cluster.hierarchy
import scipy.spatial.distance
from scipy.interpolate import UnivariateSpline

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib import cm

import plotly.graph_objects

import sklearn.metrics
import sklearn.cluster
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering
from sklearn.ensemble import IsolationForest

from DigitalCellSorter import GeneNameConverter
from DigitalCellSorter.Combat import combat

def timeMark():

    '''Print total time elapsed from the beggining of the process
        
    Args:
        None

    Returns:
        None

    Usage:
        timeMark()
    '''
        
    return print('--> Total elapsed time: %s min' % (np.round(time.thread_time() / 60., 1)), '\n')

def getStartTime():

    '''Get time (in seconds) elapsed from the epoch
    
    Args:
        None

    Returns:
        Time (in seconds)

    Usage:
        start = getStartTime()
    '''

    return time.time()
    
def getElapsedTime(start):

    '''Print total elapsed time (in minutes) elapsed from the reference point
    
    Args:
        start: reference time (in seconds)

    Returns:
        None

    Usage:
        getElapsedTime(start)
    '''

    return print('Elapsed time: ' + str(np.round((time.time() - start) / 60.,1)) + ' min' + '\n')

class DigitalCellSorter:

    '''Class of Digital Cell Sorter with methods to process single cell mRNA-seq data.
    Includes analysis and visualization tools.
    
    Usage:
        DCS = DigitalCellSorter.DigitalCellSorter()
        df_data = DCS.Clean(df_data)
    '''

    gnc = GeneNameConverter.GeneNameConverter(dictDir=os.path.join('DigitalCellSorter', 'pickledGeneConverterDict', 'ensembl_hugo_entrez_alias_dict.pythdat'))

    def __init__(self, dataName='', geneListFileName=None, mitochondrialGenes=None,
                sigmaOverMeanSigma=0.3, nClusters=5, nComponentsPCA=200, nSamplesDistribution=10000, 
                saveDir=os.path.join(''), makeMarkerSubplots=False, availableCPUsCount=os.cpu_count(), zScoreCutoff=0.3,
                subclusteringName=None, doQualityControl=True, doBatchCorrection=True, makePlots=True,
                minimumNumberOfMarkersPerCelltype=10, minimumScoreForUnknown=0.3):

        '''Initialization function. Automatically called when an instance on Digital Cell Sorter is created

        Args:
            dataName: name used in output files, Default ''
            geneListFileName: , Default None
            mitochondrialGenes: , Default None
            sigmaOverMeanSigma: threshold to consider a gene constant, Default 0.3
            nClusters: number of clusters, Default 5
            nComponentsPCA: number of pca components, Default 100
            nSamplesDistribution: , Default 10000
            saveDir: directory for output files, Default is current directory
            makeMarkerSubplots:  whether to make subplots on markers, Default True
            makePlots: whether to make all major plots, Default True
            votingScheme: voting shceme to use instead of the built-in, Default None
            availableCPUsCount: number of CPUs available, Default os.cpu_count()
            zScoreCutoff: zscore cutoff when calculating Z_mc, Default 0.3
            minimumScoreForUnknown: zscore cutoff when assigning label "Unknown", Default 0.3
            clusterName: parameter used in subclustering, Default None
            doQualityControl: whether to remove low quality cells, Default True
            doBatchCorrection: whether to correct data for batches, Default False

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter(geneListFileName='MyGeneList.xlsx')
        '''

        self.saveDir = saveDir
        self.dataName = dataName

        self.geneListFileName = geneListFileName
        self.mitochondrialGenes = mitochondrialGenes

        self.sigmaOverMeanSigma = sigmaOverMeanSigma
        self.nClusters = nClusters
        self.nComponentsPCA = nComponentsPCA
        self.zScoreCutoff = zScoreCutoff
        self.minimumNumberOfMarkersPerCelltype = minimumNumberOfMarkersPerCelltype

        self.minimumScoreForUnknown = minimumScoreForUnknown

        self.nSamplesDistribution = nSamplesDistribution
        self.availableCPUsCount = availableCPUsCount

        self.subclusteringName = dataName if subclusteringName is None else subclusteringName

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

        self.toggleMakeMarkerSubplots = makeMarkerSubplots

        return None


    # Supporting functions of class ###########################################################################################################################
    def write(self, data, fileName):
    
        '''Pickle object into a (binary) file
        
        Args:
            data: any Pyhton object, e.g. list, dictionary, file, method, variable, etc.
            fileName: path and name of the file to store binary data in

        Returns:
            None
        
        Usage:
            data = [['A', 'B', 'C'], pd.DataFrame()]
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.write(data, os.path.join('some dir 1', 'some dir 2', 'File with my data'))
        '''

        with gzip.open(fileName + '.pklz','wb') as temp_file:
            pickle.dump(data, temp_file, protocol=4)

        return

    def read(self, fileName):

        '''Unpickle object from a (binary) file

        Args:
            fileName: path and name of the file with binary data stored in

        Returns:
            Data stored in the provided file
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            data = DCS.read(os.path.join('some dir 1', 'some dir 2', 'File with my data'))
        '''

        with gzip.open(fileName + '.pklz','rb') as temp_file:
            data = pickle.load(temp_file)
            return data

        return

    def convertColormap(self, colormap):

        '''Convert colormap from the form (1.,1.,1.,1.) to 'rgba(255,255,255,1.)'

        Args:
            colormap: a dictionary, colormap to convert

        Returns:
            Converted colomap
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            colormap = DCS.convertColormap(colormap)
        '''

        return {k:'rgba' + str(tuple(list((np.array(v[:3]) * 255).astype(int)) + [v[3]])) for k, v in colormap.items()}

    def getCountsDataframe(self, se1, se2, tagForMissing='Missing'):

        '''Get a pandas.DataFrame with cross-counts (overlaps) between two pandas.Series

        Args:
            se1: pandas.Series with the first set of items
            se2: pandas.Series with the second set of items
            tagForMissing: label to assign to non-overlapping items, Default 'Missing'

        Returns:
            pandas.DataFrame with counts
        
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

    def zScoreOfSeries(self, se):
    
        '''Calculate z-score of pandas.Series and modify the Series in place
    
        Args:
            se: pandas.Series

        Returns:
            pandas.Series 

        Usage:
            se = zScoreOfSeries(se)
        '''

        se.iloc[:] = scipy.stats.zscore(se.values)

        return se
        
    def getV(self, args):

        '''Calculate the vting scores (celltypes by clusters)
    
        Args:
            args: tupple of sub-arguments
                df_M: 
                df_markers_expr: 
                cluster_index: 
                zscore_cutoff: 

        Returns:
            pandas.DataFrame containing voting scores per celltype per cluster 

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

    def getIndexOfGoodQualityCells(self, saveDir, dataName, count_depth_cutoff=0.5, number_of_genes_cutoff=0.5, mitochondrial_genes_cutoff_upper_bound=3.0):

        '''Get index of sells that satisfy the QC criteria

        Args:
            saveDir: directory for output files
            dataName: name used in output files
            count_depth_cutoff: fraction of median to take as count depth cutoff, Default 0.5
            number_of_genes_cutoff: fraction of median to take as number of genes cutoff, Default 0.5
            mitochondrial_genes_cutoff: the cutoff is median + standard_deviation * this_parameter, Default 3.0

        Returns:
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

        Args:
            se: pandas.Series with data to analyze
            cutoff: parameter for calculating the quality control cutoff
            mito: whether the analysis of mitochondrial genes fraction, Default False
            plotPathAndName: text to include in the figure title and file name
            MakeHistogramPlot: whether to make a histogram plot, Default True

        Returns:
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

    def scoreCellsByIsolationForest(self, trainingSet, testingSet, printResults=False):

        '''Function to get anomaly score of cells based on some reference

        Args:
            trainingSet: pandas.DataFrame with cells to trail isolation forest on
            testingSet: pandas.DataFrame with cells to score
            printResults: whether to print results

        Returns:
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

    def getExprOfCells(self, cells):

        '''Get expression of a set of cells.
        Run this only after function process()

        Args:
            cells: index, a set of cells of interest

        Returns:
            Pandas DataFrame with expression of the cells of interest

        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.process()
            DCS.getDfExprOfSelectedCells(cells)
        '''

        df_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r')
        df_expr.index.names = ['patient', 'cell', 'gene']
        df_expr = df_expr.unstack(level='gene', fill_value=0).T
        df_expr.index = df_expr.index.get_level_values('gene')

        return df_expr[cells]

    def getAnomalyScoresPlot(self, cells='All'):

        '''Make anomaly scores plot

        Args:
            cells: cells of interest, Default 'All'

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

        scores = self.scoreCellsByIsolationForest(df_sel, df_sel)
        print('Loaded expression data')
        timeMark()

        scores = pd.DataFrame(index=cells, data=scores).reindex(df_tsne.columns).values.T[0]

        self.makeTSNEplot(df_tsne.values, scores, self.dataName, self.saveDir, suffix='by_anomaly_score', legend=False, labels=False, colorbar=True)

        return None
    
    def getIndividualGeneExpressionPlot(self, gene, hideClusterLabels=False, outlineClusters=True):

        '''Produce individual gene expression plot on a tSNE layout

        Args:
            gene: gene of interest

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.getIndividualGeneExpressionPlot('CD4')
        '''

        df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_expr', mode='r').xs(key=gene, axis=0, level=2).T
        df_markers_expr.index = [gene]

        df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
        df_markers_expr = df_markers_expr.reindex(df_tsne.columns, axis=1).fillna(0.)

        labelled = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r').columns
        df_markers_expr.columns = labelled.to_series().reset_index().set_index(['patient', 'cell'])['cluster'].loc[df_markers_expr.columns].reset_index().set_index(['patient', 'cell', 'cluster']).index

        hugo_cd_dict = {gene: self.gnc.Convert([gene], 'hugo', 'alias', returnUnknownString=False)[0]}
        self.makeMarkerSubplots(df_markers_expr, df_tsne.values, hugo_cd_dict, self.dataName, self.saveDir, NoFrameOnFigures=True, HideClusterLabels=hideClusterLabels, outlineClusters=outlineClusters)

        return None


    # Vizualization functions of class ########################################################################################################################
    def makeMarkerExpressionPlot(self, dataName, saveDir):

        '''Produce image on marker genes and their expression on all clusters.
        Uses files generated by function DCS.Vote

        Args:
            dataName: name used in output files
            saveDir: directory for output files

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.MakeMarkerExpressionPlot(dataName, saveDir)
        '''

        df_votingResults = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='z-scores')
        votingResults = dict(zip(df_votingResults['cluster'].values, df_votingResults['Predicted cell type'].values))
        supportingMarkersList = dict(zip(df_votingResults['cluster'].values, df_votingResults['Supporting markers'].str.split(' // ')))
        allMarkersList = dict(zip(df_votingResults['cluster'].values, df_votingResults['All markers'].str.split(' // ')))
        df_markers_cluster_centroids = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='Cluster centroids', index_col=0, header=0).T

        # Y_mc.T
        X_markers_cluster_means_transpose = df_markers_cluster_centroids.values.T

        # Normalization
        for i in range(X_markers_cluster_means_transpose.shape[1]):
            X_markers_cluster_means_transpose[:,i] -= np.min(X_markers_cluster_means_transpose[:,i])
            X_markers_cluster_means_transpose[:,i] /= np.max(X_markers_cluster_means_transpose[:,i])

        ORDER = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(X_markers_cluster_means_transpose, 'ward'), no_plot=True, get_leaves=True)['leaves']
        ORDER2 = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(X_markers_cluster_means_transpose.T, 'ward'), no_plot=True, get_leaves=True)['leaves']

        X_markers_cluster_means_sorted = X_markers_cluster_means_transpose[ORDER,:][:,ORDER2]

        df_all_marker_hits = pd.DataFrame(data=np.zeros((df_markers_cluster_centroids.shape)), index=df_markers_cluster_centroids.index, columns=df_markers_cluster_centroids.columns)
        for cluster in allMarkersList:
            if not allMarkersList[cluster] is np.nan:
                for gene in allMarkersList[cluster]:
                    df_all_marker_hits.loc[gene, cluster] = 1

        df_supp_marker_hits = pd.DataFrame(data=np.zeros((df_markers_cluster_centroids.shape)), index=df_markers_cluster_centroids.index, columns=df_markers_cluster_centroids.columns)
        for cluster in supportingMarkersList:
            if not supportingMarkersList[cluster] is np.nan:
                for gene in supportingMarkersList[cluster]:
                    df_supp_marker_hits.loc[gene, cluster] = 1

        X_marker_hits = df_all_marker_hits.values.T[ORDER,:][:,ORDER2]
        X_supp_marker_hits = df_supp_marker_hits.values.T[ORDER,:][:,ORDER2]

        _figsize = np.float_(X_markers_cluster_means_transpose.shape[::-1]) / \
                    np.max(X_markers_cluster_means_transpose.shape) * 15.0 + 2.0

        _figsize[1] *= 1.5

        gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[20,1], left=0.13, right=0.99, top=0.99, bottom=0.25, wspace=0.01, hspace=0.01)

        fig = plt.figure(figsize=_figsize)
        ax = plt.subplot(gs[0])

        ax.imshow(X_markers_cluster_means_sorted,cmap='Blues',interpolation='None', aspect='auto')

        i_list,j_list = np.where(X_marker_hits.T > 0)
        color = 'w' #(1., 1., 0.7)
        ax.plot(i_list, j_list, 'k*', mec=color, mew=0.7, markersize=3)

        i_list_supp,j_list_supp = np.where(X_supp_marker_hits.T > 0)
        ax.plot(i_list_supp, j_list_supp, 'k*', mec='r', mew=0.7, markersize=3) #mec='k', alpha=0.5, markersize=6

        ax.set_xticks(range(X_markers_cluster_means_transpose.shape[1]))
        ax.set_yticks(range(X_markers_cluster_means_transpose.shape[0]))

        xtickslabels = np.array(df_markers_cluster_centroids.index[ORDER2])
        for i in range(0,len(xtickslabels),2):
            xtickslabels[i] += " ─────────"

        ax.set_xticklabels(xtickslabels,rotation=90, fontsize=5)

        ax.set_yticklabels([list(votingResults.values())[i] + ' (' + str(i) + ')' for i in ORDER], rotation=0, fontsize=8)
        #ax.set_yticklabels(['Cluster '+str(i) for i in ORDER], rotation=45,
        #fontsize=10)

        ax.set_xlim([-0.5,X_markers_cluster_means_transpose.shape[1] - 0.5])
        ax.set_ylim([-0.5,X_markers_cluster_means_transpose.shape[0] - 0.5])



        axx = plt.subplot(gs[1])
        fontsize = 5
        cells_in_clusters = df_votingResults['# cells in cluster'].values.copy()[ORDER]
        numberOfCells = cells_in_clusters.sum()

        with open(os.path.join(saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        celltypes = df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0].values.copy()[ORDER]

        axx.barh(y=range(len(cells_in_clusters)), width=cells_in_clusters, height=0.8, align='center', color=[colormap[i] for i in celltypes])
        for i in range(len(cells_in_clusters)):
            axx.text(np.max(cells_in_clusters), i, cells_in_clusters[i], ha='right',va='top', color='k', weight='bold', fontsize=fontsize)
            axx.text(0.02 * numberOfCells, i, str(round(100 * cells_in_clusters[i] / numberOfCells, 1)) + '%', ha='left',va='bottom', color='b', fontsize=fontsize)
        axx.set_xticklabels(cells_in_clusters, fontsize=fontsize)
        axx.set_yticklabels(cells_in_clusters, alpha=0)
        axx.set_xticklabels(cells_in_clusters, alpha=0)
        axx.set_xticks([])
        axx.set_yticks([])
        axx.set_xlabel('Number of\ncells in clusters', fontsize=fontsize)
        axx.set_ylim(-0.5, len(cells_in_clusters) - 0.5)

        if saveDir is not None: 
            fig.savefig(os.path.join(saveDir, dataName + '_voting.png'), dpi=600)

        return

    def makeMarkerSubplots(self, df, X_tsne, hugo_cd_dict, dataName, saveDir, NoFrameOnFigures=False, HideClusterLabels=False, outlineClusters=True):

        '''Produce subplots on each marker and its expression on all clusters

        Args:
            df: pandas.DataFrame with marker genes expression
            X_tsne: tSNE coordinates for each cell
            hugo_cd_dict: dictionary with aliases for hugo names of genes
            dataName: name used in output files
            saveDir: directory for output files
            NoFrameOnFigures: whether to include frame on the figure
            HideClusterLabels: whether to print cluster labels on the figure

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.MakeMarkerSubplot(df_markers_expr, tSNE, hugo_cd_dict, dataName, saveDir)
        '''
        
        def MarkerSubplot(counter, marker, df, votingResults, X_tsne, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, NoFrameOnFigures, HideClusterLabels, XLIM, YLIM, directory, circles):

            fig,ax = plt.subplots(figsize=(8,8))

            ax.cla()
            suffix = '(' + str(hugo_cd_dict[marker]).replace('\'', '').replace('(', '').replace(')', '').replace(' ','') + ')'
            ax.plot(np.nan,np.nan,'*',markersize=15,c=cm.seismic(1.0),label=marker + '\n' + suffix.replace(',','\n'))
            circleIndices = np.where(df.loc[marker].values == 0)[0] # cells that don't have this marker
            starIndices = np.where(df.loc[marker].values > 0)[0] # cells that have this marker
            starIndices = starIndices[np.argsort(df.loc[marker].values[starIndices])]
            args1 = [X_tsne[0,circleIndices],
                        X_tsne[1,circleIndices]]
            kwargs1 = {'marker':'o',
                        'c':'b',
                        'alpha':0.1,
                        's':6 * 3,
                        'linewidth':0,}
            args2 = [X_tsne[0,starIndices],
                        X_tsne[1,starIndices]]
            kwargs2 = {'marker':'*',
                        'c':cm.seismic(df.loc[marker].values[starIndices] / np.max(df.loc[marker].values[starIndices])),
                        's':30 * 4,
                        'linewidth':0.0,}
            ax.scatter(*args1,**kwargs1)
            ax.scatter(*args2,**kwargs2)
            for label in set(cellClusterIndexLabel):
                # cells with this label
                X_tsne2_cluster = X_tsne[:,cellClusterIndexLabel == label]
                x_mean = np.mean(X_tsne2_cluster[0,:])
                y_mean = np.mean(X_tsne2_cluster[1,:])

                _text_label = votingResults[label] if not HideClusterLabels else ''

                ax.text(x_mean,y_mean,
                        _text_label.
                            replace('-','-\n').replace(' ','\n').
                            replace('T\n','T ').replace('B\n','B ').
                            replace('\n#',' #').replace('/','/\n').
                            replace('NK\n','NK ').replace('Stem\n','Stem '),
                        fontsize=10,
                        ha='center',va='center',#alpha=0.75,
                        ).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])

                if circles:
                    radius = np.sqrt(X_tsne2_cluster.shape[1]) * 300.0
                    ax.scatter(x_mean,y_mean,s=radius * 1,facecolors='none',edgecolors='k')
            ax.set_xlim(XLIM)
            ax.set_ylim(YLIM)
            ax.legend(loc='upper right', frameon=False, fontsize=14) #loc='best',numpoints=1,fontsize=12
            ax.set_xticks([])
            ax.set_yticks([]) 
            if NoFrameOnFigures:
                #fig.patch.set_visible(False)
                ax.axis('off')
            fig.tight_layout()
            if saveDir is not None: 
                #fig.savefig('%s/marker_subplots/%s_%s_%s.png' %
                #(saveDir,dataName,marker,suffix.replace(',','_').replace('/','_')),dpi=150)
                #fig.savefig(os.path.join(saveDir, 'marker_subplots', '%s_%s_%s.pdf' % (dataName,marker,suffix.replace(',','_').replace('/','_'))))
                fig.savefig(os.path.join(saveDir, 'marker_subplots', '%s_%s_%s.png' % (dataName,marker,suffix.replace(',','_').replace('/','_'))), dpi=300)
                print(marker, end=", ", flush=True)

            plt.close(fig)

            return
        
        maxs = np.max(X_tsne,axis=1)
        mins = np.min(X_tsne,axis=1)
        maxDiffs = maxs - mins
        deltas = maxDiffs * 0.05
        XLIM = [mins[0] - deltas[0],maxs[0] + deltas[0]]
        YLIM = [mins[1] - deltas[1],maxs[1] + deltas[1]]

        directory = os.path.join(saveDir, 'marker_subplots', '')
        if not os.path.exists(directory):
            os.makedirs(directory)

        print('\nSaving marker subplots:\n')

        df_votingResults = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='z-scores')
        votingResults = dict(zip(df_votingResults['cluster'].values, df_votingResults['Predicted cell type'].values))

        cellClusterIndexLabel = df.columns.get_level_values('cluster').values

        for counter,marker in enumerate(df.index.values):
            MarkerSubplot(counter, marker, pd.DataFrame(data=np.reshape(np.array(df.loc[marker]), (1,len(df.loc[marker]))), columns=df.columns, index=[marker]), 
                            votingResults, X_tsne, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, NoFrameOnFigures, HideClusterLabels, XLIM, YLIM, directory, outlineClusters)

        print('\nDone saving marker subplots!')

        return

    def makeVotingResultsMatrixPlot(self, saveDir, dataName):

        '''Produce voting results voting matrix plot

        Args:
            dataName: name used in output files
            saveDir: directory for output files

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.MakeVotingResultsMatrixPlot(dataName, saveDir)
        '''

        df_votingResults = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='z-scores')

        cellTypes = sorted([x for x in df_votingResults.columns.values.tolist() if x not in ['cluster', 'Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers']])

        #df_votingResults['order'] =
        #scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_votingResults[cellTypes].values,
        #method='ward', metric='euclidean', optimal_ordering=True),
        #                                                               no_plot=True,get_leaves=True)['leaves']

        df_votingResults['order'] = np.argsort(np.argsort(df_votingResults['Predicted cell type']))
        df_votingResults = df_votingResults.sort_values(by='order', axis=0, ascending=False)

        numberOfCells = np.sum(df_votingResults['# cells in cluster'])
        num_of_cell_types = len(cellTypes)
        num_of_clusters = np.unique(df_votingResults['cluster']).shape[0]

        indicis_of_clusters = df_votingResults['cluster']
        assigned_names_of_clusters = df_votingResults['Predicted cell type']
        
        label_max = 0.27 * max([len(assigned_names_of_clusters[i]) for i in range(len(assigned_names_of_clusters))])

        _figsize = np.float_((num_of_cell_types,num_of_clusters)) / np.max(num_of_clusters) * 15.0
        _figsize[0] += 1.0 + label_max
        _figsize[1] += 2.0

        gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[4,1])

        fig = plt.figure(figsize=_figsize)
        ax = plt.subplot(gs[0])
        zscores = df_votingResults[cellTypes].values
        ax.imshow(zscores, aspect=1, cmap='Greens', vmin=0, vmax=0.5, interpolation='None')

        for i in range(num_of_clusters):
            for j in range(num_of_cell_types):
                if np.round(zscores[i,j],1) > 0:
                    if zscores[i,j] == np.max(zscores[i,:]):
                        ax.text(j,i,np.round(zscores[i,j],1), color='w', fontsize=70 * min([(num_of_cell_types / num_of_clusters),1]),ha='center',va='center').set_path_effects([path_effects.Stroke(linewidth=4, foreground='red'),path_effects.Normal()])
                    else:
                        ax.text(j,i,np.round(zscores[i,j],1), color='w', fontsize=60 * min([(num_of_cell_types / num_of_clusters),1]),ha='center',va='center').set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()])

        ax.set_xticks(range(num_of_cell_types))
        ax.set_yticks(range(num_of_clusters))

        ytickslabels = copy.deepcopy(assigned_names_of_clusters)
        for i in range(len(ytickslabels)):
            ytickslabels[i] = str(assigned_names_of_clusters[i]) + ' (' + str(indicis_of_clusters[i]) + ')'

        xtickslabels = np.array(copy.deepcopy(cellTypes))
        for i in range(len(xtickslabels)):
            if i % 3 == 1: xtickslabels[i] = '\n' + xtickslabels[i]
            if i % 3 == 2: xtickslabels[i] = '\n\n' + xtickslabels[i]

        ax.set_xticklabels(xtickslabels, rotation=0, fontsize=20, ha='center')
        ax.set_yticklabels(ytickslabels, fontsize=20, rotation=0) 
        ax.set_xlim([-0.5,num_of_cell_types - 0.5])
        ax.set_ylim([-0.5,num_of_clusters - 0.5])

        axx = plt.subplot(gs[1])
        cells_in_clusters = df_votingResults['# cells in cluster'].values

        with open(os.path.join(saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        celltypes = df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0].values.copy()

        axx.barh(y=range(len(cells_in_clusters)), width=cells_in_clusters, height=0.8, align='center', color=[colormap[i] for i in celltypes])
        for i in range(len(cells_in_clusters)):
            axx.text(np.max(cells_in_clusters), i, cells_in_clusters[i], ha='right',va='top', color='k', weight='bold', fontsize = 20)
            axx.text(0.02 * numberOfCells, i, str(round(100 * cells_in_clusters[i] / numberOfCells, 1)) + '%', ha='left',va='bottom', color='b', fontsize = 20)
        axx.set_xticklabels(cells_in_clusters, fontsize=20)
        axx.set_yticklabels(cells_in_clusters, alpha=0)
        axx.set_xticklabels(cells_in_clusters, alpha=0)
        axx.set_xlabel('Number of\ncells in clusters', fontsize=20)
        axx.set_ylim(-0.5, num_of_clusters - 0.5)
        
        fig.tight_layout()
        fig.savefig(os.path.join(saveDir, dataName + '_matrix_voting.png'), dpi=150)
 
        return 

    def makeHistogramNullDistributionPlot(self, saveDir, dataName):

        '''Produce histogram plot of the voting null distributions

        Args:
            dataName: name used in output files
            saveDir: directory for output files

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.MakeHistogramNullDistributionPlot(dataName, saveDir)
        '''

        df_noise_dict = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='Null distributions', index_col=0, header=[0,1], skiprows=[2])
        df_votingResultsV = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='Voting scores')
        df_votingResultsZ = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'), sheet_name='z-scores')

        predicted_cell_type_cluster = df_votingResultsZ['Predicted cell type'].values
        predicted_cell_type = df_votingResultsZ['Predicted cell type'].str.split(' #', expand=True)[0].values

        cellTypes = sorted([x for x in df_votingResultsV.columns.values.tolist() if x not in ['cluster', 'Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers']])

        df_votingResultsV = df_votingResultsV[cellTypes]
        df_votingResultsZ = df_votingResultsZ[cellTypes]

        cell_types = np.unique(df_noise_dict.columns.get_level_values(0).values)
        num_of_cell_types = cell_types.shape[0]
        clusters = np.unique(df_noise_dict.columns.get_level_values(1).values.astype(float)).astype(str)
        num_of_clusters = clusters.shape[0]

        try:
            maxx = np.round(df_noise_dict.index.values[np.where(df_noise_dict.values.max(axis=1) > 0)[0][-1]] / 100., 1) # + 0.1
            maxy = np.round(df_noise_dict.values.max(), 1) # + 0.02
        except:
            maxx = maxy = 1.0

        with open(os.path.join(saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}

        origWidth = matplotlib.rcParams['axes.linewidth']
        matplotlib.rcParams['axes.linewidth'] = 0.1

        gs = matplotlib.gridspec.GridSpec(num_of_clusters, num_of_cell_types, hspace=0.45, wspace=0.1, bottom=0.04, top=0.96, left=0.05, right=0.99)
        
        fig = plt.figure(figsize=(num_of_cell_types, num_of_clusters * 0.4))

        for i in range(num_of_cell_types):
            for j in range(num_of_clusters):

                fontsize = 3

                ax = plt.subplot(gs[i + num_of_cell_types * j])

                ax.bar(df_noise_dict.index.values / 100., df_noise_dict[(cell_types[i], eval(clusters[j]))].values, 
                       width=0.007, align='center', color=colormap[cell_types[i]])

                valueV = df_votingResultsV.loc[eval(clusters[j]), cell_types[i]]
                valueZ = df_votingResultsZ.loc[eval(clusters[j]), cell_types[i]]
                ax.axvline(x=valueV, ymin=0, ymax=1, color='k', lw=0.2)
                color = 'k' if predicted_cell_type[j] != cell_types[i] else 'r'
                ax.text(valueV + 0.02, maxy - 0.005, r'$V_{%s,%s}=$' % (i,j) + str(np.round(valueV,2)), fontsize=fontsize, va='top', ha='left', color=color)
                ax.text(valueV + 0.02, maxy - 0.045, r'$\Lambda_{%s,%s}=$' % (i,j) + str(np.round(valueZ,2)), fontsize=fontsize, va='top', ha='left', color=color)

                if j == 0:
                    ax.set_title(cell_types[i], fontdict={'color': 'b', 'size':'6'})

                if i == 0:
                    ax.text(-0.27 * maxx, maxy + 0.01, predicted_cell_type_cluster[j] + ' (Cluster %s)' % j, rotation=0, 
                            fontsize=fontsize, weight='bold', va='bottom', ha='left', color='k')
                    
                    ax.set_ylabel('Probability', fontsize=fontsize)
                    ax.set_yticklabels(np.array(range(int(20 * maxy) + 1)) / 20, fontsize=fontsize)
                else:
                    ax.set_yticklabels([], fontsize=fontsize)

                if j == num_of_clusters - 1:
                    ax.set_xlabel('Voting score', fontsize=fontsize)
                    ax.set_xticklabels(np.array(range(int(5 * maxx) + 1)) / 5, fontsize=fontsize)
                else:
                    ax.set_xticklabels([], fontsize=fontsize)

                ax.set_xticks(np.array(range(int(5 * maxx) + 1)) / 5)
                ax.set_yticks(np.array(range(int(20 * maxy) + 1)) / 20)

                ax.set_xlim(0.0, maxx)
                ax.set_ylim(0.0, maxy)

                ax.tick_params(direction='in', length=1, width=0.1, colors='k')
        
        fig.savefig(os.path.join(saveDir, dataName + '_null_distributions.png'), dpi=1200)

        matplotlib.rcParams['axes.linewidth'] = origWidth

        return

    def makeTSNEplot(self, Xtsne, cellClusterIndexLabel, dataName, saveDir, suffix, colormap=cm.jet, legend=True, labels=True, colorbar=False, fontsize=10, plotNaNs=True):

        '''Produce tSNE plot with a specified coloring scheme

        Args:
            X_tsne2: tSNE coordinates for each cell
            cellClusterIndexLabel: cluster index for each cell
            dataName: name used in output files
            saveDir: directory for output files
            suffix: text label to append to the figure name
            colormap: cell coloring sequence, can be a dictionary or cm.colormap, 
                Default matplotlib.colors.LinearSegmentedColormap.jet
            legend: whether to print legend, Default True
            labels: whether to print labels, Default True
            fontsize: labels and legend font size, Default 10
            plotNaNs: whether to plot NaN labels (in grey), Default True

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.MakeTSNEplot(tSNE, cellClusterIndexLabel, dataName, saveDir, suffix)
        '''

        def add_colorbar(fig, labels, cmap=matplotlib.colors.LinearSegmentedColormap.from_list('GR', [(0, 1, 0), (1, 0, 0)], N=100), fontsize=10):
    
            mapp=cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=np.min(labels), vmax=np.max(labels)), cmap=cmap)
            mapp.set_array(labels)
            sp = np.linspace(np.max(labels), np.min(labels), num=6, endpoint=True)

            axisColor = fig.add_axes([0.9,0.5,0.01,0.4])
            fig.colorbar(mapp, cax=axisColor, ticks=sp)

            axisColor.tick_params(labelsize=fontsize)
            axisColor.set_yticklabels(np.round(sp,2))
    
            return None

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.05,0.05,0.9,0.9])

        maxs, mins = np.max(Xtsne,axis=1), np.min(Xtsne,axis=1)

        missing = np.where(cellClusterIndexLabel!=cellClusterIndexLabel)[0]
        if len(missing)>0:
            ax.plot(Xtsne[0, missing], Xtsne[1, missing], 'o', color='grey', mew=0.5, alpha=0.2, markeredgecolor='k', label='NaN')

        nonMissing = np.where(cellClusterIndexLabel==cellClusterIndexLabel)[0]
        cellClusterIndexLabel = np.array(cellClusterIndexLabel)[nonMissing]
        Xtsne = Xtsne[:, nonMissing]

        possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel))

        if labels:
            print(possible_cluster_labels)

        for ilabel, label in enumerate(possible_cluster_labels):
            color = colormap(ilabel / len(possible_cluster_labels)) if type(colormap) is matplotlib.colors.LinearSegmentedColormap else colormap[label.split(' #')[0]]

            XtsneC = Xtsne[:,cellClusterIndexLabel == label]

            ax.plot(XtsneC[0,:], XtsneC[1,:], 'o', color=color, mew=0.5, alpha=0.3, markeredgecolor='k', label=label)

            if labels:
                ax.text(np.mean(XtsneC[0,:]), np.mean(XtsneC[1,:]), 
                    label, fontsize=fontsize, ha='center',va='center').set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])

        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.set_xlim([mins[0] - (maxs[0] - mins[0]) * 0.05, maxs[0] + (maxs[0] - mins[0]) * 0.05])
        ax.set_ylim([mins[1] - (maxs[1] - mins[1]) * 0.05, maxs[1] + (maxs[1] - mins[1]) * 0.05])

        if legend:
            plt.legend(loc='lower right', frameon=False, fontsize=fontsize)

        #fig.patch.set_visible(False)
        ax.axis('off')

        if colorbar:
            add_colorbar(fig, possible_cluster_labels, cmap=colormap, fontsize=fontsize)

        if saveDir is not None: 
            fig.savefig(os.path.join(saveDir, '%s_clusters_%s.png' % (dataName, suffix)), dpi=300)

        return

    def makeStackedBarplot(self, saveDir, dataName, clusterName):
        
        '''Produce stacked barplot on subclustering

        Args:
            saveDir: directory for output files
            dataName: name used in output files
            clusterName: label to include at the bar bottom

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.make_stacked_bar_plot_subclustering(saveDir, dataName, clusterName)
        '''

        def get_stacked_data_and_colors(saveDir):
            with open(os.path.join(saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
                colors = temp_file.readlines()
                colors = np.vstack([(color.strip('\n').split('\t')) for color in colors])
                colors = pd.DataFrame(colors.T[1], index=colors.T[0]).apply(lambda x: tuple(np.float_(x[0][1:][:-1].split(','))), axis=1)

            df = pd.read_excel(os.path.join(saveDir, dataName + '_voting.xlsx'))
            index = df['Predicted cell type']

            if not clusterName is None:
                barName = dataName + ': ' + clusterName
            else:
                barName = dataName

            index = [index[i][:len(index[i]) if index[i].find('#') - 1 == -2 else index[i].find('#') - 1].strip('*').strip('#').strip(' ') for i in range(len(index))]
            df_BM_temp = pd.DataFrame(data=df['# cells in cluster'].values, index=index, columns=[barName])
        
            df_main = pd.DataFrame(data=np.zeros((len(colors),1)), index=colors.index, columns=[barName])

            for i, item in enumerate(df_BM_temp.index):
                df_main.loc[item,barName] += df_BM_temp.iloc[i][barName] if df_BM_temp.index[i] == item else 0

            s = 'sums'
            df_main[s] = np.array(np.sum(df_main, axis=1))
            df_main.loc['Unknown',s] = 0
            df_main = df_main.apply(lambda x: 100. * x / np.sum(df_main, axis=0), axis=1).loc[np.sum(df_main, axis=1) > 0].sort_values(by=[s]).drop(columns=[s])

            return df_main, colors, clusterName

        df_Main, colors, clusterName = get_stacked_data_and_colors(saveDir)

        saveName = os.path.join(saveDir, "%s_subclustering_stacked_barplot_%s.png" % (dataName, ('All cell clusters' if clusterName == None else clusterName).replace(' ', '_').replace('*', '')))

        fig,ax = plt.subplots(figsize=(4.5,8)) #4.15

        barWidth = 0.85
        cellTypes = df_Main.index
        bottom = np.zeros((len(df_Main.columns)))

        for i in range(len(cellTypes)):
            bottom += df_Main.loc[cellTypes[i - 1]] if i > 0 else 0
            ax.bar(range(len(df_Main.columns)), list(df_Main.loc[cellTypes[i]]), bottom=list(bottom), color=colors.loc[cellTypes[i]], edgecolor='white', width=barWidth, label=cellTypes[i])
 
        plt.xticks(range(len(df_Main.columns)), list(df_Main.columns), fontsize=12)
        plt.yticks([0,20,40,60,80,100], ['0','20%','40%','60%','80%','100%'], fontsize=12)
            
        handles, labels = ax.get_legend_handles_labels()
        ms = np.max([len(item) for item in labels]) - len('cell')
        labels = [item.replace(' ','\n').replace('B\n', 'B ').replace('T\n', 'T ') if len(item) >= ms else item for item in labels[::-1]]
        ax.legend(handles[::-1], labels, loc='upper left', bbox_to_anchor=(1,1), ncol=1, frameon=False, fontsize=14, labelspacing=1, title = ''.join([' ' for _ in range(60)]))

        plt.xlim((-0.5, len(df_Main.columns) - 0.5))
        plt.ylim((0, 100))

        for spine in plt.gca().spines.values():
            spine.set_visible(False)
          
        fig.tight_layout()
        fig.savefig(saveName, dpi=300)

        print('\n=================================\nDone saving stacked bar plot: %s!\n=================================' % ('All cell clusters' if clusterName == None else clusterName))

        return
    
    def makeQualityControlHistogramPlot(self, subset, cutoff, plotPathAndName=None, N_bins=100, mito=False, displayMeasures=True, precision=4):

        '''Function to calculate QC quality cutoff and visualize it on a histogram

        Args:
            subset: data to analyze
            cutoff: cutoff to display
            plotPathAndName: text to include in the figure title and file name
            N_bins: number of bins of the histogram
            mito: whether the analysis of mitochondrial genes fraction, Default False
            displayMeasures: print vertical dashed lines along with mean, median, and standard deviation

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            cutoff = DCS.makeQualityControlHistogramPlot(subset, cutoff)
        '''

        if plotPathAndName is None:
            plotPathAndName = 'QC_Plot'

        range_min = np.min(subset)
        range_max = max(1.1*cutoff, 0.2) if mito else 5000

        hist_of_subset = scipy.stats.rv_histogram(np.histogram(subset, bins=N_bins, range=(range_min, range_max)))
        hist_data = hist_of_subset._hpdf / N_bins
        hist_bins = hist_of_subset._hbins

        fig, ax = plt.subplots(figsize=(8,8))

        bar_bin_width = range_max / N_bins
        ax.bar(hist_bins, hist_data[:-1], width=0.9 * bar_bin_width, color='b', align='center')

        ax.set_title(plotPathAndName, fontdict={'color': 'b'})
        ax.set_xlabel('Fraction' if mito else 'Count', fontsize=8)
        ax.set_ylabel('Density', fontsize=8)
        ax.set_ylim(0.,ax.get_ylim()[1])
        ax.set_xlim(range_min - 0.5 * bar_bin_width, range_max + 0.5 * bar_bin_width)

        xs = np.linspace(hist_bins[0], hist_bins[-1], 1000)
        spline_data = np.vstack((xs, UnivariateSpline(hist_bins, hist_data[:-1], k=5, s=0)(xs))).T

        sg = scipy.signal.savgol_filter(spline_data.T[1], 101, 3)
        ax.plot(spline_data.T[0], sg, 'r', lw=3, alpha=0.95)

        x, y = cutoff, sg[np.where(spline_data.T[0] >= cutoff)[0][0]]
        ax.plot([x,x], [0,y], 'k', lw=2)
        ax.plot(x, y, 'ko', ms=10, alpha=0.8)
        ax.plot(x, y, 'ro', ms=7)
        ax.text(x, -0.04 * spline_data.T[1].max(), str(np.round(cutoff, precision)), va='top', ha='center', color='r')

        fig.tight_layout()

        if displayMeasures:
            dist_std, dist_median, dist_mean = np.round(np.std(subset),precision), np.round(np.median(subset),precision), np.round(np.mean(subset),precision)
            print(plotPathAndName, '\tstd:', dist_std,  '\tmedian:', dist_median,  '\tmean:', dist_mean)

            xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
            yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
        
            ax.axvline(x=dist_mean, color='k', lw=1.0, ls='--')
            ax.text(dist_mean + 0.02 * xspan, 0.98 * yspan, r'$\mu=%s$' % (dist_mean), fontsize=10, va='top', ha='left', color='k')
        
            ax.axvline(x=dist_median, color='k', lw=1.0, ls='--')
            ax.text(dist_median + 0.02 * xspan, 0.94 * yspan, r'$M=%s$' % (dist_median), fontsize=10, va='top', ha='left', color='k')
        
            ax.axvline(x=dist_median - dist_std, color='k', lw=1.0, ls='--')
            ax.text(dist_median - dist_std + 0.02 * xspan, 0.90 * yspan, r'$M-\sigma=%s$' % (np.round(dist_median - dist_std,precision)), fontsize=10, va='top', ha='left', color='k')
        
            ax.axvline(x=dist_median + dist_std, color='k', lw=1.0, ls='--')
            ax.text(dist_median + dist_std + 0.02 * xspan, 0.90 * yspan, r'$M+\sigma=%s$' % (np.round(dist_median + dist_std,precision)), fontsize=10, va='top', ha='left', color='k')

            ax.annotate(r'$\sigma=%s$' % (dist_std), (dist_median + dist_std, 0.86 * yspan), (dist_median, 0.86 * yspan), arrowprops={'width':0, 'headwidth':0, 'headlength':0})
            ax.annotate('', (dist_median + dist_std, 0.85 * yspan), (dist_median, 0.85 * yspan), arrowprops={'arrowstyle':'<|-|>'})

        fig.savefig(plotPathAndName + '_histogram.png', dpi=300)

        return None

    def makeSankeyDiagram(self, df, dataName, colormapForIndex=None, colormapForColumns=None, linksColor='rgba(100,100,100,0.6)', title='', interactive=False, quality=4):

        '''Make a Sankey diagram, also known as 'river plot' with two groups of nodes

        Args:
            df: pandas.DataFrame with counts (overlaps)
            colormapForIndex: colors to use for nodes specified in the DataFrame index
            colormapForColumns: colors to use for nodes specified in the DataFrame columns
            linksColor: color of the non-overlapping links
            title: title to print on the diagram
            interactive: whether to launch interactive JavaScript-based graph
            quality: proportional to the resolution of the figure to save

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.sankeyDiagram(df, dataName)
        '''

        try:
            temp_index = pd.MultiIndex.from_arrays([df.index, [colormapForIndex[item] for item in df.index]], names=['label', 'color'])
            temp_columns = pd.MultiIndex.from_arrays([df.columns, [colormapForColumns[item] for item in df.columns]], names=['label', 'color'])
            df.index = temp_index
            df.columns = temp_columns
        except:
            print('Colormap error. Using default node colors')
            colormapForIndex = None
            colormapForColumns = None

        if (colormapForIndex is None) or (colormapForColumns is None):
            nodeColors = None
            nodeLabels = df.index.to_list() + df.columns.to_list()
        else:
            nodeLabels = df.index.get_level_values('label').to_list() + df.columns.get_level_values('label').to_list()
            nodeColors = df.index.get_level_values('color').to_list() + df.columns.get_level_values('color').to_list()

        sources, targets, values, labels = [], [], [], []
        for i, item in enumerate(df.index):
            sources.extend([i] * len(df.loc[item]))
            targets.extend(list(range(len(df.index), len(df.index) + len(df.loc[item]))))
            values.extend([j for j in df.loc[item].values])
            labels.extend([item[0] + ' -> ' + jtem[0] for jtem in df.loc[item].index])

        colorscales = [dict(label=label, colorscale=[[0, linksColor], [1, linksColor]]) for label in labels]

        if not nodeColors is None:
            for i in range(len(sources)):
                if nodeColors[sources[i]] == nodeColors[targets[i]]:
                    newColor = ','.join(nodeColors[sources[i]].split(',')[:3] + ['0.6)'])
                    colorscales[i] = dict(label=labels[i], colorscale=[[0, newColor], [1, newColor]])

        fig = plotly.graph_objects.Figure(data=[plotly.graph_objects.Sankey(valueformat = '', valuesuffix = '',
            node = dict(pad = 20, thickness = 40, line = dict(color = 'white', width = 0.5), label = nodeLabels, color = nodeColors,),
            link = dict(source = sources, target = targets, value = values, label = labels, colorscales = colorscales))]) #line ={'color':'rgba(255,0,0,0.8)', 'width':0.1}

        if not title is None:
            fig.update_layout(title_text=title, font_size=10)

        fig.write_image(os.path.join(dataName + '.png'), width=600, height=600, scale=quality)

        if interactive:
            fig.show()

        return

    def makePlotOfNewMarkers(self, saveDir, dataName, df_marker_cell_type, df_new_marker_cell_type):

        '''Produce plot of the new markers extracted from the annotated clusters

        Args:
            saveDir: directory for output files 
            dataName: text to append to the figure name
            df_marker_cell_type: known markers per cell types
            df_new_marker_cell_type: new markers per cell types

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.makePlotOfNewMarkers(dataName, saveDir)
        '''

        ORDERx = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_new_marker_cell_type.values.T, 'ward'), no_plot=True, get_leaves=True)['leaves']
        ORDERy = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df_new_marker_cell_type.values, 'ward'), no_plot=True, get_leaves=True)['leaves']

        genes = df_new_marker_cell_type.columns.values[ORDERx]
        celltypes = df_new_marker_cell_type.index.values[ORDERy]

        df_marker_cell_type = df_marker_cell_type[[celltype for celltype in celltypes if celltype in df_marker_cell_type.columns]]

        fig, ax = plt.subplots(figsize=(13,3))
        ax.imshow(df_new_marker_cell_type.values[ORDERy,:][:,ORDERx], cmap='Blues', interpolation='None', aspect='auto')
        ax.set_xticks([])
        ax.set_yticks(range(len(celltypes)))
        ax.set_yticklabels(celltypes, rotation=0, fontsize=8)
        ax.set_xlim([-0.5,df_new_marker_cell_type.shape[1] - 0.5])
        ax.set_ylim([-0.5,df_new_marker_cell_type.shape[0] - 0.5])

        for i, celltype in enumerate(celltypes):
            if celltype in df_marker_cell_type.columns:
                known_markers = df_marker_cell_type[celltype][df_marker_cell_type[celltype]>0].index.values
                xy = np.array([np.array([np.where(genes==marker)[0][0], i]) for marker in known_markers if marker in genes])
                print('Overlapping markers of %s: %s (%s)'%(celltype, len(xy), len(known_markers)))
                if len(xy)>0:
                    ax.plot(xy.T[0], xy.T[1], 'ro', ms=0.5)

        ax.set_title('Additional markers along with the overlapping part of the input (red)')
        fig.savefig(os.path.join(saveDir, dataName + '_new_markers.png'), dpi=600)

        return None
        

    # Main functions of class #################################################################################################################################
    def prepare(self, df_expr):

        '''Prepare pandas.DataFrame for input to function process()
        If input is pd.DataFrame validate the input whether it has correct structure.

        Args:
            df_expr: xxx
            df_expr: xxx
            df_expr: xxx

        Returns:
            Pandas DataFrame
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            df_expr = DCS.preapre(xxx)
        '''





        return df_expr

    def convertIndex(self, df_expr, nameFrom='alias', nameTo='hugo'):

        '''Convert index to hugo names, if any names in the index are
        duplicated, remove duplicates

        Args:
            df_expr: pandas.DataFrame to process

        Returns:
            Processed DataFrame

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

        Args:
            df_expr: pandas.DataFrame to clean

        Returns:
            Processed DataFrame
        
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

        Args:
            df_expr: pandas.DataFrame to normalize
            sigma_over_mean_sigma: threshold when keeping only genes with large enough standard deviation
            median: scale factor (optional), default None

        Returns:
            Processed DataFrame
        
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

        Args:
            df_expr: pandas.DataFrame to process
            PCAonly: perform Principal component analysis only, Default False
            do_fast_tsne: do FI-tSNE instead of "exact" tSNE, Default True

        Returns:
            Processed DataFrame
        
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
                    from DigitalCellSorter.FastTSNE import fast_tsne
                    X_tsne2 = fast_tsne(X_pca.T, perplexity = 30, seed=42).T
                else:
                    import fitsne
                    X_tsne2 = fitsne.FItSNE(np.array(X_pca.T, order='C')).T
            else:
                X_tsne2 = TSNE(n_components=2).fit_transform(X_pca.T).T

            return X_pca, _PCA.components_, X_tsne2

        return X_pca, _PCA.components_
    
    def cluster(self, X_pca, df_expr=None, clustering_f=AgglomerativeClustering):

        '''Cluster PCA-reduced data into a desired number of clusters

        Args:
            X_pca: PCA-reduced expression data
            df_expr: Gene expression data
            clustering_f: clustering function, if not provided then AgglomerativeClustering is used, otherwise should have .fit method and same input and output

        Returns:
            Cell cluster index labels
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            index = DCS.cluster(xPCA, 10)
        '''

        #import pynndescent
        #import networkx as nx
        #import community
        #
        #data = df_expr.values.T # X_pca.T
        #
        #n_neighbors = 40
        #print('Searching for %s nearest neighbors'%(n_neighbors))
        #knn = pynndescent.NNDescent(data, metric='euclidean', n_neighbors=n_neighbors).query(data, k=n_neighbors)
        #A = np.zeros((len(knn[0]),len(knn[0])))
        #for i in range(len(knn[0])):
        #    A[i, knn[0][i]] = knn[1][i]
        #print('Clustering the graph')
        #cellClusterIndex = pd.Series(community.best_partition(nx.from_numpy_array(A))).sort_index().values

        cellClusterIndex = clustering_f(n_clusters=self.nClusters).fit(X_pca.T).labels_

        return cellClusterIndex
    
    def vote(self, df_markers_expr, df_marker_cell_type, votingScheme):
        
        '''Produce cluster voting results

        Args:
            df_markers_expr: pandas.DataFrame with marker genes by cells expression
            df_marker_cell_type: pandas.DataFrame with marker genes by celltypes
            votingScheme: voting function

        Returns:
            Voting results, a dictionary in form of {cluster label: assigned cell type}
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            voting = DCS.vote(df_markers_expr, df_marker_cell_type)
        '''
        
        return votingScheme(df_markers_expr, df_marker_cell_type)

    def dcsVotingScheme(self, df_markers_expr, df_marker_cell_type):
        
        '''Produce cluster voting results

        Args:
            df_markers_expr: pandas.DataFrame with marker genes by cells expression
            df_marker_cell_type: pandas.DataFrame with marker genes by celltypes

        Returns:
            Voting results, a dictionary in form of {cluster label: assigned cell type}
        
        Usage:
            Function called by DCS internally
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
        df_V = self.getV((df_marker_cell_type, df_markers_expr, df_markers_expr.columns.get_level_values('cluster'), self.zScoreCutoff)).unstack()
        df_V.index.name = None

        # Generate random cluster configurations and calculate scores (Pkc) of
        # those
        print('Generating null distribution')
        startTime = getStartTime()
        pool = multiprocessing.Pool(processes = self.availableCPUsCount)
        random_df_V = pd.concat(pool.map(self.getV, [(df_marker_cell_type, df_markers_expr, randomClusterIndex[i], self.zScoreCutoff) for i in range(self.nSamplesDistribution)]), sort=False, axis=1)
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
        T[(df_L.values < minimumScoreForUnknown).all(axis=0)] = 'Unknown'
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

    def extractNewMarkerGenes(self, cluster=None, top=100, zScoreCutoff=0.3, removeUnknown=False):

        '''Extract new marker genes based on the cluster annotations

        Args:
            cluster: cluster #, if provided genes of only this culster will be returned
            top: upper bound for number of new markers per cell type
            zScoreCutoff: lower bound for a marker z-score to be significant
            removeUnknown: whether to remove type "Unknown"

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

        self.makePlotOfNewMarkers(self.saveDir, self.dataName, df_marker_cell_type, df_new_marker_genes)

        return None

    def getCells(self, celltype=None, clusterIndex=None, clusterName=None):

        '''Get cell annoations in a form of pandas.Series

        Args:
            celltype: cell type to extract, Default None

        Returns:
            Labelled cells
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.process(df_expr)
            labels = DCS.getLabelledCells()
        '''

        df = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r').columns.to_frame()

        if not clusterIndex is None:

            se = df.set_index(['patient', 'cell'])['cluster']

            condition = se==clusterIndex

            if not condition.any():
                print('No cells found')
                return None

            selectedCells = se[condition].index

            print('Selected %s cells from cluster: %s'%(len(selectedCells), clusterIndex))

            return selectedCells

        if not clusterName is None:

            se = df.set_index(['patient', 'cell'])['cluster']

            condition = se==eval(clusterName.split('#')[1])

            if not condition.any():
                print('No cells found')
                return None

            selectedCells = se[condition].index

            print('Selected %s cells from cluster: %s'%(len(selectedCells), clusterName))

            return selectedCells

        se = df.set_index(['patient', 'cell'])['label']

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

    def visualize(self):

        '''A convenient aggregate of visualization tools of this class.

        Args:
            None

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()
            DCS.process(df_expr)
            DCS.visualize()
        '''

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
        # Plot null distributions
        ##############################################################################################
        if self.toggleMakeHistogramNullDistributionPlot:
            print('Making null distributions plot')
            self.makeHistogramNullDistributionPlot(self.saveDir, self.dataName)
            timeMark()

        ##############################################################################################
        # Plot voting results matrix
        ##############################################################################################
        if self.toggleMakeVotingResultsMatrixPlot:
            print('Making voting results matrix plot')
            self.makeVotingResultsMatrixPlot(self.saveDir, self.dataName)
            timeMark()
        
        ##############################################################################################
        # Plot mean marker expression
        ##############################################################################################
        if self.toggleMakeMarkerExpressionPlot:
            print('Making marker expression plot')
            self.makeMarkerExpressionPlot(self.dataName, self.saveDir)
            timeMark()

        ##############################################################################################
        # Make tSNE plots
        ##############################################################################################
        if self.toggleMakeTSNEplotQC and self.toggleDoQualityControl:
            print('Making tSNE plots of QC')
            df_QC = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='QC', mode='r')
            df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne_pre_QC', mode='r')
            goodQUalityCells = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r').columns
            self.makeTSNEplot(df_tsne.values, df_QC['number_of_genes'].values, self.dataName, self.saveDir, 'by_number_of_genes', legend=False, labels=False, colorbar=True)
            self.makeTSNEplot(df_tsne.values, df_QC['count_depth'].values, self.dataName, self.saveDir, 'by_count_depth', legend=False, labels=False, colorbar=True)
            self.makeTSNEplot(df_tsne.values, df_QC['fraction_of_mitochondrialGenes'].values, self.dataName, self.saveDir, 'by_fraction_of_mitochondrialGenes', legend=False, labels=False, colorbar=True)
            self.makeTSNEplot(df_tsne.values, np.array([cell in goodQUalityCells for cell in df_QC.index]), self.dataName, self.saveDir, 'by_is_quality_cell', legend=False, labels=True)
            timeMark()    

        if self.toggleMakeTSNEplotClusters:
            print('Making tSNE plot by clusters')
            df_clusters = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='r')
            df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
            self.makeTSNEplot(df_tsne.values, np.array(['Cluster #%s' % (label[0]) for label in df_clusters.values]), self.dataName, self.saveDir, 'by_clusters')
            timeMark()

        if self.toggleMakeTSNEplotBatches:
            print('Making tSNE plot by patients')
            df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
            self.makeTSNEplot(df_tsne.values, df_tsne.columns.get_level_values('patient').values, self.dataName, self.saveDir, 'by_patients')
            timeMark()

        if self.toggleMakeTSNEplotAnnotatedClusters:
            print('Making tSNE plot by clusters "True" labels')
            df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r')
            df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
            df_tsne = df_tsne[df_markers_expr.groupby(level=['patient', 'cell'], sort=False, axis=1).count().columns]
            with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
                colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
            self.makeTSNEplot(df_tsne.values, df_markers_expr.columns.get_level_values('label'), self.dataName, self.saveDir, 'by_clusters_annotated',
                             colormap=colormap, legend=False)
            timeMark()

        if self.toggleAnomalyScoresTSNEplot:
            self.getAnomalyScoresPlot()
            timeMark()
           
        ##############################################################################################
        # Make stacked barplot of cell type fractions
        ##############################################################################################
        if self.toggleMakeStackedBarplot:        
            self.makeStackedBarplot(self.saveDir, self.dataName, self.subclusteringName)
            timeMark()
        
        ##############################################################################################
        # Make tSNE plots showing relative expression of different markers (one
        # for each marker)
        ##############################################################################################
        if self.toggleMakeMarkerSubplots:
            df_tsne = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_tsne', mode='r')
            df_markers_expr = pd.read_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_markers_expr', mode='r')
            df_tsne = df_tsne[pd.MultiIndex.from_arrays([df_markers_expr.columns.get_level_values('patient'), df_markers_expr.columns.get_level_values('cell')])]
            hugo_cd_dict = dict(zip(df_markers_expr.index.values.tolist(), self.gnc.Convert(list(df_markers_expr.index), 'hugo', 'alias', returnUnknownString=False)))
            self.makeMarkerSubplots(df_markers_expr, df_tsne.values, hugo_cd_dict, self.dataName, self.saveDir)
            timeMark()

        return None

    def process(self, df_expr, visualize=True):

        '''Main function

        Args:
            df_expr: gene expression data in pandas.DataFrame

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
        df_expr = self.convertIndex(df_expr)
        df_expr = self.clean(df_expr)
        timeMark()
        
        ##############################################################################################
        # Calculate QC measures
        ##############################################################################################
        if self.toggleDoQualityControl:
            if self.mitochondrialGenes is None:
                self.mitochondrialGenes = pd.read_excel(os.path.join('geneLists', 'Human.MitoCarta2.0.xls'), sheet_name='A Human MitoCarta2.0', usecols=[3], header=0)['Symbol'].values.squeeze().tolist()
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
            patients = df_expr.columns.get_level_values('patient').values.copy()

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
        cellClusterIndexLabel = self.cluster(df_xpca.values, df_expr=df_expr) #clustering_f = sklearn.cluster.KMeans
        df_clusters = pd.DataFrame(data=cellClusterIndexLabel, index=df_expr.columns)
        df_clusters.columns = ['cluster']
        df_clusters.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_clusters', mode='a', complevel=4, complib='zlib')
        df_expr = pd.concat([df_clusters, df_expr.T], sort=False, axis=1).reset_index().set_index(['patient', 'cell', 'cluster']).T
        timeMark()

        ##############################################################################################
        # Calculate average of each gene in each cluster
        ##############################################################################################
        df_gene_cluster_centroids = df_expr.groupby(level=['cluster'], sort=True, axis=1).mean()
        df_gene_cluster_centroids.to_hdf(os.path.join(self.saveDir, self.dataName + '_processed.h5'), key='df_gene_cluster_centroids', mode='a', complevel=4, complib='zlib')

        ##############################################################################################
        # Get dictionary to map from markers to cell types. Select genes from the marker list only
        ##############################################################################################
        subtypeToType = pd.read_excel(self.geneListFileName, sheet_name='CellTypesGrouped', index_col='CellType').to_dict()['CellTypeGrouped']
        df_marker_cell_type = pd.read_excel(self.geneListFileName, sheet_name='MarkerCellType', index_col='Marker')
        df_marker_cell_type.index = self.gnc.Convert(list(df_marker_cell_type.index), 'alias', 'hugo', returnUnknownString=False)
        df_marker_cell_type.columns = pd.MultiIndex.from_arrays([df_marker_cell_type.columns,[subtypeToType[subtype] for subtype in df_marker_cell_type.columns]], names=['CellType','CellTypeGrouped'])
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
        df_markers_expr = df_markers_expr.T.set_index(['patient', 'cell', 'cluster', 'label']).T.apply(pd.to_numeric)
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
        self.extractNewMarkerGenes()

        ##############################################################################################
        # Make all plots to conclude analysis
        ##############################################################################################
        if visualize:
            self.visualize()

        return None
