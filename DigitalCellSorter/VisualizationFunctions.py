import os
import copy

import numpy as np
import pandas as pd

import scipy.stats
import scipy.signal
import scipy.ndimage
import scipy.cluster.hierarchy
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import gaussian_filter

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib import cm

import plotly.express as px
import plotly.graph_objects as go
from plotly.offline import plot as plot_offline
from plotly.offline import plot_mpl
from adjustText import adjust_text

from .GenericFunctions import read, write

class VisualizationFunctions:

    '''Class of visualization functions for DigitalCellSorter'''

    def __init__(self, dataName = 'dataName', saveDir = os.path.join(''), matplotlibMode = 'Agg', safePlotting = True):

        '''Function called automatically'''

        self.saveDir = saveDir
        self.dataName = dataName
        self.matplotlibMode = matplotlibMode
        self.safePlotting = safePlotting

        return

    @property
    def matplotlibMode(self):

        return self._matplotlibMode

    @matplotlibMode.setter
    def matplotlibMode(self, value):

        self._matplotlibMode = value

        if not self.matplotlibMode is None:
            matplotlib.use(self.matplotlibMode)

        return

    def tryExcept(func):

      def internal(self, *args, **kwargs):

        if self.safePlotting:
            try:
              func(self, *args, **kwargs)

            except Exception as exception:
              print('Something went wrong while making plot: %s' % (func))
              print('\tError message: %s\n' % (exception))
        else:
            func(self, *args, **kwargs)

      return internal

    def saveFigure(self, fig, saveDir, label = 'Figure', extension = 'png', dpi = 300, close = True, attemptSavingHTML = False):

        '''Function used internally to save and close figures

        Parameters:
            saveDir: str
                Path of directories to save the object to

            label: str, Default 'Figure'
                Name of the figure to save

            extension: str, Default '.png'
                Path of directories to save the object to
            
            dpi: int, Default 300
                Figure resolution if rasterized

            close: boolean: Default True
                Whether to close the figure after saving

        Returns:
            None

        Usage:
            saveFigure(fig, saveDir, label, extension, dpi)
        '''

        if saveDir != os.path.join('') and not os.path.exists(saveDir):
            os.makedirs(saveDir)

        try:
            if not extension[0] == '.':
                extension = ''.join(['.', extension])
        except Exception as exception:
            print(exception)
            print('Figure extension/format error')
            print('Example of acceptable extension: \".png\"')

            return

        if extension in ['.png', '.jpeg', '.tiff']:
            try:
                fig.savefig(os.path.join(saveDir, label + extension), dpi=dpi)
            except Exception as exception:
                print(exception)

        elif extension in ['.svg', '.eps', '.pdf']:
            try:
                fig.savefig(os.path.join(saveDir, label + extension))
            except Exception as exception:
                print(exception)
        else:
            print('Unsupported format. Figure not saved')

        if attemptSavingHTML:
            try:
                plot_mpl(fig, filename=os.path.join(saveDir, label + '.html'), auto_open=False)
            except Exception as exception:
                print('Saving to iteractive HTML did not succeed')

        if close:
            try:
                plt.close(fig)
            except Exception as exception:
                print(exception)
                print('Error while closing figure')

        return

    # MatPlotLib-powered figures

    @tryExcept
    def makeHeatmapGeneExpressionPlot(self, df = None, genes = None, normalize = True, saveExcel = True, nameToAppend = 'heatmap', plotBy = 'cluster', figsize = (8, 4), convertGenes = False, dpi = 300, extension = 'png', **kwargs):

        '''Make heatmap gene expression plot from a provided gene expression matrix.

        Parameters:
            df: pandas.DataFrame
                Gene expression matrix

            genes: list, Default None
                List of genes to plot
    
            nameToAppend: str, Default ''
                String to append to fifure file name

            dpi: int, Default 300
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeHeatmapGeneExpressionPlot()
        '''

        if df is None:
            if self.df_expr is None:
                self.loadExpressionData()

            if self.df_expr is None:
                return
            else:
                if not genes is None:
                    if convertGenes:
                        converted = []
                        for gene in genes:
                            converted.extend(self.getHugoName(gene, printAliases=True))

                        converted = np.unique(converted)

                        common = self.df_expr.index.intersection(converted)
                    else:
                        common = self.df_expr.index.intersection(genes).drop_duplicates()

                    df = self.df_expr.loc[common].copy()
                else:
                    print('Plotting all expressed genes not supported. Provide a smaller list of genes')

                    return

        counts = df.loc[[df.index[0]]].groupby(axis=1, level=plotBy).count()
        means = df.mean(axis=1)

        df = df.groupby(axis=1, level=plotBy).mean()
        #df = df.replace(0., np.nan).groupby(axis=1, level=plotBy).mean().fillna(0.)
        df.columns = df.columns.get_level_values(plotBy)
        df.columns = list(zip(df.columns.values, counts.values[0]))

        if normalize:
            for i in range(df.shape[0]):
                df.iloc[i,:] -= np.min(df.iloc[i,:])
                df.iloc[i,:] /= np.max(df.iloc[i,:])
        
        #df = df.T.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df.T, 'ward'), no_plot=True, get_leaves=True)['leaves']].T
        df = df.iloc[scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df, 'ward'), no_plot=True, get_leaves=True)['leaves']]

        df.insert(0, ('Mean', 'All'), means)

        if saveExcel:
            df.T.to_excel(os.path.join(self.saveDir, self.dataName + '_' + nameToAppend + '_genes_by_%s' % (plotBy) + '.xlsx'))

        fig, ax = plt.subplots(figsize=figsize)

        ax.imshow(df.T.values[1:,:], cmap='Blues', interpolation='None', aspect='auto', 
                  extent=(-0.5, df.shape[0] - 0.5, df.shape[1] - 0.5, +0.5))

        ax.imshow(df.T.values[:1,:], cmap='Reds', interpolation='None', aspect='auto', 
                  extent=(-0.5, df.shape[0] - 0.5, -0.5, +0.5))

        ax.axhline(y=0.5, c='k', lw=1.5)

        ax.set_xticks(range(df.shape[0]))
        ax.set_yticks(range(df.shape[1]))

        ylabels = ['(' + str(col[1]) + ')%s#' % ('    ' if len(col[0]) <= 3 else '  ') + str(col[0]) for col in df.columns]
        ylabels[0] = 'Mean across all cells'

        ax.set_xticklabels(df.index, rotation=90, fontsize=10)
        ax.set_yticklabels(ylabels, rotation=0, fontsize=10)

        ax.set_xlim([-0.5, df.shape[0] - 0.5])
        ax.set_ylim([-0.5, df.shape[1] - 0.5])

        fig.tight_layout()

        self.saveFigure(fig, self.saveDir, self.dataName + '_' + nameToAppend + '_genes_by_%s' % (plotBy), extension=extension, dpi=dpi, **kwargs)

        return fig
    
    @tryExcept
    def makeMarkerExpressionPlot(self, fontscale = 2., dpi = 300, extension = 'png', **kwargs):

        '''Produce image on marker genes and their expression on all clusters.
        Uses files generated by function DCS.Vote

        Parameters:
            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeMarkerExpressionPlot()
        '''

        df_votingResults = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='z-scores')
        votingResults = dict(zip(df_votingResults['cluster'].values, df_votingResults['Predicted cell type'].values))
        predictedCelltypes = dict(zip(df_votingResults['cluster'].values, df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0]))
        supportingMarkersList = dict(zip(df_votingResults['cluster'].values, df_votingResults['Supporting markers'].str.split(' // ')))
        allMarkersList = dict(zip(df_votingResults['cluster'].values, df_votingResults['All markers'].str.split(' // ')))
        df_markers_cluster_centroids = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='Cluster centroids', index_col=0, header=0).T

        df_markers = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='Marker cell type weight matrix', index_col=0)
        df_markers_weighted = df_markers.copy()
        df_markers[df_markers >= 0.] = np.nan
        df_markers[df_markers < 0.] = -1.

        # Y_mc.T
        X_markers_cluster_means_transpose = df_markers_cluster_centroids.values.T

        df_means = df_markers_cluster_centroids.copy()
        
        # Normalization
        for i in range(X_markers_cluster_means_transpose.shape[1]):
            X_markers_cluster_means_transpose[:,i] -= np.min(X_markers_cluster_means_transpose[:,i])
            X_markers_cluster_means_transpose[:,i] /= np.max(X_markers_cluster_means_transpose[:,i])

        ORDER = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(X_markers_cluster_means_transpose, 'ward'), no_plot=True, get_leaves=True)['leaves']
        ORDER2 = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(X_markers_cluster_means_transpose.T, 'ward'), no_plot=True, get_leaves=True)['leaves']

        df_markers_weighted = df_markers_weighted.iloc[:, ORDER2]

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

        df_neg_supp_marker_hits = pd.DataFrame(data=np.zeros((df_markers_cluster_centroids.shape)), index=df_markers_cluster_centroids.index, columns=df_markers_cluster_centroids.columns)
        for cluster in supportingMarkersList:
            if not supportingMarkersList[cluster] is np.nan:
                for gene in allMarkersList[cluster]:
                #for gene in df_markers.columns:
                    if df_markers.loc[predictedCelltypes[cluster], gene] == -1.:
                        df_neg_supp_marker_hits.loc[gene, cluster] = 1

        X_marker_hits = df_all_marker_hits.values.T[ORDER,:][:,ORDER2]
        X_supp_marker_hits = df_supp_marker_hits.values.T[ORDER,:][:,ORDER2]
        X_neg_supp_marker_hits = df_neg_supp_marker_hits.values.T[ORDER,:][:,ORDER2]

        _figsize = np.float_(X_markers_cluster_means_transpose.shape[::-1]) / \
                    np.max(X_markers_cluster_means_transpose.shape) * 15.0 + 2.0

        _figsize[1] *= 1.5

        height_ratio = df_markers_cluster_centroids.shape[1] / (1. * df_markers_weighted.shape[0])

        gs = matplotlib.gridspec.GridSpec(2, 2, width_ratios=[20,1], height_ratios=[height_ratio,1], 
                                          left=0.13, right=0.99, top=0.99, bottom=0.25, wspace=0.01, hspace=0.04)

        fig = plt.figure(figsize=_figsize)

        if True:
            ax = plt.subplot(gs[0])

            cell_counts = df_votingResults['# cells in cluster'].values.copy()[ORDER]
            means = (df_means.iloc[ORDER2,ORDER] * cell_counts / cell_counts.sum()).sum(axis=1)

            ax.imshow(means.values[None, :], cmap='Reds', interpolation='None', aspect='auto', 
                      extent=(-0.5, means.shape[0] - 0.5, -0.5, +0.5))

            ax.imshow(X_markers_cluster_means_sorted,cmap='Blues', interpolation='None', aspect='auto',
                      extent=(-0.5, X_markers_cluster_means_sorted.shape[1] - 0.5, X_markers_cluster_means_sorted.shape[0] - 0.5 + 1.0, +0.5))

            i_list,j_list = np.where(X_marker_hits.T > 0)
            color = 'w' #(1., 1., 0.7)
            ax.plot(i_list, j_list + 1., 'k*', mec=color, mew=0.5, markersize=4)

            i_list_supp, j_list_supp = np.where(X_supp_marker_hits.T > 0)
            i_list_neg_supp, j_list_neg_supp = np.where(X_neg_supp_marker_hits.T > 0)
            ax.plot(i_list_supp, j_list_supp + 1., 'k*', mec='lime', mew=0.7, markersize=4) #mec='k', alpha=0.5, markersize=6
            ax.plot(i_list_neg_supp, j_list_neg_supp + 1., 'k*', mec='red', mew=0.7, markersize=4) #mec='k', alpha=0.5, markersize=6

            ax.set_xticks([])
            ax.set_yticks(range(X_markers_cluster_means_transpose.shape[0] + 1))

            clusterNames = list(votingResults.values())
            clusterIndices = list(votingResults.keys())

            ax.set_yticklabels(['Mean'] + [str(clusterNames[i]) + ' (' + str(clusterIndices[i]) + ')' for i in ORDER], rotation=0, fontsize=6*fontscale)
            ax.set_xlim([-0.5,X_markers_cluster_means_transpose.shape[1] - 0.5])
            ax.set_ylim([-0.5,X_markers_cluster_means_transpose.shape[0] + 1 - 0.5])

        if True:
            ax2 = plt.subplot(gs[1])
            fontsize = 5*fontscale
            cells_in_clusters = df_votingResults['# cells in cluster'].values.copy()[ORDER]
            numberOfCells = cells_in_clusters.sum()

            with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
                colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
            celltypes = df_votingResults['Predicted cell type'].str.split(' #', expand=True)[0].values.copy()[ORDER]

            ax2.barh(y=range(len(cells_in_clusters)), width=cells_in_clusters, height=0.8, align='center', color=[colormap[i] for i in celltypes])
            for i in range(len(cells_in_clusters)):
                ax2.text(np.max(cells_in_clusters), i, cells_in_clusters[i], ha='right',va='top', color='k', weight='bold', fontsize=fontsize)
                ax2.text(0.02 * numberOfCells, i, str(round(100 * cells_in_clusters[i] / numberOfCells, 1)) + '%', ha='left',va='bottom', color='b', fontsize=fontsize)
            ax2.set_xticklabels(cells_in_clusters, fontsize=fontsize)
            ax2.set_yticklabels(cells_in_clusters, alpha=0)
            ax2.set_xticklabels(cells_in_clusters, alpha=0)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlabel('Number of\ncells in clusters', fontsize=fontsize)
            ax2.set_ylim(-1.5, len(cells_in_clusters) - 0.5)

        if True:
            ax3 = plt.subplot(gs[2])

            masked = np.ma.array(df_markers_weighted.values, mask=(df_markers_weighted.values == 0.))

            cmap = plt.cm.PiYG
            #cmap = matplotlib.colors.LinearSegmentedColormap.from_list('RedGreen', [(1, 0, 0), (0,
            #1, 0)], N=100)
            cmap.set_bad('white')

            value = 0.5 * np.abs(df_markers_weighted).max().max()
            ax3.imshow(masked, cmap=cmap, vmin=-value, vmax=+value, interpolation='None', aspect='auto')

            ax3.set_xticks(range(X_markers_cluster_means_transpose.shape[1]))
            ax3.set_yticks(range(df_markers_weighted.shape[0]))

            xtickslabels = np.array(df_markers_cluster_centroids.index[ORDER2])
            for i in range(0,len(xtickslabels),2):
                xtickslabels[i] += " ─────────"

            ax3.set_xticklabels(xtickslabels, rotation=90, fontsize=5*fontscale)
            ax3.set_yticklabels(df_markers_weighted.index.values, rotation=0, fontsize=8*fontscale)

            ax3.set_xlim([-0.5, df_markers_weighted.shape[1] - 0.5])
            ax3.set_ylim([-0.5, df_markers_weighted.shape[0] - 0.5])

        self.saveFigure(fig, self.saveDir, self.dataName + '_marker_expression', extension=extension, dpi=dpi, **kwargs)

        return fig

    @tryExcept
    def internalMakeMarkerSubplots(self, df, X_projection, hugo_cd_dict, NoFrameOnFigures = False, HideClusterLabels = False, outlineClusters = True, analyzeBy = 'cluster', saveSubDir = 'marker_subplots', dpi = 300, extension = 'png', **kwargs):

        '''Produce subplots on each marker and its expression on all clusters

        Parameters:
            df: pandas.DataFrame 
                Data with marker genes expression

            X_projection: 2d numpy.array
                2D coordinates for each cell

            hugo_cd_dict: dictionary 
                With aliases for hugo names of genes

            NoFrameOnFigures: boolean, Default False
                Whether to include frame on the figure

            HideClusterLabels: boolean, Default False
                Whether to print cluster labels on the figure

            outlineClusters: boolean, Default True
                Whether to outline the clusters with circles

            analyzeBy: str, Default 'cluster'
                What level of lablels to include.
                Other possible option is 'label'

        Returns:
            None
        
        Usage:
            Function used internally

            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.internalMakeMarkerSubplots(df_markers_expr, projection, hugo_cd_dict)
        '''
        
        def MarkerSubplot(counter, marker, df, analyzeBy, X_projection, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, NoFrameOnFigures, HideClusterLabels, XLIM, YLIM, directory, circles):

            fig,ax = plt.subplots(figsize=(8,8))

            ax.cla()
            suffix = '(' + str(hugo_cd_dict[marker]).replace('\"', '').replace('\'', '').replace('(', '').replace(')', '').replace(' ','') + ')'
            ax.plot(np.nan,np.nan,'*',markersize=15,c=cm.seismic(1.0),label=marker + '\n' + suffix.replace(',','\n'))
            circleIndices = np.where(df.loc[marker].values == 0)[0] # cells that don't have this marker
            starIndices = np.where(df.loc[marker].values > 0)[0] # cells that have this marker
            starIndices = starIndices[np.argsort(df.loc[marker].values[starIndices])]
            args1 = [X_projection[0,circleIndices],
                        X_projection[1,circleIndices]]
            kwargs1 = {'marker':'o',
                        'c':'b',
                        'alpha':0.1,
                        's':6 * 3,
                        'linewidth':0,}
            args2 = [X_projection[0,starIndices],
                        X_projection[1,starIndices]]
            kwargs2 = {'marker':'*',
                        'c':cm.seismic(df.loc[marker].values[starIndices] / np.max(df.loc[marker].values[starIndices])),
                        's':30 * 4,
                        'linewidth':0.0,}
            ax.scatter(*args1,**kwargs1)
            ax.scatter(*args2,**kwargs2)
            for label in set(cellClusterIndexLabel):
                # cells with this label
                X_projection2_cluster = X_projection[:,cellClusterIndexLabel == label]
                x_mean = np.mean(X_projection2_cluster[0,:])
                y_mean = np.mean(X_projection2_cluster[1,:])

                _text_label = label if not HideClusterLabels else ''

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
                    radius = np.sqrt(X_projection2_cluster.shape[1]) * 300.0
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
            if self.saveDir is not None: 
                self.saveFigure(fig, directory, '%s_%s_%s_%s' % (self.dataName,marker,suffix.replace(',','_').replace('/','_'),analyzeBy), extension=extension, dpi=dpi, **kwargs)

                print(marker, end=" ", flush=True)

            return
        
        maxs = np.max(X_projection,axis=1)
        mins = np.min(X_projection,axis=1)
        maxDiffs = maxs - mins
        deltas = maxDiffs * 0.05
        XLIM = [mins[0] - deltas[0],maxs[0] + deltas[0]]
        YLIM = [mins[1] - deltas[1],maxs[1] + deltas[1]]

        if len(df.index) > 1:
            print('\nSaving marker expression plots:\n')
        else:
            print('Saving expression plot of:', end=' ', flush=True)

        if analyzeBy == 'celltype':
            try:
                index = df.columns.get_level_values('label').str.split(' #', expand=True).get_level_values(0).values
            except:
                index = df.columns.get_level_values('celltype').values
        else:
            index = df.columns.get_level_values(analyzeBy).values

        for counter,marker in enumerate(df.index.values):
            MarkerSubplot(counter, marker, pd.DataFrame(data=np.reshape(np.array(df.loc[marker]), (1,len(df.loc[marker]))), columns=df.columns, index=[marker]), analyzeBy, X_projection, index, hugo_cd_dict, self.dataName, self.saveDir, NoFrameOnFigures, HideClusterLabels, XLIM, YLIM, os.path.join(self.saveDir, saveSubDir, ''), outlineClusters)
        print()

        return

    @tryExcept
    def makeAnnotationResultsMatrixPlot(self, dpi = 300, extension = 'png', **kwargs):

        '''Produce voting results voting matrix plot

        Parameters:
            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeAnnotationResultsMatrixPlot()
        '''

        df_votingResults = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='z-scores')

        cellTypes = sorted([x for x in df_votingResults.columns.values.tolist() if x not in ['cluster', 'Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'Contradicting markers', 'All markers']])

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
        
        label_max = 0.5 * max([len(assigned_names_of_clusters[i]) for i in range(len(assigned_names_of_clusters))])

        _figsize = np.float_((num_of_cell_types,num_of_clusters)) / np.max(num_of_clusters) * 15.0
        _figsize[0] += 1.0 + label_max
        _figsize[1] += 2.0

        gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[3,1])

        fig = plt.figure(figsize=_figsize)
        
        ax = fig.add_axes([0.15, 0.125, 0.65, 0.85])
        axx = fig.add_axes([0.76, 0.125, 0.19, 0.85])

        zscores = df_votingResults[cellTypes].values
        ax.imshow(zscores, aspect=1, cmap='Greens', vmin=0, vmax=0.5, interpolation='None')

        for i in range(num_of_clusters):
            for j in range(num_of_cell_types):
                if np.round(zscores[i,j],1) > 0:
                    if zscores[i,j] == np.max(zscores[i,:]):
                        ax.text(j,i,np.round(zscores[i,j],1), color='w', fontsize=125 * 4 / max([num_of_cell_types, num_of_clusters]),ha='center',va='center').set_path_effects([path_effects.Stroke(linewidth=4, foreground='red'),path_effects.Normal()])
                    else:
                        ax.text(j,i,np.round(zscores[i,j],1), color='w', fontsize=125 * 3 / max([num_of_cell_types, num_of_clusters]),ha='center',va='center').set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()])

        ax.set_xticks(range(num_of_cell_types))
        ax.set_yticks(range(num_of_clusters))

        ytickslabels = copy.deepcopy(assigned_names_of_clusters)
        for i in range(len(ytickslabels)):
            ytickslabels[i] = str(assigned_names_of_clusters[i]) + ' (' + str(indicis_of_clusters[i]) + ')'

        xtickslabels = np.array(copy.deepcopy(cellTypes))
        #for i in range(len(xtickslabels)):
        #    if i % 3 == 1: xtickslabels[i] = '\n' + xtickslabels[i]
        #    if i % 3 == 2: xtickslabels[i] = '\n\n' + xtickslabels[i]

        ax.set_xticklabels(xtickslabels, rotation=30, fontsize=20, ha='right')
        ax.set_yticklabels(ytickslabels, fontsize=20, rotation=0) 
        ax.set_xlim([-0.5,num_of_cell_types - 0.5])
        ax.set_ylim([-0.5,num_of_clusters - 0.5])

        cells_in_clusters = df_votingResults['# cells in cluster'].values

        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
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
        
        self.saveFigure(fig, self.saveDir, self.dataName + '_scores_matrix', extension=extension, dpi=dpi, **kwargs)
 
        return fig

    @tryExcept
    def makeHistogramNullDistributionPlot(self, dpi = 600, extension = 'png', **kwargs):

        '''Produce histogram plot of the voting null distributions

        Parameters:
            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeHistogramNullDistributionPlot()
        '''

        try:
            df_noise_dict = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='Null distributions', index_col=0, header=[0,1], skiprows=[2])
        except Exception as exception:
            print(exception)
            print('Error loading distributions from the results file')

            return

        if len(df_noise_dict) == 0:
            print('Null distribution is empty in the results file')

            return

        df_votingResultsV = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='Voting scores', dtype={'cluster':str}).reset_index().set_index('cluster')
        df_votingResultsZ = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='z-scores', dtype={'cluster':str}).reset_index().set_index('cluster')

        predicted_cell_type_cluster = df_votingResultsZ['Predicted cell type'].values
        predicted_cell_type = df_votingResultsZ['Predicted cell type'].str.split(' #', expand=True)[0].values

        cellTypes = sorted([x for x in df_votingResultsV.columns.values.tolist() if x not in ['cluster', 'Predicted cell type', '# cells in cluster', 'Winning score', 'Supporting markers', 'All markers']])

        df_votingResultsV = df_votingResultsV[cellTypes]
        df_votingResultsZ = df_votingResultsZ[cellTypes]

        cell_types = np.unique(df_noise_dict.columns.get_level_values(0).values)
        num_of_cell_types = cell_types.shape[0]
        clusters = np.unique(df_noise_dict.columns.get_level_values(1).values).astype(str)
        num_of_clusters = clusters.shape[0]

        maxy = np.round(np.nanmax(df_noise_dict.values), 2)

        with open(os.path.join(self.saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}

        origWidth = matplotlib.rcParams['axes.linewidth']
        matplotlib.rcParams['axes.linewidth'] = 0.1

        gs = matplotlib.gridspec.GridSpec(num_of_clusters, num_of_cell_types, hspace=0.45, wspace=0.1, bottom=0.04, top=0.96, left=0.05, right=0.99)
        
        fig = plt.figure(figsize=(num_of_cell_types, num_of_clusters * 0.4))

        for i in range(num_of_cell_types):

            try:
                minx, maxx = np.round(df_noise_dict.index.values[np.where(df_noise_dict.xs(key=cell_types[i], level=0, axis=1).values.sum(axis=1) != 0.)[0][[0,-1]]], 3)
                #maxx += 0.1
                maxx = np.round(maxx, 3)
            except Exception as exception:
                print(exception)
                minx = np.round(df_noise_dict.index.values[0], 3)
                maxx = np.round(df_noise_dict.index.values[-1], 3)

            print(cell_types[i], end=':\t')

            for j in range(num_of_clusters):

                print(j, end=',', flush=True)

                fontsize = 3

                ax = plt.subplot(gs[i + num_of_cell_types * j])

                ax.bar(df_noise_dict.index.values, df_noise_dict[(cell_types[i], clusters[j])].values, 
                       width=df_noise_dict.index.values[2] - df_noise_dict.index.values[1], align='center', color=colormap[cell_types[i]])

                valueV = df_votingResultsV.loc[clusters[j], cell_types[i]]
                valueZ = df_votingResultsZ.loc[clusters[j], cell_types[i]]
                ax.axvline(x=valueV, ymin=0, ymax=1, color='k', lw=0.2)
                color = 'k' if predicted_cell_type[j] != cell_types[i] else 'r'
                xloc = 0.02 * minx + minx #valueV
                ax.text(xloc, maxy - 0.02 * maxy, r'$V_{%s,%s}=$' % (i,j) + str(np.round(valueV,2)), fontsize=fontsize, va='top', ha='left', color=color, zorder=np.inf)
                ax.text(xloc, maxy - 0.2 * maxy, r'$\Lambda_{%s,%s}=$' % (i,j) + str(np.round(valueZ,2)), fontsize=fontsize, va='top', ha='left', color=color, zorder=np.inf)

                if j == 0:
                    ax.set_title(cell_types[i], fontdict={'color': 'b', 'size':'6'})

                if i == 0:
                    ax.text(0. * maxx, maxy + 0.05 * maxy, predicted_cell_type_cluster[j] + ' (Cluster %s)' % clusters[j], rotation=0, 
                            fontsize=fontsize, weight='bold', va='bottom', ha='left', color='k')
                    
                    ax.set_ylabel('Probability', fontsize=fontsize)
                    ax.set_yticklabels([0., maxy], fontsize=fontsize)
                else:
                    ax.set_yticklabels([], fontsize=fontsize)

                if j == num_of_clusters - 1:
                    #ax.set_xlabel('Voting score', fontsize=fontsize)
                    ax.set_xticklabels([minx, maxx], fontsize=fontsize)

                    if maxx > 0.0 and minx < 0.0:
                        ax.text(0.0, -0.2 * maxy, '0.0', fontsize=fontsize, va='top', ha='center', color='k')
                else:
                    ax.set_xticklabels([], fontsize=fontsize)

                ax.set_xticks([minx, maxx])
                ax.set_yticks([0., maxy])

                ax.set_xlim(minx, maxx)
                ax.set_ylim(0., maxy)

                ax.tick_params(direction='in', length=1, width=0.1, colors='k')

            print()
        
        self.saveFigure(fig, self.saveDir, self.dataName + '_null_distributions', extension=extension, dpi=dpi, **kwargs)

        matplotlib.rcParams['axes.linewidth'] = origWidth

        return fig

    @tryExcept
    def makeProjectionPlot(self, Xprojection, cellClusterIndexLabel, suffix = '', colormap = cm.jet, legend = True, labels = True, colorbar = False, fontsize = 10, plotNaNs = True, rightShift = 0.3, dpi = 300, extension = 'png', **kwargs):

        '''Produce projection plot (2D layout) with a specified coloring scheme

        Parameters:
            Xprojection: 2D coordinates for each cell

            cellClusterIndexLabel: cluster index for each cell

            suffix: str
                Text label to append to the figure name

            colormap: cell coloring sequence, can be a dictionary or cm.colormap, 
                Default matplotlib.colors.LinearSegmentedColormap.jet

            legend: boolean, Default True
                Whether to print legend

            labels: boolean, Default True
                Whether to print labels

            colorbar: boolean, Default False
                Whether to show colorbar
                Use with non-numerical values will raise an error

            fontsize: int, Default 10
                Labels and legend font size

            plotNaNs: boolean, Default True
                Whether to plot NaN labels (in grey)
                
            rightShift: float, Default 0.3
                Fraction of space to leave on the right-hand side of the plot.
                This parameter is useful for adjusting legend overlap with data points.

            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeProjectionPlot(projection, cellClusterIndexLabel, suffix)
        '''

        def add_colorbar(fig, labels, cmap = matplotlib.colors.LinearSegmentedColormap.from_list('GR', [(0, 1, 0), (1, 0, 0)], N=100), fontsize = 10):
    
            mapp = cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=np.min(labels), vmax=np.max(labels)), cmap=cmap)
            mapp.set_array(labels)
            sp = np.linspace(np.max(labels), np.min(labels), num=6, endpoint=True)

            axisColor = fig.add_axes([0.9,0.5,0.01,0.4])
            fig.colorbar(mapp, cax=axisColor, ticks=sp)

            axisColor.tick_params(labelsize=fontsize)
            axisColor.set_yticklabels(np.round(sp,2))
    
            return None

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.05,0.05,0.9,0.9])

        maxs, mins = np.max(Xprojection,axis=1), np.min(Xprojection,axis=1)

        missing = np.where(cellClusterIndexLabel != cellClusterIndexLabel)[0]
        if len(missing) > 0:
            ax.plot(Xprojection[0, missing], Xprojection[1, missing], 'o', color='grey', mew=0.5, alpha=0.2, markeredgecolor='k', label='NaN')

        nonMissing = np.where(cellClusterIndexLabel == cellClusterIndexLabel)[0]
        cellClusterIndexLabel = np.array(cellClusterIndexLabel)[nonMissing]
        Xprojection = Xprojection[:, nonMissing]

        possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel))

        if labels:
            print(possible_cluster_labels)

        texts = []

        for ilabel, label in enumerate(possible_cluster_labels):

            if type(colormap) is matplotlib.colors.LinearSegmentedColormap:
                color = colormap(ilabel / len(possible_cluster_labels))
            elif type(colormap) is str:
                color = plt.get_cmap(colormap)(ilabel / len(possible_cluster_labels))
            else:
                color = colormap[label.split(' #')[0]]

            XprojectionC = Xprojection[:,cellClusterIndexLabel == label]

            ax.plot(XprojectionC[0,:], XprojectionC[1,:], 'o', color=color, mew=0.5, alpha=0.3, markeredgecolor='k', label=label)

            if labels:
                text = ax.text(np.median(XprojectionC[0,:]), np.median(XprojectionC[1,:]), 
                    label, fontsize=fontsize, ha='center',va='center')
                text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
                texts.append(text)

        adjust_text(texts, # arrowprops=dict(arrowstyle='-', color='k', lw=0.3, alpha=0.75),
                    expand_text=(0.9, 0.9), expand_points=(0.91, 0.9),
                    force_text=(0.01, 0.01), force_points=(0.01, 0.01))

        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.set_xlim([mins[0] - (maxs[0] - mins[0]) * 0.05, (1 + rightShift) * (maxs[0] + (maxs[0] - mins[0]) * 0.05)])
        ax.set_ylim([mins[1] - (maxs[1] - mins[1]) * 0.05, maxs[1] + (maxs[1] - mins[1]) * 0.05])

        if legend:
            plt.legend(loc='lower right', frameon=False, fontsize=fontsize)

        #fig.patch.set_visible(False)
        ax.axis('off')

        if colorbar:
            add_colorbar(fig, possible_cluster_labels, cmap=colormap, fontsize=fontsize)

        self.saveFigure(fig, self.saveDir, '%s_clusters_%s' % (self.dataName, suffix), extension=extension, dpi=dpi, **kwargs)

        return fig

    @tryExcept
    def makeStackedBarplot(self, clusterName = None, legendStyle = False, includeLowQC = True, dpi = 300, extension = 'png', **kwargs):
        
        '''Produce stacked barplot with cell fractions

        Parameters:
            clusterName: str, Deafult None
                Label to include at the bar bottom.
                If None the self.dataName value will be used

            legendStyle: boolean, Default False
                Use one out of two styles of this figure

            includeLowQC: boolean, Default True
                Wether to include low quality cells

            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeStackedBarplot(clusterName)
        '''

        def get_stacked_data_and_colors(saveDir):
            with open(os.path.join(saveDir, 'ColormapForCellTypes.txt'), 'r') as temp_file:
                colors = temp_file.readlines()
                colors = np.vstack([(color.strip('\n').split('\t')) for color in colors])
                colors = pd.DataFrame(colors.T[1], index=colors.T[0]).apply(lambda x: tuple(np.float_(x[0][1:][:-1].split(','))), axis=1)

            df = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), sheet_name='z-scores')
            index = df['Predicted cell type']

            if not clusterName is None:
                barName = self.dataName # + ': ' + clusterName
            else:
                barName = self.dataName

            index = [index[i][:len(index[i]) if index[i].find('#') - 1 == -2 else index[i].find('#') - 1].strip('*').strip('#').strip(' ') for i in range(len(index))]
            df_BM_temp = pd.DataFrame(data=df['# cells in cluster'].values, index=index, columns=[barName])
            df_BM_temp = df_BM_temp.groupby(level=0, axis=0, sort=False).sum()
        
            df_main = pd.DataFrame(data=np.zeros((len(colors),1)), index=colors.index, columns=[barName])

            for i, item in enumerate(df_BM_temp.index):
                df_main.loc[item,barName] += df_BM_temp.iloc[i][barName] if df_BM_temp.index[i] == item else 0

            s = 'sums'
            df_main[s] = np.array(np.sum(df_main, axis=1))
            df_main.loc[self.nameForUnknown, s] = 0

            if includeLowQC:
                try:
                    cells_all = pd.read_hdf(self.fileHDFpath, key='df_projection_pre_QC', mode='r').columns.get_level_values('cell')
                    cells_high = pd.read_hdf(self.fileHDFpath, key='df_projection', mode='r').columns.get_level_values('cell')
                    cells_low_count = len(cells_all.difference(cells_high))

                    del cells_all, cells_high

                    df_main.loc[self.nameForLowQC, df_main.columns[0]] = cells_low_count
                    df_main.loc[self.nameForLowQC, df_main.columns[1]] = -1

                    colors[self.nameForLowQC] = (0.6, 0.6, 0.6, 1.)
                except Exception as exception:
                    print(exception)
                    print('QC data not found')

            df_main = df_main.apply(lambda x: 100. * x / np.sum(df_main, axis=0), axis=1).loc[np.sum(df_main, axis=1) > 0].sort_values(by=[s]).drop(columns=[s])

            return df_main, colors, clusterName

        if clusterName is None:
            clusterName = self.dataName

        df_Main, colors, clusterName = get_stacked_data_and_colors(self.saveDir)

        if legendStyle:
            fig,ax = plt.subplots(figsize=(4.5,8)) #4.15
        else:
            fig = plt.figure(figsize=(4.5,8))
            ax = fig.add_axes([0.2, 0.05, 0.1, 0.9])

        barWidth = 1.0
        cellTypes = df_Main.index
        bottom = np.zeros((len(df_Main.columns)))

        centers = []
        fractions = []
        for i in range(len(cellTypes)):
            bottom += df_Main.loc[cellTypes[i - 1]].values if i > 0 else 0
            ax.bar(range(len(df_Main.columns)), list(df_Main.loc[cellTypes[i]]), bottom=list(bottom), color=colors.loc[cellTypes[i]], edgecolor='white', width=barWidth, label=cellTypes[i])

            centers.append(bottom + 0.5 * df_Main.loc[cellTypes[i]].values[0])
            fractions.append(df_Main.loc[cellTypes[i]].values[0])
 
        plt.xticks(range(len(df_Main.columns)), list(df_Main.columns), fontsize=12)
        plt.yticks([0,20,40,60,80,100], ['0','20%','40%','60%','80%','100%'], fontsize=12)
            
        handles, labels = ax.get_legend_handles_labels()
        ms = np.max([len(item) for item in labels]) - len('cell')
        labels = [item.replace(' ','\n').replace('\nCD4', ' CD4').replace('CD4\n', 'CD4 ').replace('\ncell', ' cell').replace('B\n', 'B ').replace('T\n', 'T ') if len(item) >= ms else item for item in labels[::-1]]

        if legendStyle:
            ax.legend(handles[::-1], labels, loc='upper left', bbox_to_anchor=(1,1), ncol=1, frameon=False, fontsize=14, labelspacing=1, title = ''.join([' ' for _ in range(60)]))
        else:
            fractions = np.round(np.array(fractions)[::-1], 1)
            centers = np.round(np.array(centers).T[0][::-1], 0)
            centers_orig = centers.copy()
            step = 5.

            for i in range(len(centers) - 2,0,-1):
                if (centers[i] - centers[i + 1]) < step:
                    centers[i] = centers[i + 1] + step
            
            for i in range(len(centers)):
                ax.text(1.3, centers[i], '%s%%  ' % (fractions[i]) + labels[i], fontsize=12, va='center', ha='left')
                ax.plot([0.65, 1.2], [centers_orig[i], centers[i]], c='k', lw=0.75, clip_on=False)

        plt.xlim((-0.5, len(df_Main.columns) - 0.5))
        plt.ylim((0, 100))

        for spine in plt.gca().spines.values():
            spine.set_visible(False)
          
        if legendStyle:
            fig.tight_layout()

        saveName = "%s_subclustering_stacked_barplot_%s" % (self.dataName, ('All cell clusters' if clusterName == None else clusterName).replace(' ', '_').replace('*', ''))
        self.saveFigure(fig, self.saveDir, saveName, extension=extension, dpi=dpi, **kwargs)

        print('Saved stacked bar plot: %s' % ('All cell clusters' if clusterName == None else clusterName))

        return fig
    
    @tryExcept
    def makeQualityControlHistogramPlot(self, subset, cutoff, plotPathAndName = None, N_bins = 100, mito = False, displayMeasures = True, precision = 4, quantilePlotCutoff = 0.95, dpi = 300, extension = 'png', fontScale = 1.5, includeTitle = False, **kwargs):

        '''Function to calculate QC quality cutoff and visualize it on a histogram

        Parameters:
            subset: pandas.Series
                Data to analyze

            cutoff: float
                Cutoff to display

            plotPathAndName: str, Default None
                Text to include in the figure title and file name

            N_bins: int, Default 100
                Number of bins of the histogram

            mito: boolean, Default False
                Whether the analysis of mitochondrial genes fraction

            displayMeasures: boolean, Default True
                Print vertical dashed lines along with mean, median, and standard deviation

            precision: int, Default 4
                Number of digits after decimal

            quantilePlotCutoff: float, Default 0.99
                Distributions are cut to display the range from 0 to quantilePlotCutoff

            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

            fontScale: float, Default 1.5
                Scale most of the figure fonts

            includeTitle: boolean, Default False
                Whether to include title on the figure

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            cutoff = DCS.makeQualityControlHistogramPlot(subset, cutoff)
        '''

        if plotPathAndName is None:
            plotPathAndName = 'QC_Plot'

        range_min = 0. #np.min(subset)

        if mito:
            range_max = max(1.1 * cutoff, np.quantile(subset, quantilePlotCutoff) + 0.05)
        else:
            range_max = np.quantile(subset, quantilePlotCutoff)

        hist_of_subset = scipy.stats.rv_histogram(np.histogram(subset, bins=N_bins, range=(range_min, range_max)))
        hist_data = hist_of_subset._hpdf / N_bins
        hist_bins = hist_of_subset._hbins

        fig, ax = plt.subplots(figsize=(8,8))

        bar_bin_width = range_max / N_bins
        ax.bar(hist_bins, hist_data[:-1], width=0.9 * bar_bin_width, color='b', align='center')

        try:
            title = os.path.basename(plotPathAndName)
        except Exception as exception:
            print(exception)
            title = plotPathAndName

        if includeTitle:
            ax.set_title(title, fontdict={'color': 'b'})

        ax.set_xlabel('Fraction' if mito else 'Count', fontsize=10 * fontScale)
        ax.set_ylabel('Density', fontsize=10 * fontScale)
        ax.set_ylim(0.,ax.get_ylim()[1])
        ax.set_xlim(range_min - 0.5 * bar_bin_width, range_max + 0.5 * bar_bin_width)

        ax.tick_params(labelsize=8 * fontScale)

        xs = np.linspace(hist_bins[0], hist_bins[-1], 1000)
        spline_data = np.vstack((xs, UnivariateSpline(hist_bins, hist_data[:-1], k=5, s=0)(xs))).T

        sg = scipy.signal.savgol_filter(spline_data.T[1], 101, 3)
        ax.plot(spline_data.T[0], sg, 'r', lw=3, alpha=0.95)

        try:
            x, y = cutoff, sg[np.where(spline_data.T[0] >= cutoff)[0][0]]
        except Exception as exception:
            print(exception)
            x, y = cutoff, 0.

        ax.plot([x,x], [0,y], 'k', lw=2)
        ax.plot(x, y, 'ko', ms=10, alpha=0.8)
        ax.plot(x, y, 'ro', ms=7)

        ax.text(x, -0.04 * spline_data.T[1].max(), str(np.round(cutoff, precision)), fontsize=8 * fontScale, va='top', ha='center', color='r')

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=False)

        ax.axvspan(cutoff, 1.5 * range_max if mito else -1.5 * range_min, alpha=0.1, color='red', hatch='\\', linewidth=0.1)

        fig.tight_layout()

        if displayMeasures:
            texts = []

            dist_std, dist_median, dist_mean = np.round(np.std(subset),precision), np.round(np.median(subset),precision), np.round(np.mean(subset),precision)
            print(plotPathAndName, '\tstd:', dist_std,  '\tmedian:', dist_median,  '\tmean:', dist_mean)

            xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
            yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
        
            ax.axvline(x=dist_mean, color='k', lw=1.0, ls='--')
            text = ax.text(dist_mean + 0.02 * xspan, 0.98 * yspan, r'$\mu=%s$' % (dist_mean), fontsize=fontScale * 10, va='top', ha='left', color='k')
            text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),path_effects.Normal()])
            texts.append(text)

            ax.axvline(x=dist_median, color='k', lw=1.0, ls='--')
            text = ax.text(dist_median + 0.02 * xspan, 0.94 * yspan, r'$M=%s$' % (dist_median), fontsize=fontScale * 10, va='top', ha='left', color='k')
            text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),path_effects.Normal()])
            texts.append(text)
        
            ax.axvline(x=dist_median - dist_std, color='k', lw=1.0, ls='--')
            text = ax.text(dist_median - dist_std + 0.02 * xspan, 0.90 * yspan, r'$M-\sigma=%s$' % (np.round(dist_median - dist_std,precision)), fontsize=fontScale * 10, va='top', ha='left', color='k')
            text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),path_effects.Normal()])
            texts.append(text)
        
            ax.axvline(x=dist_median + dist_std, color='k', lw=1.0, ls='--')
            text = ax.text(dist_median + dist_std + 0.02 * xspan, 0.90 * yspan, r'$M+\sigma=%s$' % (np.round(dist_median + dist_std,precision)), fontsize=fontScale * 10, va='top', ha='left', color='k')
            text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),path_effects.Normal()])
            texts.append(text)

            text = ax.text(dist_median + 0.02 * xspan, 0.76 * yspan, r'$\sigma=%s$' % (np.round(dist_std,precision)), fontsize=fontScale * 10, va='bottom', ha='left', color='k')
            text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),path_effects.Normal()])
            texts.append(text)

            ax.annotate('', (dist_median + dist_std, 0.75 * yspan), (dist_median, 0.75 * yspan), arrowprops={'arrowstyle':'<|-|>'})

            if not mito:
                text = ax.text(0.98, 0.65, 
                               '%s%%\nof distribution \nis shown' % (100. * quantilePlotCutoff), 
                               va='top', ha='right', fontsize=fontScale * 10, transform=ax.transAxes)
                text.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='white'),path_effects.Normal()])
                
            adjust_text(texts)

        self.saveFigure(fig, os.path.dirname(plotPathAndName), label=os.path.basename(plotPathAndName) + '_histogram', extension=extension, dpi=dpi, **kwargs)

        return fig

    @tryExcept
    def makePlotOfNewMarkers(self, df_marker_cell_type, df_new_marker_cell_type, dpi = 300, extension = 'png', **kwargs):

        '''Produce plot of the new markers extracted from the annotated clusters

        Parameters:
            df_marker_cell_type: pandas.DataFrame
                Known markers per cell types

            df_new_marker_cell_type: pandas.DataFrame
                New markers per cell types

            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makePlotOfNewMarkers(df_marker_cell_type, df_new_marker_cell_type)
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
                known_markers = df_marker_cell_type[celltype][df_marker_cell_type[celltype] > 0.].index.values
                xy = np.array([np.array([np.where(genes == marker)[0][0], i]) for marker in known_markers if marker in genes])
                print('Overlapping positive markers of %s: %s (%s)' % (celltype, len(xy), len(known_markers)))
                if len(xy) > 0:
                    ax.plot(xy.T[0], xy.T[1], 'go', markeredgecolor='r', ms=1.0, markeredgewidth=0.2)

                known_markers = df_marker_cell_type[celltype][df_marker_cell_type[celltype] < 0.].index.values
                xy = np.array([np.array([np.where(genes == marker)[0][0], i]) for marker in known_markers if marker in genes])
                print('Overlapping negative markers of %s: %s (%s)' % (celltype, len(xy), len(known_markers)))
                if len(xy) > 0:
                    ax.plot(xy.T[0], xy.T[1], 'ro', markeredgecolor='r', ms=1.0, markeredgewidth=0.2)

        #ax.set_title('Additional markers along with the overlapping part of the input (red)')
        self.saveFigure(fig, self.saveDir, self.dataName + '_new_markers', extension=extension, dpi=dpi, **kwargs)

        return fig

    @tryExcept
    def makeTtestPlot(self, df, dfp, label = None, reorder = True, p_value_cutoff = 0.05, dpi = 300, extension = 'png', **kwargs):

        '''Produce heatmap plot of t-test p-Values calculated gene-pair-wise
        from the annotated clusters.

        Parameters:
            df: pandas.DataFrame
                t-test statistic values
    
            dfp: pandas.DataFrame
                t-test p-Values calculated gene-pair-wise

            label: str, Default None
                Lebel to include in the plot
                
            reorder: boolean, Default True
                Reorder values to group similar

            p_value_cutoff: float, Default 0.05
                p-Value cutoff

            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeTtestPlot(df)
        '''

        if reorder:

            def metricCommonEuclidean(u,v):

                where_common = (~np.isnan(u)) * (~np.isnan(v))

                return np.sqrt(((u[where_common] - v[where_common]) ** 2).sum())

            order = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(df.values, method='average', metric=metricCommonEuclidean), no_plot=True, get_leaves=True)['leaves']

            df = df[df.columns.values[order]]
            dfp = dfp[dfp.columns.values[order]]
            df = df.loc[df.index.values[order]]
            dfp = dfp.loc[dfp.index.values[order]]

        df = df[df.columns[::-1]]
        dfp = dfp[dfp.columns[::-1]]

        fig = plt.figure(figsize=(5,5))

        ax = fig.add_axes([0.35,0.02,0.6,0.6])
        
        cmap = plt.cm.PuOr_r #BrBG #PiYG #seismic
        cmap.set_bad('grey')

        ax.imshow(df.values.astype(float), cmap=cmap, interpolation='None', aspect='auto')

        wh = np.where(dfp.values.T <= p_value_cutoff)
        ax.plot(wh[0], wh[1], '*k')

        ax.set_xticks(range(df.shape[1]))
        ax.set_yticks(range(df.shape[0]))
        ax.set_xticklabels(df.columns.values, rotation=90, fontsize=8)
        ax.set_yticklabels(df.index.values, rotation=0, fontsize=8)
        ax.set_xlim([-0.5, df.shape[1] - 0.5])
        ax.set_ylim([-0.5, df.shape[0] - 0.5])

        ax.xaxis.tick_top()

        if not label is None:
            ax.text(-0.5, 1.5, label, transform=ax.transAxes, fontsize=10, color='k', ha='left', va='top').set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='blue'),path_effects.Normal()])

        ax.set_title('Two-tailed p-Value (t-test)')

        data = df.values.flatten().astype(float)
        data = data[np.where(~np.isnan(data))]
        dataMin = np.min(data)
        dataMax = np.max(data)

        axisColor = fig.add_axes([0.22,0.75,0.08,0.02])

        norm = matplotlib.colors.Normalize(vmin=dataMin, vmax=dataMax)
        mapp = cm.ScalarMappable(norm=norm, cmap=cmap)
        mapp.set_array(data)
        fig.colorbar(mapp, cax=axisColor, ticks=[dataMax,dataMin], orientation='horizontal')
        axisColor.tick_params(labelsize=4)
        axisColor.set_xlabel('Statistic\n*p-Value < %s' % (p_value_cutoff), fontsize=5)

        axisColor.set_yticklabels([np.round(dataMax,2), np.round(dataMin,2)])

        self.saveFigure(fig, self.saveDir, self.dataName + '_ttest_%s' % (label.replace('\n', '_')), extension=extension, dpi=dpi, **kwargs)

        return fig
    
    @tryExcept
    def makeCellMarkersPiePlot(self, type1, type2, df_marker_cell_type = 'all', nameToAppend = None, listUnexpressedMarkers = True, orthogonalSectorsShift = 0.1, rotationAngle = 0, dpi = 300, extension = 'png', **kwargs):

        '''Make summary of markers comparison between two cell types.

        Parameters:
            type1: str
                Name of the first cell type to compare
                            
            type2: str
                Name of the second cell type to compare

            df_marker_cell_type: pandas.DataFrame or str, Default 'all'
                Celltypes/Markers matrix. If 'expressed', then only expressed markers will be used.
                If 'all' then all markers of the input marker list will be used.
                If an instance of a pandas.DataFrame is passed, then its all markers will be used.

            nameToAppend: str, Default None
                String to append to the figure file name.

            listUnexpressedMarkers: boolean, Default True
                List (highlight) markers that are not expressed. This option is ignored
                unless df_marker_cell_type=='all'

            orthogonalSectorsShift: float, Default 0.1
                Sectors marked as '+/-' and '-/+' are shifted off-center.
                Set this parameter to zero to have round continuous pie chart.

            rotationAngle: int or float, Default 0
                Angle in degrees that will rotate the whole pie chart counterclockwise.

            dpi: int, Default 600
                Resolution of the figure image

            extension: str, Default 'png'
                Format of the figure file

        Returns:
            Marker lists split into categories.
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeCellMarkersPiePlot('T cells', 'B cells')
        '''

        try:
            df_marker_expression = pd.read_excel(os.path.join(self.saveDir, self.dataName + '_annotation.xlsx'), 
                                                    sheet_name='Marker cell type weight matrix', index_col=0, header=0).T
        except Exception as exception:
            print(exception)
            print('Marker expression data unavailable')
            df_marker_expression = None
            listUnexpressedMarkers = False

        if type(df_marker_cell_type) is str:
            if df_marker_cell_type == 'expressed':
                listUnexpressedMarkers = False

                if not df_marker_expression is None:
                    df_marker_cell_type = df_marker_expression
                else:
                    print("Try using option 'all'")

                    return

                additional_name = 'expressed'
            elif df_marker_cell_type == 'all':
                df_marker_cell_type = self.readMarkerFile()
                additional_name = 'all'

            if nameToAppend is None:
                nameToAppend = '.'.join(os.path.basename(self.geneListFileName).split('.')[:-1])
        else:
            additional_name = 'custom'
            
            if nameToAppend is None:
                nameToAppend = ''

        df_marker_cell_type = df_marker_cell_type.fillna(0.)

        def getSet(df, celltype):

            try:
                pos = set(df.index[(df.loc[:, celltype] > 0.)].values)
                neg = set(df.index[(df.loc[:, celltype] < 0.)].values)
            except Exception as exception:
                print(exception)
                print('Cell type %s not found' % (celltype))
                print('Available celltypes are: %s' % (df.columns.values.tolist()))

                return

            return pos, neg

        def getEightAll(t1p, t2p, t1n, t2n):

            p1 = t1p.intersection(t2n)
            p3 = t1n.intersection(t2n)
            p5 = t1n.intersection(t2p)
            p7 = t1p.intersection(t2p)

            p0 = t1p - p1 - p7
            p2 = t2n - p1 - p3
            p4 = t1n - p3 - p5
            p6 = t2p - p5 - p7

            all = t1p.union(t2p).union(t1n).union(t2n)

            return [p0, p1, p2, p3, p4, p5, p6, p7], all
        
        try:
            t1p, t1n = getSet(df_marker_cell_type, type1)
            t2p, t2n = getSet(df_marker_cell_type, type2)

            sets, all = getEightAll(t1p, t2p, t1n, t2n)
        except:
            return

        if listUnexpressedMarkers:
            try:
                t1pe, t1ne = getSet(df_marker_expression, type1)
                t2pe, t2ne = getSet(df_marker_expression, type2)
                
                setsE, allE = getEightAll(t1pe, t2pe, t1ne, t2ne)
            except Exception as exception:
                print(exception)
                listUnexpressedMarkers = False

        labels = '+/*', '+/-', '*/-', '-/-', '-/*', '-/+', '*/+', '+/+'
        colors = ['limegreen', 'thistle', 'lightcoral', 'red', 'lightcoral', 'thistle', 'limegreen', 'green']
        titles = ['Positive in %s:' % (type1),
                  'Positive in %s, Negative in %s:' % (type1, type2),
                  'Negative in %s:' % (type2),
                  'Negative in both:',
                  'Negative in %s:' % (type1),
                  'Negative in %s, Positive in %s:' % (type1, type2),
                  'Positive in %s:' % (type2),
                  'Positive in both:']
        sizes = [len(item) for item in sets]
        labels = [(label if size > 0 else '') for label, size in zip(labels, sizes)]
        explode = (0.0, orthogonalSectorsShift, 0.0, 0.0, 0.0, orthogonalSectorsShift, 0.0, 0.0)

        def findAll(a, b):

            start = 0

            while True:
                start = a.find(b, start)

                if start == -1: 
                    return

                yield start

                start += len(b)

            return

        if listUnexpressedMarkers:
            str_sets = [str(sorted([i for i in list(setsE[j])])).replace("'", "").replace(']','').replace('[','').replace(' ','') for j in range(8)]
        else:
            str_sets = [str(sorted([i for i in list(sets[j])])).replace("'", "").replace(']','').replace('[','').replace(' ','') for j in range(8)]

        for i, item in enumerate(str_sets):

            all_temp = item.split(',')
            new_item = all_temp[0]
            temp_item = all_temp[0]
            limit = 75

            for gene in all_temp[1:]:
                if len(temp_item) > limit:
                    temp_item = '\n' + gene
                    new_item += '\n' + gene
                else:
                    new_item += ', ' + gene
                    temp_item += ', ' + gene

            if listUnexpressedMarkers:
                str_sets[i] = titles[i] + ' (%s):' % (len(all_temp)) + '\n' + new_item
            else:
                str_sets[i] = titles[i] + '\n' + new_item

        if listUnexpressedMarkers:
            str_setsU = [str(sorted([i for i in list(sets[j].difference(setsE[j]))])).replace("'", "").replace(']','').replace('[','').replace(' ','') for j in range(8)]
        
            for i, item in enumerate(str_setsU):

                all_temp = item.split(',')
                new_item = all_temp[0]
                temp_item = all_temp[0]
                limit = 75

                for gene in all_temp[1:]:
                    if len(temp_item) > limit:
                        temp_item = '\n' + gene
                        new_item += '\n' + gene
                    else:
                        new_item += ', ' + gene
                        temp_item += ', ' + gene

                if len(new_item) > 1:
                    str_sets[i] += '\n' + 'Not expressed (%s):' % (len(all_temp)) + '\n' + new_item

        fig = plt.figure(figsize=(8,4))
        ax = fig.add_axes([0.25,0.25,0.5,0.5])

        currentWedge = 0

        def autopctFunc(value):

            nonlocal currentWedge

            n = int(np.round((float(value) / 100. * float(np.sum(sizes))), 0))

            if listUnexpressedMarkers:
                u = len(sets[currentWedge].difference(setsE[currentWedge]))
            else:
                u = 0

            if n > 0:
                if u > 0:
                    format = "{:d}\n({:d})".format(n,u)
                else:
                    format = "{:d}".format(n)
            else:
                format = ""

            currentWedge += 1

            return format

        wedges, texts, autotexts = ax.pie(sizes, explode=explode, labels=labels, colors=colors, labeldistance=1.05, 
                                          textprops={'size':6, 'weight':'semibold', 'color':'b'},
                                          autopct=autopctFunc, wedgeprops={'linewidth': 0.5, 'edgecolor':'aqua', 'width': 0.7},
                                          shadow=False, startangle=-180 + rotationAngle, frame=False, rotatelabels=False)

        plt.setp(autotexts, size=6, weight="semibold", color='k')

        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.5)
        kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

        for i, p in enumerate(wedges):

            if len(sets[i]) == 0:
                continue

            ang = (p.theta2 - p.theta1) / 2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))

            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]

            connectionstyle = "angle,angleA=0,angleB={}".format(ang)

            kw["arrowprops"].update({"connectionstyle": connectionstyle})

            ax.annotate(str_sets[i], xy=(x, y), xytext=(1.6 * np.sign(x), 1.8 * y), fontsize=3.5, ha=horizontalalignment, **kw)

        fig.suptitle('%s & %s' % (type1, type2), fontsize=11, color='b')
        ax.axis('equal')

        ax.text(0., 0., '%s\n(%s)' % (len(all), len(all.difference(allE))) if listUnexpressedMarkers else '%s' % (len(all)), 
                color='k', fontsize=8, ha='center', va='center').set_path_effects([path_effects.Stroke(linewidth=1, foreground='blue'),path_effects.Normal()])

        saveName = 'Markers_of_%s_vs_%s_(%s)_%s' % (type1.replace('/',''), type2.replace('/',''), nameToAppend, additional_name)
        self.saveFigure(fig, self.saveDir, self.dataName + saveName, extension=extension, dpi=dpi, **kwargs)

        return dict(zip(labels, list(sets))), fig

    @tryExcept
    def makeHopfieldPCplot(self, colormap = cm.hot_r, plotTrLines = False, clusterid = 1, trID = 0, axisOff = False, fontscale = 1., trPath = None, dpi = 300, extension = 'png', **kwargs):

        '''Make radar plot of the attractors in their principal components coordinates

        Parameters:
            colormap: matplotlib.colormap or str, Default cm.hot_r
                Colormap or its string name

            plotTrLines: boolean, Default False
                Whether to plot trajectories

            clusterid: int, Default 1
                Identifier of the cluster to plot trajectories of

            trID: int, Default 0
                Identifier of the trajectories to plot

            axisOff: boolean, Default False
                Whether to hide the axes lines

            trPath: str, Default None
                Path to trajectories files

            dpi: int, Default 300
                Resolution of the figure

            extension: str, Default 'png'
                Format extension of the figure

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeHopfieldPCplot()
        '''

        if axisOff:
            ax.axis('off')

        if trPath is None:
            trPath = os.path.join(self.saveDir, 'HopfieldTrajectories')

        if not os.path.exists(trPath):
            print('Data not found', flush=True)

            return

        fig = plt.figure(figsize=(8,8))
        ax = plt.subplot(111, polar=True, theta_direction=-1, theta_offset=0.5*np.pi)

        attrs, attrs_names = read(os.path.join(trPath, 'attrs'))
        N = attrs.shape[1]
        df = pd.DataFrame(data=attrs[:N], index=attrs_names, 
                          columns=['PC%s\n%s%%'%(i+1, np.int(100.*attrs[N][i])) for i in range(N)])

        wherePC = attrs[N] > 0.001
        df = df[df.columns[wherePC]]
        N = df.shape[1]

        vmaxAt = df.max(axis=0).max()
        vminAt = df.min(axis=0).min()

        ax.set_ylim(vminAt, vmaxAt)
 
        angles = [n / float(N) * 2 * np.pi for n in range(N)] + [0.]
 
        ax.set_xticklabels([])
        ax.set_xticks(angles)

        for i, celltype in enumerate(df.index):
            values=df.loc[celltype].values.flatten().tolist()
            values.append(values[0])

            color = cm.jet(i/len(attrs_names))

            ax.plot(angles, values, color=color, linewidth=1.75, linestyle='solid', label=celltype)
            #ax.fill(angles, values, alpha=0.2, color=color, zorder=1)
            ax.fill_between(angles, 0, values, alpha=0.2, facecolor=color)

            temp_texts = ax.text(angles[np.argmax(values)], values[np.argmax(values)], celltype, color=color, fontsize=12.*fontscale, ha='center', va='center')
            temp_texts.set_path_effects([path_effects.Stroke(linewidth=1., foreground='k'), path_effects.Normal()])

        if plotTrLines:
            trajectories = read(os.path.join(trPath, 'trajectories%s')%(trID))

            initial, final, typesNames, clusterNames = read(os.path.join(trPath, 'additional'))
        
            thisTr = trajectories[:, clusterid, wherePC].T

            vmax = max(thisTr.max(axis=0).max(axis=0), vmaxAt)
            vmin = min(thisTr.min(axis=0).min(axis=0), vminAt)

            ax.set_ylim(vmin, vmax)

            inId = initial[clusterid]
            outId = final[clusterid]

            print('ClusterID: %s, Initial state: %s (%s), Final state: %s (%s)'%(clusterNames[clusterid], typesNames[inId], inId, typesNames[outId], outId))

            suffix = clusterNames[clusterid]

            fig.suptitle('Cluster %s: %s -> %s' % (clusterNames[clusterid], typesNames[inId] if inId!=-1 else 'Unknown', typesNames[outId] if outId!=-1 else 'Unknown'), fontsize=11, color='b')

            values = thisTr.T[0].tolist() + [thisTr.T[0][0]]

            timeLimit = thisTr.shape[1]

            for t in range(timeLimit):
                values = thisTr.T[t].tolist() + [thisTr.T[t][0]]

                if t==0:
                    ax.plot(angles, values, 'o-', ms=4.0, color='b', alpha=1.0, clip_on=False)

                ax.plot(angles, values, color='b', linewidth=0.5, alpha=0.04, linestyle='solid')
                ax.fill(angles, values, alpha=0.005, color='b')

                if t == timeLimit - 1:
                    ax.plot(angles, values, 'o-', ms=4.0, color='r', alpha=1.0, clip_on=False)

            ax.set_rlabel_position(0)
        else:
            vmax = vmaxAt
            vmin = vminAt
            suffix = 'attractors'

        for i, pc in enumerate(df.columns):
            temp_texts = ax.text(angles[i], 1.15 * (vmax - vmin) + vmin, pc, color='k', fontsize=14*fontscale, ha='center', va='center')
            temp_texts.set_path_effects([path_effects.Stroke(linewidth=0.5*fontscale, foreground='w'), path_effects.Normal()])

        ax.set_axisbelow(True)

        #ax.plot(angles, [0]*len(angles), color='k', linewidth=1., linestyle='-', label=celltype)

        fig.canvas.draw()
        ylabels = ax.get_yticklabels()
        ax.set_yticklabels([])

        for label in ylabels:
            ax.text(label._x, label._y, label._text, zorder=np.inf)

        self.saveFigure(fig, self.saveDir, self.dataName + '_polar_%s'%(suffix), extension=extension, dpi=dpi, **kwargs)

        return fig

    @tryExcept
    def HopfieldLandscapePlot(self, legend = False, labels = False, PCx = 0, PCy = 1, colorbar = True, fontsize = 10, plotMesh = True, plotAttractors = True, adjustText = True, axisOff = True, colorbarva = 0.75, colorbarha = 0.85, trPath = None, colormap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap', [(1, 1, 1), (0, 1, 1), (0, 0, 1), (1, 0, 0)], N=1000), dpi = 300, extension = 'png', **kwargs):

        '''Make heatmap plot of the attractors in their two principal components coordinates

        Parameters:
            legend: boolean, Default False
                Whether to add legend containing cell types names
            
            labels: boolean, Default False
                Whether to add labels 
                
            PCx: int, Default 0
                Principal component for x-coordinate of the plot
                            
            PCy: int, Default 1
                Principal component for y-coordinate of the plot
            
            colorbar: boolean, Default False
                Whether to add colorbar

            fontsize: int, Default 10
                Text labels font size
            
            plotMesh: boolean, Default False
                Whether to plot landscape heatmap
                
            plotAttractors: boolean, Default False
                Whether to plot attractor stars
                
            adjustText: boolean, Default False
                Whether to minimize text labels overlap
                
            axisOff: boolean, Default False
                Whether to hide the axes lines

            colorbarva: float, Default 0.75
                Vertical position of the bottom of the colorbar

            colorbarha: float, Default 0.85
                Horizontal position of the colorbar

            trPath: str, Default None
                Path to trajectories files

            colormap: matplotlib.colormap or str, Default matplotlib.colors.LinearSegmentedColormap.from_list('cmap', [(1, 1, 1), (0, 1, 1), (0, 0, 1), (1, 0, 0)], N=1000)
                Colormap or its string name

            dpi: int, Default 300
                Resolution of the figure

            extension: str, Default 'png'
                Format extension of the figure

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeHopfieldLandscapePlot()
        '''

        np.random.seed(0)

        colormap.set_bad('white')

        def add_colorbar(fig, labels, cmap = matplotlib.colors.LinearSegmentedColormap.from_list('GR', [(0, 1, 0), (1, 0, 0)], N=100), fontsize = 10):
    
            mapp = cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=np.min(labels), vmax=np.max(labels)), cmap=cmap)
            sp = np.linspace(np.max(labels), np.min(labels), num=6, endpoint=True)
            mapp.set_array(sp)

            axisColor = fig.add_axes([colorbarha, colorbarva, 0.01, 0.2])

            fig.colorbar(mapp, cax=axisColor, ticks=sp)

            axisColor.tick_params(labelsize=fontsize)
            axisColor.set_yticklabels(sp.astype(int))
    
            return None

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.05,0.05,0.9,0.9])

        if axisOff:
            ax.axis('off')
            
        if trPath is None:
            trPath = os.path.join(self.saveDir, 'HopfieldTrajectories')

        if not os.path.exists(trPath):
            print('Data not found', flush=True)

            return

        attrs_xpca, attrs_names = read(os.path.join(trPath, 'attrs'))
        attrs_xpca = attrs_xpca[:attrs_xpca.shape[1]]

        mesh_xpca = read(os.path.join(trPath, 'mesh'))

        mesh_energy = mesh_xpca[range(len(mesh_xpca)), -1]
        mesh_xpca = mesh_xpca[range(len(mesh_xpca)), :-1]

        data = np.vstack([attrs_xpca, mesh_xpca])

        coords = data.T[[PCx, PCy], :]
        attrs2D, mesh2D = coords[:, :attrs_xpca.shape[1]].T, coords[:, attrs_xpca.shape[1]:].T
    
        if plotMesh:
            vmin, vmax = min(mesh_energy), 0.

            vals = mesh_energy.copy()
            vals[np.where(vals > (vmax - 0.001))[0]] = vmax - 0.001
            vals[np.where(vals < (vmin + 0.001))[0]] = vmin + 0.001

            xmin, xmax = mesh2D.T[0].min(), mesh2D.T[0].max()
            ymin, ymax = mesh2D.T[1].min(), mesh2D.T[1].max()

            dx = (xmax - xmin) * 0.05
            dy = (ymax - ymin) * 0.05

            xmin -= dx
            xmax += dx
            ymin -= dy
            ymax += dy

            ngrid = 100
            grid = np.zeros((ngrid + 1, ngrid + 1))
            grid[:] = vmin

            i = (ngrid * (mesh2D.T[0] - xmin) / (xmax - xmin)).astype(int)
            j = (ngrid * (mesh2D.T[1] - ymin) / (ymax - ymin)).astype(int)

            se = pd.Series(index=zip(i,j), data=mesh_energy).groupby(axis=0, level=0).agg(np.min)
            se.index = pd.MultiIndex.from_tuples(se.index)

            grid[(se.index.get_level_values(0).values, se.index.get_level_values(1).values)] = se.values

            maskedArray = np.ma.array(grid.T, mask=np.isnan(grid.T,))

            im = ax.imshow(maskedArray[::-1], vmin=vmin, vmax=vmax, cmap=colormap, alpha=0.8,
                        extent=[xmin, xmax, ymin, ymax], interpolation='quadric', zorder=-10**8, clip_on=False)

            data = scipy.ndimage.gaussian_filter(grid.T, 1.5)
            xgrid = np.linspace(xmin, xmax, num=(ngrid+1))
            ygrid = np.linspace(ymin, ymax, num=(ngrid+1))

            tempColormap = colormap
            #tempColormap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap', [(0.75, 0.75, 0.75), (0, 1, 1), (0, 0, 1), (1, 0, 0)], N=1000)

            ax.contour(xgrid, ygrid, data, levels=10, cmap=tempColormap, linewidths=1.0, zorder=-10**8+1)
            ax.contour(xgrid, ygrid, data, levels=10, colors='k', linestyles='solid', linewidths=0.25, zorder=-10**8+2)
    
        if plotAttractors:
            texts = []

            ax.plot(attrs2D.T[0], attrs2D.T[1], '*', ms=14, color='k', alpha=1.0, zorder=-10**7, clip_on=False)
            for attr in range(attrs2D.T[0].shape[0]):
                temp_texts = ax.text(attrs2D.T[0][attr], attrs2D.T[1][attr], attrs_names[attr], fontsize=fontsize, ha='left',va='center', zorder=10 ** 10, clip_on=False)
                temp_texts.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'), path_effects.Normal()])
                texts.append(temp_texts)

            if adjustText:
                adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', lw=0.3, alpha=0.5), force_text=(0.05, 0.05))

        if plotMesh and colorbar:

            add_colorbar(fig, [vmax, vmin], cmap=colormap, fontsize=fontsize)

        self.saveFigure(fig, self.saveDir, self.dataName + '_energy_landscape_PC%s_vs_PC%s'%(PCy, PCx), extension=extension, dpi=dpi, **kwargs)

        return fig
    

    # Plotly-powered figures

    @tryExcept
    def makeSankeyDiagram(self, df, colormapForIndex = None, colormapForColumns = None, linksColor = 'rgba(100,100,100,0.6)', title = '', attemptSavingHTML = False, quality = 4, width = 400, height = 400, border = 20, nameAppend = '_Sankey_diagram'):

        '''Make a Sankey diagram, also known as 'river plot' with two groups of nodes

        Parameters:
            df: pandas.DataFrame 
                With counts (overlaps)

            colormapForIndex: dictionary, Default None
                Colors to use for nodes specified in the DataFrame index

            colormapForColumns: dictionary, Default None
                Colors to use for nodes specified in the DataFrame columns

            linksColor: str, Default 'rgba(100,100,100,0.6)'
                Color of the non-overlapping links

            title: str, Default ''
                Title to print on the diagram

            interactive: boolean , Default False
                Whether to launch interactive JavaScript-based graph

            quality: int, Default 4
                Proportional to the resolution of the figure to save

            nameAppend: str, Default '_Sankey_diagram'
                Name to append to the figure file

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeSankeyDiagram(df)
        '''

        try:
            temp_index = pd.MultiIndex.from_arrays([df.index, [colormapForIndex[item] for item in df.index]], names=['label', 'color'])
            temp_columns = pd.MultiIndex.from_arrays([df.columns, [colormapForColumns[item] for item in df.columns]], names=['label', 'color'])
            df.index = temp_index
            df.columns = temp_columns
        except Exception as exception:
            print(exception)
            print('Using default node colors')
            colormapForIndex = None
            colormapForColumns = None

        if (colormapForIndex is None) or (colormapForColumns is None):
            nodeColors = ['rgba(150,0,10,0.8)'] * len(df.index) + ['rgba(10,0,150,0.8)'] * len(df.columns)
            nodeLabels = df.index.to_list() + df.columns.to_list()
        else:
            nodeLabels = df.index.get_level_values('label').to_list() + df.columns.get_level_values('label').to_list()
            nodeColors = df.index.get_level_values('color').to_list() + df.columns.get_level_values('color').to_list()

        sources, targets, values, labels = [], [], [], []
        for i, item in enumerate(df.index):
            sources.extend([i] * len(df.loc[item]))
            targets.extend(list(range(len(df.index), len(df.index) + len(df.loc[item]))))
            values.extend([j for j in df.loc[item].values])
            if type(item) is tuple:
                labels.extend([str(item[0]) + ' -> ' + str(jtem[0]) for jtem in df.loc[item].index])
            else:
                labels.extend([str(item) + ' -> ' + str(jtem) for jtem in df.loc[item].index])

        colorscales = [dict(label=label, colorscale=[[0, linksColor], [1, linksColor]]) for label in labels]

        if not nodeColors is None:
            for i in range(len(sources)):
                if nodeColors[sources[i]] == nodeColors[targets[i]]:
                    newColor = ','.join(nodeColors[sources[i]].split(',')[:3] + ['0.6)'])
                    colorscales[i] = dict(label=labels[i], colorscale=[[0, newColor], [1, newColor]])

        fig = go.Figure(data=[go.Sankey(valueformat = '', valuesuffix = '',
            node = dict(pad = 20, thickness = 40, line = dict(color = 'white', width = 0.5), label = nodeLabels, color = nodeColors,),
            link = dict(source = sources, target = targets, value = values, label = labels, colorscales = colorscales, hoverinfo='all'))]) #line ={'color':'rgba(255,0,0,0.8)', 'width':0.1}

        if not title is None:
            fig.update_layout(title_text=title, font_size=10)

        fig.update_layout(margin=dict(l=border, r=border, t=border, b=border))

        try:
            fig.write_image(os.path.join(self.saveDir, self.dataName + nameAppend + '.png'), width=width, height=height, scale=quality)

        except Exception as exception:
            print('Cannot save static image (likely due to missing orca). Saving to interactive html')
            attemptSavingHTML = True

        if attemptSavingHTML:
            fig.update_layout(margin=dict(l=200, r=200, t=100, b=100))
            plot_offline(fig, filename=os.path.join(self.saveDir, self.dataName + nameAppend + '.html'), auto_open=False)

        return fig

    @tryExcept
    def HopfieldLandscapePlot3D(self, PCx = 0, PCy = 1, colorbar = True, fontsize = 12, plotMesh = True, plotAttractors = True, trPath = None, attemptSavingHTML=False, nameAppend = '', quality = 4, **kwargs):

        '''Make heatmap plot of the attractors in their two principal components coordinates

        Parameters:
            legend: boolean, Default False
                Whether to add legend containing cell types names
            
            labels: boolean, Default False
                Whether to add labels 
                
            PCx: int, Default 0
                Principal component for x-coordinate of the plot
                            
            PCy: int, Default 1
                Principal component for y-coordinate of the plot
            
            colorbar: boolean, Default False
                Whether to add colorbar

            fontsize: int, Default 10
                Text labels font size
            
            plotMesh: boolean, Default False
                Whether to plot landscape heatmap
                
            plotAttractors: boolean, Default False
                Whether to plot attractor stars
                
            adjustText: boolean, Default False
                Whether to minimize text labels overlap
                
            axisOff: boolean, Default False
                Whether to hide the axes lines

            colorbarva: float, Default 0.75
                Vertical position of the bottom of the colorbar

            colorbarha: float, Default 0.85
                Horizontal position of the colorbar

            trPath: str, Default None
                Path to trajectories files

            colormap: matplotlib.colormap or str, Default matplotlib.colors.LinearSegmentedColormap.from_list('cmap', [(1, 1, 1), (0, 1, 1), (0, 0, 1), (1, 0, 0)], N=1000)
                Colormap or its string name

            dpi: int, Default 300
                Resolution of the figure

            extension: str, Default 'png'
                Format extension of the figure

        Returns:
            None
        
        Usage:
            DCS = DigitalCellSorter.DigitalCellSorter()

            DCS.makeHopfieldLandscapePlot()
        '''
           
        if trPath is None:
            trPath = os.path.join(self.saveDir, 'HopfieldTrajectories')

        if not os.path.exists(trPath):
            print('Data not found', flush=True)

            return

        attrs_xpca, attrs_names = read(os.path.join(trPath, 'attrs'))
        attrs_xpca = attrs_xpca[:attrs_xpca.shape[1]]

        mesh_xpca = read(os.path.join(trPath, 'mesh'))

        mesh_energy = mesh_xpca[range(len(mesh_xpca)), -1]
        mesh_xpca = mesh_xpca[range(len(mesh_xpca)), :-1]

        data = np.vstack([attrs_xpca, mesh_xpca])

        coords = data.T[[PCx, PCy], :]
        attrs2D, mesh2D = coords[:, :attrs_xpca.shape[1]].T, coords[:, attrs_xpca.shape[1]:].T
    
        vmin, vmax = min(mesh_energy), 0.

        vals = mesh_energy.copy()
        vals[np.where(vals > (vmax - 0.001))[0]] = vmax - 0.001
        vals[np.where(vals < (vmin + 0.001))[0]] = vmin + 0.001

        xmin, xmax = mesh2D.T[0].min(), mesh2D.T[0].max()
        ymin, ymax = mesh2D.T[1].min(), mesh2D.T[1].max()

        dx = (xmax - xmin) * 0.05
        dy = (ymax - ymin) * 0.05

        xmin -= dx
        xmax += dx
        ymin -= dy
        ymax += dy

        ngrid = 100
        grid = np.zeros((ngrid + 1, ngrid + 1))
        grid[:] = vmin

        i = (ngrid * (mesh2D.T[0] - xmin) / (xmax - xmin)).astype(int)
        j = (ngrid * (mesh2D.T[1] - ymin) / (ymax - ymin)).astype(int)

        se = pd.Series(index=zip(i,j), data=mesh_energy).groupby(axis=0, level=0).agg(np.min)
        se.index = pd.MultiIndex.from_tuples(se.index)

        grid[(se.index.get_level_values(0).values, se.index.get_level_values(1).values)] = se.values

        df = se.unstack(fill_value=vmin)

        fig = go.Figure()

        if plotMesh:
            fig.add_trace(go.Surface(x=np.linspace(xmin, xmax, df.shape[0]), 
                                        y=np.linspace(ymin, ymax, df.shape[1]), 
                                        z=gaussian_filter(df.values.T, sigma=0.75), opacity=1., colorscale="blackbody_r",
                                        showscale=colorbar,
                                        hoverinfo='none',
                                        contours= {'x': {'highlight': False}, 
                                                'y': {'highlight': False}, 
                                                'z': {'highlight': False}},))

            fig.update_traces(contours_z=dict(show=True, width=3., highlightwidth=3., usecolormap=False, highlightcolor="limegreen", project=dict(x=True,y=True,z=True), highlight=True, color='grey', size=(vmax-vmin)/10.))

        annotations = []
        if plotAttractors:
            for i, point in enumerate(zip(attrs2D.T[0], attrs2D.T[1])):
                fig.add_trace(go.Scatter3d(x=[point[0], point[0]], y=[point[1], point[1]], z=[vmin, 0.5*vmin],  mode='lines', hoverinfo='none', line=dict(width=2, color='blue'), showlegend=False))

                annotations.append(dict(showarrow=False, x=point[0], y=point[1], z=0.4*vmin, text=attrs_names[i], xanchor="center", xshift=10, opacity=1, font=dict(color='black', size=fontsize)))

            fig.add_trace(go.Scatter3d(x=attrs2D.T[0], y=attrs2D.T[1], z=0. * attrs2D.T[0] + 0.5*vmin, mode='markers', 
                                        hovertext=attrs_names,
                                        hoverinfo='text',
                                        marker=dict(size=5, color='blue'),
                                        projection=dict(z=dict(show=True)),
                                        showlegend=False))

        fig.update_layout(title='Hopfield Attractors', autosize=False, width=700, height=700, margin=dict(l=75, r=75, b=75, t=90))

        fig.update_layout(scene = {'xaxis': {'title_text': 'PC1', 'nticks': 10, 'spikesides': False, 'showspikes': False, 'showbackground': False, 'showline': False, 'showticklabels': False, 'showaxeslabels': False}, 'yaxis': {'title_text': 'PC2', 'nticks': 10, 'spikesides': False, 'showspikes': False, 'showbackground': False, 'showline': False, 'showticklabels': False, 'showaxeslabels': False}, 'zaxis': {'title_text': 'Energy', 'range': (vmin, 0.), 'nticks': 10, 'showspikes': False, 'showbackground': False, 'showline': False, 'showticklabels': False, 'showaxeslabels': False}, 'aspectratio': {'x': 1, 'y': 1, 'z': 0.33}, 'annotations': annotations})

        fig.update_layout(scene_camera=dict(up=dict(x=0, y=0, z=2), center=dict(x=0, y=0, z=0), eye=dict(x=0, y=-0.25, z=1.25)))

        fileName = self.dataName + '_energy_landscape_PC%s_vs_PC%s'%(PCy, PCx) + nameAppend

        try:
            fig.write_image(os.path.join(self.saveDir, fileName + '.png'), width=700, height=700, scale=quality)

        except Exception as exception:
            print('Cannot save static image (likely due to missing orca). Saving to interactive html')
            attemptSavingHTML = True

        if attemptSavingHTML:
            plot_offline(fig, filename=os.path.join(self.saveDir, fileName + '.html'), auto_open=False)

        return fig
