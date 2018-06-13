
from importlib import reload
import sys
import os
#sys.path.insert(1, '../tools/')
import copy
cdc = copy.deepcopy
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.mixture import GaussianMixture
import sklearn.metrics
import matplotlib.patheffects as path_effects
from matplotlib import cm
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering
from tools import GeneNameConverter, MakeMarkerDict
reload(GeneNameConverter)
reload(MakeMarkerDict)
from scipy.cluster.hierarchy import linkage, dendrogram



print ('\n================\nPackages loaded!\n================')

class DigitalCellSorter:
    
    # Normalize df_expr data
    # Params - df_expr: expression data
    #          sigma_over_sigma: threshold when keeping only genes with large enough standard deviation
    #          gnc: gene name converter
    def Normalize(self, df_expr, sigma_over_mean_sigma, gnc):
        print ('Pre-filter size: %s genes, %s cells' % df_expr.shape)
        # Keep only cells with at least one expressed gene
        df_expr = df_expr[df_expr.columns[np.sum(df_expr,axis=0)>0]]
        #Keep only genes expressed in at least one cell
        df_expr = df_expr.iloc[(np.sum(df_expr,axis=1)>0).values,:]
        print ('Post-filter size: %s genes, %s cells' % df_expr.shape)

        # Scale all cells
        median = np.median(np.sum(df_expr,axis=0))*1.0 # sum of all columns/cells
        df_expr = df_expr.apply(lambda q: q*median/np.sum(q),axis=0)

        # Replace zeros with minimum value
        MIN = np.min(df_expr.values[df_expr.values>0])
        df_expr = df_expr.replace(0,MIN)

        # Take log2 of expression
        df_expr = np.log2(df_expr)
        df_expr -= np.min(df_expr.values)

        # Keep only those genes with large enough standard deviation
        df_expr = df_expr.iloc[np.where(np.std(df_expr,axis=1)/np.mean(np.std(df_expr.values))>sigma_over_mean_sigma)[0]]
        print ('Post-sigma cut size: %s genes, %s cells' % df_expr.shape)

        # Convert gene names from aliases to hugos when possible
        #df_expr.index = gnc.Convert(list(df_expr.index),'alias','hugo',returnUnknownString=False)

        # Sort rows by gene name for the heck of it
        df_expr = df_expr.sort_index()

        print ('\n=========================\nDone processing raw data!\n=========================')
        return df_expr
    
    # Project df_expr to lower dimensions
    # Params - df_expr: expression data
    #          n_components_pca: number of components for PCA dimension reduction
    def Project(self, df_expr, n_components_pca):
        print ('Performing PC projection from %s to %s features...' % (df_expr.shape[0],n_components_pca))
        X_pca = PCA(n_components=n_components_pca).fit_transform(df_expr.values.T).T
        # X_pca_2D = PCA(n_components=2).fit_transform(X_pca.T).T
        print ('Performing tSNE projection from %s to %s features...' % (n_components_pca,2))
        X_tsne = TSNE(n_components=2).fit_transform(X_pca.T).T
        
        print ('\n==================\nDone transforming!\n==================')
        return X_pca, X_tsne
    
    # Cluster df_expr
    # Params - X_pca: expression data after pca
    #          n_clusters: number of clusters
    #          clustering_f: clustering function, if not provided then AgglomerativeClustering is used, otherwise should have .fit method and same input and output
    def Cluster(self, X_pca, n_clusters, clustering_f=None):
        if clustering_f is None:
            ac = AgglomerativeClustering(n_clusters=n_clusters)
            cellClusterIndexLabel = ac.fit(X_pca.T).labels_
        else:
            cellClusterIndexLabel = clustering_f.fit(X_pca.T).labels_
        possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel))

        print ('\n================\nDone clustering!\n================')
        return cellClusterIndexLabel, possible_cluster_labels
    
    # Produce cluster voting results
    # Params - df_markers_cluster_centroids: Y_mc
    #          markerDict: markers data
    #          zscore_cutoff: zscore cutoff when calculating Z_mc
    #          votingScheme: voting function. If None then our voting function from the paper will be used. 
    #          NOTE: If you want to use your own custom voting function, make sure the output has the same format as out DCSVotingScheme function as below. You can modify inputs if needed. 
    #          dataName: name used in output files
    #          saveDir: directory for output files
    # Return - votingResults: a dictionary in form of {clusterLabel: predicted cellType}
    def Vote(self, df_markers_cluster_centroids, markerDict, zscore_cutoff, cellClusterIndexLabel, votingScheme, dataName=None, saveDir=None):
        # use the voting scheme from the paper if no other scheme is given
        if votingScheme is None:
            votingResults = self.DCSVotingScheme(df_markers_cluster_centroids, markerDict, zscore_cutoff, cellClusterIndexLabel, dataName, saveDir)
        else:
            votingResults = votingScheme(df_markers_cluster_centroids, markerDict, dataName, saveDir)
        return votingResults
    
    # Produce cluster voting results and write it to an excel file is saveDir is given
    # Params - df_markers_cluster_centroids: Y_mc
    #          markerDict: markers data
    #          zscore_cutoff: zscore cutoff when calculating Z_mc
    #          dataName: name used in output files
    #          saveDir: directory for output files
    # Return - votingResults: a dictionary in form of {clusterLabel: predicted cellType}
    def DCSVotingScheme(self, df_markers_cluster_centroids, markerDict, zscore_cutoff, cellClusterIndexLabel, dataName=None, saveDir=None):
        def f(x,zscore_cutoff):
            return (scipy.stats.zscore(x) > zscore_cutoff)*1
        
        markers = np.intersect1d(df_markers_cluster_centroids.index,list(markerDict.keys()))
        cellTypes = np.sort(np.unique([item for sublist in markerDict.values() for item in sublist]))
        
        # Build marker/cell type weight matrix
        # This is the marker/cell type matrix (M^T)_mk (the transpose of M_km as defined in the paper)
        df_marker_cellType = pd.DataFrame(np.zeros([len(markers),len(cellTypes)]),index=markers,columns=cellTypes)
        for marker in markers:
            for cellType in markerDict[marker]:
                df_marker_cellType.loc[marker,cellType] = 1

        # Normalize columns (cell types) so that the absolute number of known markers in a given cell type is irrelevant
        df_marker_cellType = df_marker_cellType.apply(lambda q: q/np.sum(q),axis=0)

        # Normalize rows (markers) by the number of cell types expressing that marker (i.e. the "specificity")
        df_marker_cellType = df_marker_cellType.apply(lambda q:q/np.sum(q>0),axis=1)

        print ('\n======================================\nDone building marker/cell type matrix!\n======================================')

        # This is Z_mc
        df_marker_hits = df_markers_cluster_centroids.apply(lambda q: f(q,zscore_cutoff),axis=1)
        df_marker_hits=pd.DataFrame.from_dict(dict(zip(df_marker_hits.index, df_marker_hits.values))).T

        fig,ax = plt.subplots(figsize=(5,3));
        ax.hist(np.sum(df_marker_hits,axis=0),10);
        ax.set_xlabel('Number of markers passing z-score threshold');
        ax.set_ylabel('Number of clusters');
        if saveDir is not None: fig.savefig(saveDir+'marker_count.pdf')


        # This is the vote matrix V_kc
        df_votes = df_marker_cellType.T.dot(df_marker_hits).apply(lambda q: q/np.sum(q),axis=0)

        # This is the cell type vector T_c
        argmaxes = np.argmax(df_votes.values,axis=0)

        # This is the cell type vector T_c (with cell type names instead of integers)
        cellTypes = pd.Series(df_votes.index[argmaxes],
                              index = df_votes.columns)
        cellTypes[df_votes.isnull().any()]='Unknown'
        
        # This is analogous to the percentage of the popular vote in a first-past-the-post political system
        scores = pd.Series([df_votes.iloc[i,j] for i,j in zip(argmaxes,range(len(argmaxes)))],
                            index = df_votes.columns)

        allMarkersList = [df_marker_hits.index[df_marker_hits[c]>0] 
                          for c in df_votes.columns] # important markers in each cluster based on Zmc
        supportingMarkersList = [np.intersect1d(allMarkers,df_marker_cellType.index[df_marker_cellType[c]>0]) 
                                 if len(allMarkers)!=0 else [] for allMarkers,c in zip(allMarkersList,cellTypes)]
        allMarkers = pd.Series([' // '.join(i) for i in allMarkersList],
                               index = df_votes.columns)
        supportingMarkers = pd.Series([' // '.join(i) for i in supportingMarkersList],
                                      index = df_votes.columns)

        # Build a table to save voting results to an excel file
        df_votingResults = pd.DataFrame(index=list(range(df_markers_cluster_centroids.shape[1])))
        df_votingResults['Cluster label'] = list(range(df_markers_cluster_centroids.shape[1]))
        df_votingResults = df_votingResults.set_index('Cluster label')
        df_votingResults['Winning score'] = scores
        df_votingResults['# cells in cluster'] = [np.sum(cellClusterIndexLabel==label) for label in df_votingResults.index]
        df_votingResults['Supporting markers'] = supportingMarkers
        df_votingResults['All markers'] = allMarkers
        df_scoresToAppend = df_votes.T
        df_scoresToAppend.columns = ['%s score'%i for i in df_scoresToAppend.columns]
        df_votingResults = pd.concat((df_votingResults,df_scoresToAppend),axis=1)


        # Append numbers if cell types are repeated (e.g. "Erythrocyte #1', "Erythrocyte #2', ...)
        for cellType in set(cellTypes):
            instances = np.where(cellTypes.values==cellType)[0]
            if len(instances)>1:
                for counter,instance in enumerate(instances,1):
                    cellTypes[instances[counter-1]] = \
                       cellTypes[instances[counter-1]] + ' #' + str(counter).zfill(len(str(len(instances))))
        df_votingResults['Predicted cell type'] = cellTypes
        
        # Build voting results dictionary with cluster labels as keys and predicted celltypes as values
        votingResults = dict(zip(range(df_markers_cluster_centroids.shape[1]), cellTypes))
    
        if saveDir is not None: df_votingResults.to_excel('%s/%s_voting.xlsx'%(saveDir,dataName))

        print ('\n============\nDone voting!\n============')
        return votingResults
    
    # Save processed data to a zip file
    # Params - df_expr: expression data
    #          votingResults: voting results dictionary
    #          cellClusterIndexLabel: cluster labels for all cells
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def SaveProcessedData(self, df_expr, votingResults, cellClusterIndexLabel, dataName, saveDir):
        df_labeled = cdc(df_expr)
        df_labeled.loc['_TRUE_LABEL'] = [list(votingResults.values())[i] for i in cellClusterIndexLabel]
        df_labeled = df_labeled.sort_values(df_labeled.last_valid_index(),axis=1)
        fileName = '%s/%s_expression_labeled.tar.gz' % (saveDir,dataName)
        print ('Saving %s...' % fileName.split('/')[-1])
        df_labeled.to_csv(fileName, compression='gzip')
        print ('Final array size written to disk: %s genes, %s cells' % df_labeled.shape)
        print ('\n=========================\nDone saving labeled data!\n=========================')
        
        
    # Produce image on genes and their expression on all clusters
    # Params - votingResults: voting results dictionary
    #          X_markers_cluster_means: Y_mc, mean expression of all markers in all clusters
    #          df_markers: markers and their expression data
    #          df_markers_cluster_centroids: Y_mc
    #          zscore_cutoff: zscore cutoff when calculating Z_mc
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def MakeMarkerExpressionPlot(self, votingResults, X_markers_cluster_means, df_markers, df_markers_cluster_centroids, zscore_cutoff, dataName, saveDir):
        def f(x,zscore_cutoff):
            return (scipy.stats.zscore(x) > zscore_cutoff)*1
        
        # Y_mc.T
        X_markers_cluster_means_transpose = X_markers_cluster_means.T

        # normalization
        for i in range(X_markers_cluster_means_transpose.shape[1]):
            X_markers_cluster_means_transpose[:,i] -= np.min(X_markers_cluster_means_transpose[:,i])
            X_markers_cluster_means_transpose[:,i] /= np.max(X_markers_cluster_means_transpose[:,i])


        Z = linkage(X_markers_cluster_means_transpose, 'ward')
        ORDER = dendrogram(Z,no_plot=True,get_leaves=True)['leaves']

        Z = linkage(X_markers_cluster_means_transpose.T, 'ward')
        ORDER2 = dendrogram(Z,no_plot=True,get_leaves=True)['leaves']

        X_markers_cluster_means_sorted = X_markers_cluster_means_transpose[ORDER,:][:,ORDER2]

        df_marker_hits = df_markers_cluster_centroids.apply(lambda q: f(q,zscore_cutoff),axis=1)
        df_marker_hits=pd.DataFrame.from_dict(dict(zip(df_marker_hits.index, df_marker_hits.values))).T
        
        X_marker_hits = df_marker_hits.values.T[ORDER,:][:,ORDER2]

        fig,ax = plt.subplots(
            figsize=\
              np.float_(X_markers_cluster_means_transpose.shape[::-1])/\
                        np.max(X_markers_cluster_means_transpose.shape)*20.0+2.0
        )
        ax.imshow(X_markers_cluster_means_sorted,cmap='Blues',interpolation='None')
        ax.set_aspect('auto')

        i_list,j_list = np.where(X_marker_hits.T>0)
        ax.plot(i_list,j_list,'k*',mec='w')

        ax.set_xticks(range(X_markers_cluster_means_transpose.shape[1]))
        ax.set_yticks(range(X_markers_cluster_means_transpose.shape[0]))

        xtickslabels = np.array(df_markers.index[ORDER2])

        ax.set_xticklabels(xtickslabels,rotation=90)
        ax.set_yticklabels([list(votingResults.values())[i]+' ('+str(i)+')' for i in ORDER])#,rotation=24,ha='right',va='top')

        ax.set_xlim([-0.5,X_markers_cluster_means_transpose.shape[1]-0.5])
        ax.set_ylim([-0.5,X_markers_cluster_means_transpose.shape[0]-0.5])

        fig.tight_layout()
        if saveDir is not None: fig.savefig('%s/%s_voting.png'%(saveDir,dataName),dpi=100)

        print ('\n================================\nDone plotting marker expression!\n================================')

    
    # Produce subplots on each marker and its expression on all clusters
    # Params - df_expr: expression data
    #          votingResults: voting results dictionary
    #          X_tsne: expression data after tSNE projection
    #          markers: list of markers
    #          cellClusterIndexLabel: cluster labels for all cells
    #          hugo_cd_dict: dictionary that converts gene names from aliases to hugos
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def MakeMarkerSubplot(self, df_expr, votingResults, X_tsne, markers, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir):
        maxs = np.max(X_tsne,axis=1)
        mins = np.min(X_tsne,axis=1)
        maxDiffs = maxs - mins
        deltas = maxDiffs*0.05
        XLIM = [mins[0]-deltas[0],maxs[0]+deltas[0]]
        YLIM = [mins[1]-deltas[1],maxs[1]+deltas[1]]

        directory = '%s/marker_subplots/' % saveDir
        if not os.path.exists(directory):
            os.makedirs(directory)

        q = int(np.ceil(np.sqrt(len(markers))))

        fig,ax = plt.subplots(figsize=(8,8))

        for counter,marker in enumerate(markers):
                ax.cla()
                if hugo_cd_dict[marker] == marker:
                    label = marker
                else:
                    label = '%s (%s)' % (hugo_cd_dict[marker],marker)
                ax.plot(np.nan,np.nan,'*',markersize=15,c=cm.seismic(1.0),label=label)
                circleIndices = np.where(df_expr.loc[marker].values==0)[0] # cells that don't have this marker
                starIndices = np.where(df_expr.loc[marker].values>0)[0] # cells that have this marker
                starIndices = starIndices[np.argsort(df_expr.loc[marker].values[starIndices])]
                args1 = [X_tsne[0,circleIndices],
                         X_tsne[1,circleIndices]]
                kwargs1 = {'marker':'o',
                           'c':'b',
                           'alpha':0.1,
                           's':6*3,
                           'linewidth':0,}
                args2 = [X_tsne[0,starIndices],
                         X_tsne[1,starIndices]]
                kwargs2 = {'marker':'*',
                           'c':cm.seismic(df_expr.loc[marker].values[starIndices]/np.max(df_expr.loc[marker].values[starIndices])),
                           's':30*4,
                           'linewidth':0.0,}
                ax.scatter(*args1,**kwargs1)
                ax.scatter(*args2,**kwargs2)
                for label in set(cellClusterIndexLabel):
                    # cells with this label
                    X_tsne_cluster = X_tsne[:,cellClusterIndexLabel==label]
                    x_mean = np.mean(X_tsne_cluster[0,:])
                    y_mean = np.mean(X_tsne_cluster[1,:])
                    ax.text(x_mean,y_mean,
                            (list(votingResults.values())[label]).
                                replace('-','-\n').replace(' ','\n').
                                replace('T\n','T ').replace('B\n','B ').
                                replace('\n#',' #').replace('/','/\n').
                                replace('NK\n','NK ').replace('Stem\n','Stem '),
                            fontsize=10,
                            ha='center',va='center',#alpha=0.75,
                            ).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
                    radius = np.sqrt(X_tsne_cluster.shape[1])*300.0
                    ax.scatter(x_mean,y_mean,s=radius*1,facecolors='none',edgecolors='k')
                ax.set_xlim(XLIM)
                ax.set_ylim(YLIM)
                ax.legend(loc='best',numpoints=1,fontsize=12)
                ax.set_xticks([]); ax.set_yticks([]); fig.tight_layout()
                if saveDir is not None: 
                    fig.savefig('%s/marker_subplots/%s_%s_%s.png' % (saveDir,dataName,hugo_cd_dict[marker],marker),dpi=300)
                    print ('\n=========================\nDone Saving marker subplot %s_%s labeled data!\n=========================' % (hugo_cd_dict[marker],marker))
                    
                    
    # Produce plot on clusters
    # Params - df_votingResults: voting results dictionary
    #          X_tsne: expression data after tSNE projection
    #          cellClusterIndexLabel: cluster labels for all cells
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def MakeClusterPlot(self, votingResults, X_tsne, cellClusterIndexLabel, dataName, saveDir):
        fig,ax = plt.subplots(figsize=(8,8))
        for label in set(cellClusterIndexLabel):
            # cells in that cluster, 2*number of cells
            X_tsne_cluster = X_tsne[:,cellClusterIndexLabel==label]
            # X_tsne_cluster = X_pca_2D[:,cellClusterIndexLabel==label]

            ax.plot(X_tsne_cluster[0,:],X_tsne_cluster[1,:],'o',
                    color=cm.jet(label*1.0/len(set(cellClusterIndexLabel))),mew=0.5,alpha=0.3,markeredgecolor='k')
            x_mean = np.mean(X_tsne_cluster[0,:])
            y_mean = np.mean(X_tsne_cluster[1,:])
            # draw the line to the centroid
            for i in range(X_tsne_cluster.shape[1]):
                ax.plot([x_mean,X_tsne_cluster[0,i]],[y_mean,X_tsne_cluster[1,i]],'k-',alpha=0.03,zorder=-np.inf)
            ax.text(x_mean,y_mean,
                            list(votingResults.values())[label].
                                replace('-','-\n').replace(' ','\n').
                                replace('T\n','T ').replace('B\n','B ').
                                replace('\n#',' #').replace('/','/\n').
                                replace('NK\n','NK ').replace('Stem\n','Stem '),
                    fontsize=10,
                    ha='center',va='center',
                    ).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])

        ax.set_xticks([])
        ax.set_yticks([])


        maxs = np.max(X_tsne,axis=1)
        mins = np.min(X_tsne,axis=1)
        maxDiffs = maxs - mins
        deltas = maxDiffs*0.05
        XLIM = [mins[0]-deltas[0],maxs[0]+deltas[0]]
        YLIM = [mins[1]-deltas[1],maxs[1]+deltas[1]]
        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM)

        fig.tight_layout() 
        if saveDir is not None: fig.savefig('%s/%s_clusters.png'%(saveDir,dataName),dpi=300)
        print ('\n================================\nDone plotting cluster image!\n================================')
    
        
    # Main function
    # Params - df_expr: expression data
    #          dataName: name used in output files
    #          sigma_over_mean_sigma: threshold when keeping only genes with large enough standard deviation
    #          n_clusters: number of clusters
    #          n_components_pca: number of pca components
    #          zscore_cutoff: zscore cutoff when calculating Z_mc
    #          saveDir: directory for output files, if None then output not saved
    #          marker_expression_plot: whether to produce marker expression plot
    #          tSNE_cluster_plot: whether to produce cluster plot
    #          save_processed_data: whether to save processed data
    #          marker_subplot: whether to make subplots on markers
    def Process(self, df_expr, dataName, sigma_over_mean_sigma=0.3, n_clusters=11, n_components_pca=100, zscore_cutoff=0.3, saveDir=None, marker_expression_plot=True, tSNE_cluster_plot=True, save_processed_data=True, marker_subplot=True, votingScheme=None):
        print ("MARK")
        np.random.seed(0)
        gnc = GeneNameConverter.GeneNameConverter(dictDir='tools/pickledGeneConverterDict/ensembl_hugo_entrez_alias_dict.pythdat')
        
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
            
        ##############################################################################################
        # Normalize
        ##############################################################################################
        df_expr = self.Normalize(df_expr, sigma_over_mean_sigma, gnc)
        
        try:
            if savedName != dataName: asdf
            print ('\n====================\nAlready transformed!\n====================')
        except NameError:
            X_pca, X_tsne = self.Project(df_expr, n_components_pca)
            savedName = cdc(dataName)
  

        ##############################################################################################
        # Cluster
        ##############################################################################################
        cellClusterIndexLabel, possible_cluster_labels = self.Cluster(X_pca, n_clusters)
        
        
        ##############################################################################################
        # Get dictionary to map from markers to cell types
        ##############################################################################################
        markerDict,hugo_cd_dict = MakeMarkerDict.MakeMarkerDict('geneLists/cd_marker_handbook2.xlsx',gnc=gnc)

        print ('\n=========================\nDone loading marker data!\n=========================')
        
        
        ##############################################################################################
        # Compute mean expression of each marker in each cluster
        ##############################################################################################
        markers = np.intersect1d(df_expr.index,list(markerDict.keys()))
        df_markers = df_expr.loc[markers]
        df_markers=df_markers.groupby(level=0).mean()
        X_markers = df_markers.values

        means = []
        for cluster in possible_cluster_labels:
            means.append( np.mean(X_markers[:,cellClusterIndexLabel==cluster],axis=1) )

        X_markers_cluster_means = np.vstack(means).T

        # This is the marker/centroid matrix Y_mc
        df_markers_cluster_centroids = pd.DataFrame(X_markers_cluster_means,
                                                    index=markers,
                                                    columns=possible_cluster_labels)

        print ('\n=================================\nDone identifying cluster markers!\n=================================')
        
        
        ##############################################################################################
        # Vote on cell type
        ##############################################################################################
        votingResults = self.Vote(df_markers_cluster_centroids, markerDict, zscore_cutoff, cellClusterIndexLabel, votingScheme, dataName, saveDir)
        
        
        ##############################################################################################
        # Plot mean marker expression for each cluster
        ##############################################################################################
        if marker_expression_plot:
            self.MakeMarkerExpressionPlot(votingResults, X_markers_cluster_means, df_markers, df_markers_cluster_centroids, zscore_cutoff, dataName, saveDir)
        
        
        ##############################################################################################
        # tSNE picture of final clustering and cell types
        ##############################################################################################
        if tSNE_cluster_plot:
            self.MakeClusterPlot(votingResults, X_tsne, cellClusterIndexLabel, dataName, saveDir)
           
        
        ##############################################################################################
        # Save labelled expression data to disk
        ##############################################################################################
        if save_processed_data:
            self.SaveProcessedData(df_expr, votingResults, cellClusterIndexLabel, dataName, saveDir)
        
        
        ##############################################################################################
        # Make a bunch of tSNE plots showing relative expression of different markers
        ##############################################################################################
        if marker_subplot:
            self.MakeMarkerSubplot(df_expr, votingResults, X_tsne, markers, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir)
            
        return cellClusterIndexLabel