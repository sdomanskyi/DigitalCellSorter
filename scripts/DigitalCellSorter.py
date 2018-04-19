
% matplotlib inline
from importlib import reload
import sys
sys.path.insert(1, '../tools/')
import os
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
import GeneNameConverter; reload(GeneNameConverter)
import MakeMarkerDict; reload(MakeMarkerDict)
from scipy.cluster.hierarchy import linkage, dendrogram



print ('\n================\nPackages loaded!\n================')

class DigitalCellSorter:
    # Load raw data
    def __init__(self, dataName):
        dir_expr = '/'.join(['../data',dataName,'matrix.mtx'])
        dir_geneNames = '/'.join(['../data',dataName,'genes.tsv'])
        with open(dir_expr) as myfile:
            ijv = myfile.readlines()
        header = ijv[:3]
        ijv = np.vstack([np.int_(i.strip('\n').split(' ')) for i in ijv[3:]])
        ijv[:,:2] -= 1 # -1 for first two columns
        imax,jmax = np.int_(header[-1].split(' ')[:2])
        self.df_expr = np.zeros([imax,jmax])
        self.df_expr[ijv[:,0],ijv[:,1]] = ijv[:,2]
        self.df_expr = pd.DataFrame(self.df_expr,index=pd.read_csv(dir_geneNames,delimiter='\t',header=None).values[:,1])
        print ('\n======================\nDone loading raw data!\n======================')
        
        
    def Process(sigma_over_mean_sigma=0.3, n_clusters=11, n_components_pca=100, zscore_cutoff=0.3, saveDir=None):
        #sigma_over_mean_sigma = 0.3 # theta from the paper
        #n_clusters = 11 # number of clusters
        #n_components_pca = 100 # number of principal components to keep for clustering
        #zscore_cutoff = 0.3 # zeta from the paper
        np.random.seed(0)
        gnc = GeneNameConverter.GeneNameConverter(dictDir='../tools/pickledGeneConverterDict/ensembl_hugo_entrez_alias_dict.pythdat')
        
        print ('Pre-filter size: %s genes, %s cells' % self.df_expr.shape)
        #Keep only cells with at least one expressed gene
        self.df_expr = self.df_expr[self.df_expr.columns[np.sum(self.df_expr,axis=0)>0]]
        #Keep only genes expressed in at least one cell
        self.df_expr = self.df_expr.iloc[(np.sum(self.df_expr,axis=1)>0).values,:]
        print ('Post-filter size: %s genes, %s cells' % self.df_expr.shape)

        #Scale all cells
        median = np.median(np.sum(self.df_expr,axis=0))*1.0 # sum of all columns/cells
        self.df_expr = self.df_expr.apply(lambda q: q*median/np.sum(q),axis=0)

        #Replace zeros with minimum value
        MIN = np.min(self.df_expr.values[self.df_expr.values>0])
        self.df_expr = self.df_expr.replace(0,MIN)

        #Take log2 of expression
        self.df_expr = np.log2(self.df_expr)
        self.df_expr -= np.min(self.df_expr.values)

        #Keep only those genes with large enough standard deviation
        self.df_expr = self.df_expr.iloc[np.where(np.std(self.df_expr,axis=1)/np.mean(np.std(self.df_expr.values))>sigma_over_mean_sigma)[0]]
        print ('Post-sigma cut size: %s genes, %s cells' % self.df_expr.shape)

        #Convert gene names from aliases to hugos when possible
        self.df_expr.index = gnc.Convert(list(self.df_expr.index),'alias','hugo',returnUnknownString=False)

        #Sort rows by gene name for the heck of it
        self.df_expr = self.df_expr.sort_index()

        print ('\n=========================\nDone processing raw data!\n=========================')
        
        
        try:
            if savedName != dataName: asdf
            print ('\n====================\nAlready transformed!\n====================')
        except NameError:
            
            print ('Performing PC projection from %s to %s features...' % (self.df_expr.shape[0],n_components_pca))
            X_pca = PCA(n_components=n_components_pca).fit_transform(self.df_expr.values.T).T
            # X_pca_2D = PCA(n_components=2).fit_transform(X_pca.T).T
            print ('Performing tSNE projection from %s to %s features...' % (n_components_pca,2))
            X_tsne = TSNE(n_components=2).fit_transform(X_pca.T).T
            savedName = cdc( dataName )
            print ('\n==================\nDone transforming!\n==================')
        
        
        ##############################################################################################
        # Cluster
        ##############################################################################################

        ac = AgglomerativeClustering(n_clusters=n_clusters)
        cellClusterIndexLabel = ac.fit(X_pca.T).labels_
        possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel))

        print ('\n================\nDone clustering!\n================')
        
        
        ##############################################################################################
        # Get dictionary to map from markers to cell types
        ##############################################################################################
        markerDict,hugo_cd_dict = MakeMarkerDict.MakeMarkerDict('../geneLists/cd_marker_handbook.xlsx',gnc=gnc)

        print ('\n=========================\nDone loading marker data!\n=========================')
        
        
        ##############################################################################################
        # Build marker/cell type weight matrix
        ##############################################################################################

        markers = np.intersect1d(self.df_expr.index,list(markerDict.keys()))
        cellTypes = np.sort(np.unique([item for sublist in markerDict.values() for item in sublist]))

        #This is the marker/cell type matrix (M^T)_mk (the transpose of M_km as defined in the paper)
        df_marker_cellType = pd.DataFrame(np.zeros([len(markers),len(cellTypes)]),index=markers,columns=cellTypes)
        for marker in markers:
            for cellType in markerDict[marker]:
                df_marker_cellType.loc[marker,cellType] = 1

        #Normalize columns (cell types) so that the absolute number of known markers in a given cell type is irrelevant
        df_marker_cellType = df_marker_cellType.apply(lambda q: q/np.sum(q),axis=0)

        #Normalize rows (markers) by the number of cell types expressing that marker (i.e. the "specificity")
        df_marker_cellType = df_marker_cellType.apply(lambda q:q/np.sum(q>0),axis=1)

        print ('\n======================================\nDone building marker/cell type matrix!\n======================================')

        
        ##############################################################################################
        # Compute mean expression of each marker in each cluster
        ##############################################################################################

        df_markers = self.df_expr.loc[markers]
        X_markers = df_markers.values

        means = []
        for cluster in possible_cluster_labels:
            means.append( np.mean(X_markers[:,cellClusterIndexLabel==cluster],axis=1) )

        X_markers_cluster_means = np.vstack(means).T

        #This is the marker/centroid matrix Y_mc
        df_markers_cluster_centroids = pd.DataFrame(X_markers_cluster_means,
                                                    index=markers,
                                                    columns=possible_cluster_labels)

        print ('\n=================================\nDone identifying cluster markers!\n=================================')
        
        
        ##############################################################################################
        # Vote on cell type
        ##############################################################################################

        def f(x,zscore_cutoff):
            return (scipy.stats.zscore(x) > zscore_cutoff)*1

        #This is Z_mc
        df_marker_hits = df_markers_cluster_centroids.apply(lambda q: f(q,zscore_cutoff),axis=1)

        fig,ax = plt.subplots(figsize=(5,3));
        ax.hist(np.sum(df_marker_hits,axis=0),10);
        ax.set_xlabel('Number of markers passing z-score threshold');
        ax.set_ylabel('Number of clusters');
        # if saveDir is not None: fig.savefig(saveDir+'marker_count.pdf')


        #This is the vote matrix V_kc
        df_votes = df_marker_cellType.T.dot(df_marker_hits).apply(lambda q: q/np.sum(q),axis=0)

        #This is the cell type vector T_c
        argmaxes = np.argmax(df_votes.values,axis=0)

        #This is the cell type vector T_c (with cell type names instead of integers)
        cellTypes = pd.Series(df_votes.index[argmaxes],
                              index = df_votes.columns)

        #This is analogous to the percentage of the popular vote in a first-past-the-post political system
        scores = pd.Series([df_votes.iloc[i,j] for i,j in zip(argmaxes,range(len(argmaxes)))],
                            index = df_votes.columns)

        allMarkersList = [df_marker_hits.index[df_marker_hits[c]>0] 
                          for c in df_votes.columns] # important markers in each cluster based on Zmc
        supportingMarkersList = [np.intersect1d(allMarkers,df_marker_cellType.index[df_marker_cellType[c]>0]) 
                                 for allMarkers,c in zip(allMarkersList,cellTypes)]
        allMarkers = pd.Series([' // '.join(i) for i in allMarkersList],
                               index = df_votes.columns)
        supportingMarkers = pd.Series([' // '.join(i) for i in supportingMarkersList],
                                      index = df_votes.columns)

        df_votingResults = pd.DataFrame(index=possible_cluster_labels)
        df_votingResults['Cluster label'] = possible_cluster_labels
        df_votingResults = df_votingResults.set_index('Cluster label')
        df_votingResults['Predicted cell type'] = cellTypes
        df_votingResults['Winning score'] = scores
        df_votingResults['# cells in cluster'] = [np.sum(cellClusterIndexLabel==label) for label in df_votingResults.index]
        df_votingResults['Supporting markers'] = supportingMarkers
        df_votingResults['All markers'] = allMarkers
        df_scoresToAppend = df_votes.T
        df_scoresToAppend.columns = ['%s score'%i for i in df_scoresToAppend.columns]
        df_votingResults = pd.concat((df_votingResults,df_scoresToAppend),axis=1)


        #Append numbers if cell types are repeated (e.g. "Erythrocyte #1', "Erythrocyte #2', ...)
        for cellType in set(df_votingResults['Predicted cell type'].values):
            instances = np.where(df_votingResults['Predicted cell type'].values==cellType)[0]
            if len(instances)>1:
                for counter,instance in enumerate(instances,1):
                    df_votingResults.loc[instance,'Predicted cell type'] = \
                       df_votingResults.loc[instance,'Predicted cell type'] + ' #' + str(counter).zfill(len(str(len(instances))))


        if saveDir is not None: df_votingResults.to_excel('%s/%s_voting.xlsx'%(saveDir,dataName))

        print ('\n============\nDone voting!\n============')

        df_votingResults[df_votingResults.columns[:np.where(df_votingResults.columns=='Supporting markers')[0][0]]]
        
        
        ##############################################################################################
        # Plot mean marker expression for each cluster
        ##############################################################################################
        # Y_mc.T
        X_markers_cluster_means_transpose = X_markers_cluster_means.T

        # i: marker
        # normalization
        for i in range(X_markers_cluster_means_transpose.shape[1]):
            X_markers_cluster_means_transpose[:,i] -= np.min(X_markers_cluster_means_transpose[:,i])
            X_markers_cluster_means_transpose[:,i] /= np.max(X_markers_cluster_means_transpose[:,i])

        
        Z = linkage(X_markers_cluster_means_transpose, 'ward')
        ORDER = dendrogram(Z,no_plot=True,get_leaves=True)['leaves']

        Z = linkage(X_markers_cluster_means_transpose.T, 'ward')
        ORDER2 = dendrogram(Z,no_plot=True,get_leaves=True)['leaves']

        X_markers_cluster_means_sorted = X_markers_cluster_means_transpose[ORDER,:][:,ORDER2]

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
        ax.set_yticklabels([df_votingResults['Predicted cell type'][i]+' ('+str(i)+')' for i in ORDER])#,rotation=24,ha='right',va='top')

        ax.set_xlim([-0.5,X_markers_cluster_means_transpose.shape[1]-0.5])
        ax.set_ylim([-0.5,X_markers_cluster_means_transpose.shape[0]-0.5])

        fig.tight_layout()
        if saveDir is not None: fig.savefig('%s/%s_voting.png'%(saveDir,dataName),dpi=100)

        print ('\n================================\nDone plotting marker expression!\n================================')
        
        
        ##############################################################################################
        # tSNE picture of final clustering and cell types
        ##############################################################################################

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
                            df_votingResults['Predicted cell type'][label].
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
            
            
            
        ##############################################################################################
        # Save labelled expression data to disk
        ##############################################################################################

        df_labeled = cdc(self.df_expr)
        df_labeled.loc['_TRUE_LABEL'] = [df_votingResults['Predicted cell type'][i] for i in cellClusterIndexLabel]
        df_labeled = df_labeled.sort_values(df_labeled.last_valid_index(),axis=1)
        fileName = '%s/%s_expression_labeled.tar.gz' % (saveDir,dataName)
        print ('Saving %s...' % fileName.split('/')[-1])
        df_labeled.to_csv(fileName, compression='gzip')
        print ('Final array size written to disk: %s genes, %s cells' % df_labeled.shape)
        print ('\n=========================\nDone saving labeled data!\n=========================')
        
        
        
        ##############################################################################################
        # Make a bunch of tSNE plots showing relative expression of different markers
        ##############################################################################################
        # max column is where max is
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
                circleIndices = np.where(self.df_expr.loc[marker].values==0)[0] # cells that don't have this marker
                starIndices = np.where(self.df_expr.loc[marker].values>0)[0] # cells that have this marker
                print (starIndices)
                starIndices = starIndices[np.argsort(self.df_expr.loc[marker].values[starIndices])]
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
                           'c':cm.seismic(self.df_expr.loc[marker].values[starIndices]/np.max(self.df_expr.loc[marker].values[starIndices])),
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
                            (df_votingResults['Predicted cell type'][label]).
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
                if saveDir is not None: fig.savefig('%s/marker_subplots/%s_%s_%s.png' % (saveDir,dataName,hugo_cd_dict[marker],marker),dpi=300)