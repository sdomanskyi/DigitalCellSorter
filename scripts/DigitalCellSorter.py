from importlib import reload
import random
import sys
import os
#sys.path.insert(1, '../tools/')
import copy
cdc = copy.deepcopy
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
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
import pickle
import gzip

import asdf
import multiprocessing

def dumm(i):
    return i

#print ('\n================\nPackages loaded!\n================')

import time

def getStartTime():
    result = time.time()
    return result
    
def getElapsedTime(_start):
    print('Elapsed time: ' + str(np.round(time.time() - _start,2)) + ' sec' +'\n')
    return None

def f(x,zscore_cutoff):
    return 1*(scipy.stats.zscore(x) > zscore_cutoff)

def Get_S_kc(_arguments):

    cellClusterIndexLabel_random, X_markers, markers, zscore_cutoff, df_marker_cellType = _arguments
           
    np.random.seed()
    np.random.shuffle(cellClusterIndexLabel_random)

    possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel_random))

    means = []
    for cluster in possible_cluster_labels:
        means.append( np.mean(X_markers[:,cellClusterIndexLabel_random==cluster],axis=1) )

    X_markers_cluster_means = np.vstack(means).T

    # This is the marker/centroid matrix Y_mc
    df_markers_cluster_centroids_random = pd.DataFrame(X_markers_cluster_means, index=markers, columns=possible_cluster_labels)

    # Construct Zr_mc
    df_marker_hits_random = df_markers_cluster_centroids_random.apply(lambda q: f(q,zscore_cutoff),axis=1)
    
    df_marker_hits_random=pd.DataFrame.from_dict(dict(zip(df_marker_hits_random.index, df_marker_hits_random.values))).T

    #Construct the vote matrix S_kc
    data = df_marker_cellType.T.values.dot(df_marker_hits_random.values)
    df_votes_random = pd.DataFrame(data=data, index=df_marker_cellType.T.index).apply(lambda q: q if np.sum(q)==0 else q/np.sum(q),axis=0)

    return np.ndarray.flatten(df_votes_random.values)


class DigitalCellSorter:
 
    # Normalize df_expr data
    # Params - df_expr: expression data
    #          sigma_over_sigma: threshold when keeping only genes with large enough standard deviation
    #          gnc: gene name converter
    def Normalize(self, df_expr, sigma_over_mean_sigma, gnc, saveDir, dataName, clusterIndex):
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

        if not clusterIndex==None:

            def get_index_of_cluster(clusterIndex):
                temp_file = gzip.open(saveDir + "/X_pca/INDEX_%s_cluster_%s.csv" %(dataName, clusterIndex) + '.pklz')
                INDEX = pickle.load(temp_file)
                temp_file.close()
                return INDEX

            print('\n')
            print(clusterIndex)
            print('\n')

            cells = []
            for cluster in clusterIndex:
                 for cell in get_index_of_cluster(cluster):
                     cells.append(cell)

            #remove all cells from the Data but selected cluster
            columns_to_drop = []
            for item in list(df_expr.columns):
                if not item in cells:
                    columns_to_drop.append(item)

            df_expr = df_expr.drop(columns=columns_to_drop)

            #Remove genes unexpressed in the selected cells
            df_expr = df_expr.iloc[(np.sum(df_expr,axis=1)>0).values,:]

            print('\nRemoved all cells but cluster(s) #%s. Sample: %s. %s cells remaining!\n'%(clusterIndex, dataName, df_expr.shape[1]))

        print ('\n=========================\nDone processing raw data!\n=========================')
        return df_expr
    
    # Project df_expr to lower dimensions
    # Params - df_expr: expression data
    #          n_components_pca: number of components for PCA dimension reduction
    def Project(self, df_expr, n_components_pca):
        startTime = getStartTime()
        print ('Performing PC projection from %s to %s features...' % (df_expr.shape[0],n_components_pca))
        X_pca = PCA(n_components=n_components_pca).fit_transform(df_expr.values.T).T
        getElapsedTime(startTime)

        startTime = getStartTime()
        make_n_components=2
        print ('Performing tSNE projection from %s to %s features...' % (n_components_pca,make_n_components))
        X_tsne2 = TSNE(n_components=make_n_components).fit_transform(X_pca.T).T
        getElapsedTime(startTime)

        print ('\n==================\nDone transforming!\n==================')
        return X_pca, X_tsne2 
    
    # Cluster df_expr
    # Params - X_pca: expression data after pca
    #          n_clusters: number of clusters
    #          clustering_f: clustering function, if not provided then AgglomerativeClustering is used, otherwise should have .fit method and same input and output
    def Cluster(self, X_pca, n_clusters, clustering_f=None):
        startTime = getStartTime()

        if clustering_f is None:
            ac = AgglomerativeClustering(n_clusters=n_clusters)
            cellClusterIndexLabel = ac.fit(X_pca.T).labels_
        else:
            cellClusterIndexLabel = clustering_f.fit(X_pca.T).labels_
        possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel))

        getElapsedTime(startTime)

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
    def Vote(self, df_markers_cluster_centroids, X_markers, markerDict, markerDictnegative, geneListFileName, zscore_cutoff, cellClusterIndexLabel, N_samples_for_distribution, votingScheme, dataName=None, saveDir=None, AvailableCPUsCount = 1):
        # use the voting scheme from the paper if no other scheme is given
        if votingScheme is None:
            votingResults = self.DCSVotingScheme(df_markers_cluster_centroids, X_markers, markerDict, markerDictnegative, geneListFileName, zscore_cutoff, cellClusterIndexLabel, N_samples_for_distribution, dataName, saveDir, AvailableCPUsCount=AvailableCPUsCount)
        else:
            votingResults = votingScheme(df_markers_cluster_centroids, markerDict, markerDictnegative, dataName, saveDir)
        return votingResults
    
    # Produce cluster voting results and write it to an excel file is saveDir is given
    # Params - df_markers_cluster_centroids: Y_mc
    #          markerDict: markers data
    #          zscore_cutoff: zscore cutoff when calculating Z_mc
    #          dataName: name used in output files
    #          saveDir: directory for output files
    # Return - votingResults: a dictionary in form of {clusterLabel: predicted cellType}
    def DCSVotingScheme(self, df_markers_cluster_centroids, X_markers, markerDict, markerDictnegative, geneListFileName, zscore_cutoff, cellClusterIndexLabel, N_samples_for_distribution,
                        dataName=None, 
                        saveDir=None, 
                        MakeVotingResultsMatrixPlot=True, 
                        MakeMarkerHitsCountPlot=True, 
                        MakeHistogramNoisePlot = True,
                        AvailableCPUsCount = 1):

        markers = np.intersect1d(df_markers_cluster_centroids.index,list(markerDict.keys()))
        cellTypes = np.sort(np.unique([item for sublist in markerDict.values() for item in sublist]))

        with open(saveDir + 'ColormapForClusters.txt', 'w') as temp_file:
            for cellType_index in range(len(cellTypes)):
                temp_file.write(cellTypes[cellType_index] + '\t' + str(cm.jet(cellType_index/len(cellTypes))) + '\n')
            temp_file.write('Unknown' + '\t' + str(cm.jet(1.0)) + '\n')

        # Build marker/cell type weight matrix
        # This is the marker/cell type matrix (M^T)_mk (the transpose of M_km as defined in the paper)
        df_marker_cellType = pd.DataFrame(np.zeros([len(markers),len(cellTypes)]),index=markers,columns=cellTypes)
        for marker in markers:
            if marker in markerDict:
                for cellType in markerDict[marker]:
                    df_marker_cellType.loc[marker,cellType] = +1

        ## Normalize columns (cell types) so that the absolute number of known markers in a given cell type is irrelevant
        #df_marker_cellType = df_marker_cellType.apply(lambda q: q/np.sum(q),axis=0)

        ## Normalize rows (markers) by the number of cell types expressing that marker (i.e. the "specificity")
        #df_marker_cellType = df_marker_cellType.apply(lambda q:q/np.sum(q>0),axis=1)

        df_marker_cellType.to_excel('%s/%s_marker-cell_type_weight_matrix.xlsx'%(saveDir,dataName))

        print ('\n======================================\nDone building marker/cell type matrix!\n======================================')


        # Generate distributions for calculating voting p-values, and z-scores
        def Get_S_kc_Distribution(N_trials, RecordHistData = True, N_bins = 100, runParallel = True):

            if runParallel:
                PoolOfWorkersSize = AvailableCPUsCount
                _Pool = multiprocessing.Pool(processes = PoolOfWorkersSize)
                return_values = _Pool.map(Get_S_kc, [(copy.deepcopy(cellClusterIndexLabel), X_markers, markers, zscore_cutoff, df_marker_cellType) for i in range(N_trials)])
                _Pool.close()
                _Pool.join()
                S_kc = np.vstack(return_values)
            else:
                S_kc = [Get_S_kc((copy.deepcopy(cellClusterIndexLabel), X_markers, markers, zscore_cutoff, df_marker_cellType))]
                for i in range(N_trials-1):
                    S_kc = np.concatenate((S_kc, [Get_S_kc((copy.deepcopy(cellClusterIndexLabel), X_markers, markers, zscore_cutoff, df_marker_cellType))]), axis=0)

            if RecordHistData:
                cellTypes = np.sort(np.unique([item for sublist in markerDict.values() for item in sublist]))
                lenCellTypes = len(cellTypes)
                possible_cluster_labels = np.sort(np.unique(cellClusterIndexLabel))
                n_clusters = len(possible_cluster_labels)

                def get_hdf(data):
                    return scipy.stats.rv_histogram(np.histogram(data, bins=N_bins, range=(0,1)))._hpdf

                hist_noise = [get_hdf(S_kc.T[0])]
                for i in range(1,n_clusters*lenCellTypes):
                    hist_noise = np.concatenate((hist_noise, [get_hdf(S_kc.T[i])]), axis=0)

                hist_noise /= N_bins
            
                index_for_dict_noise = []

                for i in range(lenCellTypes):
                    for j in range(n_clusters):
                        index_for_dict_noise = np.append(index_for_dict_noise, cellTypes[i] + ': Cluster #' + str(possible_cluster_labels[j]))

                hist_noise = np.hstack((index_for_dict_noise.reshape((lenCellTypes*n_clusters,1)), hist_noise)).T
                hist_noise = np.hstack((np.vstack(['Bins'] + [str(i/N_bins) for i in range(N_bins + 2)]), hist_noise))

                np.savetxt('%s/%s_dict_noise.csv'%(saveDir,dataName), hist_noise[:-1], delimiter=',', fmt='%s')

            return S_kc.T

        print ('\nGenerating noise distributions...')

        startTime = getStartTime()
        dict_noise = Get_S_kc_Distribution(N_samples_for_distribution)
        getElapsedTime(startTime)

        print ('\n=================================\nDone generating niose distributions!\n=================================')

        # This is Z_mc
        df_marker_hits = df_markers_cluster_centroids.apply(lambda q: f(q,zscore_cutoff),axis=1)

        df_marker_hits=pd.DataFrame.from_dict(dict(zip(df_marker_hits.index, df_marker_hits.values))).T

        if MakeMarkerHitsCountPlot:
            fig,ax = plt.subplots(figsize=(5,3));
            ax.hist(np.sum(df_marker_hits,axis=0),10);
            ax.set_xlabel('Number of markers passing z-score threshold');
            ax.set_ylabel('Number of clusters');
            if saveDir is not None: 
                fig.savefig(saveDir+'marker_count.pdf')

        # Construct the vote matrix V_kc
        df_votes = pd.DataFrame(data=df_marker_cellType.T.values.dot(df_marker_hits.values), 
                                index=df_marker_cellType.T.index).apply(lambda q: q if np.sum(q)==0 else q/np.sum(q),axis=0)

        n_clusters = len(np.sort(np.unique(cellClusterIndexLabel)))

        df_pvalue, df_Z_score = copy.deepcopy(df_votes), copy.deepcopy(df_votes)

        for i in range(df_pvalue.shape[0]):
            for j in range(df_pvalue.shape[1]):
                df_pvalue.values[i][j] = np.count_nonzero(np.nan_to_num(dict_noise[i * n_clusters + j]) >= np.nan_to_num(df_votes.values[i][j])) / dict_noise.shape[1]

                slice = dict_noise[i * n_clusters + j][dict_noise[i * n_clusters + j]>0]
                std = np.std(slice)
                df_Z_score.values[i][j] = -3 if std==0 or std!=std else (np.nan_to_num(df_votes.values[i][j]) - np.mean(slice)) / std

        # This is the cell type vector T_c
        argmaxes = np.argmax(df_Z_score.values,axis=0)
        # This is the cell type vector T_c (with cell type names instead of integers)
        cellTypes = pd.Series(df_Z_score.index[argmaxes], index = df_Z_score.columns)
        cellTypes[(df_Z_score<0).all()]='Unknown'

        # This is analogous to the percentage of the popular vote in a first-past-the-post political system
        scores = pd.Series([df_votes.iloc[i,j] for i,j in zip(argmaxes,range(len(argmaxes)))],
                            index = df_votes.columns)

        allMarkersList = [df_marker_hits.index[df_marker_hits[c]>0] for c in df_votes.columns] # important markers in each cluster based on Zmc
        supportingMarkersList = [np.array([]) if len(allMarkers)==0 else np.intersect1d(allMarkers,df_marker_cellType.index[c in df_marker_cellType.columns and df_marker_cellType[c]>0]) 
                                 for allMarkers,c in zip(allMarkersList,cellTypes)]
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

        df_pvalueToAppend = df_pvalue.T
        df_pvalueToAppend.columns = ['%s p-value'%i for i in df_pvalueToAppend.columns]
        df_votingResults = pd.concat((df_votingResults,df_pvalueToAppend),axis=1)

        df_Z_scoreToAppend = df_Z_score.T
        df_Z_scoreToAppend.columns = ['%s z-score'%i for i in df_Z_scoreToAppend.columns]
        df_votingResults = pd.concat((df_votingResults,df_Z_scoreToAppend),axis=1)

        # Append numbers if cell types are repeated (e.g. "Erythrocyte #1', "Erythrocyte #2', ...)
        for cellType in set(cellTypes):
            instances = np.where(cellTypes.values==cellType)[0]
            if len(instances)>1:
                for counter,instance in enumerate(instances,1):
                    cellTypes[instances[counter-1]] = \
                       cellTypes[instances[counter-1]] + ' #' + str(counter).zfill(len(str(len(instances))))

        for i, sublist in enumerate(list(supportingMarkersList)):
            if len(sublist)<=3 and len(sublist)>=1:
                cellTypes[i] +='*'

        df_votingResults['Predicted cell type'] = cellTypes

        # Build voting results dictionary with cluster labels as keys and predicted celltypes as values
        votingResults = dict(zip(range(df_markers_cluster_centroids.shape[1]), cellTypes))
    
        if saveDir is not None: df_votingResults.to_excel('%s/%s_voting.xlsx'%(saveDir,dataName))

        print ('\n============\nDone voting!\n============')

        if saveDir is not None: 
            writer = pd.ExcelWriter('%s/%s_voting.xlsx'%(saveDir,dataName))
            header_columns = list(df_votingResults.columns[0:4])
            label_column = list(['Predicted cell type'])
            df_votingResults.to_excel(writer, 'Scores', columns=header_columns + list(df_scoresToAppend.columns) + label_column)
            df_votingResults.to_excel(writer, 'p-values', columns=header_columns + list(df_pvalueToAppend.columns) + label_column)
            df_votingResults.to_excel(writer, 'z-scores', columns=header_columns + list(df_Z_scoreToAppend.columns) + label_column)
            writer.save()

        num_of_cell_types = len(np.sort(np.unique([item for sublist in markerDict.values() for item in sublist])))
        num_of_clusters = len(df_votingResults.index)

        if MakeHistogramNoisePlot:
            startTime = getStartTime()
            df_noise_dict = pd.read_csv('%s/%s_dict_noise.csv'%(saveDir,dataName))
        
            gs = matplotlib.gridspec.GridSpec(5, 4)
            fig = plt.figure(figsize=(15,15))

            for i in range(min(num_of_cell_types,20)):
                ax = plt.subplot(gs[i])

                for j in range(num_of_clusters):
                    ax.bar(df_noise_dict['Bins'] + 0.001*j, df_noise_dict[df_noise_dict.columns[1+i*num_of_clusters + j]], width=0.001, align='center', color=cm.tab10(j/num_of_clusters))

                ax.set_title(df_noise_dict.columns[1+i*num_of_clusters][:-12], fontdict={'color': 'b'})

                ax.set_xlabel('Score', fontsize=8)
                ax.set_ylabel('Probability', fontsize=8)

                ax.set_xlim(0.0, 1.0) #0.3
                ax.set_ylim(0, 0.6)

        
            fig.tight_layout()
            fig.savefig(saveDir + dataName + '_noise_dict.pdf',dpi=300)

            print('\n================================\nDone plotting NOISE histograms!\n================================')
            getElapsedTime(startTime)

        if MakeVotingResultsMatrixPlot:
            startTime = getStartTime()
            df_votingResults = pd.read_excel('%s/%s_voting.xlsx'%(saveDir,dataName), sheet_name='z-scores')

            numberOfCells = np.sum(df_votingResults.T.loc['# cells in cluster'])
            num_of_cell_types = df_votingResults.values.shape[1]-5-1
            num_of_clusters = df_votingResults.values.shape[0]

            scores = np.flip(np.array(df_votingResults.values.T[5:5+num_of_cell_types], dtype=float),1)
            winning_scores = np.flip(df_votingResults.values.T[1],0)
            indicis_of_clusters = np.flip(df_votingResults.values.T[0],0)
            assigned_names_of_clusters = np.flip(df_votingResults.values.T[5+num_of_cell_types],0)
            cells_in_clusters = np.flip(df_votingResults.values.T[2],0)
            
            score_columns = df_votingResults.T.index[5:5+num_of_cell_types]

            label_max = 0.27 * max([len(assigned_names_of_clusters[i]) for i in range(len(assigned_names_of_clusters))])

            _figsize = np.float_((num_of_cell_types,num_of_clusters)) / np.max(num_of_clusters) * 15.0
            _figsize[0] += 1.0 + label_max + 3.0 + 1.0
            _figsize[1] += 2.0

            gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[6,1])

            fig = plt.figure(figsize=_figsize)
            ax = plt.subplot(gs[0])
            axx = plt.subplot(gs[1])

            ax.imshow(np.transpose(scores), aspect=1, cmap='Greens', vmin=0, vmax=0.5, interpolation='None')#num_of_clusters/num_of_cell_types

            for i in range(num_of_cell_types):
                for j in range(num_of_clusters):
                    if np.round(scores[i,j],1) > 0:
                        if scores[i,j] == np.max(scores[:,j]):
                            ax.text(i,j,np.round(scores[i,j],1), color='w', fontsize=40 * min([(num_of_cell_types / num_of_clusters),1]),ha='center',va='center').set_path_effects(
                                [path_effects.Stroke(linewidth=10, foreground='red'),path_effects.Normal()])
                        else:
                            ax.text(i,j,np.round(scores[i,j],1), color='w', fontsize=34 * min([(num_of_cell_types / num_of_clusters),1]),ha='center',va='center').set_path_effects(
                                [path_effects.Stroke(linewidth=3, foreground='black'),path_effects.Normal()])

            ax.set_xticks(range(num_of_cell_types))
            ax.set_yticks(range(num_of_clusters))

            ytickslabels = copy.deepcopy(assigned_names_of_clusters)
            for i in range(len(ytickslabels)):
                ytickslabels[i] = str(assigned_names_of_clusters[i]) + ' (' + str(indicis_of_clusters[i]) + ')'

            xtickslabels = np.array(copy.deepcopy(score_columns))
            for i in range(len(xtickslabels)):
                xtickslabels[i] = xtickslabels[i][0:len(xtickslabels[i]) - 6 - 2]
                if i%3 == 1: xtickslabels[i] = '\n' + xtickslabels[i]
                if i%3 == 2: xtickslabels[i] = '\n\n' + xtickslabels[i]

            ax.set_xticklabels(xtickslabels, rotation=0, fontsize=30, ha='center')
            ax.set_yticklabels(ytickslabels, fontsize=30, rotation=0) 
            ax.set_xlim([-0.5,num_of_cell_types - 0.5])
            ax.set_ylim([-0.5,num_of_clusters - 0.5])

            indicis_of_clusters = np.flip(indicis_of_clusters, 0)
            axx.barh(indicis_of_clusters, cells_in_clusters, align='center', color='m')

            for i in range(len(cells_in_clusters)):
                axx.text(np.max(cells_in_clusters), indicis_of_clusters[i], cells_in_clusters[i], ha='right',va='center', color='k', weight='bold', fontsize = 30)
            
            for i in range(len(cells_in_clusters)):
                axx.text(0.02*numberOfCells, indicis_of_clusters[i], str(round(100 * cells_in_clusters[i] / numberOfCells, 1)) + '%', ha='left',va='center', color='b', fontsize = 30)

            axx.set_xticklabels(cells_in_clusters, fontsize=30)
            axx.set_yticklabels(indicis_of_clusters,alpha=0)
            axx.set_xticklabels(indicis_of_clusters,alpha=0)
            axx.set_xlabel('Number of\ncells in clusters', fontsize=30)
            axx.set_ylim(-0.5, num_of_clusters - 0.5)
            fig.tight_layout()
            fig.savefig(saveDir + dataName + '_matrix_voting.pdf',dpi=300)
 
            print('\n================================\nDone plotting MATRIX voting plot!\n================================')
            getElapsedTime(startTime)

        return votingResults, supportingMarkersList
    
    # Save processed data to a zip file
    # Params - df_expr: expression data
    #          votingResults: voting results dictionary
    #          cellClusterIndexLabel: cluster labels for all cells
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def SaveProcessedData(self, df_expr, X_tsne2, votingResults, cellClusterIndexLabel, dataName, saveDir):
        df_labeled = cdc(df_expr)
        df_labeled.loc['_X_coordinate_tSNE'] = [str(X_tsne2[0,i]) for i in range(X_tsne2[0,:].shape[0])]
        df_labeled.loc['_Y_coordinate_tSNE'] = [str(X_tsne2[1,i]) for i in range(X_tsne2[1,:].shape[0])]
        df_labeled.loc['_TRUE_LABEL'] = [list(votingResults.values())[i] for i in cellClusterIndexLabel]
        df_labeled = df_labeled.sort_values(df_labeled.last_valid_index(),axis=1)
        fileName = '%s/%s_expression_labeled.tar.gz' % (saveDir,dataName)
        print ('Saving %s...' % fileName.split('/')[-1])

        temp_file = gzip.open(fileName + '.pklz','wb')
        pickle.dump(df_labeled,temp_file)
        temp_file.close()

        ##To open the DataFrame later for analysis use the following:
        #temp_file = gzip.open(fileName + '.pklz','wb')
        #df_labeled = pickle.load(temp_file)
        #temp_file.close()

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
    def MakeMarkerExpressionPlot(self, votingResults, supportingMarkersList, X_markers_cluster_means, df_markers, df_markers_cluster_centroids, zscore_cutoff, dataName, saveDir):
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

        df_supp_marker_hits = pd.Series(cdc(df_marker_hits.values), index=cdc(df_marker_hits.index))

        for i in range(len(df_supp_marker_hits)):
            df_supp_marker_hits[i].fill(0)
            curr_gene_name = df_supp_marker_hits.index[i]
            for j in range(len(supportingMarkersList)):
                if curr_gene_name in supportingMarkersList[j]:
                    df_supp_marker_hits[i][j] = 1

        df_marker_hits=pd.DataFrame.from_dict(dict(zip(df_marker_hits.index, df_marker_hits.values))).T
        df_supp_marker_hits=pd.DataFrame.from_dict(dict(zip(df_supp_marker_hits.index, df_supp_marker_hits.values))).T
        
        X_marker_hits = df_marker_hits.values.T[ORDER,:][:,ORDER2]
        X_supp_marker_hits = df_supp_marker_hits.values.T[ORDER,:][:,ORDER2]

        _figsize =\
            np.float_(X_markers_cluster_means_transpose.shape[::-1])/\
                    np.max(X_markers_cluster_means_transpose.shape)*15.0+2.0

        _figsize[1] *= 1.5

        fig,ax = plt.subplots(figsize=_figsize)
        ax.imshow(X_markers_cluster_means_sorted,cmap='Blues',interpolation='None')
        ax.set_aspect('auto')


        i_list,j_list = np.where(X_marker_hits.T>0)
        ax.plot(i_list,j_list,'k*',mec='w', markersize=4)

        i_list_supp,j_list_supp = np.where(X_supp_marker_hits.T>0)
        ax.plot(i_list_supp,j_list_supp,'r*', mfc='None', markersize=4) #mec='k', alpha=0.5, markersize=6

        ax.set_xticks(range(X_markers_cluster_means_transpose.shape[1]))
        ax.set_yticks(range(X_markers_cluster_means_transpose.shape[0]))

        xtickslabels = np.array(df_markers.index[ORDER2])
        for i in range(0,len(xtickslabels),2):
            xtickslabels[i] += " ─────────"

        ax.set_xticklabels(xtickslabels,rotation=90, fontsize=10)
        ax.set_yticklabels([list(votingResults.values())[i]+' ('+str(i)+')' for i in ORDER], rotation=45, fontsize=10) #,rotation=24,ha='right',va='top')
        #ax.set_yticklabels(['Cluster '+str(i) for i in ORDER], rotation=45, fontsize=10)

        ax.set_xlim([-0.5,X_markers_cluster_means_transpose.shape[1]-0.5])
        ax.set_ylim([-0.5,X_markers_cluster_means_transpose.shape[0]-0.5])

        fig.tight_layout()
        if saveDir is not None: 
            fig.savefig('%s/%s_voting.png'%(saveDir,dataName),dpi=300)
            fig.savefig('%s/%s_voting.pdf'%(saveDir,dataName),dpi=300)

        print ('\n================================\nDone plotting marker expression!\n================================')

    
    # Produce subplots on each marker and its expression on all clusters
    # Params - df_expr: expression data
    #          votingResults: voting results dictionary
    #          X_tsne2: expression data after tSNE projection
    #          markers: list of markers
    #          cellClusterIndexLabel: cluster labels for all cells
    #          hugo_cd_dict: dictionary that converts gene names from aliases to hugos
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def MakeMarkerSubplot(self, df_expr, votingResults, X_tsne2, markers, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, AvailableCPUsCount=1, NoFrameOnFigures=False, HideClusterLabels=False):
        
        def MarkerSubplot(_arguments):

            counter, marker, df_expr, votingResults, X_tsne2, markers, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, NoFrameOnFigures, HideClusterLabels, XLIM, YLIM, directory = _arguments

            fig,ax = plt.subplots(figsize=(8,8))

            ax.cla()
            if hugo_cd_dict[marker] == marker:
                label = marker
            else:
                label = '%s (%s)' % (hugo_cd_dict[marker],marker)
            ax.plot(np.nan,np.nan,'*',markersize=15,c=cm.seismic(1.0),label=label)
            circleIndices = np.where(df_expr.loc[marker].values==0)[0] # cells that don't have this marker
            starIndices = np.where(df_expr.loc[marker].values>0)[0] # cells that have this marker
            starIndices = starIndices[np.argsort(df_expr.loc[marker].values[starIndices])]
            args1 = [X_tsne2[0,circleIndices],
                        X_tsne2[1,circleIndices]]
            kwargs1 = {'marker':'o',
                        'c':'b',
                        'alpha':0.1,
                        's':6*3,
                        'linewidth':0,}
            args2 = [X_tsne2[0,starIndices],
                        X_tsne2[1,starIndices]]
            kwargs2 = {'marker':'*',
                        'c':cm.seismic(df_expr.loc[marker].values[starIndices]/np.max(df_expr.loc[marker].values[starIndices])),
                        's':30*4,
                        'linewidth':0.0,}
            ax.scatter(*args1,**kwargs1)
            ax.scatter(*args2,**kwargs2)
            for label in set(cellClusterIndexLabel):
                # cells with this label
                X_tsne2_cluster = X_tsne2[:,cellClusterIndexLabel==label]
                x_mean = np.mean(X_tsne2_cluster[0,:])
                y_mean = np.mean(X_tsne2_cluster[1,:])

                _text_label = (list(votingResults.values())[label])
                if HideClusterLabels: _text_label = ''

                ax.text(x_mean,y_mean,
                        _text_label.
                            replace('-','-\n').replace(' ','\n').
                            replace('T\n','T ').replace('B\n','B ').
                            replace('\n#',' #').replace('/','/\n').
                            replace('NK\n','NK ').replace('Stem\n','Stem '),
                        fontsize=10,
                        ha='center',va='center',#alpha=0.75,
                        ).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
                radius = np.sqrt(X_tsne2_cluster.shape[1])*300.0
                ax.scatter(x_mean,y_mean,s=radius*1,facecolors='none',edgecolors='k')
            ax.set_xlim(XLIM)
            ax.set_ylim(YLIM)
            ax.legend(loc='upper right', frameon=False, fontsize=14) #loc='best',numpoints=1,fontsize=12
            ax.set_xticks([]); ax.set_yticks([]); 
            if NoFrameOnFigures:
                fig.patch.set_visible(False)
                ax.axis('off')
            fig.tight_layout()
            if saveDir is not None: 
                fig.savefig('%s/marker_subplots/%s_%s_%s.pdf' % (saveDir,dataName,hugo_cd_dict[marker],marker),dpi=300)
                print ('%s_%s' % (hugo_cd_dict[marker],marker), end=", ", flush=True)

            plt.close(fig)

            return
        
        maxs = np.max(X_tsne2,axis=1)
        mins = np.min(X_tsne2,axis=1)
        maxDiffs = maxs - mins
        deltas = maxDiffs*0.05
        XLIM = [mins[0]-deltas[0],maxs[0]+deltas[0]]
        YLIM = [mins[1]-deltas[1],maxs[1]+deltas[1]]

        directory = '%s/marker_subplots/' % saveDir
        if not os.path.exists(directory):
            os.makedirs(directory)

        print('\nSaving marker subplots:\n')

        for counter,marker in enumerate(markers):
            MarkerSubplot((counter, marker, pd.DataFrame(data=np.reshape(np.array(df_expr.loc[marker]), (1,len(df_expr.loc[marker]))), columns=df_expr.columns, index=[marker]), 
                            votingResults, X_tsne2, markers, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, NoFrameOnFigures, HideClusterLabels, XLIM, YLIM, directory))

        print('\nDone saving marker subplots!')

        return
                    
                    
    # Produce plot on clusters
    # Params - df_votingResults: voting results dictionary
    #          X_tsne2: expression data after tSNE projection
    #          cellClusterIndexLabel: cluster labels for all cells
    #          dataName: name used in output files
    #          saveDir: directory for output files
    def MakeClusterPlot(self, votingResults, X_tsne2, cellClusterIndexLabel, dataName, saveDir):

        NoFrameOnFigures=False
        HideClusterLabels=False

        fig,ax = plt.subplots(figsize=(8,8))

        with open(saveDir + 'ColormapForClusters.txt', 'r') as temp_file:
            colors = temp_file.readlines()
            
        colors = np.vstack([(color.strip('\n').split('\t')) for color in colors])
        colors = pd.DataFrame(colors.T[1], index=colors.T[0])
        colors = colors.apply(lambda x: tuple(np.float_(x[0][1:][:-1].split(','))), axis=1)

        for label in set(cellClusterIndexLabel):
            # cells in that cluster, 2*number of cells
            X_tsne2_cluster = X_tsne2[:,cellClusterIndexLabel==label]

            cellType = list(votingResults.values())[label]

            for i in range(len(cellType)):
                if cellType[i]=="#": 
                    cellType = cellType[:i-len(cellType)-1]
                    break
                    
            cellType = cellType.strip('*')

            ax.plot(X_tsne2_cluster[0,:],X_tsne2_cluster[1,:],'o',
                    color=colors.loc[cellType],mew=0.5,alpha=0.3,markeredgecolor='k')
            x_mean = np.mean(X_tsne2_cluster[0,:])
            y_mean = np.mean(X_tsne2_cluster[1,:])
            # draw the line to the centroid
            #for i in range(X_tsne2_cluster.shape[1]):
            #    ax.plot([x_mean,X_tsne2_cluster[0,i]],[y_mean,X_tsne2_cluster[1,i]],'k-',alpha=0.03,zorder=-np.inf)

            actual_label = list(votingResults.values())[label]
            general_unknown_label = 'Cluster #' + str(label + 1)

            if HideClusterLabels:
                text_label = ''
            else:
                text_label = actual_label

            ax.text(x_mean,y_mean,
                            text_label.
                                replace('-','-\n').replace(' ','\n').
                                replace('T\n','T ').replace('B\n','B ').
                                replace('\n#',' #').replace('/','/\n').
                                replace('NK\n','NK ').replace('Stem\n','Stem '),
                    fontsize=10,
                    ha='center',va='center',
                    ).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])

        ax.set_xticks([])
        ax.set_yticks([])


        maxs = np.max(X_tsne2,axis=1)
        mins = np.min(X_tsne2,axis=1)
        maxDiffs = maxs - mins
        deltas = maxDiffs*0.05
        XLIM = [mins[0]-deltas[0],maxs[0]+deltas[0]]
        YLIM = [mins[1]-deltas[1],maxs[1]+deltas[1]]
        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM)

        if NoFrameOnFigures:
            fig.patch.set_visible(False)
            ax.axis('off')

        fig.tight_layout() 

        _par_Unlabelled = '_Unlabelled' if HideClusterLabels else ''
        _par_NoFrame = '_NoFrame' if NoFrameOnFigures else ''

        if saveDir is not None: fig.savefig('%s/%s%s%s_clusters2D.pdf'%(saveDir,dataName,_par_NoFrame,_par_Unlabelled),dpi=300)
        print ('\n================================\nDone plotting cluster image!\n================================')

        return
     

    # Produce stacked barplot on subclustering
    # Params - clusterIndex: selected cluster for processing
    #          dataName: name used in output files
    #          saveDir: directory for output files     
    def make_stacked_bar_plot_subclustering(self, saveDir, dataName, clusterIndex, clusterName):
        
        def get_stacked_data_and_colors(saveDir):
            with open(saveDir + 'ColormapForClusters.txt', 'r') as temp_file:
                colors = temp_file.readlines()
                colors = np.vstack([(color.strip('\n').split('\t')) for color in colors])
                colors = pd.DataFrame(colors.T[1], index=colors.T[0]).apply(lambda x: tuple(np.float_(x[0][1:][:-1].split(','))), axis=1)

            df = pd.read_excel(saveDir + dataName + '_voting.xlsx')
            index = df['Predicted cell type']

            if not clusterIndex==None:
                barName = dataName[:-5] + ': ' + clusterName
            else:
                barName = dataName[:-5]

            index = [index[i][:len(index[i]) if index[i].find('#')-1==-2 else index[i].find('#')-1].strip('*').strip('#').strip(' ') for i in range(len(index))]
            df_BM_temp = pd.DataFrame(data=df['# cells in cluster'].values, index=index, columns=[barName])
        
            df_main = pd.DataFrame(data=np.zeros((len(colors),1)), index=colors.index, columns=[barName])

            for i, item in enumerate(df_BM_temp.index):
                df_main.loc[item,barName] += df_BM_temp.iloc[i][barName] if df_BM_temp.index[i]==item else 0

            s = 'sums'
            df_main[s] = np.array(np.sum(df_main, axis=1))
            df_main.loc['Unknown',s] = 0
            df_main = df_main.apply(lambda x: 100.*x/np.sum(df_main, axis=0), axis=1).loc[np.sum(df_main, axis=1)>0].sort_values(by=[s]).drop(columns=[s])

            return df_main, colors, clusterName

        df_Main, colors, clusterName = get_stacked_data_and_colors(saveDir)

        saveName = saveDir + "%s_subclustering_stacked_barplot_%s.png"%(dataName, ('All cell clusters' if clusterName==None else clusterName).replace(' ', '_').replace('*', ''))

        fig,ax = plt.subplots(figsize=(4.5,8)) #4.15

        barWidth = 0.85
        cellTypes = df_Main.index
        bottom = np.zeros((len(df_Main.columns)))

        for i in range(len(cellTypes)):
            bottom += df_Main.loc[cellTypes[i-1]] if i>0 else 0
            ax.bar(range(len(df_Main.columns)), list(df_Main.loc[cellTypes[i]]), bottom=list(bottom), color=colors.loc[cellTypes[i]], edgecolor='white', width=barWidth, label=cellTypes[i])
 
        plt.xticks(range(len(df_Main.columns)), list(df_Main.columns), fontsize=12)
        plt.yticks([0,20,40,60,80,100], ['0','20%','40%','60%','80%','100%'], fontsize=12)
            
        handles, labels = ax.get_legend_handles_labels()
        ms = np.max([len(item) for item in labels]) - len('cell')
        labels = [item.replace(' ','\n').replace('B\n', 'B ').replace('T\n', 'T ') if len(item)>=ms else item for item in labels[::-1]]
        ax.legend(handles[::-1], labels, loc='upper left', bbox_to_anchor=(1,1), ncol=1, frameon=False, fontsize=14, labelspacing=1, title = ''.join([' ' for _ in range(60)]))

        plt.xlim((-0.5, len(df_Main.columns) -0.5))
        plt.ylim((0, 100))

        for spine in plt.gca().spines.values():
            spine.set_visible(False)
          
        fig.tight_layout()
        fig.savefig(saveName,dpi=300)

        print ('\n=================================\nDone saving subclustering plot: %s!\n================================='%('All cell clusters' if clusterName==None else clusterName))

        return


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
    def Process(self, df_expr, dataName, geneListFileName, 
                sigma_over_mean_sigma=0.3, 
                n_clusters=5, n_components_pca=100, zscore_cutoff=0.3, 
                saveDir=None, marker_expression_plot=True, tSNE_cluster_plot=True, 
                save_processed_data=True, marker_subplot=True, votingScheme=None,
                N_samples_for_distribution = 10000,
                SaveTransformedData = True,
                attemptLoadingSavedTransformedData = True,
                SaveXpcaDataCSV = True,
                AvailableCPUsCount = 20,
                clusterIndex=None,
                clusterName=None):

        startTime = getStartTime()

        print ("MARK")
        np.random.seed(0)
        gnc = GeneNameConverter.GeneNameConverter(dictDir='tools/pickledGeneConverterDict/ensembl_hugo_entrez_alias_dict.pythdat')
        
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
            
        ##############################################################################################
        # Normalize
        ##############################################################################################
        startTime = getStartTime()
        df_expr = self.Normalize(df_expr, sigma_over_mean_sigma, gnc, saveDir, dataName, clusterIndex)

        getElapsedTime(startTime)

        savedDataExists = False
        
        if attemptLoadingSavedTransformedData: 
            try:
                asdf_file = asdf.open(saveDir + dataName + '.asdf') #all_array_compression='zlib'
                savedDataExists = True
                print('\nLoading transformed data...\n')
            except OSError:
                print('\nWARNING: No saved data file could be openned...\n')

        if savedDataExists: 
            X_pca = copy.deepcopy(asdf_file.tree['X_pca'])
            X_tsne2 = copy.deepcopy(asdf_file.tree['X_tsne2'])
            cellClusterIndexLabel = copy.deepcopy(asdf_file.tree['cellClusterIndexLabel'])
            possible_cluster_labels = copy.deepcopy(asdf_file.tree['possible_cluster_labels'])

            print ('\n====================\nLoaded transformed data!\n====================')

            #override saved data
            ##############################################################################################
            # Cluster
            ##############################################################################################
            cellClusterIndexLabel, possible_cluster_labels = self.Cluster(X_pca, n_clusters)
        else:

            ##############################################################################################
            # Transform data
            ##############################################################################################
            X_pca, X_tsne2 = self.Project(df_expr, n_components_pca)

            
            ##############################################################################################
            # Cluster
            ##############################################################################################
            cellClusterIndexLabel, possible_cluster_labels = self.Cluster(X_pca, n_clusters)


            ##############################################################################################
            # Save transformed data
            ##############################################################################################
            if SaveTransformedData:
                tree = {
                        'X_pca': X_pca,
                        'X_tsne2': X_tsne2,
                        'cellClusterIndexLabel': cellClusterIndexLabel,
                        'possible_cluster_labels': possible_cluster_labels
                    }

                asdf_file = asdf.AsdfFile(tree)
                asdf_file.write_to(saveDir + dataName + '.asdf')


        if SaveXpcaDataCSV and (clusterIndex==None):
            def p_save(fileName, data):
                if not os.path.exists(saveDir + '/X_pca/'):
                    os.makedirs(saveDir + '/X_pca/')

                temp_file = gzip.open(fileName,'wb')
                pickle.dump(data,temp_file)
                temp_file.close()
                return

            p_save(saveDir + "/X_pca/X_pca_%s.csv" %(dataName) + '.pklz', X_pca.T)

            for cluster in list(possible_cluster_labels):
                p_save(saveDir + "/X_pca/X_pca_%s_cluster_%s.csv" %(dataName, cluster) + '.pklz', X_pca.T[cellClusterIndexLabel==cluster])

            for cluster in list(possible_cluster_labels):
                df = pd.DataFrame(data=cellClusterIndexLabel)
                df = df.iloc[np.where(df==cluster)]

                p_save(saveDir + "/X_pca/INDEX_%s_cluster_%s.csv" %(dataName, cluster) + '.pklz', list(df.index))
                
            print ('\n=========================\nDone saving X_pca (%s-component) data!\n=========================' %(n_components_pca))

        ##############################################################################################
        # Get dictionary to map from markers to cell types
        ##############################################################################################
        
        def makeMergedList(ListFileName):
            df_CellTypesGrouped = pd.read_excel(ListFileName, sheet_name='CellTypesGrouped', index_col='CellType')
            Groups = np.sort(np.unique(df_CellTypesGrouped.values))
            AllTypes = np.sort(np.unique(df_CellTypesGrouped.index))

            df_MarkerCellType = pd.read_excel(ListFileName, sheet_name='MarkerCellType', index_col='Marker')
            df_MarkerCellType_NEW = pd.DataFrame(np.zeros((len(df_MarkerCellType.index),len(Groups))), index=df_MarkerCellType.index, columns=Groups)

            for j in range(len(Groups)):
                setTypes = []
                for k in range(len(AllTypes)):
                    _type = df_CellTypesGrouped.loc[AllTypes[k]][0]
                    if Groups[j] == _type:
                        setTypes.append(AllTypes[k])

                for i in range(len(df_MarkerCellType_NEW.index)):
                    df_MarkerCellType_NEW.values[i][j] = 1*(df_MarkerCellType.filter(items=setTypes).values[i].any())
        
            ##drop markers overlapping between cell types
            #df_MarkerCellType_NEW = df_MarkerCellType_NEW.iloc[(np.sum(df_MarkerCellType_NEW,axis=1)<2).values,:]

            dummy,dummy2,temp_hugo_cd_dict = MakeMarkerDict.MakeMarkerDict(geneListFileName,gnc=gnc)

            #drop non-expressed genes
            expressed_genes = np.intersect1d(df_expr.index,list(temp_hugo_cd_dict.keys())) 
            expressed_markers = [temp_hugo_cd_dict[gene] for gene in expressed_genes]
            df_MarkerCellType_NEW = df_MarkerCellType_NEW.iloc[[marker in expressed_markers for marker in df_MarkerCellType_NEW.index],:]

            #drop cell types with <=3 markers only
            columns = df_MarkerCellType_NEW.columns
            columnsToDrop = []
            for i in range(len(columns)):
                if df_MarkerCellType_NEW[columns[i]].sum() <= 3: #3 #0
                    columnsToDrop.append(columns[i])
            df_MarkerCellType_NEW = df_MarkerCellType_NEW.drop(columns=columnsToDrop)

            columns = df_MarkerCellType_NEW.columns
            subsets = []
            for cellType in list(columns):
                subsets.append(list(df_MarkerCellType_NEW[cellType].index[df_MarkerCellType_NEW[cellType]>0]))

            expression_per_gene_sorted = pd.Series(data=np.sum(df_expr, axis=1), index=df_expr.index).iloc[[gene in expressed_genes for gene in df_expr.index]].sort_values(ascending=False)
            
            limit_of_markers_per_cellType = 20
            markers_to_keep = []

            chosen_subsets = []
            for i, cellType in enumerate(list(columns)):
                temp_count = 0
                
                gene_expression_in_subset = expression_per_gene_sorted.iloc[[marker in subsets[i] for marker in [temp_hugo_cd_dict[gene] for gene in expression_per_gene_sorted.index]]].index
                markers_to_genes_subset = [temp_hugo_cd_dict[gene] for gene in gene_expression_in_subset]

                chosen_subsets.append(markers_to_genes_subset[:limit_of_markers_per_cellType if len(markers_to_genes_subset)>limit_of_markers_per_cellType else len(markers_to_genes_subset)])

                for j, marker in enumerate(markers_to_genes_subset):
                    if temp_count < limit_of_markers_per_cellType:
                        markers_to_keep.append(marker)
                        temp_count += 1

            df_MarkerCellType_NEW = df_MarkerCellType_NEW.iloc[[marker in markers_to_keep for marker in df_MarkerCellType_NEW.index],:]

            for i, cellType in enumerate(list(df_MarkerCellType_NEW.columns)):
                for marker in list(df_MarkerCellType_NEW.index):
                    df_MarkerCellType_NEW.loc[marker,cellType] = 0

                for marker in list(df_MarkerCellType_NEW.index):
                    if marker in chosen_subsets[i]:
                        df_MarkerCellType_NEW.loc[marker,cellType] = 1

            df_MarkerCellType_NEW.index.name = 'Marker'
            newName = ListFileName[:-5] + '_merged.xlsx'
            df_MarkerCellType_NEW.to_excel(newName)

            return newName
        markerDict,markerDictnegative,hugo_cd_dict = MakeMarkerDict.MakeMarkerDict(makeMergedList(geneListFileName),gnc=gnc)
        
        print ('\n=========================\nDone loading marker data!\n=========================')
        
        ##############################################################################################
        # Compute mean expression of each marker in each cluster
        ##############################################################################################
        markers = np.intersect1d(df_expr.index,list(markerDict.keys())) ###
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
        votingResults, supportingMarkersList = self.Vote(df_markers_cluster_centroids, X_markers, markerDict, 
                                                         markerDictnegative, geneListFileName, zscore_cutoff, cellClusterIndexLabel, 
                                                         N_samples_for_distribution, votingScheme, dataName, saveDir, AvailableCPUsCount)
        
        
        ##############################################################################################
        # Plot mean marker expression for each cluster
        ##############################################################################################
        if marker_expression_plot:
            startTime = getStartTime()
            self.MakeMarkerExpressionPlot(votingResults, supportingMarkersList, X_markers_cluster_means, df_markers, df_markers_cluster_centroids, zscore_cutoff, dataName, saveDir)
            getElapsedTime(startTime)
        
        ##############################################################################################
        # tSNE picture of final clustering and cell types
        ##############################################################################################
        if tSNE_cluster_plot:
            startTime = getStartTime()
            self.MakeClusterPlot(votingResults, X_tsne2, cellClusterIndexLabel, dataName, saveDir)
            getElapsedTime(startTime)
        
        ##############################################################################################
        # Save labelled expression data to disk
        ##############################################################################################
        if save_processed_data:
            startTime = getStartTime()
            self.SaveProcessedData(df_expr, X_tsne2, votingResults, cellClusterIndexLabel, dataName, saveDir)
            getElapsedTime(startTime)
        
        ##############################################################################################
        # Make a bunch of tSNE plots showing relative expression of different markers
        ##############################################################################################
        if marker_subplot:
            startTime = getStartTime()
            self.MakeMarkerSubplot(df_expr, votingResults, X_tsne2, markers, cellClusterIndexLabel, hugo_cd_dict, dataName, saveDir, AvailableCPUsCount)
            getElapsedTime(startTime)
           
        ##############################################################################################
        # Make stacked barplot for a requested cluster
        ##############################################################################################  
        self.make_stacked_bar_plot_subclustering(saveDir, dataName, clusterIndex, clusterName)

        return cellClusterIndexLabel