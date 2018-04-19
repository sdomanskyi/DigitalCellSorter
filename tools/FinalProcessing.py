# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:45:03 2018

@author: szedlak1
"""


import copy
cdc = copy.deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import scipy.sparse
#from scipy.optimize import minimize
import pandas as pd
import TTest





def FinalProcessing(df_labeled,numGenesToKeepPerCellTypeFromTTest):
    
    ##############################################################################################
    #Tag duplicate genes, e.g. if gene X shows up thrice, they become X_0, X_1, X_2
    ##############################################################################################
    index_array = np.array(df_labeled.index)
    q = np.unique([i for i in index_array if np.sum(index_array==i)>1])
    if len(q)>0:
        for j in q:
            indices = np.where(index_array==j)[0]
            for c,i in enumerate(indices):
                index_array[i] = '%s_%s' % (index_array[i],c)
        df_labeled.index = index_array
    
    print 'Done tagging duplicate genes...'
    
    ##############################################################################################
    #Replace all zeros with small random numbers to balance between spin-up and spin-down
    ##############################################################################################
    
    df_labeled2 = cdc(df_labeled)
    for i in xrange(df_labeled2.shape[0]-1):
        MIN = np.min([j for j in df_labeled2.iloc[i].values if j!=0.0])
        df_labeled2.iloc[i] = df_labeled2.iloc[i] + np.random.rand(df_labeled2.shape[1])*MIN/100000.0
    
    print 'Done replacing zeros...'
    
    ##############################################################################################
    # Create attractors from clusters
    ##############################################################################################
    
    np.random.seed(2)
    
    df_attrsCont = df_labeled2.iloc[:-1,[]]
    
    uniqueLabels = np.unique(df_labeled2.loc['_TRUE_LABEL'])
    allTestCols = []
    allTrainCols = []
    for label in uniqueLabels:
        candidateColumns = [i for i in df_labeled2.columns if df_labeled2.loc['_TRUE_LABEL'][i]==label]
        trainCols,testCols = np.hsplit(np.random.permutation(candidateColumns),[len(candidateColumns)/2])
        trainCols = np.sort(trainCols)
        testCols = np.sort(testCols)
        df_attrsCont.loc[:,label] = np.mean(df_labeled2[trainCols].iloc[:-1],axis=1)
        allTrainCols.extend(trainCols)
        allTestCols.extend(testCols)
    
    np.random.seed()
    
    
    
    
    genesToKeep = TTest.TTest(df_labeled2[allTrainCols],numGenesToKeepPerCellTypeFromTTest)
    df_labeled2 = df_labeled2.loc[genesToKeep]
    df_attrsCont = df_attrsCont.loc[genesToKeep[:-1]]
    
    
    
    df_logreg_train = df_labeled2[allTrainCols]
    
    df_attrs = df_attrsCont.apply(lambda q: (q>np.median(q))*2.0-1.0,axis=1)
    x_train = df_attrs.values.T
    y_train_labels = df_attrs.columns
    
    df_test = df_labeled2[allTestCols]
    df_test.iloc[:-1] = df_test.iloc[:-1].apply(lambda q: (q>np.median(q))*2.0-1.0,axis=1)
    
    # for i in df_test.index[:-1]:
    #     df_test.loc[i,:] = ( df_test.loc[i,:]>np.median(df_logreg_train.loc[i]) )*2.0-1.0
    
    x_test = df_test.iloc[:-1].values.T
    y_test_labels = df_test.iloc[-1].values
    
    print 'Done transforming data...'
    
    return (x_train,y_train_labels,
            x_test,y_test_labels,
            df_logreg_train,)















