# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:45:03 2018

@author: szedlak1
"""


import copy
cdc = copy.deepcopy
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.cm as cmap
#import scipy.sparse
#from scipy.optimize import minimize
#import pandas as pd
import TTest





def FinalProcessing(df_expr,numGenesToKeepPerCellTypeFromTTest):
    
    df = cdc(df_expr)
    
    ##############################################################################################
    #Tag duplicate genes; e.g. if gene X shows up thrice, they become X_0, X_1, X_2
    ##############################################################################################
    
    print 'Tagging duplicate genes...'
    index_array = np.array(df.index)
    q = np.unique([i for i in index_array if np.sum(index_array==i)>1])
    if len(q)>0:
        for j in q:
            indices = np.where(index_array==j)[0]
            for c,i in enumerate(indices):
                index_array[i] = '%s_%s' % (index_array[i],c)
        df.index = index_array
    
    
    ##############################################################################################
    #Replace all zeros with small random numbers to balance between spin-up and spin-down
    ##############################################################################################
    
    print 'Replacing zeros...'
    for i in xrange(df.shape[0]-1):
        MIN = np.min([j for j in df.iloc[i].values if j!=0.0])
        df.iloc[i] = df.iloc[i] + np.random.rand(df.shape[1])*MIN/100000.0
    
    
    ##############################################################################################
    # Create attractors from clusters
    ##############################################################################################
    
    print 'Doing train/test split and building attractors...'
    np.random.seed(2)
    df_attrsCont = df.iloc[:,[]]
    uniqueLabels = np.unique(df.loc['_TRUE_LABEL'])
    allTestCols = []
    allTrainCols = []
    for label in uniqueLabels:
        candidateColumns = [i for i in df.columns if df.loc['_TRUE_LABEL'][i]==label]
        trainCols,testCols = np.hsplit(np.random.permutation(candidateColumns),[len(candidateColumns)/2])
        trainCols = np.sort(trainCols)
        testCols = np.sort(testCols)
        #The "continuous attractors" are the centroids of each cell type
        df_attrsCont.loc[:,label] = np.mean(df[trainCols].iloc[:-1],axis=1).tolist()+[label]
        allTrainCols.extend(trainCols)
        allTestCols.extend(testCols)
    np.random.seed()
    
    
    ##############################################################################################
    # Perform t-tests to identify which genes should be kept
    ##############################################################################################
    
    print 'Performing t-tests...'
    genesToKeep = TTest.TTest(df[allTrainCols],
                              numGenesToKeepPerCellTypeFromTTest)
    df = df.loc[genesToKeep]
    df_attrsCont = df_attrsCont.loc[genesToKeep]
    
    
    ##############################################################################################
    # Final reformatting
    ##############################################################################################
    
    print 'Reformatting...'
    
    df_continuous_train = df[allTrainCols]
    df_continuous_test = df[allTestCols]
    
    
    
    df_bool_train = cdc(df_attrsCont)
    df_bool_train.iloc[:-1] = df_bool_train.iloc[:-1].apply(lambda q: (q>np.median(q))*2.0-1.0,axis=1)
    
    df_bool_test = cdc(df_continuous_test)
    df_bool_test.iloc[:-1] = df_bool_test.iloc[:-1].apply(lambda q: (q>np.median(q))*2.0-1.0,axis=1)
    
    
    print 'Done processing data.'
    
    return (df_continuous_train,
            df_continuous_test,
            df_bool_train,
            df_bool_test)















