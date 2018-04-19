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
import TTest2





def FinalProcessing(df_expr,targetNumGenes):
    
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
    MINS = np.min(df.replace(0.0,np.inf).values[:-1,:],axis=1)
    df.iloc[:-1] += np.expand_dims(MINS,axis=1) * np.random.rand(df.shape[0]-1,df.shape[1])
    
    ##############################################################################################
    # Create attractors from clusters
    ##############################################################################################
    
    print 'Building attractors...'
    df_attrsCont = df.iloc[:-1,[]]
    for label in np.unique(df.loc['_TRUE_LABEL']):
        cols = [i for i in df.columns if df.loc['_TRUE_LABEL'][i]==label]
        df_attrsCont.loc[:,label] = np.mean(df[cols].iloc[:-1],axis=1)
    
    ##############################################################################################
    # Perform t-tests to identify which genes should be kept
    ##############################################################################################
    
    print 'Performing t-tests...'
    keepIndices = TTest2.TTest(df,targetNumGenes)
    genesToKeep = df.index[keepIndices]
    df = df.loc[genesToKeep]
    df_attrsCont = df_attrsCont.loc[genesToKeep]
    
    
    ##############################################################################################
    # Final reformatting
    ##############################################################################################
    
    print 'Reformatting...'
    
    df_attrsBool = cdc(df_attrsCont)
    df_attrsBool = df_attrsBool.apply(lambda q: (q>np.median(q))*2.0-1.0,axis=1)
    
    print 'Done processing data.'
    
    return (df_attrsCont,
            df_attrsBool)















