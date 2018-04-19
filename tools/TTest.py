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
from scipy.stats import ttest_ind




def TTest(df,numGenesToKeepPerCellType):
    uniqueCellTypes = np.unique(df.loc['_TRUE_LABEL'])
    keepIndices = []
    for cellType in uniqueCellTypes:
        x_cellType = np.float_( df[df.columns[df.loc['_TRUE_LABEL'].values==cellType]].values[:-1,:] )
        x_notCellType = np.float_( df[df.columns[df.loc['_TRUE_LABEL'].values!=cellType]].values[:-1,:] )
        pList = ttest_ind(x_cellType,x_notCellType,axis=1,equal_var=False)[1]
#         _, pList = zip(*[ttest_ind(x_cellType[i,:],x_notCellType[i,:],equal_var=False) for i in xrange(x_cellType.shape[0])])
        keepIndices.extend( np.sort(np.argsort(pList)[:numGenesToKeepPerCellType]) )
    keepIndices = np.unique(keepIndices).tolist()+[df.shape[0]-1]
    return list(df.index[keepIndices])

#def TTest(df,numGenesToKeepPerCellType):
#    uniqueCellTypes = np.unique(df.loc['_TRUE_LABEL'])
#    keepIndices = []
#    for cellType in uniqueCellTypes:
#        x_cellType = np.float_( df[df.columns[df.loc['_TRUE_LABEL'].values==cellType]].values )
#        x_notCellType = np.float_( df[df.columns[df.loc['_TRUE_LABEL'].values!=cellType]].values )
#        pList = ttest_ind(x_cellType,x_notCellType,axis=1,equal_var=False)[1]
#        keepIndices.extend( np.sort(np.argsort(pList)[:numGenesToKeepPerCellType]) )
#    keepIndices = np.unique(keepIndices).tolist()
#    return list(df.index[keepIndices])














