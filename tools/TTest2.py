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


import statsmodels.sandbox.stats.multicomp
FDR = statsmodels.sandbox.stats.multicomp.multipletests



def TTest(df,targetNumGenes):
    uniqueCellTypes = np.unique(df.loc['_TRUE_LABEL'])
    
    all_pvals = []
    keepIndices = set()
    
    for c,cellType1 in enumerate(uniqueCellTypes[:-1]):
        for cellType2 in uniqueCellTypes[c+1:]:
            x_cellType1 = np.float_( df[df.columns[df.loc['_TRUE_LABEL'].values==cellType1]].values[:-1,:] )
            x_cellType2 = np.float_( df[df.columns[df.loc['_TRUE_LABEL'].values==cellType2]].values[:-1,:] )
            pList = ttest_ind(x_cellType1,x_cellType2,axis=1,equal_var=False)[1]
            all_pvals.append(pList)
    
#    for qqq in xrange(2):
    
    all_pvals = np.hsplit(FDR(np.hstack(all_pvals),method='fdr_bh')[1],np.cumsum([len(i) for i in all_pvals[:-1]]))
    
    xvals = np.unique(np.hstack(all_pvals))
    yvals = []
    for cutoff in xvals:
        yvals.append( len(set(np.hstack([np.where(pLst<=cutoff)[0] for pLst in all_pvals]))) )
        if yvals[-1] >= targetNumGenes:
            xvals = xvals[:len(yvals)]
            break
    
    fig,ax = plt.subplots(figsize=(12,6))
    ax.plot(xvals,yvals)
    ax.set_xscale('log')
    
#    all_pvals = np.hsplit(FDR(np.hstack(all_pvals),method='fdr_bh')[1],np.cumsum([len(i) for i in all_pvals[:-1]]))
    
#    asdf
    
    keepIndices = np.hstack([np.where(pLst<=cutoff)[0] for pLst in all_pvals])
    
    return np.sort(list(set(keepIndices)))

#    for cutoff in 
#    np.sort(list(set(np.hstack([np.where(pList<cutoff)[0] for pList in all_pvals]))))
    
#    return np.sort(list(keepIndices))
#                keepIndices.extend( np.sort(np.argsort(pList)[:numGenesToKeepPerCellType]) )
#    keepIndices = np.unique(keepIndices).tolist()+[df.shape[0]-1]
#    return list(df.index[keepIndices])















