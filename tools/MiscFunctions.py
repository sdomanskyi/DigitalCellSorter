# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:45:03 2018

@author: szedlak1
"""

from scipy.stats import ttest_ind
import numpy as np

#Get indices to keep by computing high/low risk p-vals from t-scores for each gene
def TTestIndices(x,y,y_bool,numGenesToKeep):
    highRisk = x[y_bool==True,  :]
    lowRisk  = x[y_bool==False, :]
    #Compute t-test p-values
    _, pList = ttest_ind(highRisk,lowRisk,equal_var=False,axis=0)
    #Ensure the special features are kept
#    pList[-numSpecialFeatures:] = -np.inf
    pList = np.float_(pList)
    pList[np.isnan(pList)] = 1.0
    keepIndices = np.sort(np.argsort(pList)[:numGenesToKeep])
    return keepIndices

#Entropy weights
def EntropyWeights(x,
                   progressionThreshold_days,
                   weight_window,
                   ):
    sigmoid = np.vectorize(lambda x: (1+ np.exp(-(progressionThreshold_days-x)/weight_window))**(-1))
    entropy = np.vectorize(lambda x: -x*np.log2(x) - (1-x)*np.log2(1-x) if 0<x<1 else 0.)
    sampleWeights = 1.-entropy(sigmoid(x))
    return sampleWeights