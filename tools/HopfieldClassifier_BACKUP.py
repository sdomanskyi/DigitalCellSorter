# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:45:03 2018

@author: szedlak1
"""


import sys
sys.path.insert(1, '../tools/')
import copy
cdc = copy.deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import scipy.sparse
#from scipy.optimize import minimize
import pandas as pd













class HopfieldClassifier:
    def __init__(self,
                 tmax = 500,
                 T = 0.0001,
                 fractionToUpdate = 0.5,
                 numRuns = 100,
                 ):
        self.tmax = tmax
        self.T = T
        self.fractionToUpdate = fractionToUpdate
        self.numRuns = numRuns
    
    def fit(self,
            x_train,
            y_train_labels,):
        self.x_train = x_train.T
        self.y_train_labels = cdc(y_train_labels)
        self.attrs = np.matrix( np.apply_along_axis(lambda q: (q>np.median(q))*2.0-1.0,1,self.x_train) )
    
    def predict_proba(self,
                      x_test,
                      T = None,
                      numRuns = None,
                      verbose = False,
                      ):
        if T is None:
            T = self.T
        if numRuns is None:
            numRuns = self.numRuns
        sigma0 = np.matrix( x_test.T )
        assert sigma0.shape[0] == self.attrs.shape[0]
        # Convert expression data from continuous to Boolean using the training
        # distributions as a reference
        sigma0 = np.vstack( [sigma0[i,:]>np.median(self.x_train[i,:]) 
                             for i in xrange(sigma0.shape[0])] )*2.0-1.0
        if verbose:
            import progressbar
            bar = progressbar.ProgressBar()
            iterable = bar(xrange(numRuns))
        else:
            iterable = xrange(numRuns)
        predList = []
        for sample in iterable:
            ovs = Run(self.attrs,self.tmax,self.T,sigma0,fractionToUpdate=self.fractionToUpdate,)
            predList.append( np.array(np.argmax(ovs,axis=0)).flatten() )
        pred_array = np.int_(np.vstack(predList))
        preds = []
        for i in xrange(pred_array.shape[1]):
            preds.append( [pred_array[:,i].tolist().count(j) for j in xrange(len(self.y_train_labels))] )
        df_pred_freq = pd.DataFrame(preds,columns=self.y_train_labels)
        return df_pred_freq,ovs
































def Run(attrs,
        tmax,
        T,
        sigma,
        fractionToUpdate=1.0,
        randomSeed = None,
        ):
    
    sigma_local = cdc(sigma)
    if randomSeed is not None: np.random.seed(randomSeed)
#    Q_inv = scipy.sparse.eye(attrs.shape[1])/attrs.shape[0]
    Q_inv = np.linalg.inv(attrs.T*attrs)
    kwargs = {'attrs':attrs,
              'sigma':sigma_local,
              'T':T,
              'Q_inv':Q_inv,
              'fractionToUpdate':fractionToUpdate,}
    for t in xrange(tmax):
        sigma_local,overlap = Update(**kwargs)
    return overlap






def Update(attrs,sigma,T,Q_inv=None,fractionToUpdate=1.0):
    (N,numSamples) = np.shape(sigma)
    if Q_inv is None:
        Q_inv = scipy.sparse.eye(attrs.shape[1])/N
    h = attrs*(Q_inv*(attrs.T * sigma)) / N
    if T==0.0:
        sigma_new = np.matrix( np.sign(h+(np.random.rand(N,numSamples)-0.5)*1e-100) )
    else:
        #The following is equivalent to this:
        #    probSpinUp = 1.0 / (1.0+np.exp(-2.0*h/T))
        #    diceRolls = np.random.rand(N,numSamples)
        #    sigma_new = np.matrix(
        #                          (probSpinUp>diceRolls)*2.0 - 1.0
        #                          )
        sigma_new = np.matrix( np.float_(1.0/(1.0+np.exp(-2.0*h/T))>np.random.rand(N,numSamples))*2.0 - 1.0 )
    for mu in xrange(numSamples):
        whichToUpdate = np.random.rand(N)<=fractionToUpdate
        sigma[whichToUpdate,mu] = np.array(
                                           sigma_new[whichToUpdate,mu]
                                           ).flatten()
    overlap = Overlap(attrs,sigma,Q_inv)
    return sigma,overlap

def Overlap(attrs,sigma,Q_inv): return attrs.T*sigma/sigma.shape[0]
#def Overlap(attrs,sigma,Q_inv): return Q_inv*attrs.T*sigma

















if __name__=='__main__':
    numAttr = 50
    numFeat = 3000
    numSamp_test = 3
    x_train = np.random.rand(numAttr,numFeat)
    y_train_labels = np.array( 
            [''.join(np.random.choice(list('abcdefghijklmnopqrstuvwxyz'.upper()),3)) for i in xrange(numAttr)]
                              )
    x_test = cdc(x_train[:numSamp_test,:])
    x_test += np.random.randn(*x_test.shape)/100.0
    y_test_labels = cdc(y_train_labels[:numSamp_test])
    kwargs = {
    'T'       : 0.0000,
    'numRuns' : 10,
    }
    model = HopfieldClassifier(**kwargs)
    print ''
    print 'Fitting...'
    x = model.fit(x_train,y_train_labels)
    print ''
    print 'Predicting...'
    df_pred = model.predict_proba(x_test,verbose=True)
    
    predictions = df_pred.columns[np.argmax(df_pred.values,axis=1)]
    
    print ''
    print 'True, Predicted:'
    for i in zip(y_test_labels,predictions):
        print i
    
#    print ''
#    print np.hstack( (pred,np.matrix(y_train_labels).T) )















