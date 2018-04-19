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














class HopfieldClassifier:
    def __init__(self,
                 tmax = 500,
                 T = 0.0001,
                 fractionToUpdate = 0.5,
                 numRuns = 100,
                 CV = False,
                 ):
        self.tmax = tmax
        self.T = T
        self.fractionToUpdate = fractionToUpdate
        self.numRuns = numRuns
        self.CV = CV
    
    def fit(self,
            x_train,
            y_train_bool,
            sample_weight = None,
            ):
        
        self.x_train = x_train.T
        self.y_train_bool = cdc(y_train_bool)
        self.attrs = np.matrix( np.apply_along_axis(lambda q: (q>np.median(q))*2.0-1.0,1,self.x_train) )
        
        if sample_weight is not None: 
            assert len(sample_weight)==len(y_train_bool)
            self.sample_weight = cdc(sample_weight)
        else: self.sample_weight = np.ones_like(self.y_train_bool)*1.0

        if self.CV: return self.CrossValidate()
    
    def CrossValidate(self,
                      verbose = True,
                      ):
        
        from noisyopt import minimizeCompass
        
        options = {}
        if verbose:
            options.update( {'disp': True} )
        
        
        print 'Cross-validating...\n'
        res = minimizeCompass(self.CVObjFn, 
                              [self.T], 
                              bounds = ((0.0,100.0),), 
                              options = options, 
                              tol = 1e-3,
                              paired = False,
                              errorcontrol = False,
                              )
        self.T = res.x[0]
        
        return res
        
        print "Optimal temperature: " % self.T
    
    def CVObjFn(self,
                T, #has to be an array instead of a scalar
                k=3, #num CV slices
                ):
#        return T**2.0
        numRuns = 20
        runList = []
        allIndices = range(self.attrs.shape[1])
        np.random.shuffle(allIndices)
        for testIndices in np.array_split(allIndices,k):
            trainIndices = np.setdiff1d(allIndices,testIndices)
            runList.append(
                   -np.sum(
                      self.predict_proba(self.attrs[:,testIndices].T,
                                         T = T[0],
                                         numRuns=numRuns,
                                         attrs = self.attrs[:,trainIndices],
                                         )[:,1] * \
                      (self.y_train_bool[testIndices]*2.0+1.0)
                           )
                           )
#            runList.append(
#                   -np.sum(
#                      (self.predict_proba(self.attrs[:,testIndices].T,
#                                          T = T[0],
#                                          numRuns=numRuns,
#                                          attrs = self.attrs[:,trainIndices],
#                                          )[:,1]*2.0+1.0 ) * \
#                      (self.y_train_bool[testIndices]*2.0+1.0)
#                           )
#                           )
        
        output = -(np.mean( runList )+1.0)/2.0/self.attrs.shape[0]
        print output
        
        return output
    
    def predict_proba(self,
                      x_test,
                      T = None,
                      numRuns = None,
                      attrs = 'default',
                      verbose = False,
                      ):
        if type(attrs)==str:
            attrs = self.attrs
        if T is None:
            T = self.T
        if numRuns is None:
            numRuns = self.numRuns
        sigma0 = np.matrix( x_test.T )
        assert sigma0.shape[0] == attrs.shape[0]
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
            ovs = Run(attrs,
                      self.tmax,
                      self.T,
                      fractionToUpdate=self.fractionToUpdate,
                      sigma=sigma0,
                     )
            
            mu_max = np.array(np.argmax(ovs[-1],axis=0)).flatten()
            m_max = [ np.array(ovs[-1])[i,j] for (i,j) in zip(mu_max,range(ovs[-1].shape[1])) ]
            label = self.y_train_bool[mu_max].flatten() * 2.0 - 1.0
            weight = self.sample_weight[mu_max]
            pred = m_max * label * weight
            
            predList.append( pred )
        
        predList = np.mean(np.vstack(predList),axis=0)
        prob_high_risk = (predList+1.0)/2.0
        return np.float_( zip( 1.0-prob_high_risk,
                               prob_high_risk ) )




































def Run(
        attrs,
        tmax,
        T,
        sigma='random',
        numSamples=None,
        fractionToUpdate=1.0,
        returnSigmas=False,
        randomSeed = None,
        updateFrequency = None,
        plotFrequency = None,
        boolLabels = None,
        trueLabel = 'none provided',
        ):
    if randomSeed is not None:
        np.random.seed(randomSeed)
    N = attrs.shape[0]
    if type(sigma)==str:
        numSamples = 1
        sigma = np.matrix( (np.random.rand(N,numSamples)>0.5)*2.0-1.0 )
    else: numSamples = sigma.shape[1]
    
#    Q_inv = scipy.sparse.eye(attrs.shape[1])/N
    Q_inv = np.linalg.inv(attrs.T*attrs)
    
    randomVectorForPlottingStuff = np.random.rand(attrs.shape[1])
    overlaps = [Overlap(attrs,sigma,Q_inv)]
    sigmas = [sigma]
    if plotFrequency is not None:
        fig,ax = plt.subplots(figsize=(10,4))
    for t in xrange(tmax):
        if updateFrequency is not None:
            if np.mod(t,updateFrequency)==0:
                print '%s/%s' % (t,tmax)
        sigma,overlap = Update(attrs,
                               sigma,
                               T,
                               Q_inv=Q_inv,
                               fractionToUpdate=fractionToUpdate,
                               )
        if returnSigmas: sigmas.append(sigma)
        else: 
            overlaps.append(overlap)
            if plotFrequency is not None:
                if np.mod(t,plotFrequency)==0 or t==tmax-1:
                    overlaps_temp = np.array( np.hstack(overlaps) )
                    ax.cla()
                    ax.set_xlabel('time')
                    ax.set_ylabel('overlap')
                    for i in xrange(overlaps_temp.shape[0]):
                        if boolLabels is not None:
                            if boolLabels[i]: c='r'
                            else: c='b'
                        else:
                            c = cmap.nipy_spectral(np.random.rand())
                        ax.plot(overlaps_temp[i,:],
                                c,
                                zorder=randomVectorForPlottingStuff[i],
                                linewidth=2,
                                alpha=0.75
                                )
                    ax.set_ylim([-1.03,1.03])
                    ax.set_title('High risk: %s'%trueLabel)
                    for _ in xrange(10):
                        try: fig.savefig('../plots/trajectory.png',dpi=50); break
                        except IOError: pass
#    plt.close(fig)
    if returnSigmas:
        return sigmas
    return overlaps






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
#    overlap = attrs.T*sigma/N
    return sigma,overlap

def Overlap(attrs,sigma,Q_inv): return attrs.T*sigma/sigma.shape[0]
#def Overlap(attrs,sigma,Q_inv): return Q_inv*attrs.T*sigma

















if __name__=='__main__':
    
    
    numSamp = 20
    numFeat = 100
    
    numSamp_test = 3
    
    x_train = np.random.rand(numSamp,numFeat)
    y_train_bool = np.random.rand(numSamp)>0.5
    
#    x_test = np.random.rand(numSamp_test,numFeat)
    x_test = cdc(x_train)
    
    kwargs = {
    'T'       : 0.0001,
    'numRuns' : 3,
    }
    
    model = HopfieldClassifier(**kwargs)
    
    print ''
    print 'Fitting...'
    x = model.fit(x_train,y_train_bool)
#    asdf
    
    print ''
    print 'Predicting...'
    pred = model.predict_proba(x_test)
    
    print ''
    print np.hstack( (pred,np.matrix(y_train_bool).T) )















