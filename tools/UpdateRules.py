# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 19:29:59 2016

@author: Tony
"""


import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt
import matplotlib.cm as cmap




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
    h = attrs*(Q_inv*(attrs.T * sigma))
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










