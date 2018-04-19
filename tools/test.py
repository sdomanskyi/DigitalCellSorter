# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 12:21:52 2018

@author: szedlak1
"""

import sys
sys.path.insert(1, '../scripts/')

import numpy as np
from noisyopt import minimizeCompass

import HopfieldRegression

import copy
cdc = copy.deepcopy



#def F(x):
#    return (x**2).sum() + 0.1*np.random.randn()
#
#res = minimizeCompass(F, x0=[1.0, 2.0], deltatol=0.1, paired=False)

#def F(x):
#    return x*x + 0.1*np.random.randn()



attrs = np.matrix( (np.random.rand(100,13)>0.5)*2.0-1.0 )



F = lambda T: np.argmax(
                HopfieldRegression.Run(attrs,10,T,sigma=cdc(attrs[:,0]),numSamples=10)[-1]
                        )
#asdf
res = minimizeCompass(F, x0=[1.0], bounds = ((0.0,10.0),), paired=False)

print res.x