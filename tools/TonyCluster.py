# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:26:14 2018

@author: szedlak1
"""




from sklearn.cluster import KMeans
import numpy as np

def GetCentroids(X,k,returnDistribution=False):
    q = KMeans(n_clusters=k,).fit(X)
    if returnDistribution: return (q.cluster_centers_,np.sort([np.sum(q.labels_==i) for i in set(q.labels_)])[::-1])
    else: return q.cluster_centers_






#X = np.random.rand(10,3)
#k = 4
#
#a,b = GetCentroids(X,k,returnDistribution=True)
#
#print a,b