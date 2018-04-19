# -*- coding: utf-8 -*-
"""
Created on Mon Apr 07 21:02:03 2014

@author: Tony
"""

import pickle

#Save data to "pickled" file.
#Usage:
#   FileDump.Save(x,"data/someFile.pythdat")
def Save(x,pathToFile, protocol = 'default'):
    kwargs = {}
    if type(protocol)!=str:
        kwargs['protocol'] = protocol
    with open(pathToFile, 'wb') as f:
        pickle.dump(x, f, **kwargs)

#Load "pickled" data.
#Usage:
#    x=FileDump.Load("data/someFile.pythdat")
def Load(pathToFile):
    with open(pathToFile, 'rb') as f:
        my_list = pickle.load(f)
        return my_list