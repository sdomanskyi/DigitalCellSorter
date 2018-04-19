# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import warnings
from imblearn.over_sampling import SMOTE
from imblearn.combine import SMOTETomek, SMOTEENN
from imblearn.under_sampling import EditedNearestNeighbours






#Wrapper for oversampling method that includes ENN edge trimming
def Oversample_SMOTEENN(X,            #2D feature array (num patients x num features)
                        y,            #1D continuous output array (num patients)
                        y_bool,       #1D Boolean output array (num patients) (0=low risk, 1=high risk)
                        newNumSamples, #Desired final number of patients+pseudopatients (balanced so that high/low is 50%/50%)
                        ):
    with warnings.catch_warnings():
        # SMOTE gives silly warnings that I know 
        # are irrelevant for our task
        warnings.simplefilter("ignore")
        my_enn=EditedNearestNeighbours(n_neighbors=10)
        
        if newNumSamples is None:
            newNumSamples = np.max(np.sum(y_bool),np.sum(1.0-y_bool))*2
        
        thing = SMOTEENN(ratio={0:newNumSamples/2, 1:newNumSamples/2},
                         enn=my_enn,
                         )
        #Combine X and y arrays so that D_PFS is interpolated too
        X_for_fitting = np.hstack((X,np.expand_dims(y,axis=1)))
        #Get oversampled array and Boolean labels
        X_oversampled,y_bool_oversampled = thing.fit_sample(X_for_fitting,y_bool)
        #Extract continuous D_PFS
        y_oversampled = X_oversampled[:,-1]
        #Remove D_PFS from the feature array
        X_oversampled = X_oversampled[:,:-1]
    return (X_oversampled,
            y_oversampled,
            y_bool_oversampled)






#Wrapper for oversampling method that includes Tomek edge trimming
def Oversample_SMOTETomek(X,            #2D feature array (num patients x num features)
                          y,            #1D continuous output array (num patients)
                          y_bool,       #1D Boolean output array (num patients) (0=low risk, 1=high risk)
                          newNumSamples, #Desired final number of patients+pseudopatients (balanced so that high/low is 50%/50%)
                          ):
    with warnings.catch_warnings():
        # SMOTE gives silly warnings that I know 
        # are irrelevant for our task
        warnings.simplefilter("ignore")
        
        if newNumSamples is None:
            newNumSamples = np.max(np.sum(y_bool),np.sum(1.0-y_bool))*2
        
        thing = SMOTETomek(
                           ratio={0:newNumSamples/2,   #half of the samples are FALSEs...
                                  1:newNumSamples/2},  #... and half are TRUEs
                           )
        #Combine X and y arrays so that D_PFS is interpolated too
        X_for_fitting = np.hstack((X,np.expand_dims(y,axis=1)))
        #Get oversampled array and Boolean labels
        X_oversampled,y_bool_oversampled = thing.fit_sample(X_for_fitting,y_bool)
        #Extract continuous D_PFS
        y_oversampled = X_oversampled[:,-1]
        #Remove D_PFS from the feature array
        X_oversampled = X_oversampled[:,:-1]
    return (X_oversampled,
            y_oversampled,
            y_bool_oversampled)



#Wrapper for oversampling method
def Oversample_SMOTE(X,            #2D feature array (num patients x num features)
                     y,            #1D continuous output array (num patients)
                     y_bool,       #1D Boolean output array (num patients) (0=low risk, 1=high risk)
                     newNumSamples, #Desired final number of patients+pseudopatients (balanced so that high/low is 50%/50%)
                     ):
    with warnings.catch_warnings():
        # SMOTE gives silly warnings that I know 
        # are irrelevant for our task
        warnings.simplefilter("ignore")
        
        if newNumSamples is None:
            newNumSamples = np.max(np.sum(y_bool),np.sum(1.0-y_bool))*2
        
        thing = SMOTE(
                      ratio={0:newNumSamples/2,   #half of the samples are FALSEs...
                             1:newNumSamples/2},  #... and half are TRUEs
                      )
        #Combine X and y arrays so that D_PFS is interpolated too
        X_for_fitting = np.hstack((X,np.expand_dims(y,axis=1)))
        #Get oversampled array and Boolean labels
        X_oversampled,y_bool_oversampled = thing.fit_sample(X_for_fitting,y_bool)
        #Extract continuous D_PFS
        y_oversampled = X_oversampled[:,-1]
        #Remove D_PFS from the feature array
        X_oversampled = X_oversampled[:,:-1]
    return (X_oversampled,
            y_oversampled,
            y_bool_oversampled)



##Function that takes list of training feature arrays, a matching list of 
##continuous outputs, and a matching list of Boolean outputs and returns 
##balanced data with pseudopatients using the SMOTE oversampling method
#def OversampleLists(X_list,
#                    y_list,
#                    y_bool_list,
#                    newNumSamplesPerDataSet=None,
#                    verbose=False,
#                    tomek=False,
#                    enn=False,
#                    ):
#    newNumSamplesList = [newNumSamplesPerDataSet]*len(y_list)
#    if tomek:
#        OversampleFunction = Oversample_SMOTETomek
#    elif enn:
#        OversampleFunction = Oversample_SMOTEENN
#    else:
#        OversampleFunction = Oversample_SMOTE
#    (X_over_list,
#     y_over_list,
#     y_bool_over_list) = zip(*map(OversampleFunction,
#                                  X_list,
#                                  y_list,
#                                  y_bool_list,
#                                  newNumSamplesList,
#                                  ))
#    if verbose:
#        print ''
#        print '================================================'
#        print '=========     OVERSAMPLING SUMMARY     ========='
#        print '================================================'
#        print 'Old X shapes: ', [X.shape for X in X_list]
#        print 'New X shapes: ', [X.shape for X in X_over_list]
#        print 'Old y shapes: ', [y.shape[0] for y in y_list]
#        print 'New y shapes: ', [y.shape[0] for y in y_over_list]
#        print 'Old fraction TRUEs: ', np.round([np.sum(y)*1.0/len(y) for y in y_bool_list],4)
#        print 'New fraction TRUEs: ', np.round([np.sum(y)*1.0/len(y) for y in y_bool_over_list],4)
#        print '================================================'
#        print ''
#    return (X_over_list,
#            y_over_list,
#            y_bool_over_list)













#Function that takes list of training feature arrays, a matching list of 
#continuous outputs, and a matching list of Boolean outputs and returns 
#balanced data with pseudopatients using the SMOTE oversampling method
def OversampleLists(X_list,
                    y_list,
                    y_bool_list,
                    newNumSamplesPerDataSet=None, #if None, this function 
                                                  #determines the minimum 
                                                  #number of pseudopatients to 
                                                  #generate so that all classes
                                                  #are balanced and all data 
                                                  #sets have the same number of
                                                  #pseudopatients
                    verbose=False,
                    tomek=False,
                    enn=False,
                    ):
    
    if newNumSamplesPerDataSet is None:
        #Count number of members of majority class in each data set
        majorityClassCount = [ np.max((np.sum(y_bool),
                                       np.sum(1.0-y_bool))) for y_bool in y_bool_list ]
        #Find max majority class size across all data sets
        maxSizeMajority = np.max(majorityClassCount)
        #Set the new number of samples to return for each data set (2*number of
        #patients in largest majority class across all data sets)
        newNumSamplesPerDataSet = int(np.round(2*maxSizeMajority))
    
    newNumSamplesList = [newNumSamplesPerDataSet for i in xrange(len(y_list))]
    
#    print newNumSamplesList
    
    if tomek:
        OversampleFunction = Oversample_SMOTETomek
    elif enn:
        OversampleFunction = Oversample_SMOTEENN
    else:
        OversampleFunction = Oversample_SMOTE
    (X_over_list,
     y_over_list,
     y_bool_over_list) = zip(*map(OversampleFunction,
                                  X_list,
                                  y_list,
                                  y_bool_list,
                                  newNumSamplesList,
                                  ))
    if verbose:
        print ''
        print '================================================'
        print '=========     OVERSAMPLING SUMMARY     ========='
        print '================================================'
        print 'Old X shapes: ', [X.shape for X in X_list]
        print 'New X shapes: ', [X.shape for X in X_over_list]
        print 'Old y shapes: ', [y.shape[0] for y in y_list]
        print 'New y shapes: ', [y.shape[0] for y in y_over_list]
        print 'Old fraction TRUEs: ', np.round([np.sum(y)*1.0/len(y) for y in y_bool_list],4)
        print 'New fraction TRUEs: ', np.round([np.sum(y)*1.0/len(y) for y in y_bool_over_list],4)
        print '================================================'
        print ''
    return (X_over_list,
            y_over_list,
            y_bool_over_list)









if __name__=='__main__':
#    #Example oversampling
#    ns = 100 #number of original samples
#    nf = 200 #number of features
#    nt = 15  #number of TRUEs
#    nb = 3   #number of batches (data sets)
#    ns_final = 500 #total number of samples desired (TRUEs plus FALSEs)
#    X_list = [np.random.rand(ns,nf) for i in range(nb)]      #list of random feature arrays
#    y_list = [np.arange(0.0,ns*1.0,1.0) for i in range(nb)]           #list of imbalanced Boolean output arrays
#    y_bool_list = [i<nt for i in y_list]           #list of imbalanced Boolean output arrays
#    (X_over_list,
#     y_over_list,
#     y_bool_over_list) = OversampleLists(X_list,
#                                              y_list,
#                                              y_bool_list,
#                                              500,
#                                              verbose=True
#                                              ) #get balanced data from raw data

    #Example oversampling
    ns = 100 #number of original samples
    nf = 200 #number of features
    nt = 15  #number of TRUEs
    nb = 3   #number of batches (data sets)
    ns_final = 500 #total number of samples desired (TRUEs plus FALSEs)
    X_list = [np.random.rand(ns,nf) for i in range(nb)]      #list of random feature arrays
    y_list = [np.arange(0.0,ns*1.0,1.0) for i in range(nb)]           #list of imbalanced Boolean output arrays
    y_bool_list = [i<nt for i in y_list]           #list of imbalanced Boolean output arrays
    (X_over_list,
     y_over_list,
     y_bool_over_list) = OversampleLists(X_list,
                                              y_list,
                                              y_bool_list,
                                              500,
                                              verbose=True,
                                              tomek=False,
                                              ) #get balanced data from raw data



