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
import matplotlib.patheffects as path_effects


def ConfusionMatrix(preds,y_train_labels,y_test_labels,title=None,figsize=(12,9),
                    xlabel=None,
                    ylabel=None,
                    source=None,
                    target=None,
                    makePlot = True,
                    ):
    if makePlot:
        fig,ax = plt.subplots(figsize=figsize)
    realPredPairs = zip(y_test_labels,[preds.columns[i] for i in np.argmax(preds.values,axis=1)])
    label2index = {y_train_labels[i]:i for i in xrange(len(y_train_labels))}
    x = np.zeros([len(y_train_labels)]*2)
    for real,pred in realPredPairs:
        x[ label2index[real] , label2index[pred] ] += 1
    
    if source is not None and target is not None and makePlot:
        source_y = np.where(source==y_train_labels)[0][0]
        target_x = np.where(target==y_train_labels)[0][0]
        head_length = 0.15
        head_width = 0.15
        c='r'
        ax.arrow(-0.5, source_y, target_x-head_length, 0,
                 head_width=head_width, head_length=head_length, fc=c, ec=c)
        ax.arrow(target_x, source_y-0.5, 
                 0, head_length-source_y,
                 head_width=head_width, head_length=head_length, fc=c, ec=c)
        import matplotlib.patches as patches
        ax.add_patch(
            patches.Rectangle(
                (target_x-0.5, source_y-0.5),   # (x,y)
                1,          # width
                1,          # height
                facecolor='None',edgecolor=c,lw=3.0,
            )
        )
    
    if makePlot:
        for i in xrange(x.shape[0]):
            for j in xrange(x.shape[1]):
                ax.text(j,i,str(int(x[i,j])),zorder=np.inf,ha='center',va='center',fontsize=22,
                        ).set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),path_effects.Normal()])
    
    rowCounts = np.sum(x,axis=1)
    colCounts = np.sum(x,axis=0)
    
    for i in xrange(len(y_train_labels)):
        x[i,:] /= np.sum(x[i,:])
    
    if makePlot:
        cax = ax.imshow(x,cmap='Purples',interpolation='None',vmin=0,vmax=1)
        ax.xaxis.tick_top()
        ax.set_xticks(range(len(y_train_labels)))
        ax.set_yticks(range(len(y_train_labels)))
        fs_ticks = 12
    #     ax.set_xticklabels(y_train_labels,rotation=20,ha='left',fontsize=fs_ticks)
        xticklabels = [lab+(' (%s)'%int(colCounts[c])) for c,lab in enumerate(y_train_labels)]
        yticklabels = [lab+(' (%s)'%int(rowCounts[c])) for c,lab in enumerate(y_train_labels)]
        yticklabels = [i.replace('/','/\n') for i in yticklabels]
        ax.set_xticklabels(xticklabels,fontsize=fs_ticks,rotation=20,ha='left')
        ax.set_yticklabels(yticklabels,fontsize=fs_ticks)
    #     ax.set_yticklabels([i.replace(' ','\n').replace('T\n','T ').replace('B\n','B ').replace('\nand',' and').replace('multiple\nmyelomas','multiple myelomas').replace('\n#',' #') for i in y_train_labels],fontsize=fs_ticks)
        ax.xaxis.set_label_position('top')
        
        if xlabel is None: xlabel = 'Predicted label'
        if ylabel is None: xlabel = 'Real label'
        
        ax.set_ylabel(ylabel,fontsize=30)
        if title is not None:
            xlabel = xlabel + (' (%s)'%title)
        ax.set_xlabel(xlabel,fontsize=30)
        fig.colorbar(cax)
        XLIM = ax.get_xlim()
        YLIM = ax.get_ylim()[::-1]
        ax.plot(XLIM,YLIM,'w',linewidth=4)
        ax.plot(XLIM,YLIM,'k',linewidth=2)
        
        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM[::-1])
        
        fig.tight_layout()
    
    if makePlot:
        return x,fig,ax
    else:
        return x




































