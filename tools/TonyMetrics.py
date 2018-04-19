# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 12:18:42 2017

@author: Tony
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc, precision_recall_curve,\
                            average_precision_score, precision_recall_fscore_support,\
                            matthews_corrcoef, classification_report, confusion_matrix






#Return metrics in a pandas dataframe and make scatter plots of real and predicted outputs
def AllMetrics(predictionList,
               realList,
               shortNames,
               progressionThreshold_days,
               window_c,
               simpleOutput = False,
               makePlots = True,
               ):
    
    #progressionThreshold_days = 30.5*18.0
    #window_c = 9.0 * progressionThreshold_days / 18.0
    sig = np.vectorize(lambda x: (1+ np.exp(-(progressionThreshold_days-x)/window_c))**(-1))
    
    if makePlots:
        fig,axList = plt.subplots(2,2,figsize=(10,10))
        axList = [ax for row in axList for ax in row]
    else:
        axList = [None]*4
    
    columnNames = ['iAUC','Bal','AUC','F1','Matt','Prec','Recall']
    dummyArray = [[i]+np.zeros([len(columnNames)]).tolist() for i in shortNames+['Wtd avg']]
    df = pd.DataFrame(dummyArray, columns = ['Data set']+columnNames).set_index('Data set')
    
    numPositiveList = []
    
    for counter,(name,prediction,real,ax) in enumerate(zip(shortNames,
                                                           predictionList,
                                                           realList,
                                                           axList)):
        
        prediction_bool = np.int_(prediction>0.5)
        
        real_bool = np.int_(real<progressionThreshold_days)
        real_bool_i_list = [np.int_(real < (
                                            progressionThreshold_days * (18.0+i)/18.0
                                            )
                                    ) for i in range(-6,7,2)]
        
        numPositiveList.append(np.sum(real_bool_i_list))
        
        fpr, tpr, _ = roc_curve(real_bool, prediction)
        roc_auc = auc(fpr, tpr)
        
        roc_auc_i = []
        for i, real_bool_i in enumerate(real_bool_i_list):
            fpr_i, tpr_i, _ = roc_curve(real_bool_i, prediction)
            roc_auc_i.append( auc(fpr_i, tpr_i) )
        
#        print zip([progressionThreshold_days * (18.0+i)/18.0 for i in range(-6,7,2)],roc_auc_i)
#        print ''
        
        cnf_matrix = confusion_matrix(real_bool, prediction_bool)
        recall=cnf_matrix[1,1]/float(cnf_matrix[1,1]+cnf_matrix[1,0])
        prec=cnf_matrix[1,1]/float(cnf_matrix[1,1]+cnf_matrix[0,1])
        f1 = 2.*(prec*recall)/(prec+recall)
        b_acc = (cnf_matrix[1,1]/float(cnf_matrix[1,1]+cnf_matrix[1,0])+cnf_matrix[0,0]/float(cnf_matrix[0,0]+cnf_matrix[0,1]))/2.
        matt = matthews_corrcoef(real_bool, prediction_bool)
        
        df['AUC'][name] = roc_auc
        df['iAUC'][name] = np.mean(roc_auc_i)
        df['Bal'][name] = b_acc
        df['F1'][name] = f1
        df['Matt'][name] = matt
        df['Prec'][name] = prec
        df['Recall'][name] = recall
        
        if makePlots:
            real_transformed = sig(real)
            ax.plot(prediction,real_transformed,'o',color='navy',markersize=10,alpha=0.2)
            ax.axvline(0.5, linestyle='--')
            ax.axhline(0.5, linestyle='--')
            [ax.axhline(sig(progressionThreshold_days+30.5*i), linestyle='-',alpha=0.2) for i in range(-6,7,2)]
    #         ax.set_title(name)
            fs = 12
            if counter>=2:
                ax.set_xlabel('prediction',fontsize=fs)
            else:
                ax.set_xticklabels(['' for i in ax.get_xticklabels()])
            if counter%2==0:
                ax.set_ylabel('sig(real)',fontsize=fs)
            else:
                ax.set_yticklabels(['' for i in ax.get_yticklabels()])
            delta = 0.2
            ax.set_xlim([-delta,1.0+delta])
            ax.set_ylim([-delta,1.0+delta])
            text_fudge_factor = 0.75
            
            TP = np.sum( (prediction>0.5)*(real_transformed>0.5) )*1.0
            FP = np.sum( (prediction>0.5)*(real_transformed<=0.5) )*1.0
            TN = np.sum( (prediction<=0.5)*(real_transformed<=0.5) )*1.0
            FN = np.sum( (prediction<=0.5)*(real_transformed>0.5) )*1.0
            TOT_T = np.sum(real_transformed>0.5)*1.0
            TOT_F = np.sum(real_transformed<=0.5)*1.0
            TP_percent = np.round( 100.0*TP/TOT_T , 2 )
            FN_percent = np.round( 100.0*FN/TOT_T , 2 )
            TN_percent = np.round( 100.0*TN/TOT_F , 2 )
            FP_percent = np.round( 100.0*FP/TOT_F , 2 )
            TP_message = str(int(TP))+' ('+str(TP_percent)+'%)'
            FP_message = str(int(FP))+' ('+str(FP_percent)+'%)'
            TN_message = str(int(TN))+' ('+str(TN_percent)+'%)'
            FN_message = str(int(FN))+' ('+str(FN_percent)+'%)'
            
            ax.text(-delta*text_fudge_factor,-delta*text_fudge_factor,'true negative\n'+TN_message,ha='left',va='bottom',fontsize=fs)
            ax.text(-delta*text_fudge_factor,1.0+delta*text_fudge_factor,'false negative\n'+FN_message,ha='left',va='top',fontsize=fs)
            ax.text(1.0+delta*text_fudge_factor,-delta*text_fudge_factor,'false positive\n'+FP_message,ha='right',va='bottom',fontsize=fs)
            ax.text(1.0+delta*text_fudge_factor,1.0+delta*text_fudge_factor,'true positive\n'+TP_message,ha='right',va='top',fontsize=fs)
            
            ax.text(0.5,1.0+delta*text_fudge_factor,name,fontsize=fs*1.4,color='g',ha='center',va='top',
                    bbox=dict(facecolor='white', edgecolor='green', boxstyle='round'))
            
            ax.set_title('iAUC = %s'%np.mean(roc_auc_i),fontsize=fs)
    
    
    for columnName in columnNames:
#        print numPositiveList
#        print df[columnName].values[:-1]
#        asdf
        df[columnName]['Wtd avg'] = np.dot(numPositiveList,df[columnName].values[:-1])/np.sum(numPositiveList)
    
    if makePlots:
        fig.tight_layout()
    
    if simpleOutput:
        df = df[df.columns[:3]]
    
    if makePlots:
        return df,fig,ax
    
    return df









#Return metrics in a pandas dataframe and make scatter plots of real and predicted outputs
def AllMetrics_Simple(prediction,
                      real,
                      progressionThreshold_days,
                      window_c,
                      shortName = None,
                      simpleOutput = False,
                      makePlot = True,
#                      savePlot = False,
                      ):
    
    #progressionThreshold_days = 30.5*18.0
    #window_c = 9.0 * progressionThreshold_days / 18.0
    sig = np.vectorize(lambda x: (1+ np.exp(-(progressionThreshold_days-x)/window_c))**(-1))
    
    if makePlot:
        fig,ax = plt.subplots(figsize=(6,6))
    
    if '.png' in shortName:
        savePlot = True
    
    if shortName is not None:
        simpleName = shortName.split('/')[-1].strip('.png')
    else:
        simpleName = 'no_name'
    
    columnNames = ['iAUC','Bal','AUC','F1','Matt','Prec','Recall']
    dummyArray = [simpleName]+np.zeros([len(columnNames)]).tolist()
    df = pd.Series(dummyArray, index = ['Data set']+columnNames)
    prediction_bool = np.int_(prediction>0.5)
    
    real_bool = np.int_(real<progressionThreshold_days)
    real_bool_i_list = [np.int_(real < (
                                        progressionThreshold_days * (18.0+i)/18.0
                                        )
                                ) for i in range(-6,7,2)]
    
    fpr, tpr, _ = roc_curve(real_bool, prediction)
    roc_auc = auc(fpr, tpr)
    
    roc_auc_i = []
    for i, real_bool_i in enumerate(real_bool_i_list):
        fpr_i, tpr_i, _ = roc_curve(real_bool_i, prediction)
        roc_auc_i.append( auc(fpr_i, tpr_i) )
    
#        print zip([progressionThreshold_days * (18.0+i)/18.0 for i in range(-6,7,2)],roc_auc_i)
#        print ''
    
    cnf_matrix = confusion_matrix(real_bool, prediction_bool)
    recall=cnf_matrix[1,1]/float(cnf_matrix[1,1]+cnf_matrix[1,0])
    prec=cnf_matrix[1,1]/float(cnf_matrix[1,1]+cnf_matrix[0,1])
    f1 = 2.*(prec*recall)/(prec+recall)
    b_acc = (cnf_matrix[1,1]/float(cnf_matrix[1,1]+cnf_matrix[1,0])+cnf_matrix[0,0]/float(cnf_matrix[0,0]+cnf_matrix[0,1]))/2.
    matt = matthews_corrcoef(real_bool, prediction_bool)
    
    df['AUC'] = roc_auc
    df['iAUC'] = np.mean(roc_auc_i)
    df['Bal'] = b_acc
    df['F1'] = f1
    df['Matt'] = matt
    df['Prec'] = prec
    df['Recall'] = recall
    
    if makePlot:
        real_transformed = sig(real)
        ax.plot(prediction,real_transformed,'o',color='navy',markersize=10,alpha=0.2)
        ax.axvline(0.5, linestyle='--')
        ax.axhline(0.5, linestyle='--')
        [ax.axhline(sig(progressionThreshold_days+30.5*i), linestyle='-',alpha=0.2) for i in range(-6,7,2)]
#         ax.set_title(name)
        fs = 12
        ax.set_xlabel('prediction',fontsize=fs)
        ax.set_ylabel('sig(real)',fontsize=fs)
        delta = 0.2
        ax.set_xlim([-delta,1.0+delta])
        ax.set_ylim([-delta,1.0+delta])
        text_fudge_factor = 0.75
        
        TP = np.sum( (prediction>0.5)*(real_transformed>0.5) )*1.0
        FP = np.sum( (prediction>0.5)*(real_transformed<=0.5) )*1.0
        TN = np.sum( (prediction<=0.5)*(real_transformed<=0.5) )*1.0
        FN = np.sum( (prediction<=0.5)*(real_transformed>0.5) )*1.0
        TOT_T = np.sum(real_transformed>0.5)*1.0
        TOT_F = np.sum(real_transformed<=0.5)*1.0
        TP_percent = np.round( 100.0*TP/TOT_T , 2 )
        FN_percent = np.round( 100.0*FN/TOT_T , 2 )
        TN_percent = np.round( 100.0*TN/TOT_F , 2 )
        FP_percent = np.round( 100.0*FP/TOT_F , 2 )
        TP_message = str(int(TP))+' ('+str(TP_percent)+'%)'
        FP_message = str(int(FP))+' ('+str(FP_percent)+'%)'
        TN_message = str(int(TN))+' ('+str(TN_percent)+'%)'
        FN_message = str(int(FN))+' ('+str(FN_percent)+'%)'
        
        ax.text(-delta*text_fudge_factor,-delta*text_fudge_factor,'true negative\n'+TN_message,ha='left',va='bottom',fontsize=fs)
        ax.text(-delta*text_fudge_factor,1.0+delta*text_fudge_factor,'false negative\n'+FN_message,ha='left',va='top',fontsize=fs)
        ax.text(1.0+delta*text_fudge_factor,-delta*text_fudge_factor,'false positive\n'+FP_message,ha='right',va='bottom',fontsize=fs)
        ax.text(1.0+delta*text_fudge_factor,1.0+delta*text_fudge_factor,'true positive\n'+TP_message,ha='right',va='top',fontsize=fs)
        
        if shortName is not None:
            ax.text(0.5,1.0+delta*text_fudge_factor,simpleName,fontsize=fs*1.4,color='g',ha='center',va='top',
                    bbox=dict(facecolor='white', edgecolor='green', boxstyle='round'))
        
        if savePlot:
            fig.savefig(shortName,dpi=100)
    
    if makePlot:
        fig.tight_layout()
    
    if simpleOutput:
        df = df[df.columns[:3]]
    
    if makePlot:
        return df,fig,ax
    
    return df























