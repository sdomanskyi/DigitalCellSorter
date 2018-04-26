
#import scripts.DigitalCellSorter as DCS
import numpy as np
import pandas as pd
from scripts.DigitalCellSorter import DigitalCellSorter as DCS

if __name__ == '__main__':
    
    dataName = 'aml035pre'

    dir_expr = '/'.join(['data',dataName,'matrix.mtx'])
    dir_geneNames = '/'.join(['data',dataName,'genes.tsv'])
    with open(dir_expr) as myfile:
        ijv = myfile.readlines()
    header = ijv[:3]
    ijv = np.vstack([np.int_(i.strip('\n').split(' ')) for i in ijv[3:]])
    ijv[:,:2] -= 1 # -1 for first two columns
    imax,jmax = np.int_(header[-1].split(' ')[:2])
    df_expr = np.zeros([imax,jmax])
    df_expr[ijv[:,0],ijv[:,1]] = ijv[:,2]
    df_expr = pd.DataFrame(df_expr,index=pd.read_csv(dir_geneNames,delimiter='\t',header=None).values[:,1])
    print ('\n======================\nDone loading raw data!\n======================')
    
    dcs = DCS()
    dcs.Process(df_expr, dataName, saveDir = 'demo_output/')