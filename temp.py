import os
import pandas as pd
import DigitalCellSorter.DigitalCellSorter as DigitalCellSorter
    
if __name__ == '__main__':

    DCS = DigitalCellSorter.DigitalCellSorter(dataName='BM1', saveDir=os.path.join('output', 'BM1', ''))
    DCS.getAnomalyScoresPlot(cells=DCS.getCells(clusterName='Cluster #2'))