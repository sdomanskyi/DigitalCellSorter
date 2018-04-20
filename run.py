#import scripts.DigitalCellSorter as DCS
from scripts.DigitalCellSorter import DigitalCellSorter as DCS

if __name__ == '__main__':
   dcs = DCS('aml035pre')
   dcs.Process(saveDir = 'dcs_output/')
