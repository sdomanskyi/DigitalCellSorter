'''Ggeneral functions for conveniece of use
'''

import os
import pickle
import gzip
import time
import numpy as np

def write(data, fileName):
    
    '''Pickle object into a (binary) file
        
    Parameters:
        data: any Pyhton object, e.g. list, dictionary, file, method, variable, etc.
        fileName: path and name of the file to store binary data in

    Returns:
        None
        
    Usage:
        data = [['A', 'B', 'C'], pd.DataFrame()]
        write(data, os.path.join('some dir 1', 'some dir 2', 'File with my data'))
    '''

    with gzip.open(fileName + '.pklz','wb') as temp_file:
        pickle.dump(data, temp_file, protocol=4)

    return None

def read(fileName):

    '''Unpickle object from a (binary) file

    Parameters:
        fileName: path and name of the file with binary data stored in

    Returns:
        Data stored in the provided file
        
    Usage:
        read(os.path.join('some dir 1', 'some dir 2', 'File with my data'))
    '''

    with gzip.open(fileName + '.pklz','rb') as temp_file:
        data = pickle.load(temp_file)
        return data

    return None

def timeMark():

    '''Print total time elapsed from the beggining of the process
    from which the function is called
        
    Parameters:
        None

    Returns:
        None

    Usage:
        timeMark()
    '''
        
    return print('--> Total elapsed time: %s min' % (np.round(time.thread_time() / 60., 1)), '\n')

def getStartTime():

    '''Get time (in seconds) elapsed from the epoch
    
    Parameters:
        None

    Returns:
        Time (in seconds)

    Usage:
        start = getStartTime()
    '''

    return time.time()
    
def getElapsedTime(start):

    '''Print total elapsed time (in minutes) elapsed from the reference point
    
    Parameters:
        start: float or int 
            Reference time (in seconds)

    Returns:
        None

    Usage:
        getElapsedTime(start)
    '''

    return print('Elapsed time: ' + str(np.round((time.time() - start) / 60.,1)) + ' min' + '\n')

