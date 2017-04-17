'''
#####################################################################
#-------------------DEV TOOLS----------------------------------------
#####################################################################
'''

from time import time
from utils import spec2norm
from shutil import copyfile

def resetmas():
    '''
    Rename xdtest.mas and xdtest.inp to xd.mas and xd.inp.
    '''
    try:
        copyfile('xdtest.mas', 'xd.mas')
    except Exception:
        pass
    try:
        copyfile('xdtest.inp', 'xd.inp')
    except Exception:
        pass

def timeDec(f):
    '''
    Decorator to print total running time of a function.
    '''
    def timeFunc(*args, **kwargs):
        tzero = time()
        rtn = f(*args, **kwargs)
        tfin = time()
        print('{0:40}{1:.7f} s'.format(f.__name__, tfin-tzero))
        return rtn
    return timeFunc

def atomsInPair(atoms, pair):
    inPair = True
    pair = [spec2norm(item) for item in pair]
    for atom in atoms:
        if atom not in pair:
            inPair = False
            
    return inPair
