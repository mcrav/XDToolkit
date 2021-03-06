import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.stats import probplot
import matplotlib.pyplot as plt
from quickplot import makeResMap
from ast import literal_eval
from traceback import print_exception, format_exception
import subprocess
import hashlib
import time
import copy
import os
import sys
import webbrowser
from shutil import copyfile, rmtree
from collections import Counter
import multiprocessing as mp
from functools import partial
import datetime
import csv

from xdcool import Ui_MainWindow
from pref import Ui_pref
from about import Ui_aboutBox
from wizard import Ui_wizard
from resmap import Ui_resmap
from sendBug import Ui_sendBug
from sendSugg import Ui_sendSugg
from checkneebs import Ui_checkneebs
from autoTOPXD import Ui_autoTOPXD
from PyQt5 import QtOpenGL
from PyQt5.QtWidgets import (
        QWidget, QMessageBox, QLabel, QDialogButtonBox, QSplashScreen,
        QPushButton, QApplication, QDialog, QFileDialog, QMainWindow, QGridLayout,
        QScrollArea, QSizePolicy, QSpacerItem, QFileSystemModel, QTreeView, QOpenGLWidget)
from PyQt5.QtCore import QSettings, QThread, pyqtSignal, Qt, QMetaObject, QDir, pyqtSlot, QSize, QModelIndex, QRectF
from PyQt5.QtGui import QPixmap, QFont, QStandardItem, QStandardItemModel, QPalette, QColor, QBrush, QIcon, QRegion, QPainterPath
from devtools import resetmas, timeDec
from backup import backup, loadBackup
from databank import helpTexts
from emailfuncs import sendEmail
from xderrfix import check4errors, fixLsmCif, removePhantomAtoms, fixBrokenLabels, addNCST, fixCuNobleGasError
from xdfiletools import addSnlCutoff, setupmas, resetKeyTable, multipoleKeyTable, addDUM, addCustomLocCoords
from resetbond import armRBs, disarmRBs, check4RBHs, check4RB, resetBond, autoResetBond, delResetBond
from wizardfuncs import (seqMultRef, wizAddResetBond, wizAddLocCoords, wizAddCHEMCON, wizAddMultipoles,
                         wizAddCustomLCS)

from initfuncs import ins2all
from results import (FFTDetective, FOUcell, FOU3atoms, grd2values, setupPROPDpops, getDorbs, readSUs,
                     getRF2, getKrauseParam, getConvergence, getDMSDA, setup3AtomLapmap, setupmasTOPXD,
                     geom2Angles, tabPropSetupMas)

from utils import (convert2XDLabel, lab2type, spec2norm, rawInput2labels, labels2list, formatLabels, isfloat,
                   getCellParams, findElements, getNumAtoms, getEleNum, res2inp, findMasCHEMCON, totalEstTime,
                   coords2tuple, listjoin, getAtomList, atomTableBegins, atomTableEnds, seconds2TimeStr,
                   sevenSpacedListjoin, printExc, getProgramFont)

from chemcon import (getEnvSig, removeCHEMCON, check4CHEMCON, writeCHEMCON, findCHEMCONbyInputElement,
                     findCHEMCONbyInputAtoms, CPPgetEnvSig)

'''
This file contains the GUI and functions that rely on global variables.
Some functions that do not rely on global variables are included here as they
are related to functions that heavily rely on global variables (i.e. refinement setup functions).
Likewise, some functions that rely on global variables are contained in other files
to maintain ordered grouping (i.e. chemcon functions).
'''
'''
#####################################################################
#-------------------COMPOUND INITIALIZATION--------------------------
#####################################################################
'''

def makeHklHash(file):
    '''
    Create sha256 hash of shelx.hkl file. Return hash as string.
    '''
    hklHash = ''

    if os.path.isfile(file):
        with open(file,'r') as hkl:
            hklTxt = hkl.read()
            hklHash = hashlib.sha256(bytes(hklTxt, 'utf-8')).hexdigest()

    return hklHash

def inCache(file):
    '''
    Find out if project is in cache. Return result as bool and md5 hash of shelx.hkl.
    '''
    global cachePath
    global cacheHashPath
    inCache = False

    hklHash = makeHklHash(file)
    cacheHashPath = '{}/{}'.format(cachePath, hklHash)
    if os.path.isdir(cacheHashPath):
        inCache = True

    return (inCache, hklHash)

def clearCache():
    '''
    Delete everything from cache.
    '''
    print('Clearing cache')
    rmtree(cachePath)
    os.makedirs(cachePath)


def initialiseGlobVars():
    '''
    Initialise global variables of nearest neighbour, bond angle and atomic position dictionaries.
    '''
    global globAtomLabs
    global globAtomTypes
    global globAtomAngles
    global globAtomPos
    global cachePath

    if os.path.isfile('shelx.hkl'):
        cacheRes = inCache('shelx.hkl')
        hklHash = cacheRes[1]
        isInCache = cacheRes[0]

        #If global dicts are in the cache folder, load them from there instead of making them again.
        try:
            if not isInCache:
                print('Structure not in cache. Initializing compound.')
                raise Exception
            else:
                print('Structure found in cache.')
                with open('{}/{}/atomLabs.buckfast'.format(cachePath, hklHash),'r') as atomLabs:
                    globAtomLabs = literal_eval(atomLabs.read())

                with open('{}/{}/atomTypes.buckfast'.format(cachePath, hklHash),'r') as atomTypes:
                    globAtomTypes = literal_eval(atomTypes.read())

                with open('{}/{}/atomAngles.buckfast'.format(cachePath, hklHash),'r') as atomAngles:
                    globAtomAngles = literal_eval(atomAngles.read())

                with open('{}/{}/atomPos.buckfast'.format(cachePath, hklHash),'r') as atomPos:
                    globAtomPosRaw = literal_eval(atomPos.read())
                    globAtomPos = {}
                    for item, value in globAtomPosRaw.items():
                        globAtomPos[item] = (np.array(value[0]), value[1])

        except Exception as e:
            if os.path.isfile('shelx.ins'):
                print('starting ins2all')
                x = ins2all()
                print(str(e))
                globAtomLabs = x[0]
                globAtomTypes = x[1]
                globAtomAngles = x[2]
                globAtomPos = x[4]

                #Write dicts to cache folder
                newCacheFolderPath = '{}/{}'.format(cachePath, hklHash)
                if os.path.isdir(newCacheFolderPath):
                    rmtree(newCacheFolderPath)

                os.makedirs(newCacheFolderPath)

                with open('{}/atomLabs.buckfast'.format(newCacheFolderPath),'w') as atomLabs:
                    atomLabs.write(str(globAtomLabs))

                with open('{}/atomTypes.buckfast'.format(newCacheFolderPath),'w') as atomTypes:
                    atomTypes.write(str(globAtomTypes))

                with open('{}/atomAngles.buckfast'.format(newCacheFolderPath),'w') as atomAngles:
                    atomAngles.write(str(globAtomAngles))

                with open('{}/atomPos.buckfast'.format(newCacheFolderPath),'w') as atomPos:
                    newPos = {}
                    for item, value in x[4].items():
                        newPos[item] = (tuple(value[0]), value[1])

                    atomPos.write(str(newPos))

            else:
                globAtomLabs = {}
                globAtomTypes = {}
                globAtomAngles = {}
                globAtomPos = {}

    else:
        globAtomLabs = {}
        globAtomTypes = {}
        globAtomAngles = {}
        globAtomPos = {}


'''
#####################################################################
#-------------------CHEMCON------------------------------------------
#####################################################################
'''

@timeDec
def findCHEMCON():
    '''
    Find chemical equivalency in structure. Return CHEMCON dictionary.
    '''
    global globAtomEnv          #Global dictionary of atoms and their chemical environment hash values i.e. {'C(1)':'f11390f0a9cadcbb4f234c8e8ea8d236'}
    global cacheHashPath
    globAtomEnv = {}
    atomLabs = copy.copy(globAtomLabs)
    atoms = [atom + ',asym' for atom in getAtomList()]
    chemconFilePath = '{}/chemcon.buckfast'.format(cacheHashPath)
    atomEnvFilePath = '{}/atomEnv.buckfast'.format(cacheHashPath)

    if os.path.isfile(chemconFilePath) and os.path.isfile(atomEnvFilePath):
        with open(chemconFilePath, 'r') as chemconCache:
            CHEMCON = literal_eval(chemconCache.read())

        with open(atomEnvFilePath,'r') as envCache:
            globAtomEnv = literal_eval(envCache.read())

    else:

        print('CHEMCON not in cache. Finding chemical equivalency...')
        envs = {}
        CHEMCON = {}
        pool = mp.Pool(5)
        atomEnv = pool.map_async(partial(getEnvSig, atomLabs), atoms).get()
        pool.close()
        for item in atomEnv:
            globAtomEnv[spec2norm(item[0])] = item[1]
            envs.setdefault(item[1],[]).append(item[0])
        #Organise hash value dictionary into dictionary of parent atoms that appear first in ATOM table,
        #and children that appear further down.
        hashParents = {}
        with open('xd.mas','r') as mas:
            atomTab = False
            for line in mas:
                if atomTableEnds(line):
                    atomTab = False

                elif atomTab:
                    row = line.upper().split()
                    atomHash = globAtomEnv[row[0]]
                    if atomHash not in hashParents:
                        hashParents[atomHash] = row[0]
                        CHEMCON[row[0]] = []
                    else:
                        CHEMCON[hashParents[atomHash]].append(row[0])

                elif atomTableBegins(line):
                    atomTab = True

        with open(chemconFilePath, 'w') as chemconCache:
            chemconCache.write(str(CHEMCON))

        with open(atomEnvFilePath, 'w') as envCache:
            envCache.write(str(globAtomEnv))
    print(CHEMCON)
    print('CHEMCON\n~~~~~~~')
    return CHEMCON


'''
########################################################################
---------------------REFINEMENTS----------------------------------------
########################################################################
'''

def scaleFacRef():
    '''
    Setup xd.mas to refine scale factors.
    '''
    setupmas()
    resetKeyTable()

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('SELECT cycle'):
                row = line.split()
                row[2] = '-10'
                rowStr = ' '.join(row)
                newmas.write(rowStr + '\n')

            else:
                newmas.write(line)

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def highAngleRef(sinthlMin, sinthlMax):
    '''
    Setup xd.mas to refine high angle non-H positions and ADPs.
    '''
    setupmas()
    resetKeyTable()

    keyTab = False

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('KAPPA'):
                keyTab = False

            if line.startswith('SKIP'):
                row = line.split()
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOn}  {snlMin:<5.3f} {snlMax:<5.3f}'.format(*row, snlOn = '*sinthl', snlMin = sinthlMin, snlMax = sinthlMax)
                newmas.write(rowStr + '\n')

            elif keyTab:
                row = line.split()

                if line[:1] != 'H':
                    rowStr = '{0:8}{1}'.format(row[0], '111 111111 0000000000 000000000000000 00 000 00000 0000000 000000000')
                    newmas.write(rowStr + '\n')
                else:
                    rowStr = '{0:8}{1}'.format(row[0], '000 000000 0000000000 000000000000000 00 000 00000 0000000 000000000')
                    newmas.write(rowStr + '\n')

            elif line.startswith('!RESET BOND'):
                newmas.write(line[1:])
            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def lowAngleRef(sinthlMin,sinthlMax):
    '''
    Setup xd.mas to refine low angle H-positions and isotropic ADPs.
    '''
    setupmas()
    resetKeyTable()

    keyTab = False

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('KAPPA'):
                keyTab = False

            if keyTab:
                row = line.split()

                if line[:1] != 'H':
                    rowStr = '{0:8}{1}'.format(row[0], '000 000000 0000000000 000000000000000 00 000 00000 0000000 000000000')
                    newmas.write(rowStr + '\n')
                else:
                    rowStr = '{0:8}{1}'.format(row[0], '111 100000 0000000000 000000000000000 00 000 00000 0000000 000000000')
                    newmas.write(rowStr + '\n')

            elif line.startswith('RESET BOND'):
                newmas.write('!' + line)
#                newmas.write(line)

            elif line.startswith('SKIP'):
                row = line.split()
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOn}  {snlMin:<5.2f} {snlMax:<5.2f}'.format(*row, snlOn = '*sinthl', snlMin = sinthlMin, snlMax = sinthlMax)
                newmas.write(rowStr + '\n')

            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def setupKappa(mode = 'default'):
    '''
    Setup number of kappa parameters in xd.mas and xd.inp based on CHEMCON.
    '''
    kappas = {}
    kappaWritten = False

    atomTab = False
    i = 0
    x = 0
    Hpresent = False
    inpTable = []
    eleList = []

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if atomTableEnds(line):
                atomTab = False

            if atomTab:

                row = line.split()
                #Find parent atoms of CHEMCON
                if len(row) == 12:
                    if not(line.startswith('H(') and Hpresent):
                        i+=1
                    row[9] = str(i)
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                    newmas.write(rowStr + '\n')
                    kappas[row[0].upper()]= str(i)

                    if line[0:2].isalpha():
                        if line[0:2].upper() not in eleList:
                            x+=1
                            inpTable.append(x)
                            eleList.append(line[0:2].upper())
                        else:
                            inpTable.append(x)

                    else:
                        if line[0:1].upper() not in eleList:
                            x+=1
                            inpTable.append(x)
                            eleList.append(line[0:1].upper())
                            if line[0:1] == 'H':
                                Hpresent = True
                        else:
                            inpTable.append(x)

                else:
                    row[9] = kappas[row[12]]
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)
                    newmas.write(rowStr + '\n')
            #Find KEEP KAPPA
            elif line.startswith('KEEP') or line.startswith('!KEEP'):

                row = line.split()
                if row[1].upper() == 'KAPPA':
                    j = 1
                    rowStr = 'KEEP     KAPPA'

                    while j <= i:
                        rowStr = rowStr + '  ' + str(j)
                        j+=1

                    newmas.write(rowStr + '\n')
                else:
                    newmas.write(line)

            #Add write number of KAPPA lines below key table
            elif line.startswith('KAPPA'):

                if kappaWritten == False:
                    k = 1

                    rowStr = 'KAPPA   000000'
                    while k<= i:
                        newmas.write(rowStr + '\n')
                        k+=1
                    kappaWritten = True

            else:
                newmas.write(line)

            if atomTableBegins(line):
                atomTab = True

        kapInpRes('xd.res', inpTable, i)
        kapInpRes('xd.inp', inpTable, i)

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

    return i


def kapInpRes(fileName, inpTable, i):
    '''
    Setup xd.inp or xd.res with correctly formatted kappa parameters.
    '''
    inpTableBool = False
    Hpresent = False

    #If user hasn't pressed res2inp also update res file so it is copied over correctly to inp
    if os.path.exists(fileName):
        with open(fileName, 'r') as res, open('xdnew.buckfast','w') as newres:

            for line in res:

                if line.startswith('H('):
                    Hpresent = True


                if line.startswith('USAGE'):
                    row = line.split()
                    row[4] = str(i)
                    rowStr = '{0:10}{1:6}{2:3}{3:4}{4:4}{5:4}{6:4}{7:4}{8:4}{9:4}{10:4}{11:6}{12:3}{13:4}{14}'.format(*row)
                    newres.write(rowStr + '\n')


                elif line.startswith('  1  '):
                    #z count to end of inpTable
                    z = 0


                    if Hpresent == True:
                        while z <= (len(inpTable)-2):
                            newres.write('  ' + str(inpTable[z]) + '  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000' + '\n')
                            z+=1
                        newres.write('  ' + str(inpTable[z]) +  '  1.200000  1.200000  1.200000  1.200000  1.200000  1.200000' +'\n')

                    else:
                         while z <= (len(inpTable)-1):
                            newres.write('  ' + str(inpTable[z]) + '  1.000000  1.000000  1.000000  1.000000  1.000000  1.000000' + '\n')
                            z+=1

                    inpTableBool = True

                else:
                    if not inpTableBool:
                        newres.write(line)
                    else:
                        if not line.startswith('  '):
                            inpTableBool = False
                            newres.write(line)

        os.remove(fileName)
        os.rename('xdnew.buckfast',fileName)


def kapMonRef():
    '''
    Setup xd.mas to refine kappa and monopoles.
    '''
    setupmas()
    resetKeyTable()
    i = getEleNum()

    kapInpRes('xd.res', range(1, i+1), i)
    kapInpRes('xd.inp', range(1, i+1), i)

    j=1

    keyTab = False

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('EXTCN'):
                keyTab = False

            if keyTab:
                row = line.split()

                if line[:5] != 'KAPPA':
                    rowStr = '{0:8}{1}'.format(row[0], '000 000000 0000000000 000000000000000 10 000 00000 0000000 000000000')
                    newmas.write(rowStr + '\n')
                else:
                    if j<i:
                        newmas.write('KAPPA   100000' + '\n')
                        j+=1
                    else:
                        newmas.write('KAPPA   000000' + '\n')

            elif line.startswith('SKIP'):
                row = line.split()
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOff}  {8} {9}'.format(*row, snlOff = 'sinthl')
                newmas.write(rowStr + '\n')

            elif line.startswith('!RESET BOND'):
                newLine = line[1:]
                newmas.write(newLine)

            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def multipoleMagician(trackAtom = None):
    '''
    Setup xd.mas to refine multipoles. Return list of atoms for which problems were encountered.
    '''
    setupmas()
    resetKeyTable()
    kapMonRef()
    writeSITESYM(findSITESYM(trackAtom))
    warnings = makeLCS()
    multipoleKeyTable()
    return warnings


def nonHPosADPKey():
    '''
    Setup key table in xd.mas for refinement of multipoles, and non-H positions and ADPs.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        keyTab = False

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('KAPPA'):
                    keyTab = False

            if keyTab:
                row = line.split()

                if line[:1] != 'H':
                    rowStr = '{0:8}{1} {2}'.format(row[0], '111 111111 0000000000 000000000000000', ' '.join(row[5:10]))
                    newmas.write(rowStr + '\n')
                else:
                    rowStr = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', ' '.join(row[5:10]))
                    newmas.write(rowStr + '\n')
            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def posADPMultRef():
    '''
    Setup xd.mas to refine multipoles, and non-H positions and ADPs.
    '''
    #For an unknown reason calling KapMonRef messes up
    warnings = multipoleMagician()
    nonHPosADPKey()

    return warnings

'''
########################################################################
#--------------------LOWER SYMMETRY-------------------------------------
########################################################################
'''

def mm2Tom():
    '''
    Convert mm2 and cyl to m in atom table.
    '''
    XAtoms = ('F(','CL','BR','I(', 'O(', 'N(')
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        atomTab = False

        for line in mas:

            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.split()

                if row[11] == 'mm2':
                    row[11] = 'm'

                if row[0][:2] in XAtoms and row[11] == 'cyl':
                    row[11] = 'm'

                if len(row) == 12:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                else:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)

                newmas.write(rowStr + '\n')

            else:
                newmas.write(line)

            if atomTableBegins(line):
                atomTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def lowerSym():
    '''
    Setup xd.mas to refine multipoles of lowered symmetry, and non-H positions and ADPs.
    '''
    writeLocalCoordSys(getLowLocalCoordSys())
    mm2Tom()
    multipoleKeyTable()
    nonHPosADPKey()

def lowerSymTo1():
    '''
    Lower all symmetry to 1 and remove all CHEMCON from xd.mas.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        keyTab = False
        atomTab = False

        for line in mas:

            if line.startswith('KAPPA'):
                keyTab = False

            elif line.startswith('END ATOM'):
                atomTab = False

            if keyTab:
                row = line.split()

                if line[:2] != 'H(':
                    if line.startswith('C('):
                        rowStr = '{0:8}{1} {2}'.format(row[0], ' '.join(row[1:5]), '10 111 11111 1111111 000000000' )
                        newmas.write(rowStr + '\n')
                    else:
                        rowStr = '{0:8}{1} {2}'.format(row[0], ' '.join(row[1:5]), '10 111 11111 1111111 111111111' )
                        newmas.write(rowStr + '\n')
                else:
                    rowStr = '{0:8}{1}'.format(row[0], ' '.join(row[1:10]))
                    newmas.write(rowStr + '\n')

            elif atomTab:
                newline = line[:69] + '\n'
                newmas.write(newline)

            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

            elif atomTableBegins(line):
                atomTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def getLowLocalCoordSys():
    '''
    Create local coordinate system for each atom based on low symmetry.
    Return local coordinate systems.
    '''
    atomNeebsRaw = copy.copy(globAtomLabs)   #Get dictionary of atoms and their nearest neighbour labels i.e. {C(2):['C(3)','C(4)','H(2)','H(3)']}
    atomSymsRaw = findSITESYM()         #Get dictionary of atoms and their local symmetry labels and their configuration of different/identical atoms i.e. ('mm2',2+2) could mean a C bonded to 2Hs and 2Cs
    atomLocCoords = {}                      #Initialise dictionary that will be returned and passed to writeLocalCoordSys format {'C(2)':['DUM0','Z','C(3),'Y']}
    i = 1	                                  #Default value for dumI in addDUM(atom1,atom2,dumI). addDUM will find the lowest unused index itself.

    atomNeebs = {}
    atomSyms = {}
    #Convert ',asym' to xd atom label and create dummys of other atoms
    for atom, sym in atomSymsRaw.items():

        atomSyms[atom] = list(sym[:2])
        for i,item in enumerate(sym[2:]):
            i+=2                                        #correct index to correspond with sym
            if isinstance(item, str):
                splitLab = item.split(',')
                if splitLab[1] != 'asym':
                    atomSyms[atom].append(spec2masLab(item))
                else:
                    atomSyms[atom].append(splitLab[0])
            else:
                atomSyms[atom].append([])
                for item2 in item:
                    splitLab = item2.split(',')
                    if splitLab[1] != 'asym':
                        atomSyms[atom][i].append(spec2masLab(item2))
                    else:
                        atomSyms[atom][i].append(splitLab[0])


    for atom, neebs in atomNeebsRaw.items():
        atomNeebs[atom] = []
        for neeb in neebs:
            splitLab = neeb.split(',')
            if splitLab[1] == 'asym':
                atomNeebs[atom].append(splitLab[0])
            else:
                atomNeebs[atom].append(spec2masLab(neeb))

    for atom,sym in atomSyms.items():       #Go through atomSym dict and do different things depending on nearest neighbour configuration

        symTag = sym[1]
        atom = atom.upper()
        #For all atoms all different configurations up to 4-coordinate are checked for, and a different local coordinate system is setup for each possible configuration.
        #2 identical neighbours mm2
        if symTag[:2] == '1 ':
            if symTag == '1 m':
                atomLocCoords[atom] = atomNeebs[atom][0]

                yAtom = ''
                for atom2 in atomNeebs[atomNeebs[atom][0]]:
                    if atom2 != atom:
                        yAtom = atom2
                        break
                atomLocCoords[atom] = [atomNeebs[atom][0], 'X', yAtom, 'Y']


        elif symTag == '2':
            #Mirror plane along 2 neighbour atoms
            atomLocCoords[atom]=[atomNeebs[atom][1], 'X', atomNeebs[atom][0],'Y']

        #2 different neighbours m
        elif symTag == '11':
            atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

        #3 identical neighbours m
        #If 3 neighbours and 3 angles all ~120, assign m (could be 3m but how would axes be setup)
        #If 3 neighbours and 1+1+1 angles and sum ~360, assign m
        #If 3 neighbours and 2+1 angles and sum ~360, assign mm2
        #If 3 neighbours and 3 angles but sum < 360, assign m (could be 3m but how would axes be setup?)
        #If 3 neighbours and 1+1+1 angles and sum < 360, assign 1
        #If 3 neighbours and 2+1 angles and sum < 360, assign m
        elif symTag[:2] == '3 ':

            if symTag[-1:] == 'p':                                 #If 3 neighbours in plane assign m in plane
                atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

            elif symTag == '3 3n':                                              #Think NH3
                x = addDUM(atomNeebs[atom][0],atomNeebs[atom][1],i)            #Make X-axis in centre of arbritary atoms
                atomLocCoords[atom] = ['DUM' + str(x), 'X', atomNeebs[atom][2],'Y'] #X-axis in centre of two atoms then Y-axis on opposite atoms so Z is perpendicular to m

            elif symTag == '3 21n':
                x = addDUM(sym[3][0],sym[3][1],i)            #Make Z-axis in centre of 1 angle.
                atomLocCoords[atom] = ['DUM' + str(x), 'X', sym[2],'Y'] #Opposite atom to dummy atom needed for m perpendicular to Z

            #1+1+1 neighbours: '111 p' gives m, '111 n' gives 1 symmetry
        elif symTag == '111 p':
            atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

            #2+1 different neighbours mm2
                #If 2+1 neighbours and 3 angles and sum ~360, assign mm2
                #If 2+1 neighbours and 1+1+1 angles and sum ~360, assign m
                #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle is opposite odd neighbour, assign mm2
                #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle isn't opposite odd neighbour, assign m
                #If 2+1 neighbours and 3 angles and sum < 360, assign m
                #If 2+1 neighbours and 1+1+1 angles and sum < 360, assign 1
                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle is opposite odd neighbour, assign m
                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle isn't opposite odd neighbour, assign 1

        elif symTag[:3] == '21 ':
            if symTag[-1:] == 'p' or symTag == '21 21pmm2':                                 #If 3 neighbours in plane assign m in plane
                atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

            elif symTag == '21 3n' or symTag == '21 21nm':         #If 2+1 neighbours and 3 angles and sum < 360, assign m
                                                                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle is opposite odd neighbour, assign m
                x = addDUM(sym[2][0],sym[2][1],i)            #Make X-axis between 2 atoms
                atomLocCoords[atom] = ['DUM' + str(x), 'X', sym[3],'Y']

        if symTag[:2] == '4 ':
            if symTag == '4 t':
                #assign m through arbitrary 2 atoms
                atomLocCoords[atom]=[atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

            elif symTag == '4 sp':
                atomLocCoords[atom] = [atomNeebs[atom][0], 'X', atomNeebs[atom][1], 'Y']

            elif symTag == '4 t/sp':
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

            elif symTag == '4 et':
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

        #2+1+1
        elif symTag[:4] == '211 ':
            if symTag in ('211 t','211 t oa'):
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

            elif symTag == '211 spmm2':
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

            elif symTag == '211 spm':
                atomLocCoords[atom] = [atomNeebs[0], 'X', atomNeebs[1], 'Y']

            elif symTag == '211 t/sp':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

        elif symTag[:3] == '31 ':
            if symTag == '31 t':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 sp':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 t/sp':

                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 epy':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 et':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 toa':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

        #2+2
        elif symTag[:3] == '22 ':

            if symTag == '22 t':
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

            elif symTag == '22 sp':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '22 t/spmm2':
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

            elif symTag == '22 etmm2':
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

        #1+1+1+1
        elif symTag == '1111 sp':
            atomLocCoords[atom] = [atomNeebs[atom][0], 'X', atomNeebs[atom][1], 'Y']

        elif sym[0] == '1':                      #Define arbritrary system for '1' symmetry so it is not returned as not being defined
            atomLocCoords[atom] = [atomNeebs[atom][0], 'Z', atomNeebs[atom][1], 'Y']

    return (atomLocCoords)

def getMultRes():
    '''
    Get populations of all multipoles from xd.res. Return populations.
    '''
    atomLabs = copy.copy(globAtomLabs)
    i = 0
    multPops = {}
    atomFound = False
    try:
        if os.path.isfile(os.getcwd() + '\\xd.res'):
            res = open('xd.res', 'r')
        else:
            res = open('xd.inp', 'r')

        for line in res:
            row = line.split()

            if row[0].upper() in atomLabs.keys():
                i = 0
                multPops[row[0].upper()] = []
                atomFound = True
                parentAtom = row[0].upper()

            if 1 < i < 5 and atomFound:
                popFloats = [float(item) for item in row]
                multPops[parentAtom].extend(popFloats)

            i += 1

    finally:
        res.close()

    return multPops

def checkMultRes():
    '''
    Setup xd.mas to refine only multipoles that have significant populations.
    '''
    multPops = getMultRes()
    multI = {}

    for atom, pops in multPops.items():
        multI[atom] = []
        for pop in pops:
            if abs(pop) > 0.02:
                multI[atom].append('1')
            else:
                multI[atom].append('0')

    with open('xd.mas', 'r') as mas, open('xdnew.mas', 'w') as newmas:

        keyTab = False
        j = 0
        for line in mas:
            j = 0
            if line.startswith('KAPPA'):
                keyTab = False

            if keyTab:
                row = line.split()

                if line[:1] != 'H':
                    rowStr = line[:46]
                    for char in line[46:-1]:                #-2 there to not include \n at end of line

                        if char != ' ':
                            rowStr += str(multI[row[0].upper()][j])
                            j+=1
                        else:
                            rowStr += char

                    newmas.write(rowStr + '\n')
                else:
                    rowStr = '{0:8}{1} {2}\n'.format(row[0], ' '.join(row[1:5]), '10 001 00000 0000000 000000000' )
                    newmas.write(rowStr)
            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


'''
########################################################################
#--------------------SITESYM--------------------------------------------
########################################################################
'''

def findSITESYM(trackAtom = 'C(01A)'):
    '''
    Find local symmetry for all atoms based on nearest neighbours and bond angles.
    Return local symmetry for every atom.
    '''
    #Get nearest neighbours from findNeighbours()
    atoms = copy.copy(globAtomTypes)
    atomLabs = copy.copy(globAtomLabs)
    bondAngles = copy.copy(globAtomAngles)
    atomSyms = {}

    #Tag system: first number is neighbour config, second is angle config, p = planar, n = not planar

    #Work out symmetry and add it to list based on neighbour atoms
    for atom,neighbours in atoms.items():

        atom = atom.upper()

        #Atom bonded to 1 neighbour = cyl
        if len(neighbours) == 1:
            if atom[:2] == 'H(':
                atomSyms[atom] = ('cyl','1 c')
            else:
                atomSyms[atom] = ('cyl','1 m')


        #Atom bonded to 2 identical neighbours = mm2, 2 different neighbours = m
        elif len(neighbours) == 2:
            if neighbours[0] == neighbours[1]:
                atomSyms[atom]= ['mm2','2']
            else:
                atomSyms[atom] = ['m','11']

        #Atom bonded to 3 identical neighbours = m, 2 identical neighbours = mm2, all different neighbours = m
        elif len(neighbours) == 3:

            angles = sorted([float(entry[0]) for entry in bondAngles[atom]])

            #3 neighbours
            #If 3 neighbours and 3 angles all ~120, assign m (could be 3m but how would axes be setup)
            #If 3 neighbours and 1+1+1 angles and sum ~360, assign m
            #If 3 neighbours and 2+1 angles and sum ~360, assign mm2
            #If 3 neighbours and 3 angles but sum < 360, assign m (could be 3m but how would axes be setup?)
            #If 3 neighbours and 1+1+1 angles and sum < 360, assign 1
            #If 3 neighbours and 2+1 angles and sum < 360, assign m
            if neighbours[1:] == neighbours[:-1]:

                #If 3 atoms in plane
                if 356 < sum(angles):

                    #If all angles are approximately equal, assign mm2 symmetry (could be 3m but how would axes be setup).
                    if max(angles) - min(angles) < 4:

                        atomSyms[atom] = ['mm2','3 3p']

                    #Check for 2+1 angles
                    elif check421angles(angles)[0]:
                        oddAngle = check421angles(angles)[1]
                        oddAngleAtoms = [(entry[1][0],entry[1][2]) for entry in bondAngles[atom] if entry[0] == oddAngle]

                        for atomLab in atomLabs[atom]:
                            if atomLab not in oddAngleAtoms[0]:
                                ZAtom = atomLab
                                break

                        atomSyms[atom] = ('mm2', '3 21p', ZAtom, oddAngleAtoms[0][0])

                    #If 1+1+1 angles, assign m symmetry
                    else:
                        atomSyms[atom] = ['m','3 111p']

                #If 3 atoms not in plane
                elif sum(angles) < 356:

                    #If all angles equal assign 'm' (could be 3m but how would axes work?).
                    if max(angles) - min(angles) < 7:
                        atomSyms[atom] = ['m', '3 3n']

                    #Check for 2+1 angles
                    elif check421angles(angles)[0]:
                        oddAngle = check421angles(angles)[2]
                        oddAngleAtoms = [(entry[1][0],entry[1][2]) for entry in bondAngles[atom] if entry[0] == oddAngle]

                        for atomLab in atomLabs[atom]:
                            if atomLab not in oddAngleAtoms:
                                XAtom = atomLab
                                break

                        atomSyms[atom] = ('m', '3 21n', XAtom, (oddAngleAtoms[0],oddAngleAtoms[2]))

                    #1+1+1 angles detected and '1' assigned.
                    else:
                        atomSyms[atom] = ['1', '3 111n']

            #1+1+1 neighbours
            elif len(set(neighbours)) == 3:
                #If 1+1+1 and bond angles are all different and sum ~360, assign m symmetry (atoms in plane)
                #If 1+1+1 and 2+1 bond angles and sum ~360, assign m symmetry (atoms in plane)
                #If 1+1+1 and all bond angles are between 118-122, m symmetry in horizontal plane
                #Simplification of 3 situations: If sum of bond angles ~360, atoms are in plane so m symmetry.
                if 356 < sum(angles):
                    atomSyms[atom]= ['m','111 p']
                #If 1+1+1 and bond angles are all the same but not ~120, assign 1 symmetry (atoms not in plane)
                #If 1+1+1 and bond angles are all different and sum < 360, assign 1 symmetry (atoms not in plane)
                #If 1+1+1 neighbours and 2+1 bond angles and sum < 360, assign 1 symmetry (atoms not in plane)
                #Simplification of 3 situations: If sum of bond angles < 360, atoms not in plane so 1 symmetry
                elif sum(angles) < 356:
                    atomSyms[atom] = ['1','111 n']

            #2+1 neighbours
                #If 2+1 neighbours and 3 angles and sum ~360, assign mm2
                #If 2+1 neighbours and 1+1+1 angles and sum ~360, assign m
                #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle is opposite odd neighbour, assign mm2
                #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle isn't opposite odd neighbour, assign m
                #If 2+1 neighbours and 3 angles and sum < 360, assign m
                #If 2+1 neighbours and 1+1+1 angles and sum < 360, assign 1
                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle is opposite odd neighbour, assign m
                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle isn't opposite odd neighbour, assign 1
            else:
                #Make set of all angles rounded to 5
                oddAtom = findOddAtom(neighbours, atomLabs[atom])

                #If atoms in plane
                if 356 < sum(angles):
                    #3 angles so assign mm2
                    if max(angles) - min(angles) < 4:
                        atomSyms[atom] = ['mm2','21 3p', oddAtom]

                    #Check for 2+1 angles
                    elif check421angles(angles)[0]:

                        oddAngle = check421angles(angles)[1]
                        oddAngleAtoms = [(entry[1][0],entry[1][2]) for entry in bondAngles[atom] if entry[0] == oddAngle]

                        if oddAtom not in oddAngleAtoms[0]:
                            atomSyms[atom] = ('mm2', '21 21pmm2', oddAtom, oddAngleAtoms[0])
                        else:
                            atomSyms[atom] = ('m', '21 21pm')

                    #If 1+1+1 angles
                    else:

                        atomSyms[atom] = ['m', '21 111p']

                #If atoms not in plane
                else:
                    #3 angles so assign m
                    if max(angles) - min(angles) < 4:

                        XDumAtoms = []

                        for atomLab in atomLabs[atom]:
                            if atomLab != oddAtom:
                                XDumAtoms.append(atomLab)

                        atomSyms[atom] = ['m','21 3n', XDumAtoms, oddAtom]

                    #Check for 2+1 angles
                    elif check421angles(angles)[0]:

                        oddAngle = check421angles(angles)[1]

                        oddAngleAtoms = [(entry[1][0],entry[1][2]) for entry in bondAngles[atom] if entry[0] == oddAngle]

                        if oddAtom not in oddAngleAtoms[0]:
                            atomSyms[atom] = ('m', '21 21nm', [oddAngleAtoms[0][0], oddAngleAtoms[0][1]], oddAtom)
                        else:
                            atomSyms[atom] = ('1', '21 21n1')

                    else:
                        atomSyms[atom] = ['1', '21 111n']

#--------------------------------4 NEIGHBOURS------------------------------------------------------
        #t = tetrahedral, sp = square planar, t/sp = inbetween tetrahedral and square planar, et = elongated tetrahedral, epy = elongated pyramid
        #Atom bonded to 4 identical neighbours = mm2, 1+1+2 neighbours = m, 2+2 neighbours = mm2, 1+3 neighbours = m, all different neigbours = 1
        elif len(neighbours) == 4:
            angles = [float(entry[0]) for entry in bondAngles[atom]]
            i = 0
            squPla = False
            sortedAngles = sorted(angles)

            #Check for square planar
            for angle in angles:            #If two angles are ~180 then coordination is square planar
                if 177 < angle:
                    i+=1
            if i == 2:
                squPla = True

            #4 identical neighbours
            if neighbours[1:] == neighbours[:-1]:

                if max(angles) - min(angles) < 7:
                    atomSyms[atom] = ['mm2','4 t']    #6 angles ~109.5, assign mm2, tetrahedral

                elif squPla:
                    atomSyms[atom] = ['mm2','4 sp']     #square planar, assign mm2 (could use 4m but how would you define 4-axis)

                elif max(sortedAngles[:4]) - min(sortedAngles[:4]) < 8 and max(sortedAngles[4:]) - min(sortedAngles[4:]) < 8:   #Get distortion of tetrahedral towards square planar

                    ZDums = ()

                    for angle in bondAngles[atom]:

                        if float(angle[0]) in sortedAngles[4:]:
                            ZDums = (angle[1][0], angle[1][2])
                            break

                    for atomLab in atomLabs[atom]:
                        if atomLab not in ZDums:
                            YAtom = atomLab
                            break

                    atomSyms[atom] = ['mm2', '4 t/sp', ZDums, YAtom]

                elif max(sortedAngles[:2]) - min(sortedAngles[:2]) < 8 and max(sortedAngles[2:]) - min(sortedAngles[2:]) < 8:  #Get elongated tetrahedral geometry

                    ZDums = ()

                    for angle in bondAngles[atom]:

                        if float(angle[0]) in sortedAngles[:2]:
                            ZDums = (angle[1][0], angle[1][2])
                            break

                    for atomLab in atomLabs[atom]:

                        if atomLab not in ZDums:
                            YAtom = atomLab
                            break

                    atomSyms[atom] = ['mm2', '4 et', ZDums, YAtom]




            #All different
            elif len(set(neighbours)) == 4:
                if not squPla:
                    atomSyms[atom] = ['1','1111']            #Otherwise assign 1
                else:
                    atomSyms[atom] = ['m','1111 sp']          #If square planar assign m

            #2+1+1 neighbours
            #If square planar, assign mm2 if 1+1 are opposite each other, otherwise m
            #If tetrahedral, assign m   ???????
            #If square based pyramid, assign  ???????
            #All angles different???
            elif len(set(neighbours)) == 3:

                twoatoms = findIdentAtoms(neighbours, atomLabs[atom], 2)[0]  #Find 2atoms in 211neighbours
                xyatoms = []                                          #Get 11 atoms for xy axes
                angleList = []

                for atomLab in atomLabs[atom]:
                    if atomLab not in twoatoms:
                        xyatoms.append(atomLab)

                for angle in bondAngles[atom]:
                    if xyatoms[0] in angle[1] and xyatoms[1] in angle[1]:
                        pass
                    else:
                        angleList.append(float(angle[0]))

                if max(angles) - min(angles) < 7:       #If tetrahedral assign m

                      atomSyms[atom] = ['m', '211 t', xyatoms]

                elif max(angleList) - min(angleList) < 7:               #All angles the same except 1+1 angle different
                    atomSyms[atom] = ['m', '211 t oa', xyatoms]

                elif squPla:                                                #If square planar
                    oppAtoms = []                                           #List to store atoms opposite each other
                    twoatoms = findIdentAtoms(neighbours,atomLabs[atom],2)[0]    #Get 2 atoms in 211 neighbours
                    opp11 = False                                           #Bool to store if 11 atoms are opposite each other

                    for angle in bondAngles[atom]:                      #Go through angles and check for 180 i.e. opposites atoms
                        if 177 < angle[0]:
                            oppAtoms.append((angle[1][0],angle[1][2]))  #Store list of tuples of opposite atoms

                    for pair in oppAtoms:
                        if frozenset(pair) == frozenset(twoatoms):    #If 2 atoms are opposite then 11opp = True
                            opp11 == True
                            break

                    if opp11:
                        for atomLab in atomLabs[atom]:         #If 11atoms are opposite assign mm2
                            if atomLab not in twoatoms:           #And use 11atom for Z-axis and 2atom for Y axis
                                ZAtom = atomLab
                                break
                        atomSyms[atom] = ['mm2', '211 spmm2', (ZAtom, twoatoms[0])]

                    else:                                   #If 11atoms are next to each other assign m
                        atomSyms[atom] = ['m', '211 spm']

                elif max(sortedAngles[:4]) - min(sortedAngles[:4]) < 7 and max(sortedAngles[4:]) - min(sortedAngles[4:]) < 7:   #Get distortion of tetrahedral towards square planar
                    xyatoms = findIdentAtoms(neighbours,atomLabs[atom],1)

                    atomSyms[atom] = ['m', '211 t/sp', xyatoms[0][0], xyatoms[1][0]]

                elif max(sortedAngles[:2]) - min(sortedAngles[:2]) < 8 and max(sortedAngles[2:]) - min(sortedAngles[2:]) < 8:  #Get elongated tetrahedral geometry

                    xyatoms = findIdentAtoms(neighbours, atomLabs[atom], 1)
                    atomSyms[atom] = ['m', '211 et', xyatoms[0], xyatoms[1]]

            #2+2 and 3+1 neighbours
            elif len(set(neighbours)) == 2:
                #Find duplicity in neighbours
                x = dict(Counter(neighbours))
                #3+1 neighbours
                if 3 in x.values():
                    if max(angles) - min(angles) < 7:
                        XAtom = findOddAtom(neighbours, atomLabs[atom])
                        threeAtoms = findIdentAtoms(neighbours, atomLabs[atom], 3)
                        atomSyms[atom] = ['m','31 t', XAtom, threeAtoms[0][0]]

                    elif squPla:
                        ZAtom = findOddAtom(neighbours, atomLabs[atom])

                        for angle in bondAngles[atom]:
                            if ZAtom in angle[1] and angle[0] < 177:
                                for atomLab in angle[1]:
                                    if atomLab not in [atom, ZAtom]:
                                        YAtom = atomLab
                                        break

                        atomSyms[atom] = ['mm2', '31 sp', ZAtom, YAtom]

                    elif max(sortedAngles[:4]) - min(sortedAngles[:4]) < 7 and max(sortedAngles[4:]) - min(sortedAngles[4:]) < 7:   #Get distortion of tetrahedral towards square planar
                        xAtom = findIdentAtoms(neighbours,atomLabs[atom],1)[0]
                        for bondAngle in bondAngles[atom]:
                            if bondAngle[0] >= min(sortedAngles[4:]):
                                for atomLab in bondAngle[1]:
                                    if atomLab not in [atom,xAtom[0]]:
                                        yAtom = atomLab
                                        break

                        atomSyms[atom] = ['m', '31 t/sp', xAtom[0], yAtom]

                    elif max(sortedAngles[:3]) - min(sortedAngles[:3]) < 8 and max(sortedAngles[3:]) - min(sortedAngles[3:]) < 8:  #Get Y3-C-X with different Y-X and Y-Y angles

                            xAtom = findIdentAtoms(neighbours,atomLabs[atom],1)[0]
                            yAtom = findIdentAtoms(neighbours,atomLabs[atom],3)
                            atomSyms[atom] = ['m', '31 epy', xAtom[0], yAtom[0][0]]

                    elif max(sortedAngles[:2]) - min(sortedAngles[:2]) < 8 and max(sortedAngles[2:]) - min(sortedAngles[2:]) < 8:  #Get elongated tetrahedral geometry

                        xAtom = findIdentAtoms(neighbours, atomLabs[atom], 1)[0]

                        for bondAngle in bondAngles[atom]:
                            if bondAngle[0] <= min(sortedAngles[:2]) and xAtom in bondAngle[1]:
                                for atomLab in bondAngle[1]:
                                    if atomLab not in [atom,xAtom[0]]:
                                        yAtom = atomLab
                                        break
                        atomSyms[atom] = ['m', '31 et', xAtom[0][0], yAtom]

                    elif (max(sortedAngles[:5]) - min(sortedAngles[:5]) < 8):
                        yAtom = findIdentAtoms(neighbours,atomLabs[atom],1)[0][0]

                        for bondAngle in bondAngles[atom]:
                            if bondAngle[0] > sortedAngles[4]:
                                oaAtoms = (bondAngle[1][0],bondAngle[1][2])

                                break
                        for atomLab in atomLabs[atom]:
                            if atomLab not in [yAtom,oaAtoms]:
                                xAtom = atomLab
                                break

                        atomSyms[atom] = ['mm2', '31 toa', xAtom, yAtom]

                    elif (max(sortedAngles[1:]) - min(sortedAngles[1:]) < 8):
                        yAtom = findIdentAtoms(neighbours,atomLabs[atom],1)[0][0]

                        for bondAngle in bondAngles[atom]:
                            if bondAngle[0] < sortedAngles[1]:
                                oaAtoms = (bondAngle[1][0],bondAngle[1][2])
                                break
                        for atomLab in atomLabs[atom]:
                            if atomLab not in [yAtom,oaAtoms[0],oaAtoms[1]]:
                                xAtom = atomLab
                                break

                        atomSyms[atom] = ['mm2', '31 toa', xAtom, yAtom]


                #2+2 neighbours
                else:
                    #2+2 tetrahedral
                    if max(angles) - min(angles) < 7:
                        atoms22 = findIdentAtoms(neighbours, atomLabs[atom], 2)   #Returns [[C(1),C(2)],[H(1),H(2)]]
                        atomSyms[atom] = ['mm2', '22 t', atoms22[0], atoms22[1][0]]

                    #2+2 square planar
                    elif squPla:
                        atoms22 = findIdentAtoms(neighbours, atomLabs[atom], 2)         #Returns [[C(1),C(2)],[H(1),H(2)]]
                        atomSyms[atom] = ['mm2', '22 sp', atoms22[0], atoms22[1][0]]    #Should be correct wherever 2+2 atoms are

                    elif max(sortedAngles[:4]) - min(sortedAngles[:4]) < 7 and max(sortedAngles[4:]) - min(sortedAngles[4:]) < 7:   #Get distortion of tetrahedral towards square planar

                        atoms22 = findIdentAtoms(neighbours, atomLabs[atom], 2)
                        oppAtoms = False

                        for bondAngle in bondAngles[atom]:
                            if bondAngle[0] >= min(sortedAngles[4:]):
                                    if (atoms22[0][0] in bondAngle[1] and atoms22[0][1] in bondAngle[1]) or (atoms22[1][0] in bondAngle[1] and atoms22[1][1] in bondAngle[1]):
                                        oppAtoms = True
                                        break
                        if oppAtoms:
                            atomSyms[atom] = ['mm2', '22 t/spmm2', atoms22[0], atoms22[1][0]]
                        else:
                            atomSyms[atom] = ['1', '22 t/sp1']

                    elif max(sortedAngles[:2]) - min(sortedAngles[:2]) < 8 and max(sortedAngles[2:]) - min(sortedAngles[2:]) < 8:  #Get elongated tetrahedral geometry

                        oppAtoms = False
                        atoms22 = findIdentAtoms(neighbours, atomLabs[atom], 2)

                        for bondAngle in bondAngles[atom]:
                            if bondAngle[0] <= min(sortedAngles[:2]):
                                if (atoms22[0][0] in bondAngle[1] and atoms22[0][1] in bondAngle[1]) or (atoms22[1][0] in bondAngle[1] and atoms22[1][1] in bondAngle[1]):
                                    oppAtoms = True
                                    break

                        if oppAtoms:
                            ZDums = (atoms22[0][0],atoms22[0][1])

                            atomSyms[atom] = ['mm2', '22 etmm2', ZDums, atoms22[1][0]]                   #If 2+2atoms form bite angle assign mm2

                        else:
                            atomSyms[atom] = ['1', '22 et1', atomLabs[atom][0], atomLabs[atom][1]]    #Otherwise assign 1

#-----------------------------6 NEIGHBOURS---------------------------------------------------------
#        elif len(neighbours) == 6:
#
#            angles = bondAngles[atom]
#            sortedAngles = sorted(angles)
#
#            #6 neighbours
#            if len(set(neighbours)) == 1:
#                #Approximate octahedron
#                if min(sortedAngles[11:]) > 170 and max(sortedAngles[11:] < 190):
#                    if min(sortedAngles[:11]) > 80 and max(sortedAngles[11:]) < 100:
#
#                        atomSyms[atom] = ['4m', '6 o']

        if atom == trackAtom:
            print(atom)
            print(neighbours)
            print(bondAngles[atom])
            print(atomSyms[atom])

    return atomSyms

def check421angles(angles):
    '''
    Check 3 angles to see if 2 are very similar and 1 is different. Return result as bool.
    '''
    angles21 = False
    oddAngle = None
    #Check for 2+1 angles
    if abs(angles[0] - angles[1]) < 4 and not abs(angles[0] - angles[2]) < abs(angles[0] - angles[1]) and not abs(angles[1] - angles[2]) < abs(angles[0] - angles[1]):
        angles21 = True
        oddAngle = angles[2]
    elif abs(angles[0] - angles[2]) < 4 and not abs(angles[0] - angles[1]) < abs(angles[0] - angles[2]) and not abs(angles[1] - angles[2]) < abs(angles[0] - angles[2]):
        angles21 = True
        oddAngle = angles[1]
    elif abs(angles[1] - angles[2]) < 4 and not abs(angles[0] - angles[1]) < abs(angles[1] - angles[2]) and not abs(angles[0] - angles[2]) < abs(angles[1] - angles[2]):
        angles21 = True
        oddAngle = angles[0]

    return (angles21, oddAngle)

def findOddAtom(neighbourTypes, neighbourLabs):
    '''
    Find an atom which is the only one of its kind in a list of nearest neighbours. Return atom label.
    '''
    oddElement = ''.join([item for item in neighbourTypes if neighbourTypes.count(item) == 1])

    oddAtom = ''
    if len(oddElement) == 1:
        #Use odd element to find odd atom
        for atomLab in neighbourLabs:
            if atomLab[:1] == oddElement:
                oddAtom = atomLab
                break
    else:
        for atomLab in neighbourLabs:
            if atomLab[:2] == oddElement:
                oddAtom = atomLab
                break

    return oddAtom

def findIdentAtoms(neighbourTypes, neighbourLabs, countNum):
    '''
    Find atoms in a list of nearest neighbours with a given count in this list. Return atoms.
    '''
    oddElements = []

    for item in neighbourTypes:
        if neighbourTypes.count(item) == countNum:
            oddElements.append(item)

    atoms = []

    for i,element in enumerate(oddElements):
        atoms.append([])
        if len(element) == 1:
        #Use odd element to find odd atom
            for atomLab in neighbourLabs:
                if atomLab[:1] == element:
                    atoms[i].append(atomLab)

        else:
            for atomLab in neighbourLabs:
               if atomLab[:2] == element:
                    atoms[i].append(atomLab)

    return atoms

def writeSITESYM(atomSymDict):
    '''
    Adds SITESYM column to atom table in xd.mas. Returns any atoms for which SITESYM wasn't added.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        atomTab = False
        SITESYMs = atomSymDict
        unaddedSym = []

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.split()

                #If atom is in dictionary it is added with sym label
                if row[0].upper() in SITESYMs.keys():

                    row[11] = SITESYMs[row[0].upper()][0]

                    if len(row) == 12:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                    else:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)
                    newmas.write(rowStr + '\n')

                #Add cyl symmetry for Hs
                elif row[0][0]=='H':

                    row[11] = 'cyl'

                    if len(row)==12:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                    else:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)
                    newmas.write(rowStr + '\n')

                else:
                    newmas.write(line)
                    unaddedSym.append(row[0])

            else:
                newmas.write(line)

            if atomTableBegins(line):
                atomTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

    return unaddedSym

def findUnaddedSym():
    '''
    Find any atoms with 'NO' in the SITESYM column of the atom table. Return atoms.
    '''
    with open('xd.mas', 'r') as mas:

        atomTab = False
        noSymAtoms = []

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.split()

                #If sym label is NO add to list of atoms with unadded SITESYM
                if len(row)==11 or row[11] == 'NO':
                    noSymAtoms.append(row[0])

            if atomTableBegins(line):
                atomTab = True

    return noSymAtoms


'''
########################################################################
#--------------------LOCAL COORDINATE SYSTEM----------------------------
########################################################################
'''

def spec2masLab(specLab):
    '''
    Convert 'C(1),asym' type label to C(1).
    Convert 'C(1),dum42' type label to dummy atom.
    Return new label of atom.
    '''
    splitLab = specLab.split(',')

    if splitLab[1] == 'asym':
        newLab = splitLab[0]
    else:
        newLab = addDUMsym(specLab)

    return newLab

def addDUMsym(atomLab, dumI = 1):
    '''
    Add dummy atom to xd.mas on atomic position of given atom. Return dummy atom label.
    '''
    global globAtomPos
    #Open "xd.inp" to read
    with open('xd.mas','r') as mas:

        addedDUMs = {}                          #Variable to know what dummy atoms are already added

        #Go through mas file line by line
        #Find any dummy atoms already there
        #[:-1] removes \n from end, save every dummy atom coords as key with index as value
        for line in mas:
            if line.startswith('DUM'):

                row = line.split()
                addedDUMs[tuple(row[1:])] = row[0][3:]

            elif line.startswith('END ATOM'):
                break

        dumPos = globAtomPos[atomLab][0]
        #Find the average of the x,y and z fractional coordinates of the 2 atoms
        dumX = dumPos[0]
        dumY = dumPos[1]
        dumZ = dumPos[2]

    #Open xd.mas and create new mas file xdnew.mas
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Reprint xd.mas and add dummy atom in correct place to xdnew.mas
        for line in mas:
            if 'END ATOM' in line:                      #Find correct place in xd.mas to write dummy atoms

                newCoords = [dumX,dumY,dumZ]
                newCoords = tuple(('{0:.4f}'.format(item) for item in newCoords))

                if newCoords in addedDUMs:      #See if new coordinates belong to an existing dummy atom

                    dumI = addedDUMs[newCoords] #If they do, just return the index of that dummy atom

                else:                                               #If they don't, create a new dummy atom
                    while str(dumI) in addedDUMs.values():          #While dumI is already a dummy atom in the mas file, increment it by 1
                        dumI = int(dumI) + 1
                    #When a dummy atom index not already in the mas file is found, write the new dummy atom in the correct format
                    dumAtom = '!DUM{0} is on a {4} position\nDUM{0:<5}{1:8.4f}{2:8.4f}{3:8.4f}\n'.format(dumI,dumX,dumY,dumZ,atomLab.split(',')[0])
                    newmas.write(dumAtom)

            newmas.write(line)                          #Write every other line unchanged

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

    return ('DUM' + str(dumI))                                   #Return dumI so that getLocalCoordSys knows what dummy atom label to use

def neebDict2masLabs(atomNeebsRaw):
    '''
    Take raw neighbour dictionary and convert every label to either a normal mas label or a dummy atom label.
    Also update global environment dictionary with new labels.
    Return dictionary.
    '''
    atomNeebs = {}

    for atom, neebs in atomNeebsRaw.items():

        atomNeebs[atom] = []
        for neeb in neebs:
            splitLab = neeb.split(',')
            if splitLab[1] == 'asym':
                atomNeebs[atom].append(splitLab[0])
            else:
                newLab = spec2masLab(neeb)
                atomNeebs[atom].append(newLab)
                globAtomEnv[newLab] = globAtomEnv[splitLab[0]]

    return atomNeebs

def getLocalCoordSys():
    '''
    Create local coordinate system around every atom based on local symmetry.
    Return local coordinate systems.
    '''
    atomNeebsRaw = copy.copy(globAtomLabs)   #Get dictionary of atoms and their nearest neighbour labels i.e. {C(2):['C(3)','C(4)','H(2)','H(3)']}
    atomSymsRaw = findSITESYM()         #Get dictionary of atoms and their local symmetry labels and their configuration of different/identical atoms i.e. ('mm2',2+2) could mean a C bonded to 2Hs and 2Cs
    atomLocCoords = {}                      #Initialise dictionary that will be returned and passed to writeLocalCoordSys format {'C(2)':['DUM0','Z','C(3),'Y']}

    i = 1	                                  #Default value for dumI in addDUM(atom1,atom2,dumI). addDUM will find the lowest unused index itself.

    atomSyms = {}
    #Convert ',asym' to xd atom label and create dummys of other atoms
    for atom, sym in atomSymsRaw.items():

        atomSyms[atom] = list(sym[:2])
        for i,item in enumerate(sym[2:]):
            i+=2                                        #correct index to correspond with sym
            if isinstance(item, str):
                splitLab = item.split(',')
                if splitLab[1] != 'asym':
                    atomSyms[atom].append(spec2masLab(item))
                else:
                    atomSyms[atom].append(splitLab[0])
            else:
                atomSyms[atom].append([])
                for item2 in item:
                    splitLab = item2.split(',')
                    if splitLab[1] != 'asym':
                        atomSyms[atom][i].append(spec2masLab(item2))
                    else:
                        atomSyms[atom][i].append(splitLab[0])

    atomNeebs = neebDict2masLabs(atomNeebsRaw)

    for atom,sym in atomSyms.items():       #Go through atomSym dict and do different things depending on nearest neighbour configuration

        symTag = sym[1]

        atom = atom.upper()
        #For all atoms all different configurations up to 4-coordinate are checked for, and a different local coordinate system is setup for each possible configuration.
        #2 identical neighbours mm2
        if symTag[:2] == '1 ':
            if symTag == '1 m':
                atomLab = atomNeebs[atom][0]

                yAtom = ''
                for atom2 in atomNeebs[atomNeebs[atom][0]]:
                    if atom2 != atom:
                        yAtom = atom2
                        break
                atomLocCoords[atom] = [atomNeebs[atom][0], 'Z', yAtom, 'Y']

        elif symTag == '2':
            #Add dummy atom for Z-axis between two connected atoms
            x = addDUM(atomNeebs[atom][0],atomNeebs[atom][1],i)             #x stores dumI value returned by addDUM i.e. DUM3, x = 3
            atomLocCoords[atom]=['DUM'+str(x),'Z',atomNeebs[atom][0],'Y']

        #2 different neighbours m
        elif symTag == '11':
            atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

        #3 identical neighbours m
        #If 3 neighbours and 3 angles all ~120, assign m (could be 3m but how would axes be setup)
        #If 3 neighbours and 1+1+1 angles and sum ~360, assign m
        #If 3 neighbours and 2+1 angles and sum ~360, assign mm2
        #If 3 neighbours and 3 angles but sum < 360, assign m (could be 3m but how would axes be setup?)
        #If 3 neighbours and 1+1+1 angles and sum < 360, assign 1
        #If 3 neighbours and 2+1 angles and sum < 360, assign m
        elif symTag[:2] == '3 ':

            if symTag == '3 3p':
                atomLocCoords[atom] = [atomNeebs[atom][0], 'Z', atomNeebs[atom][1], 'Y']

            elif symTag == '3 111p':                                  #If 3 neighbours and 3 angles all ~120, assign m (could be 3m but how would axes be setup)
                atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']   #If 3 neighbours and 1+1+1 angles and sum ~360, assign m

            elif symTag == '3 21p':         #If 3 neighbours and 2+1 angles and sum ~360, assign mm2
                atomLocCoords[atom] = [sym[2],'Z', sym[3],'Y']     #Random atom chosen for Y axis, should work

            elif symTag == '3 3n':                                              #Think NH3
                x = addDUM(atomNeebs[atom][0],atomNeebs[atom][1],i)            #Make X-axis in centre of arbritary atoms
                atomLocCoords[atom] = ['DUM' + str(x), 'X', atomNeebs[atom][2],'Y'] #X-axis in centre of two atoms then Y-axis on opposite atoms so Z is perpendicular to m

            elif symTag == '3 21n':
                x = addDUM(sym[3][0],sym[3][1],i)            #Make Z-axis in centre of 1 angle.
                atomLocCoords[atom] = ['DUM' + str(x), 'X', sym[2],'Y'] #Opposite atom to dummy atom needed for m perpendicular to Z

            #1+1+1 neighbours: '111 p' gives m, '111 n' gives 1 symmetry
        elif symTag == '111 p':
            atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

            #2+1 different neighbours mm2
                #If 2+1 neighbours and 3 angles and sum ~360, assign mm2
                #If 2+1 neighbours and 1+1+1 angles and sum ~360, assign m
                #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle is opposite odd neighbour, assign mm2
                #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle isn't opposite odd neighbour, assign m
                #If 2+1 neighbours and 3 angles and sum < 360, assign m
                #If 2+1 neighbours and 1+1+1 angles and sum < 360, assign 1
                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle is opposite odd neighbour, assign m
                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle isn't opposite odd neighbour, assign 1

        elif symTag[:3] == '21 ':
            if symTag == '21 3p' or symTag == '21 21pmm2':           #If 2+1 neighbours and 3 angles and sum ~360, assign mm2
                for atomLab in atomNeebs[atom]:                                                   #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle is opposite odd neighbour, assign mm2
                    if atomLab != sym[2]:
                        YAtom = atomLab
                atomLocCoords[atom] = [sym[2], 'Z', YAtom,'Y']     #Y assigned to arbitrary neighbour which should be okay?

            elif symTag == '21 111p' or symTag == '21 21pm':           #If 2+1 neighbours and 1+1+1 angles and sum ~360, assign m
                                                                    #If 2+1 neighbours and 2+1 angles and sum ~360 and odd angle isn't opposite odd neighbour, assign m
                atomLocCoords[atom] = [atomNeebs[atom][0],'X',atomNeebs[atom][1],'Y']

            elif symTag == '21 3n' or symTag == '21 21nm':         #If 2+1 neighbours and 3 angles and sum < 360, assign m
                                                                #If 2+1 neighbours and 2+1 angles and sum < 360, and odd angle is opposite odd neighbour, assign m
                x = addDUM(sym[2][0],sym[2][1],i)            #Make X-axis between 2 atoms
                atomLocCoords[atom] = ['DUM' + str(x), 'X', sym[3],'Y']

        if symTag[:2] == '4 ':
            #4 identical neighbours mm2, 6 identical angles
            if symTag == '4 t':
                #Add dummy atom for Z-axis between two of the atoms
                x = addDUM(atomNeebs[atom][0],atomNeebs[atom][1],i)                 #x stores dumI value returned by addDUM i.e. DUM3, x = 3
                atomLocCoords[atom]=['DUM'+str(x),'Z',atomNeebs[atom][2],'Y']

            elif symTag == '4 sp':
                atomLocCoords[atom] = [atomNeebs[atom][0], 'Z', atomNeebs[atom][1], 'Y']

            elif symTag == '4 t/sp':
                x = addDUM(sym[2][0],sym[2][1],i)
                atomLocCoords[atom] = ['DUM' + str(x), 'Z', sym[3], 'Y']

            elif symTag == '4 et':

                x = addDUM(sym[2][0], sym[2][1], i)
                atomLocCoords[atom] = ['DUM' + str(x), 'Z', sym[3], 'Y']

        #2+1+1
        elif symTag[:4] == '211 ':
            if symTag in ('211 t','211 t oa'):
                atomLocCoords[atom] = [sym[2][0], 'X', sym[2][1], 'Y']

            elif symTag == '211 spmm2':
                atomLocCoords[atom] = [sym[2][0], 'Z', sym[2][1], 'Y']

            elif symTag == '211 spm':
                atomLocCoords[atom] = [atomNeebs[0], 'X', atomNeebs[1], 'Y']

            elif symTag == '211 t/sp':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

        elif symTag[:3] == '31 ':
            if symTag == '31 t':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 sp':
                atomLocCoords[atom] = [sym[2], 'Z', sym[3], 'Y']

            elif symTag == '31 t/sp':

                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 epy':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 et':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

            elif symTag == '31 toa':
                atomLocCoords[atom] = [sym[2], 'X', sym[3], 'Y']

        #2+2
        elif symTag[:3] == '22 ':

            if symTag == '22 t':
                x = addDUM(sym[2][0],sym[2][1],i)
                atomLocCoords[atom] = ['DUM' + str(x), 'Z', sym[3], 'Y']

            elif symTag == '22 sp':
                atomLocCoords[atom] = [sym[2], 'Z', sym[3], 'Y']

            elif symTag == '22 t/spmm2':
                x = addDUM(sym[2][0], sym[2][1], i)
                atomLocCoords[atom] = ['DUM' + str(x), 'Z', sym[3], 'Y']

            elif symTag == '22 etmm2':
                x = addDUM(sym[2][0], sym[2][1], i)
                atomLocCoords[atom] = ['DUM' + str(x), 'Z', sym[3], 'Y']

        #1+1+1+1
        elif symTag == '1111 sp':
            atomLocCoords[atom] = [atomNeebs[atom][0], 'X', atomNeebs[atom][1], 'Y']

        elif sym[0] == '1':                      #Define arbritrary system for '1' symmetry so it is not returned as not being defined
            atomLocCoords[atom] = [atomNeebs[atom][0], 'Z', atomNeebs[atom][1], 'Y']

    return (atomLocCoords, atomNeebs)

def makeLCS():
    '''
    Write local coordinate systems to mas file for all atoms. Return atoms for which local coordinate system could not be found as list.
    '''
    global globAtomEnv

    x = getLocalCoordSys()

    if globAtomEnv:
        try:
            y = writeLocalCoordSys(getCHEMCONLocalCoordSys(x[0], x[1]))
        except KeyError:
            y = writeLocalCoordSys(x[0])
    else:
        y = writeLocalCoordSys(x[0])

    return y

def getCHEMCONLocalCoordSys(atomLocCoordsDict, atomNeebs):
    '''
    Update local coordinate systems to make sure all chemically equivalent atoms have equivalent local coordinate systems.
    Return new local coordinate system dictionary.
    '''
    envs = copy.copy(globAtomEnv)
    chemcon = findMasCHEMCON()
    children = chemcon.keys()
    newLocCoords = {}

    for atom, locCoords in atomLocCoordsDict.items():

        if atom in children:

            pLocCoordSys = atomLocCoordsDict[chemcon[atom]]

            cLocCoordSys = ['', pLocCoordSys[1], '' ,pLocCoordSys[3]]

            neebs = atomNeebs[atom]

            pHash1 = envs[pLocCoordSys[0]]
            pHash2 = envs[pLocCoordSys[2]]

            for neeb in neebs:

                if envs[neeb] == pHash1:

                    cLocCoordSys[0] = neeb

                elif envs[neeb] == pHash2:
                    cLocCoordSys[2] = neeb

            if cLocCoordSys[0] and cLocCoordSys[2]:
                newLocCoords[atom] = cLocCoordSys

        else:
            newLocCoords[atom] = locCoords

    return newLocCoords

def writeLocalCoordSys(atomLocCoordsDict):
    '''
    Write local coordinate systems to atom table.
    Return any atoms for which a local coordinate system was not added.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        atomTab = False
        coordSystems = atomLocCoordsDict

        unaddedLocCoords = []

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.split()
                atom = row[0].upper()
                #If atom is in dictionary it is added with local coordinate system
                if atom in coordSystems.keys():

                    row[1] = coordSystems[atom][0]
                    row[2] = coordSystems[atom][1]
                    row[4] = coordSystems[atom][2]
                    row[5] = coordSystems[atom][3]

                    if len(row)==12:

                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                    else:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)
                    newmas.write(rowStr + '\n')

                else:
                    newmas.write(line)
                    if line[:2].upper() != 'H(':
                        unaddedLocCoords.append(row[0])
            else:
                newmas.write(line)

            if atomTableBegins(line):
                atomTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

    return unaddedLocCoords                 #List of atoms for which local coordinates haven't been added

def check4CustomLCS():

    masAtoms = {}
    customAtoms = {}
    atomTab = False

    with open('xd.mas','r') as mas:

        for line in mas:
            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                line = line.upper()
                row = line.split()
                masAtoms[row[0]] = line

            if atomTableBegins(line):
                atomTab = True

    copyfile('xd.mas','xdcopy.mas')
    multipoleMagician()
    atomTab = False

    with open('xd.mas') as mas:

        for line in mas:
            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                line = line.upper()
                row = line.split()
                if line != masAtoms[row[0]]:
                    customAtoms[row[0]] = masAtoms[row[0]]

            if atomTableBegins(line):
                atomTab = True

    os.remove('xd.mas')
    os.rename('xdcopy.mas', 'xd.mas')

    print(customAtoms)
    print("CUSTOM ATOMS\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    return customAtoms


def getMonoFont(size = 10):
    '''
    Choose font based on OS. Return QFont.
    '''
    font = QFont()
    font.setPointSize(size)
    if sys.platform.startswith('linux'):
        font.setFamily("Bitstream Vera Sans Mono")

    elif sys.platform=='win32':
        font.setFamily("Consolas")

    return font

def getIcon():
    icon = QIcon()
    icon.addPixmap(QPixmap(iconAbsPath), QIcon.Normal, QIcon.Off)
    return icon

'''
#########################################################################
#--------------------GUI-------------------------------------------------
#########################################################################
'''

class XDLSM(QThread):
    '''
    XDLSM
    '''
    startSignal = pyqtSignal()
    finishedSignal = pyqtSignal()
    warningSignal = pyqtSignal()

    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
        '''
        Run XDLSM and add unit cell parameters to xd_lsm.cif afterwards.
        '''

        try:
            self.xdlsmRunning = subprocess.Popen(xdlsmAbsPath, shell = False, cwd = os.getcwd(), stdout=subprocess.PIPE)
            self.startSignal.emit()
            #self.xdlsmRunning.wait()
            print(self.xdlsmRunning.communicate())
            try:
                fixLsmCif()
            except Exception:
                pass
            self.finishedSignal.emit()

        except Exception as e:
            print()
            self.warningSignal.emit()

class XDProg(QThread):
    '''
    Run XD Program in QThread.
    '''
    startSignal = pyqtSignal()
    finishedSignal = pyqtSignal()
    warningSignal = pyqtSignal()

    def __init__(self, prog, args=None):
        QThread.__init__(self)
        self.xdProgName = prog
        self.args = ''
        if args:
            self.args = args

        global xdProgAbsPaths
        try:
            self.xdProg = xdProgAbsPaths[prog]

        except Exception:
            self.xdProg = ''

    def __del__(self):
        self.wait()

    def run(self):
        '''
        Run XD program
        '''
        if self.xdProg:
            try:
                if prog == 'topxd':
                    self.xdProgRunning = subprocess.Popen([self.xdProg, 'topxd.out'], shell = False, cwd = os.getcwd(), stdout = subprocess.PIPE)
                else:
                    self.xdProgRunning = subprocess.Popen([self.xdProg], shell = False, cwd = os.getcwd(), stdout = subprocess.PIPE)

                self.startSignal.emit()
                #self.xdProgRunning.wait()
                print(self.xdProgRunning.communicate())
                self.finishedSignal.emit()

            except Exception:
                pass
        else:
            self.warningSignal.emit()

class XDINI(QThread):
    '''
    XDINI
    '''
    startSignal = pyqtSignal()
    finishedSignal = pyqtSignal()
    warningSignal = pyqtSignal()

    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
        '''
        Run XDINI
        '''
        print('starting XDINI')
        fixBrokenLabels(copy.copy(globAtomLabs))
        print('Fixed broken labels XDINI')
        self.xdiniRunning = subprocess.Popen([xdiniAbsPath, ''.join(compoundID4XDINI.split()), 'shelx'], shell = False, cwd = os.getcwd())
        print('Running XDINI')
        self.startSignal.emit()
        print('Waiting for XDINI')
        self.xdiniRunning.wait()
        #print(self.xdiniRunning.communicate())
        removePhantomAtoms()
        print('removed phantoms')
        try:
            initialiseGlobVars()
            print('initialized globs')
        except Exception as e:
            print(e)
            pass

        self.finishedSignal.emit()

class FindCHEMCONThread(QThread):
    '''
    Find CHEMCON in a separate thread so it doesn't freeze the whole program.
    '''
    finishedSignal = pyqtSignal()
    chemcon = {}
    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
        self.chemcon = findCHEMCON()
        self.finishedSignal.emit()


class aboutBox(QWidget, Ui_aboutBox):
    '''
    About XD Toolkit window.
    '''
    def __init__(self, parent=None):
        super(aboutBox, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.setStyleSheet('''QWidget{background-color: #e8e8e8;
                                                   border-style: outset;
                                                   border-width: 2px;
                                                   border-radius: 10px;
                                                   border-color: transparent;
                                                   padding: 6px;}''')

class sendBug(QDialog, Ui_sendBug):
    '''
    Send bug report window.
    '''
    def __init__(self, parent=None):
        super(sendBug, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.buttonBox.button(QDialogButtonBox.Ok).setText("Send")

class sendSugg(QDialog, Ui_sendSugg):
    '''
    Send suggestion window.
    '''
    def __init__(self, parent=None):
        super(sendSugg, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.buttonBox.button(QDialogButtonBox.Ok).setText("Send")

class prefGui(QDialog, Ui_pref):
    '''
    Preferences window
    '''
    def __init__(self, parent=None):
        '''
        Initialise all buttons in preferences window.
        '''
        super(prefGui, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.chooseMoleCoolPath.clicked.connect(self.chooseMCQt)
        self.chooseMercPathBut.clicked.connect(self.chooseMerc)
        self.chooseTextPathBut.clicked.connect(self.chooseTextEd)
        self.chooseXDPathBut.clicked.connect(self.chooseXD)

    def chooseXD(self):
        '''
        Prompt user to choose folder containing XD executables. Set variable to folder path.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Select XD Folder"))
        xdFolder = ''

        if folder:
            if sys.platform == 'win32':
                if 'xdlsm.exe' not in os.listdir(folder):
                    if 'bin' in os.listdir(folder) and 'xdlsm.exe' in os.listdir(folder + '/bin'):
                        xdFolder = folder + '/bin'
                else:
                    xdFolder = folder

                if not os.environ.get('XD_DATADIR'):
                    self.msg = 'Make sure you have the XD_DATADIR environment variable setup. It should be:\n\nXD_DATADIR={}lib/xd\n\n'.format(xdFolder[:-3])
                    self.envMsg = QMessageBox()
                    self.envMsg.setWindowTitle('XD_DATADIR environment variable')
                    self.envMsg.setText(self.msg)
                    self.envMsg.show()

            elif sys.platform.startswith('linux'):
                if 'xdlsm' not in os.listdir(folder):
                    if 'bin' in os.listdir(folder) and 'xdlsm' in os.listdir(folder + '/bin'):
                        xdFolder = folder + '/bin'
                else:
                    xdFolder = folder

                if not os.environ.get('XD_DATADIR'):
                    self.msg = 'Make sure you have the XD_DATADIR environment variable setup. To do this go to the terminal and enter:\n\nsudo gedit /etc/environment\n\nNow add the following line to the environment file and save it:\n\nXD_DATADIR={}lib/xd\n\nLogout and log back in for the changes to take effect.'.format(xdFolder[:-3])
                    self.envMsg = QMessageBox()
                    self.envMsg.setWindowTitle('XD_DATADIR environment variable')
                    self.envMsg.setText(self.msg)
                    self.envMsg.show()

            self.xdpath = xdFolder
            self.settingsXDLab.setText('Current path: ' + self.xdpath)

    def chooseMCQt(self):
        '''
        Prompt user to choose MoleCoolQt executable file. Set variable to file path.
        '''
        file = QFileDialog.getOpenFileName(None, "Select MoleCoolQt executable", (os.path.expanduser('~')))

        if file[0]:
            self.molecoolpath = file[0]

        try:
            self.settingsMoleCoolLab.setText('Current path: ' + self.molecoolpath)
        except Exception:
            self.settingsMoleCoolLab.setText('Current path: ')

    def chooseMerc(self):
        '''
        Prompt user to choose Mercury executable file. Set variable to file path.
        '''
        file = QFileDialog.getOpenFileName(None, "Select mercury executable", (os.path.expanduser('~')))

        if file[0]:
            self.mercurypath = file[0]
        try:
            self.chooseMercPathLab.setText('Current path: ' + self.mercurypath)
        except Exception:
            self.settingsMercPathLab.setText('Current path: ')

    def chooseTextEd(self):
        '''
        Prompt user to choose text editor executable file. Set variable to file path.
        '''
        file = QFileDialog.getOpenFileName(None, "Select text editor executable", (os.path.expanduser('~')))

        if file:
            self.textedpath = file[0]
        try:
            self.chooseTextPathLab.setText('Current path: ' + self.textedpath)
        except Exception:
            self.chooseTextPathLab.setText('Current path: ')

class mercury(QThread):
    '''
    Mercury
    '''
    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
        '''
        Open shelx.ins in Mercury.
        '''
        try:
            if os.path.exists(os.getcwd() + '/shelx.ins'):
                self.mercuryOpen = subprocess.Popen([globMercAbsPath, 'shelx.ins'], shell = False, cwd = os.getcwd())
            else:
                self.mercuryOpen = subprocess.Popen([globMercAbsPath], shell = False, cwd = os.getcwd())
        except Exception:
            pass

class molecool(QThread):
    '''
    MoleCoolQt
    '''
    def __init__(self):
        QThread.__init__(self)

    def __del__(self):
        self.wait()

    def run(self):
        '''
        Open xd.res/xd.inp in MoleCoolQt.
        '''
        try:
            if os.path.exists(os.getcwd() + '/xd.res'):
                self.molecoolOpen = subprocess.Popen([molecoolQtAbsPath, 'xd.res'], shell = False, cwd = os.getcwd())
            elif os.path.exists(os.getcwd() + '/xd.inp'):
                self.molecoolOpen = subprocess.Popen([molecoolQtAbsPath, 'xd.inp'], shell = False, cwd = os.getcwd())
            else:
                self.molecoolOpen = subprocess.Popen([molecoolQtAbsPath], shell = False, cwd = os.getcwd())
        except Exception:
            pass

class resmap(QWidget, Ui_resmap):
    '''
    Residual density map window.
    '''
    def __init__(self, grdFile, parent=None):
        super(resmap, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.grdFile = grdFile
        self.setup()

    def setup(self):
        '''
        Create residual density map and display it with save button.
        '''
        i=0
        self.tempFileName = 'temp_quickplot{}.png'.format(i)
        while os.path.isfile(self.tempFileName):
            i+=1
            self.tempFileName = 'temp_quickplot{}.png'.format(i)
        fig = makeResMap(self.grdFile, self.tempFileName)
        fig.set_dpi(90)
        canvas = FigureCanvas(fig)
        canvas.setParent(self)
        saveBut = QPushButton('Save hi-res PNG file')
        saveBut.clicked.connect(self.savePng)
        saveBut.setFixedWidth(200)
        self.saveLab = QLabel()
        self.resmapLayout.addWidget(canvas, 0)
        self.resmapLayout.addWidget(saveBut,1)
        self.resmapLayout.addWidget(self.saveLab,2)
        self.setLayout(self.resmapLayout)
        # prevent the canvas to shrink beyond a point
        # original size looks like a good minimum size
        canvas.setMinimumSize(canvas.size())

    def savePng(self):
        '''
        Save residual density map as png file.
        '''
        filename = QFileDialog.getSaveFileName(self,"PNG of residual density map",os.getcwd(), "PNG Files (*.png)")[0]
        if filename:
            if filename[:-4] != '.png':
                filename += '.png'
            copyfile(self.tempFileName, filename)
            self.saveLab.setText('Residual map saved to <i>"{}"</i>'.format(filename))

    def closeEvent(self, event):
        for file in os.listdir(os.getcwd()):
            if file.startswith('temp_quickplot'):
                os.remove(file)
        event.accept()


class NPP(QWidget, Ui_resmap):
    '''
    Show normal probability plot in new window.
    '''
    def __init__(self, parent=None):

        super(NPP, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.setWindowTitle('Normal probability plot')
        self.setup()

    def setup(self):
        '''
        Setup normal probability plot window.
        '''
        self.values = grd2values()
        self.maxValue = max(self.values)
        self.minValue = min(self.values)
        print(self.maxValue)
        fig = plt.figure(facecolor = '#dddddd', figsize=(7,5), dpi=80)
        ax = fig.add_subplot(1,1,1)
        probplot(self.values, plot = plt)
        plt.title('Normal probability plot')
        ax.get_lines()[0].set_marker('.')
        ax.get_lines()[0].set_markersize(0.8)
        canvas = FigureCanvas(fig)
        saveBut = QPushButton('Save PNG file')
        saveBut.clicked.connect(self.savePng)
        saveBut.setFixedWidth(150)
        self.infoLab = QLabel()
        self.infoLab.setText('Minimum residual = {0: .3f}\nMaximum residual = {1: .3f}'.format(self.minValue, self.maxValue))
        self.saveLab = QLabel()
        self.resmapLayout.addWidget(canvas, 0)
        self.resmapLayout.addWidget(saveBut,1)
        self.resmapLayout.addWidget(self.saveLab,2)
        self.resmapLayout.addWidget(self.infoLab,3)
        self.setLayout(self.resmapLayout)

        # prevent the canvas to shrink beyond a point
        # original size looks like a good minimum size
        canvas.setMinimumSize(canvas.size())

    def savePng(self):
        '''
        Save residual density map as png file.
        '''
        filename = QFileDialog.getSaveFileName(self,"PNG of normal probability plot",os.getcwd(), "PNG Files (*.png)")
        if filename[0]:
            fig = plt.figure(facecolor = '#dddddd', figsize=(11,9), dpi=80)
            ax = fig.add_subplot(1,1,1)
            probplot(self.values, plot = plt)
            plt.title('Normal probability plot')
            ax.get_lines()[0].set_marker('.')
            ax.get_lines()[0].set_markersize(0.8)
            plt.savefig(filename[0])
            self.saveLab.setText('Normal probability plot saved to <i>"{}"</i>'.format(filename[0]))

class checkNeebs(QWidget, Ui_checkneebs):
    '''
    Window where user can check that automatically generated neighbours are correct.
    '''
    def __init__(self, parent=None):
        super(checkNeebs, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        atomList = getAtomList()
        if len(atomList) < 30:
            self.tableLab.setText(self.writeNeebTable(atomList))

        else:
            midPt = int(len(atomList)/2)
            self.tableLab.setText(self.writeNeebTable(atomList[:midPt]))
            self.tableLab2.setText(self.writeNeebTable(atomList[midPt:]))

        self.scrollArea.setWidget(self.scrollAreaWidgetContents_2)
        self.gridLayout.addWidget(self.scrollArea, 0, 0, 1, 1)

        self.retranslateUi(self)
        QMetaObject.connectSlotsByName(self)

    def writeNeebTable(self, atomList):
        tableStr = '\n{0:10}{1}\n'.format('ATOM', 'NEIGHBOURS')
        for i, atom in enumerate(atomList):
            neebs = globAtomLabs[atom]
            tableStr += '\n<b>{1:10}</b><i>{2}</i>'.format(i, atom, sevenSpacedListjoin([spec2norm(item) for item in neebs], ' '))
        tableStr = '<pre>' + tableStr + '</pre>'

        return tableStr

class wizardRunning(QDialog, Ui_wizard):
    '''
    XD Wizard execution window.
    '''
    finishedSignal = pyqtSignal()

    def __init__(self, parent=None):
        '''
        Intialise xdini, xdlsm and various wizard variables.
        '''
        super(wizardRunning, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.xdini = XDINI()
        self.xdlsm = XDLSM()
        self.xdlsm.finishedSignal.connect(self.xdWizRef)
        self.i = 0
        self.xxxDanger = 0          #Stop runaway recursive function calling if error can't be fixed that should be fixed.
        self.errorFixed = False     #Stop error check if error has been fixed and lsm needs to run again
        self.startingTime = 0
        self.finishingTime = 0
        self.refList = ['0 - XDINI', '1 - Scale factors', '2 - High angle non-H positions and ADPs',
                        '3 - Low angle H positions and isotropic ADPs', '4 - Kappa and monopoles', '5 - Multipoles',
                        '6 - Multipoles, and non-H positions and ADPs', '7 - Lower symmetry', '8 - Final refinement']

    def closeEvent(self, event):
        '''
        Stop xdlsm if window closed and emit finished signal.
        '''
        self.xdlsm.xdlsmRunning.terminate()
        self.finishedSignal.emit()
        event.accept()

    def xdWiz(self, backupFolder):
        '''
        Start XD Wizard.
        '''
        if self.collectDat:
            self.tzero = time.time()
        print('wiz called')
        print(backupFolder)
        self.i = 0
        resStr = '<br>{0:47}{1:23}{2:14}{3:16}{4}<br>'.format('Refinement', 'RF<sup>2</sup>', 'Convergence', 'Average DMSDA', 'Max DMSDA')
        resStr = '<pre>' + resStr + '</pre>'
        self.wizResLab.setText(resStr)
        print(resStr)
        self.folder = backupFolder
        print(self.folder)
        try:
            os.makedirs('Backup/' + self.folder)
        except Exception as e:
            printExc('xdWiz',e)
            print("Invalid backup folder name. Can't continue.")
        if os.path.isfile('shelx.ins') and os.path.isfile('shelx.hkl'):
            self.wizStatusLab.setText('Initializing compound...')
            self.xdini.finishedSignal.connect(self.xdWizRef)
            self.xdini.start()

        else:
            self.xdWizRef()


    def xdWizRef(self):
        '''
        Backup, print resutls and run next refinement.
        Call itself repeatedly until all refinements are done.
        Finally emit finished signal.
        '''
        print('xdWizRef called')

#        if self.collectDat:
#            global timeFileAbsPath
#            if self.i > 0 and self.i < len(self.refList):
#                try:
#                    with open(timeFileAbsPath, 'a') as timeFile:
#
#                        self.finishingTime = time.time()
#
#                        runningTime = self.finishingTime - self.startingTime
#                        numAtoms = getNumAtoms()
#                        ref = self.refList[self.i]
#                        refSplit = ref.split()
#                        refName = ''.join(refSplit[2:])
#
#                        try:
#                            refName = refName[:13]
#                        except Exception:
#                            pass
#
#                        timeStr = '{0:21}{1:10}{2:<13.2f}{3}\n'.format(refName, str(numAtoms), runningTime, sys.platform)
#                        timeFile.write(timeStr)
#
#                except Exception:
#                    pass

        global wizUniSnlMax
        self.i += 1

        if self.i > 1 and not self.errorFixed:
            error = check4errors()

            if not error[0]:
                self.xxxDanger += 1
                print('error found')
                print(error[1])
                if error[1].startswith('  * * * Parameter T and/or NCST should be increased'):
                    addNCST()
                    self.i = self.i-2
                    if self.xxxDanger < 4:
                        self.errorFixed = True
                        self.xdWizRef()
                        return
                    else:
                        self.wizStatusLab.setText('The following error occurred:<br><br>' + error[1] + '<br><br>XD wizard cannot proceed.')
                        return

                elif error[1].startswith('Noble gas configuration not recognized for element Cu'):
                    fixCuNobleGasError()
                    self.i = self.i-2
                    if self.xxxDanger < 4:
                        self.errorFixed = True
                        self.xdWizRef()
                        return
                    else:
                        self.wizStatusLab.setText('The following error occurred:<br><br>' + error[1] + '<br><br>XD wizard cannot proceed.')
                        return

                else:
                    self.wizStatusLab.setText('The following error occurred:<br><br>' + error[1] + '<br><br>XD wizard cannot proceed.')
                    return
            else:
                print('no errors')
                if self.i > 1:
                    backup('{}/{}'.format(self.folder, self.refList[self.i-1]))
                    print('backed up')
                if self.i > 1:
                    c = ''
                    resStr = str(self.wizResLab.text())
                    rf2 = getRF2()
                    conv = getConvergence('xd_lsm.out')
                    if conv:
                        c = 'Yes'
                    else:
                        c = 'No'
                    try:
                        dmsda = getDMSDA('xd_lsm.out')
                    except Exception:
                        dmsda = ['','']

                    ref = self.refList[self.i-1]
                    refSplit = ref.split()
                    refName = ' '.join(refSplit[2:])
                    refNum = refSplit[0]
                    newRes = '<br>{0:>2} - {1:40}{2:6.2f}{3:8}{4:14}{5:16}{6}'.format(refNum, refName, rf2, ' %', c, dmsda[0], dmsda[1])
                    resStr = resStr[:-6] + newRes + '</pre>'
                    self.wizResLab.setText(resStr)
                    self.finalRes = resStr
                    if self.i < len(self.refList):
                        self.wizStatusLab.setText('Setting up xd.mas for ' + refName[:1].lower() + refName[1:])
                        res2inp()
                        print('res2inp')

        try:
            if self.i == len(self.refList):

                self.wizStatusLab.setText('Making normal probability plot...')
                self.wizStatusLab.repaint()
                QApplication.processEvents()
                krauseParam = getKrauseParam()
                if krauseParam:
                    self.finalRes += '\nKrause parameter = {0:.2f}'.format(krauseParam)
                self.xdlsm.finishedSignal.disconnect(self.xdWizRef)
                self.finishedSignal.emit()
                return

            ref = self.refList[self.i]
            refSplit = ref.split()
            refNum = refSplit[0]
            refName = ' '.join(refSplit[2:])

            print(refName)

            noUniSnl = False

            if refName.upper().startswith('SCALE'):
                scaleFacRef()
                print('scalefacref')

            elif refName.upper().startswith('HIGH ANGLE'):
                sinthlMin = wizHighSnlMin
                sinthlMax = wizHighSnlMax
                highAngleRef(sinthlMin,sinthlMax)
                wizAddResetBond()
                noUniSnl = True

            elif refName.upper().startswith('LOW ANGLE'):
                print('LOW ANGLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                sinthlMin = wizLowSnlMin
                sinthlMax = wizLowSnlMax
                lowAngleRef(sinthlMin,sinthlMax)
                noUniSnl = True

            elif refName.upper().startswith('KAPPA'):
                wizAddCHEMCON()
                kapMonRef()

            elif refName.upper() == ('MULTIPOLES'):
#                multipoleMagician()
                multipoleKeyTable()

            elif refName.upper() == ('DIPOLES'):
                #multipoleMagician()
                wizAddMultipoles()
                seqMultRef(1)

            elif refName.upper().startswith('QUADRUPOLES'):
#                multipoleMagician()
                wizAddMultipoles()
                seqMultRef(2)

            elif refName.upper().startswith('OCTUPOLES'):
#                multipoleMagician()
                wizAddMultipoles()
                seqMultRef(3)

            elif refName.upper().startswith('HEXADECAPOLES'):
 #               multipoleMagician()
                wizAddMultipoles()
                seqMultRef(4)

            elif refName.upper().startswith('MULTIPOLES,'):
                #posADPMultRef()
                nonHPosADPKey()

            elif refName.upper().startswith('LOWER'):
                lowerSym()
                print(self.customLCS)
                wizAddCustomLCS(self.customLCS)
                print('sym lowered')

            elif refName.upper().startswith('FINAL'):
                checkMultRes()
                removeCHEMCON()

            if wizUniSnlMax and not noUniSnl:
                addSnlCutoff(snlmin = 0.0, snlmax = wizUniSnlMax)

            if self.collectDat:
                self.startingTime = time.time()
            self.xdlsm.start()
            self.errorFixed = False

            print('xdlsm started')
            self.wizStatusLab.setText('{0} {1}/{2} - {3}'.format('Running XDLSM', refNum, len(self.refList)-1, refName))

        except Exception as e:
            print(e)
            self.wizStatusLab.setText("Couldn't setup xd.mas for {}".format(refName))

class AutoTOPXDThread(QThread):
    finishedSignal = pyqtSignal()
    def __init__(self, atom):
        QThread.__init__(self)
        self.atom = atom

    def __del__(self):
        self.wait()

    def run(self):
        if sys.platform.startswith('linux'):
                with open('topxd_' + self.atom + '.out','w') as outputFile:
                    self.topxd = subprocess.Popen([topxdAbsPath], shell=False, cwd = os.getcwd(), stdout=outputFile)
                    self.topxd.communicate()

        elif sys.platform=='win32':
            self.topxd = subprocess.Popen([topxdAbsPath, 'topxd_' + self.atom + '.out'], shell=False, cwd=os.getcwd())
            self.topxd.wait()

        self.finishedSignal.emit()

class AutoTOPXD(QWidget, Ui_autoTOPXD):
    '''
    Modal window that shows while autoTOPXD is running.
    '''
    terminatedSignal = pyqtSignal()
    finishedSignal = pyqtSignal()

    def __init__(self, atomList, phi, theta, parent=None):
        super(AutoTOPXD, self).__init__(parent)
        self.setupUi(self)
        self.setFont(getMonoFont())
        self.setWindowIcon(getIcon())
        self.i=0
        self.phi = phi
        self.theta = theta
        self.cancelBut.clicked.connect(lambda: self.terminatedSignal.emit())
        self.atomList = atomList
        self.atomsDone = 0
        self.runTimes = []
        self.totalNumAtoms = len(atomList)
        self.progressBar.setFormat('%p / {0}'.format(self.totalNumAtoms))
        self.progressBar.setMaximum(self.totalNumAtoms)

    def run(self):
        self.atom = self.atomList[self.i]
        setupmasTOPXD(self.atom, self.phi, self.theta)
        self.statusLab.setText('TOPXD running on {}...'.format(self.atom))
        self.startTime = time.time()
        self.topxd = AutoTOPXDThread(self.atom)
        self.topxd.finishedSignal.connect(self.topxdFinished)
        self.topxd.start()

    def topxdFinished(self):
        self.endTime = time.time()
        self.topxd.terminate()
        self.atomsDone += 1
        self.progressBar.setValue(self.atomsDone)
        self.i+=1
        self.runTimes.append(self.endTime-self.startTime)
        self.avgRunTime = sum(self.runTimes)/len(self.runTimes)
        self.estRemTime = len(self.atomList[self.i:])*self.avgRunTime
        self.now = datetime.datetime.now()
        self.estFinTime = self.now + datetime.timedelta(seconds = self.estRemTime)

        self.statusStr = 'Estimated finishing time: {}'.format(self.estFinTime.strftime('%Y-%m-%d %H:%M:%S'))

        os.rename('topxd_' + self.atom + '.out', 'topxd/topxd_' + self.atom + '.out')

        if self.i<len(self.atomList):
            self.atom = self.atomList[self.i]
            setupmasTOPXD(self.atom, self.phi, self.theta)
            self.statusStr += '\n\nTOPXD running on {}...'.format(self.atom)
            self.statusLab.setText(self.statusStr)
            self.topxd = AutoTOPXDThread(self.atom)
            self.topxd.finishedSignal.connect(self.topxdFinished)
            self.topxd.start()
            self.startTime = time.time()
        else:
            self.finishedSignal.emit()

class XDToolGui(QMainWindow, Ui_MainWindow):
    '''
    Main window.
    '''
    def __init__(self, parent=None):
        '''
        Initialise everything in main window.
        '''
        super(XDToolGui, self).__init__(parent)
        self.setupUi(self)

        self.versionNum = '0.8.0'

        self.helpTexts = helpTexts
        #Set font depending on OS
        self.setFont(getProgramFont())

        self.cwdStatusLab = QLabel()
        self.settings = QSettings('prefs')
        self.initialiseSettings()	#Load user preferences

        initialiseGlobVars()
        self.checkNeebsButs = (self.checkNeebsBut, self.wizCheckNeebsBut)
        self.enableCheckNeebsButs()
        self.changeUserIns()

        #Check for XD files and if they are not there prompt user to find directory.
        if not self.settings.value('xdpath'):
            self.findXD()

        self.treeModel = QFileSystemModel()
        self.treeModel.setRootPath(QDir.rootPath())
        self.folderTree.setModel(self.treeModel)
        self.folderTree.setRootIndex(self.treeModel.index(QDir.currentPath()))
        self.folderTree.setColumnHidden(1, True)
        self.folderTree.setColumnHidden(2, True)
        self.folderTree.setColumnHidden(3, True)
        self.folderTree.doubleClicked.connect(self.folderTreeDoubleClicked)
        #StyleSheet should be moved to QtDesigner
        self.folderTree.setStyleSheet('''QTreeView{background-color: #ffffff;
                                                   alternate-background-color: #eeeeff;
                                                   border-style: outset;
                                                   border-width: 0px;
                                                   border-radius: 10px;
                                                   border-color: transparent;
                                                   padding: 6px;}''')
        
        self.addedLocCoords = {}
        self.ins = ''
        self.forbiddenChars = ['*', '?', '"', '/', '\\', '<', '>', ':', '|']
        self.invalidFolderNameStr = ('Invalid folder name given.\n'
                                     'Remove the following characters: * ? " / \ < > : |')
        self.permErrorMsg = 'Close MoleCoolQT and try again.'
        self.lstMissingErrorMsg = 'Select add shelx.ins to project folder and try again.'
        self.atomNoExistErrorMsg = 'Atom not in structure.'
        self.backupConfirmStr = 'Current files backed up to: '
        self.labList = [self.wizBackupInput, self.xdWizINILab, self.wizTestStatusLab, self.xdWizardStatusLab, self.pkgXDLab, self.setupFOURStatusLab, self.getDpopsStatusLab, self.XDINILab, self.manRefSetupLab, self.manRefBackupInput, self.manRefSnlMin, self.manRefSnlMax, self.manRefResLab, self.editRBLab, self.editRBLab, self.autoResetBondStatusLab, self.resetBondStatusLab, self.editRBLab, self.CHEMCONStatusLab, self.inputElementCHEMCON, self.inputAtomCHEMCON, self.addDUMLab, self.alcsStatusLab, self.multKeyStatusLab, self.showResLab, self.resNPPLab, self.loadBackupLab]
        self.tabWidget.currentChanged.connect(self.updateHelpText)
        self.tabWidget.setCurrentIndex(0)
        toolboxes = [self.resToolbox, self.toolsToolbox]
        self.toolsToolbox.currentChanged.connect(self.updateHelpText)
        for item in toolboxes:
            item.setCurrentIndex(0)
        self.lsmResLine.setVisible(False)
        #Display current working directory on startup.
        self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())
        self.statusbar.addWidget(self.cwdStatusLab)
        self.toolbarRefreshCwd.triggered.connect(self.refreshFolder)
        self.toolbarSetFolder.triggered.connect(self.setFolder)
        self.toolbarOpenmas.triggered.connect(self.openMas)
        self.menuSendBug.triggered.connect(self.bugRepDialog)
        self.menuSendSugg.triggered.connect(self.sendSuggDialog)
        self.menuClearCache.triggered.connect(lambda: clearCache())
        #Run addDUM() with user input when button is pressed.
        self.dumBut.clicked.connect(self.addDUMPress)
        self.DUMnum.returnPressed.connect(self.addDUMPress)
        self.alcsBut.clicked.connect(self.alcs)
        self.alcsLocSymInput.returnPressed.connect(self.alcs)
        #Run multipoleKeyTable() when button is pressed.
        self.multKeyBut.clicked.connect(self.multKeyPress)
        #Select working directory when button is pressed.
        self.dirBut.triggered.connect(self.setFolder)
        #Run resetBond() with user input when button is pressed.
        self.resetBondBut.clicked.connect(self.resetBondInput)
        self.resetLengthInput.returnPressed.connect(self.resetBondInput)
        self.resetHInput.returnPressed.connect(self.resetBondInput)
        self.armRBBut.clicked.connect(self.armRBPress)
        self.disarmRBBut.clicked.connect(self.disarmRBPress)
        #Disable RESET BOND textbox when 'All' is checked.
        self.allResetBox.stateChanged.connect(self.disableResetText)
        self.inputElementCHEMCON.returnPressed.connect(self.runCHEMCON)
        self.inputAtomCHEMCON.returnPressed.connect(self.runCHEMCON)
        self.inputElementCHEMCON.textChanged.connect(lambda: self.chooseElementCHEMCON.setChecked(True))
        self.inputAtomCHEMCON.textChanged.connect(lambda: self.chooseAtomCHEMCON.setChecked(True))
        self.autoCHEMCONBut.clicked.connect(self.autoCHEMCONPress)
        #Run writeCHEMCON() with argument based on user input
        self.addCHEMCONBut.clicked.connect(self.runCHEMCON)
        #Run res2inp when button is pressed
        self.toolbarRes2Inp.triggered.connect(self.res2inpPress)
        #Backup buttons
        self.customBackupInput.returnPressed.connect(lambda: self.genBackupPress(str(self.customBackupInput.text()), self.customBackupLab))
        self.customBackupBut.clicked.connect(lambda: self.genBackupPress(str(self.customBackupInput.text()), self.customBackupLab))
        #Load backup Button
        self.loadBackupBut.clicked.connect(self.loadBackupPress)
        #Get results from backup
        self.resBackupBut.clicked.connect(self.showResBackup)
        self.resBackupSumBut.clicked.connect(self.resBackupSum)
        #Open bond lengths paper
        self.brammerBondLengthsBut.clicked.connect(lambda: webbrowser.open('http://pubs.rsc.org/en/content/articlepdf/1987/p2/p298700000s1'))
        #Run autoResetBond
        self.autoResetBondBut.clicked.connect(self.autoResetBondPress)
        #Add XDFOUR instructions when button pressed
        self.setupFOURBut.clicked.connect(self.addFOURIns)
        self.setupFOURInputText.returnPressed.connect(self.addFOURIns)
        self.setupFOURInputText.textChanged.connect(lambda: self.setupFOURInputBox.setChecked(True))
        #Laplacian tools
        self.lapGrdBut.clicked.connect(lambda: self.showLapmap(self.grdLapmapLab))
        self.taLapmapBut.clicked.connect(self.make3AtomLapmap)

        #Tabulate properties
        self.tabPropBut.clicked.connect(self.tabPropRun)

        #Auto TOPXD
        self.autoTopxdBut.clicked.connect(self.autoTopxd)
        self.autoTopxdInput.textChanged.connect(self.autoTopxdInputHandler)

        self.delResetBondBut.clicked.connect(self.delResetBondPress)
        self.getResBut.clicked.connect(self.getResPress)
        #Dorb pops
        self.getDpopsBut.clicked.connect(self.getDpops)
        self.resNPPBut.clicked.connect(self.makeNPP)
        #Open project folder
        self.menuOpenCwd.triggered.connect(self.openCwd)
        #Open settings dialog
        self.menuLoadIAM.triggered.connect(self.loadIAM)
        self.menuPref.triggered.connect(self.openPrefs)
        self.toolbarSettings.triggered.connect(self.openPrefs)
        self.menuManual.triggered.connect(self.openManual)
        self.menuAbout.triggered.connect(self.openAbout)
        #Run xdlsm.exe
        self.xdProgButs = [self.resNPPBut, self.pkgTOPXDBut, self.pkgXDFFTBut, self.pkgXDFOURBut, self.pkgXDGEOMBut,
                          self.pkgXDPDFBut, self.pkgXDLSMBut, self.pkgXDPROPBut, self.xdWizINIBut, self.xdWizardBut, self.runXDLSMBut, self.runXDINIBut, self.setupFOURBut, self.getDpopsBut]
        self.xdlsm = XDLSM()
        self.runXDLSMBut.clicked.connect(lambda: self.xdlsm.start())
        self.xdlsm.startSignal.connect(self.startXDLSM)
        self.xdlsm.finishedSignal.connect(self.finishedXDLSM)
        self.xdlsm.warningSignal.connect(self.findXD)
        self.xdfour = XDProg('xdfour')
        self.xdfour.startSignal.connect(self.startXDFOUR)
        self.xdfour.finishedSignal.connect(self.finishedXDFOUR)

        self.xdfft = XDProg('xdfft')
        self.xdfft.startSignal.connect(self.startXDFFT)
        self.xdfft.finishedSignal.connect(self.finishedXDFFT)
        self.xdprop = XDProg('xdprop')
        self.xdprop.startSignal.connect(self.startXDPROP)
        self.xdprop.finishedSignal.connect(self.finishedXDPROP)
        self.xdini = XDINI()
        self.xdini.startSignal.connect(self.startXDINI)
        self.xdini.finishedSignal.connect(self.finishedXDINI)
        self.xdini.warningSignal.connect(self.badLabelWarning)
        self.runXDINIBut.clicked.connect(self.runXDINI)
        self.manRefIDInput.returnPressed.connect(self.runXDINI)
        self.molecoolqt = molecool()
        self.toolbarMolecool.triggered.connect(self.openMolecool)
        self.mercury = mercury()
        self.toolbarMercury.triggered.connect(self.openMerc)
        self.toolbarRes2Inp.triggered.connect(self.res2inpPress)

        self.wiz2ndStageObjects = [self.wizAdvOptBut, self.wizBackupInput, self.xdWizardBut, self.wizReuseMasBut, self.wizUniSnlMax, self.wizHighSnlMin, self.wizHighSnlMax, self.wizLowSnlMin, self.wizLowSnlMax]
        self.wizSnlInput = [self.wizHighSnlMin, self.wizHighSnlMax, self.wizLowSnlMin, self.wizLowSnlMax]
        for item in self.wizSnlInput:
            item.editingFinished.connect(self.wizTest)

        self.xdWizINIBut.clicked.connect(self.xdWizINI)
        self.xdWizardBut.clicked.connect(self.xdWizRun)
        self.xdWizCmpID.returnPressed.connect(self.xdWizINI)
        self.wizReuseMasBut.clicked.connect(self.wizReuseMas)
        self.changeWizBackup()

        #Manual refinement
        self.manRefDropMenu.currentIndexChanged.connect(self.refChosen)
        self.manRefSetupBut.clicked.connect(self.setupRef)
        self.manRefResBut.clicked.connect(self.showResRef)
        self.manRefBackupInput.returnPressed.connect(lambda: self.genBackupPress(str(self.manRefBackupInput.text()), self.manRefBackupLab))
        self.manRefBackupBut.clicked.connect(lambda: self.genBackupPress(str(self.manRefBackupInput.text()), self.manRefBackupLab))
        self.manRefSnlMax.returnPressed.connect(self.setupRef)

        #XD package
        self.XDLSMLabs = [self.manRefSetupLab, self.pkgXDLab]
        self.XDLSMButs = [self.runXDLSMBut, self.pkgXDLSMBut]
        self.pkgXDButs = [self.pkgTOPXDBut, self.pkgXDFFTBut, self.pkgXDFOURBut, self.pkgXDGEOMBut,
                          self.pkgXDLSMBut, self.pkgXDPROPBut]
        self.pkgXDLSMBut.clicked.connect(lambda: self.xdlsm.start())

        self.XDFOURButs = [self.setupFOURBut, self.pkgXDFOURBut, self.resNPPBut]
        self.XDFOURLabs = [self.setupFOURStatusLab, self.pkgXDLab, self.resNPPLab]
        self.pkgXDFOURBut.clicked.connect(lambda: self.xdfour.start())

        self.XDFFTButs = [self.setupFOURBut, self.pkgXDFFTBut]
        self.XDFFTLabs = [self.setupFOURStatusLab, self.pkgXDLab]
        self.pkgXDFFTBut.clicked.connect(lambda: self.xdfft.start())

        self.XDPROPButs = [self.getDpopsBut, self.pkgXDPROPBut]
        self.XDPROPLabs = [self.getDpopsStatusLab, self.pkgXDLab]
        self.pkgXDPROPBut.clicked.connect(lambda: self.xdprop.start())

        self.xdgeom = XDProg('xdgeom')
        self.xdgeom.startSignal.connect(self.startXDGEOM)
        self.xdgeom.finishedSignal.connect(self.finishedXDGEOM)
        self.XDGEOMLabs = [self.pkgXDLab]
        self.XDGEOMButs = [self.pkgXDGEOMBut]
        self.pkgXDGEOMBut.clicked.connect(lambda: self.xdgeom.start())

        self.xdpdf = XDProg('xdpdf')
        self.xdpdf.startSignal.connect(self.startXDPDF)
        self.xdpdf.finishedSignal.connect(self.finishedXDPDF)
        self.XDPDFLabs = [self.pkgXDLab]
        self.XDPDFButs = [self.pkgXDPDFBut]
        self.pkgXDPDFBut.clicked.connect(lambda: self.xdpdf.start())

        self.topxd = XDProg('topxd')
        self.topxd.startSignal.connect(self.startTOPXD)
        self.topxd.finishedSignal.connect(self.finishedTOPXD)
        self.TOPXDLabs = [self.pkgXDLab]
        self.TOPXDButs = [self.pkgTOPXDBut]
        self.pkgTOPXDBut.clicked.connect(lambda: self.topxd.start())

        self.xdProgsRunning = [self.xdgeom, self.xdfour, self.xdfft,
                               self.xdprop, self.xdpdf, self.topxd]

        for but in self.checkNeebsButs:
            but.clicked.connect(self.checkNeebsPress)

        #If XD program is run and XD path hasn't been set, program prompts user to set it.
        for prog in self.xdProgsRunning:
            prog.warningSignal.connect(self.findXD)

        self.model = QStandardItemModel()
        self.wizSeqMultBox = QStandardItem("Add multipoles one by one.")
        self.wizSeqMultBox.setCheckState(Qt.Checked)
        self.wizSeqMultBox.setCheckable(True)
        self.model.appendRow(self.wizSeqMultBox)

        self.wizLowerSymBox = QStandardItem("Try lowering symmetry at end of refinement.")
        self.wizLowerSymBox.setCheckState(Qt.Checked)
        self.wizLowerSymBox.setCheckable(True)
        self.model.appendRow(self.wizLowerSymBox)
        self.wizAdvList.setModel(self.model)

        self.wizAdvList.setStyleSheet('''QListView{background-color: #dddddd;
                                                   border-style: outset;
                                                   border-width: 2px;
                                                   border-radius: 10px;
                                                   border-color: transparent;
                                                   padding: 6px;}
                                        QListView:item{background-color: #dddddd;
                                                   border-style: outset;
                                                   border-width: 2px;
                                                   border-radius: 10px;
                                                   border-color: transparent;
                                                   padding: 3px;}
                                                 ''')

        self.wizAdvOptBut.clicked.connect(self.wizAdvOptToggle)
        self.wizAdvOptOpen = False

        self.showMaximized()

        self.menubar.setNativeMenuBar(True)

    def folderTreeDoubleClicked(self, index):
        '''
        Handle item in folder tree being double clicked.
        '''
        filePath = self.treeModel.filePath(index)
       #fileName = self.treeModel.fileName(index)

        global textEdAbsPath

        #Don't open if a folder has been clicked/
        if os.path.isdir(filePath):
            pass

        #Opens everything in chosen text editor on windows, on OSX and Linux opens in default opener.
        #Would be nice to tailor it more to the files that windows doesn't open as text as default.
        else:
            try:
                if sys.platform == 'win32':
                    if not filePath[-6:]=='xd.mas':
                        os.startfile(filePath)
                    else:
                        subprocess.call([textEdAbsPath, filePath])
                else:
                    opener ="open" if sys.platform == "darwin" else "xdg-open"
                    subprocess.call([opener, filePath])

            except Exception as e:
                printExc('folderTreeDoubleClicked',e)
                pass


    def wizAdvOptToggle(self):
        '''
        Handler to expand/collapse advanced wizard options.
        '''
        if not self.wizAdvOptOpen:
            self.wizAdvList.setMaximumSize(2000, 70)
            self.wizAdvOptOpen = True
        else:
            self.wizAdvList.setMaximumSize(0,0)
            self.wizAdvOptOpen = False
#########################################################################
#--------------------EMAIL-----------------------------------------------
#########################################################################

    def bugRepDialog(self):
        '''
        Handle bug report dialog.
        '''
        self.sendDialog = sendBug()
        self.sendDialog.show()
        self.sendDialog.accepted.connect(self.sendBugReport)

    def sendBugReport(self):
        '''
        Send bug report with input from dialog.
        '''
        self.sendDialog.accepted.disconnect(self.sendBugReport)
        try:
            copyfile('shelx.ins','shelx.buckfast')
            msg = str(self.sendDialog.bugReportInput.toPlainText())
            emailStr = str(self.sendDialog.emailAdInput.text())

            if self.sendDialog.incFilesBox.isChecked():
                sendEmail(body = msg, email = emailStr, subject = 'Bug Report: XD Toolkit', attachments = ['shelx.buckfast','shelx.hkl'])
            else:
                sendEmail(body = msg, email = emailStr, subject = 'Bug Report: XD Toolkit')
            if os.path.isfile('shelx.buckfast'):
                os.remove('shelx.buckfast')
            self.sendDialog.accept()
            self.thanks = QMessageBox()
            self.thanks.setText('Thanks for your feedback!')
            self.thanks.setWindowTitle('Submit Bug Report')
            self.thanks.show()
        except Exception:
            self.msg = QMessageBox()
            self.msg.setText('An error occurred, report not sent. You can instead email "mcrav@chem.au.dk" with your message.')
            self.msg.setWindowTitle('Submit Bug Report Failed')
            self.msg.show()

    def sendSuggDialog(self):
        '''
        Handle send suggestion dialog.
        '''
        self.sendDialog = sendSugg()
        self.sendDialog.show()
        self.sendDialog.accepted.connect(self.sendSuggestion)

    def sendSuggestion(self):
        '''
        Send suggestion with input from dialog.
        '''
        self.sendDialog.accepted.disconnect(self.sendSuggestion)
        try:
            msg = str(self.sendDialog.suggInput.toPlainText())
            emailStr = str(self.sendDialog.emailAdInput.text())

            sendEmail(body = msg, email = emailStr, subject = 'Suggestion: XD Toolkit')

            self.sendDialog.accept()
            self.thanks = QMessageBox()
            self.thanks.setText('Thanks for your feedback!')
            self.thanks.setWindowTitle('Submit Suggestion')
            self.thanks.show()

        except Exception:
            self.msg = QMessageBox()
            self.msg.setText('An error occurred, report not sent. You can instead email "mcrav@chem.au.dk" with your message.')
            self.msg.setWindowTitle('Submit Bug Report Failed')
            self.msg.show()


#########################################################################
#--------------------WIZARD----------------------------------------------
#########################################################################

    def changeWizBackup(self):

        k = 1
        while os.path.isdir('Backup/XD Wizard Run ' + str(k)):
            k += 1
        self.wizBackupFolder = 'XD Wizard Run ' + str(k)
        self.wizBackupInput.setText(self.wizBackupFolder)

    def xdWizINI(self):
        '''
        Run XDINI with compound ID from wiz screen.
        '''
        if os.path.isfile('shelx.ins'):
            global compoundID4XDINI
            compoundID4XDINI = str(self.xdWizCmpID.text())

            if compoundID4XDINI.strip() != '':
                self.xdini.start()
                self.addedLocCoords = {}
                self.disableXDButs()
                self.xdini.finishedSignal.connect(self.enableXDButs)
                self.xdini.finishedSignal.connect(self.wizCheckIni)
                self.xdWizINILab.setText('Initializing compound...')

            else:
                self.xdWizINILab.setText('Type compound ID and run again.')

        else:
            self.xdWizINILab.setText('Invalid project folder.')


    #Check that XDINI has created xd.mas, xd.hkl and xd.inp.
    def wizCheckIni(self):
        '''
        Check that XDINI has created all files and test them.
        Unlock second stage of XD Wizard if all files are present, regardless of test result.
        '''
        self.xdini.finishedSignal.disconnect(self.wizCheckIni)
        if os.path.isfile('xd.mas') and os.path.isfile('xd.inp') and os.path.isfile('xd.hkl'):
            self.xdWizINILab.setText('Compound initialized successfully.')
            self.manRefIDInput.setText(str(self.xdWizCmpID.text()))
            try:
                writeCHEMCON(findCHEMCON())
                autoResetBond(copy.copy(globAtomLabs), copy.copy(globAtomTypes))
                multipoleMagician()
                self.wizTest()
            except Exception as e:
                printExc('wizCheckIni',e)
                pass

            for item in self.wiz2ndStageObjects:
                item.setEnabled(True)

            self.changeWizBackup()

        else:
            self.xdWizINILab.setText('An error occurred. Check xd_ini.out in project folder.')

    def wizReuseMas(self):
        '''
        Rename xdwiz.mas to xd.mas.
        '''
        try:
            if os.path.isfile('xdwiz.mas'):
                os.remove('xd.mas')
                copyfile('xdwiz.mas','xd.mas')
                self.wizTest()
            else:
                self.wizTestStatusLab.setText('''Couldn't find xdwiz.mas file from last run.''')

        except Exception:
            self.wizTestStatusLab.setText('An error occurred.')


    def wizTest(self):
        '''
        Check CHEMCON, reset bond and sin(theta/lambda) have all been inputted appropriately.
        Update GUI labels and return bools of results.
        '''
        c = check4CHEMCON()
        r = check4RBHs()
        s = True
        f = os.path.isfile('xd.mas') and os.path.isfile('xd.inp') and os.path.isfile('xd.hkl')

        missingSym = findUnaddedSym()

        if missingSym:
            m = False

        else:
            m = True

#        try:
#            for item in self.wizSnlInput:
#                a = float(str(item.text()))     #checks if input can be converted to float
#
#                if str(self.wizUniSnlMax.text()).strip():
#                    a = float(str(item.text()))
#
#            s = True
#
#
#        except ValueError:
#            s = False
        for item in self.wizSnlInput:
            if not isfloat(str(item.text())):
                s = False

        uniSnlStr = str(self.wizUniSnlMax.text()).strip()
        if uniSnlStr:
            if not isfloat(uniSnlStr):
                s = False

        if c and not r and s and f and m:

            if self.wizLowSnlMin.isEnabled():
                self.xdWizardBut.setEnabled(True)
#            estTime = totalEstTime()
#            minutes, seconds = divmod(estTime, 60)
#            hours, minutes = divmod(minutes, 60)
#            timeStr = '{0:.0f}hrs {1:.0f}mins {2:.0f}secs'.format(hours, minutes, seconds)
            #wizStr = '{0} Estimated running time is {1}'.format('Ready to run XD Wizard.', timeStr)
            wizStr = 'Ready to run XD Wizard.'
            self.wizTestStatusLab.setText(wizStr)

        else:
            wizStr = ''

            if not c:
                wizStr += '-Add chemical constraints.<br>'

            if r:
                wizStr += '-Add bond constraints for atoms:\n\n{}<br>'.format(listjoin(r, ', '))

            if not m:
                wizStr += '{0}{1}'.format('-Add local coordinate system and local symmetry manually for ', ', '.join(missingSym).strip(', '))

            if not s:
                wizStr += '-Invalid input for sin(&theta;/&lambda;) cutoffs.<br>'
                self.xdWizardBut.setEnabled(False)

            if not f:
                wizStr += '-Starting files not found. Initialize compound.'

            self.wizTestStatusLab.setText('<p>{}</p>'.format(wizStr.strip('<br>')))

        return (c,r,(m,missingSym),s,f)


    def wizFinalTest(self):
        '''
        If Wizard test criteria have not been met, ask user if they want to proceed anyway.
        '''
        testRes = self.wizTest()
        procMsg = '\nProceed anyway?'
        rWarnMsg = 'No bond restraints detected for {}.\n'.format(listjoin(testRes[1], ', '))
        cWarnMsg = 'No chemical constraints detected.\n'
        mWarnMsg = 'No local coordinate system added for {}'.format(listjoin(testRes[2][1], ', '))
        testPassed = True

        if testRes[1] or not testRes[0] or testRes[2][1] or not testRes[3] or not testRes[4]:
            testPassed = False

        if testPassed:
            return True

        else:
            if testRes[3] and testRes[4]:

                msg = procMsg
                if not testRes[0]:
                    msg = cWarnMsg + msg
                if testRes[1]:
                    msg = rWarnMsg + msg
                if not testRes[2][0]:
                    msg = mWarnMsg + msg
                warningMsg = QMessageBox.question(self, 'Warning', msg, QMessageBox.Yes, QMessageBox.No)

                if warningMsg == QMessageBox.Yes:
                    return True
                else:
                    return False
            else:
                return False


    def xdWizRun(self):
        '''
        Start wizard with given settings.
        '''
        try:
            if self.wizFinalTest():
                self.xdWizRunning = wizardRunning()

                try:
                    if self.settings.value('senddata') == 'yes':
                        self.xdWizRunning.collectDat = True
                        print(self.xdWizRunning.collectDat)
                    else:
                        self.xdWizRunning.collectDat = False
                except Exception:
                    self.xdWizRunning.collectDat = False

                self.xdWizRunning.rejected.connect(self.xdWizFinished)
                self.xdWizRunning.finishedSignal.connect(self.xdWizFinished)
                if self.wizSeqMultBox.checkState():
                    self.xdWizRunning.refList = ['0 - XDINI', '1 - Scale factors', '2 - High angle non-H positions and ADPs',
                                '3 - Low angle H positions and isotropic ADPs', '4 - Kappa and monopoles', '5 - Dipoles',
                                '6 - Quadrupoles', '7 - Octupoles', '8 - Hexadecapoles',
                                '9 - Multipoles, and non-H positions and ADPs', '10 - Lower symmetry',
                                '11 - Final refinement']

                if not self.wizLowerSymBox.checkState():
                    self.xdWizRunning.refList = self.xdWizRunning.refList[:-2]

                global wizHighSnlMin
                wizHighSnlMin = float(str(self.wizHighSnlMin.text()))
                global wizHighSnlMax
                wizHighSnlMax = float(str(self.wizHighSnlMax.text()))
                global wizLowSnlMin
                wizLowSnlMin = float(str(self.wizLowSnlMin.text()))
                global wizLowSnlMax
                wizLowSnlMax = float(str(self.wizLowSnlMax.text()))
                global wizUniSnlMax
                if str(self.wizUniSnlMax.text()):
                    wizUniSnlMax = float(str(self.wizUniSnlMax.text()))
                else:
                    wizUniSnlMax = 2.0

                self.xdWizRunning.customLCS = check4CustomLCS()

                copyfile('xd.mas','xdwiz.mas')

                self.xdWizRunning.show()
                self.xdWizRunning.xdWiz(str(self.wizBackupInput.text()))
        except Exception as e:
            print(e)
            self.wizTestStatusLab.setText('An error occured, cannot run wizard.')


    def xdWizFinished(self):
        '''
        When wizard has finished, create normal probability plot and show all results.
        '''
        global timeFileAbsPath

        try:
            self.xdWizRunning.xdlsm.xdlsmRunning.terminate()
        except Exception:
            pass

        try:
            results = self.xdWizRunning.finalRes
            FOUcell()
            self.xdfour.start()
            self.xdfour.wait()
            values = grd2values()
            try:
                self.nppCanvas.setParent(None)
            except Exception:
                print('''couldn't remove canvas''')
            self.fig = plt.figure(facecolor = '#dddddd', figsize=(5,3), dpi=80)
            self.ax = self.fig.add_subplot(1,1,1)

            probplot(values, plot = plt)
            plt.title('Normal probability plot')
            plt.tight_layout()
            self.ax.get_lines()[0].set_marker('.')
            self.ax.get_lines()[0].set_markersize(0.8)
            self.nppCanvas = FigureCanvas(self.fig)
            self.nppWizLayout.addWidget(self.nppCanvas, 0)
            self.setLayout(self.nppWizLayout)
            self.nppCanvas.setMinimumSize(self.nppCanvas.size())
            self.nppCanvas.setMaximumSize(700,200000)
            self.xdWizardStatusLab.setText(results)

            if self.settings.value('senddata') == 'yes':
                self.xdWizRunning.tfin = time.time()
                runningTime = self.xdWizRunning.tfin - self.xdWizRunning.tzero

                with open(timeFileAbsPath,'a') as lsmTimes:
                    lsmTimes.write('{0:10}{1:<15}{2:<13.2f}{3:<13}{4}\n'.format('WIZ', getNumAtoms(), runningTime,  len(self.xdWizRunning.refList), sys.platform))

        except Exception as e:   #Handle user cancelling halfway through
            print(e)
            pass

        finally:

            self.xdWizRunning.accept()
            self.resetWizInput()

            try:
                self.wizScrollArea.verticalScrollBar().setValue(self.wizScrollArea.verticalScrollBar().maximum())
            except Exception:
                print('auto scroll error')

    def resetWizInput(self):
        '''
        Disable wizard 2nd stage.
        '''
        self.wizTestStatusLab.setText('')

        for item in self.wiz2ndStageObjects:
            item.setEnabled(False)

        self.xdWizardBut.setEnabled(False)


#########################################################################
#--------------------XD PROGRAMS-----------------------------------------
#########################################################################

    def disableXDButs(self, survivors = []):
        '''
        Disable all buttons that run XD programs, except given survivors.
        '''
        for but in self.xdProgButs:
            if but not in survivors:
                but.setEnabled(False)
        self.xdWizCmpID.returnPressed.disconnect(self.xdWizINI)

    def enableXDButs(self):
        '''
        Enable all buttons that run XD programs.
        '''
        try:
            self.xdlsm.finishedSignal.disconnect(self.enableXDButs)
        except Exception:
            pass

        try:
            self.xdini.finishedSignal.disconnect(self.enableXDButs)
        except Exception:
            pass

        try:
            self.xdfour.finishedSignal.disconnect(self.enableXDButs)
        except Exception:
            pass

        for but in self.xdProgButs:
            if but != self.xdWizardBut:
                but.setEnabled(True)

        self.xdWizCmpID.returnPressed.connect(self.xdWizINI)


    def startXDLSM(self):
        '''
        Handle XDLSM starting.
        '''
        try:
            self.xdlsm.finishedSignal.disconnect(self.finishedXDLSM)
        except Exception:
            pass
        self.xdlsm.finishedSignal.connect(self.finishedXDLSM)   #Sets up the finishedSignal in case xdlsm.exe was killed the last time and the finishedSignal was disconnected
        for lab in self.XDLSMLabs:
            lab.setText('XDLSM running...')

        for but in self.XDLSMButs:
            but.setText('Cancel')
            but.disconnect()
            but.clicked.connect(self.killXDLSM)

        self.runXDLSMBut.setText(' Cancel')

        self.disableXDButs(survivors = self.XDLSMButs)
        self.xdlsm.startTime = time.time()


    #Called when XDLSM finishes
    def finishedXDLSM(self):
        '''
        Handle XDLSM finishing.
        '''
        for but in self.XDLSMButs:
            try:
                but.disconnect()
            except Exception:
                pass
            but.clicked.connect(lambda: self.xdlsm.start())    #Sets up button again to start XDLSM
            but.setText('Run XDLSM')

        self.runXDLSMBut.setText(' Run XDLSM')

        for lab in self.XDLSMLabs:
            lab.setText('XDLSM finished')                         #Sets status label to 'XDLSM finished'
        try:
            self.check4res()

            self.xdlsm.time = time.time() - self.xdlsm.startTime

            if self.xdlsm.time < 15:
                error = check4errors()
                if not error[0]:
                    if error[1].startswith('  * * * Parameter T and/or NCST should be increased'):

                        addNCST()
                        self.xdlsm.start()

                    elif error[1].startswith('Noble gas configuration not recognized for element Cu'):

                        fixCuNobleGasError()
                        self.xdlsm.start()
        except Exception:
            pass

        try:
            self.xdWizardBut.disconnect()
        except Exception:
            pass
        self.xdWizardBut.clicked.connect(self.xdWizRun)
        self.enableXDButs()


    def killXDLSM(self):
        '''
        Kill XDLSM.
        '''
        try:
            self.enableXDButs()
            self.xdlsm.finishedSignal.disconnect() #Disconnects finished signal so that QThread doesn't changes status label to 'XDLSM finished'. Finished signal reconnected in startXDLSM()
            self.xdlsm.xdlsmRunning.terminate()             #Kills xdlsm.exe
            for but in self.XDLSMButs:
                but.setText('Run XDLSM')
                but.disconnect()
                but.clicked.connect(lambda: self.xdlsm.start())

            self.runXDLSMBut.setText(' Run XDLSM')

            for lab in self.XDLSMLabs:
                lab.setText('XDLSM terminated')       #Sets status label to 'XDLSM terminated'

            self.xdWizardStatusLab.setText('XD Wizard terminated.')
            self.xdWizardBut.disconnect()
            self.xdWizardBut.clicked.connect(self.xdWizRun)

        except Exception:
            pass


    def startXDFOUR(self):
        '''
        Handle XDFOUR starting.
        '''
        for lab in self.XDFOURLabs:
            lab.setText('XDFOUR running...')

        for but in self.XDFOURButs:
            but.setText('Cancel')
            but.disconnect()
            but.clicked.connect(self.killXDFOUR)

        self.xdfour.finishedSignal.connect(self.finishedXDFOUR)   #Sets up the finishedSignal in case XDFOUR.exe was killed the last time and the finishedSignal was disconnected
        self.disableXDButs(survivors = self.XDFOURButs)

    def resetXDFOURButs(self):
        '''
        Reset all XDFOUR buts to their default text and connections.
        '''
        for but in self.XDFOURButs:
            but.disconnect()

        self.setupFOURBut.setText('Make residual density map')
        self.resNPPBut.setText('Make normal probability plot')
        self.pkgXDFOURBut.setText('Run XDFOUR')

        self.setupFOURBut.clicked.connect(self.addFOURIns) #Sets button back to starting XDFOUR
        self.pkgXDFOURBut.clicked.connect(lambda: self.xdfour.start())
        self.resNPPBut.clicked.connect(self.makeNPP)

    def finishedXDFOUR(self):
        '''
        Handle XDFOUR finishing.
        '''
        self.enableXDButs()
        for lab in self.XDFOURLabs:
            lab.setText('XDFOUR finished')                         #Sets status label to 'XDFOUR finished'

        self.resetXDFOURButs()

        self.pkgXDLab.setText('XDFOUR finished')

    def killXDFOUR(self):
        '''
        Kill XDFOUR.
        '''
        try:
            self.xdfour.finishedSignal.disconnect() #Disconnects finished signal so that QThread doesn't changes status label to 'XDFOUR finished'. Finished signal reconnected in startXDFOUR()
            self.xdfour.xdProgRunning.terminate()             #Kills XDFOUR.exe

            for lab in self.XDFOURLabs:
                lab.setText('XDFOUR terminated')       #Sets status label to 'XDFOUR terminated'
            self.enableXDButs()
            self.resetXDFOURButs()
        except Exception:
            pass


    def startXDFFT(self):
        '''
        Handle XDFFT starting.
        '''
        for lab in self.XDFFTLabs:
            lab.setText('XDFFT running')                  #Sets the status label to 'XDFFT running'

        for but in self.XDFFTButs:
            but.setText('Cancel')                      #Sets the button text to 'Cancel'
            but.disconnect()
            but.clicked.connect(self.killXDFFT)        #Makes the button kill XDFFT.exe if clicked
        self.xdfft.finishedSignal.connect(self.finishedXDFFT)   #Sets up the finishedSignal in case XDFFT.exe was killed the last time and the finishedSignal was disconnected
        self.disableXDButs(survivors = self.XDFFTButs)

    def finishedXDFFT(self):
        '''
        Handle XDFFT finishing.
        '''
        self.enableXDButs()
        for lab in self.XDFFTLabs:
            lab.setText('XDFFT finished.')                         #Sets status label to 'XDFFT finished'

        for but in self.XDFFTButs:
            but.disconnect()

        self.pkgXDFFTBut.setText('Run XDFFT')
        self.setupFOURBut.setText('Run XDFOUR')
        self.setupFOURBut.clicked.connect(self.addFOURIns)    #Sets up button again to start XDFOUR
        self.pkgXDFFTBut.clicked.connect(lambda: self.xdfft.start())
        self.pkgXDLab.setText('XDFFT finished')

    def killXDFFT(self):
        '''
        Kill XDFFT.
        '''
        try:
            self.xdfft.finishedSignal.disconnect() #Disconnects finished signal so that QThread doesn't changes status label to 'XDFFT finished'. Finished signal reconnected in startXDFFT()
            self.xdfft.xdProgRunning.terminate()             #Kills XDFFT.exe
            self.pkgXDFFTBut.setText('Run XDFFT')
            self.setupFOURBut.setText('Run XDFOUR')           #Sets XDFFT button back to 'Run XDFFT'
            self.setupFOURStatusLab.setText('XDFFT terminated')       #Sets status label to 'XDFFT terminated'
            for but in self.XDFFTButs:
                but.disconnect()
            self.setupFOURBut.clicked.connect(self.addFOURIns) #Sets button back to starting XDFFT
            self.pkgXDFFTBut.clicked.connect(lambda: self.xdfft.start())
            self.enableXDButs()
        except Exception:
            pass


    def startXDPROP(self):
        '''
        Handle XDPROP starting.
        '''
        for lab in self.XDPROPLabs:
            lab.setText('XDPROP running')

        for but in self.XDPROPButs:
            but.setText('Cancel')
            but.disconnect()
            but.clicked.connect(self.killXDPROP)

        self.disableXDButs(self.XDPROPButs)
        self.xdprop.finishedSignal.disconnect(self.finishedXDPROP)
        self.xdprop.finishedSignal.connect(self.finishedXDPROP)   #Sets up the finishedSignal in case XDPROP.exe was killed the last time and the finishedSignal was disconnected

    def finishedXDPROP(self):
        '''
        Handle XDPROP finishing.
        '''
        for but in self.XDPROPButs:
            but.setText('Run XDPROP')                           #Sets button back to 'Run XDPROP'
            but.disconnect()
        self.getDpopsBut.setText('Get d-orbital populations')
        self.getDpopsBut.clicked.connect(self.getDpops)    #Sets up button again to make XDPROP instructions and run XDPROP
        self.pkgXDPROPBut.clicked.connect(lambda: self.xdprop.start())
        self.pkgXDLab.setText('XDPROP finished')

        self.enableXDButs()

    def killXDPROP(self):
        '''
        Kill XDPROP.
        '''
        try:
            self.xdprop.finishedSignal.disconnect()
            self.xdprop.xdProgRunning.terminate()             #Kills XDPROP.exe
            self.getDpopsBut.setText('Run XDPROP')           #Sets XDPROP button back to 'Run XDPROP'
            self.getDpopsStatusLab.setText('XDPROP terminated')       #Sets status label to 'XDPROP terminated'
            self.getDpopsBut.disconnect()
            self.getDpopsBut.clicked.connect(self.getDpops) #Sets button back to starting XDPROP
             #Disconnects finished signal so that QThread doesn't changes status label to 'XDPROP finished'. Finished signal reconnected in startXDPROP()
            self.pkgXDPROPBut.setText('Run XDPROP')
            self.pkgXDPROPBut.disconnect()
            self.pkgXDPROPBut.clicked.connect(lambda: self.xdprop.start())
        except Exception:
            pass
        finally:
            self.enableXDButs()


    def startXDGEOM(self):
        '''
        Handle XDGEOM starting.
        '''
        for lab in self.XDGEOMLabs:
            lab.setText('XDGEOM running')

        for but in self.XDGEOMButs:
            but.setText('Cancel')
            but.disconnect()
            but.clicked.connect(self.killXDGEOM)

        self.disableXDButs(survivors = self.XDGEOMButs)

        self.xdgeom.finishedSignal.disconnect()
        self.xdgeom.finishedSignal.connect(self.finishedXDGEOM)

    def finishedXDGEOM(self):
        '''
        Handle XDGEOM finishing.
        '''
        for but in self.XDGEOMButs:
            but.setText('Run XDGEOM')
            but.disconnect()
            but.clicked.connect(lambda: self.xdgeom.start())

        self.enableXDButs()

        for lab in self.XDGEOMLabs:
            lab.setText('XDGEOM finished')

    def killXDGEOM(self):
        '''
        Kill XDGEOM.
        '''
        try:
            self.xdgeom.finishedSignal.disconnect()
            self.xdgeom.xdProgRunning.terminate()
            for but in self.XDGEOMButs:
                but.setText('Run XDGEOM')
                but.disconnect()
                but.clicked.connect(lambda: self.xdgeom.start())
            for lab in self.XDGEOMLabs:
                lab.setText('XDGEOM terminated')
        except Exception:
            pass
        finally:
            self.enableXDButs()


    def startXDPDF(self):
        '''
        Handle XDPDF starting.
        '''
        for lab in self.XDPDFLabs:
            lab.setText('XDPDF running')

        for but in self.XDPDFButs:
            but.setText('Cancel')
            but.disconnect()
            but.clicked.connect(self.killXDPDF)

        self.disableXDButs(survivors = self.XDPDFButs)

        self.xdpdf.finishedSignal.disconnect()
        self.xdpdf.finishedSignal.connect(self.finishedXDPDF)

    def finishedXDPDF(self):
        '''
        Handle XDPDF finishing.
        '''
        for but in self.XDPDFButs:
            but.setText('Run XDPDF')
            but.disconnect()
            but.clicked.connect(lambda: self.xdpdf.start())

        self.enableXDButs()

        for lab in self.XDPDFLabs:
            lab.setText('XDPDF finished')

    def killXDPDF(self):
        '''
        Kill XDPDF.
        '''
        try:
            self.xdpdf.finishedSignal.disconnect()
            self.xdpdf.xdProgRunning.terminate()
            for but in self.XDPDFButs:
                but.setText('Run XDPDF')
                but.disconnect()
                but.clicked.connect(lambda: self.xdpdf.start())
            for lab in self.XDPDFLabs:
                lab.setText('XDPDF terminated')
        except Exception:
            pass
        finally:
            self.enableXDButs()


    def startTOPXD(self):
        '''
        Handle TOPXD starting.
        '''
        for lab in self.TOPXDLabs:
            lab.setText('TOPXD running')

        for but in self.TOPXDButs:
            but.setText('Cancel')
            but.disconnect()
            but.clicked.connect(self.killTOPXD)

        self.disableXDButs(survivors = self.TOPXDButs)

        self.topxd.finishedSignal.disconnect()
        self.topxd.finishedSignal.connect(self.finishedTOPXD)

    def finishedTOPXD(self):
        '''
        Handle TOPXD finishing.
        '''
        for but in self.TOPXDButs:
            but.setText('Run TOPXD')
            but.disconnect()
            but.clicked.connect(lambda: self.topxd.start())

        self.enableXDButs()

        for lab in self.TOPXDLabs:
            lab.setText('TOPXD finished')

    def killTOPXD(self):
        '''
        Kill TOPXD.
        '''
        try:
            self.topxd.finishedSignal.disconnect()
            self.topxd.xdProgRunning.terminate()
            for but in self.TOPXDButs:
                but.setText('Run TOPXD')
                but.disconnect()
                but.clicked.connect(lambda: self.topxd.start())
            for lab in self.TOPXDLabs:
                lab.setText('TOPXD terminated')
        except Exception:
            pass
        finally:
            self.enableXDButs()


    def runXDINI(self):
        '''
        Run XDINI with compound ID in manual refinement.
        '''
        global compoundID4XDINI
        compoundID4XDINI = str(self.manRefIDInput.text())
        if compoundID4XDINI.strip() != '':
            self.disableXDButs()
            if os.path.isfile('shelx.ins') and os.path.isfile('shelx.hkl'):
                self.xdini.finishedSignal.connect(self.enableXDButs)
                self.xdini.start()
            else:
                self.XDINILab.setText('Invalid project folder.')
        else:
            self.XDINILab.setText('Type compound ID and run again.')

    def startXDINI(self):
        '''
        Handle XDINI starting.
        '''
        self.XDINILab.setText('Initializing compound...')                  #Sets the status label to 'XDINI running'
        self.runXDINIBut.setText('Cancel')                      #Sets the button text to 'Cancel'
        self.runXDINIBut.clicked.connect(self.killXDINI)        #Makes the button kill XDINI.exe if clicked
        self.xdini.finishedSignal.disconnect(self.finishedXDINI)
        self.xdini.finishedSignal.connect(self.finishedXDINI)   #Sets up the finishedSignal in case XDINI.exe was killed the last time and the finishedSignal was disconnected

    def finishedXDINI(self):
        '''
        Handle XDINI finishing.
        '''
        print('finishedINI')
        self.XDINILab.setText('Compound initialized successfully. Ready to begin refinement.')                         #Sets status label to 'XDINI finished'
        self.xdWizINILab.setText('Compound initialized successfully. Follow instructions below and click "Test".')
        self.xdWizCmpID.setText(str(self.manRefIDInput.text()))
        self.refChosen()
        self.runXDINIBut.setText('Run XDINI')
        self.check4res()
        self.enableXDButs()
        self.changeUserIns()
        self.runXDINIBut.disconnect()                         #Sets button back to 'Run XDINI'
        self.runXDINIBut.clicked.connect(self.runXDINI)    #Sets up button again to start XDINI

        for item in self.wiz2ndStageObjects:
            item.setEnabled(True)

    def killXDINI(self):
        '''
        Kill XDINI.
        '''
        try:
            self.xdini.xdiniRunning.terminate()             #Kills XDINI.exe
            self.runXDINIBut.disconnect()
            self.runXDINIBut.setText('Run XDINI')           #Sets XDINI button back to 'Run XDINI'
            self.XDINILab.setText('XDINI terminated')       #Sets status label to 'XDINI terminated'
            self.runXDINIBut.clicked.connect(self.runXDINI) #Sets button back to starting XDINI
            self.xdini.finishedSignal.disconnect() #Disconnects finished signal so that QThread doesn't changes status label to 'XDINI finished'. Finished signal reconnected in startXDINI()

        except Exception:
            pass


#########################################################################
#--------------------UTILITIES-------------------------------------------
#########################################################################

    def openMerc(self):

        if self.settings.value('mercurypath'):
            self.mercury.start()
        else:
            msg = '''Can't find Mercury. Download or select executable file below.'''
            mercMsg = QMessageBox()
            mercMsg.setText(msg)

            mercMsg.addButton(QPushButton('Download'), QMessageBox.YesRole)
            mercMsg.addButton(QPushButton('Select Mercury executable file'), QMessageBox.NoRole)
            result = mercMsg.exec_()

            if result == 0:
                webbrowser.open('https://www.ccdc.cam.ac.uk/support-and-resources/Downloads/')
            elif result == 1:
                self.openPrefs()

    def openMolecool(self):
        try:
            if self.settings.value('molecoolpath'):
                self.molecoolqt.start()
            else:
                msg = '''Can't find MoleCoolQt. Download or select executable file below.'''
                moleMsg = QMessageBox()
                moleMsg.setText(msg)

                moleMsg.addButton(QPushButton('Download'), QMessageBox.YesRole)
                moleMsg.addButton(QPushButton('Select MoleCoolQt executable file'), QMessageBox.NoRole)
                result = moleMsg.exec_()

                if result == 0:
                    webbrowser.open('http://www.molecoolqt.de/dl.php')
                elif result == 1:
                    self.openPrefs()
        except Exception:
            pass


    def findXD(self):

        msg = '''Couldn't find folder containing the XD programs. Select folder now?'''
        findXD = QMessageBox.question(self, 'Find XD programs', msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)

        if findXD == QMessageBox.Yes:

            folder = str(QFileDialog.getExistingDirectory(None, "Select XD Folder"))
            xdFolder = ''

            if sys.platform == 'win32':
                if 'xdlsm.exe' not in os.listdir(folder):
                    if 'bin' in os.listdir(folder) and 'xdlsm.exe' in os.listdir(folder + '/bin'):
                        xdFolder = folder + '/bin'
                    self.warningMsg = QMessageBox.warning(self, 'Warning', 'Invalid folder chosen.<br><br>Change folder in preferences.')
                    self.warningMsg.show()
                else:
                    xdFolder = folder

                if not os.environ.get['XD_DATADIR']:
                    self.msg = 'Make sure you have the XD_DATADIR environment variable setup. It should be:\n\nXD_DATADIR={}lib/xd\n\n'.format(xdFolder[:-3])
                    self.envMsg = QMessageBox()
                    self.envMsg.setWindowTitle('XD_DATADIR environment variable')
                    self.envMsg.setText(self.msg)
                    self.envMsg.show()

            elif sys.platform.startswith('linux'):
                if 'xdlsm' not in os.listdir(folder):
                    if 'bin' in os.listdir(folder) and 'xdlsm' in os.listdir(folder + '/bin'):
                        xdFolder = folder + '/bin'
                    else:
                        self.warningMsg = QMessageBox.warning(self, 'Warning', 'Invalid folder chosen.<br><br>Change folder in preferences.')
                        self.warningMsg.show()
                else:
                    xdFolder = folder

                if not os.environ.get('XD_DATADIR'):
                    self.msg = 'Make sure you have the XD_DATADIR environment variable setup. To do this go to the terminal and enter:\n\nsudo gedit /etc/environment\n\nNow add the following line to the environment file and save it:\n\nXD_DATADIR={}lib/xd\n\nLogout and log back in for the changes to take effect.'.format(xdFolder[:-3])
                    self.envMsg = QMessageBox()
                    self.envMsg.setWindowTitle('XD_DATADIR environment variable')
                    self.envMsg.setText(self.msg)
                    self.envMsg.show()

            self.settings.setValue('xdpath', xdFolder)
            self.initialiseSettings()


    def badLabelWarning(self):
        '''
        DOESN'T WORK: Warn user if labels in shelx.ins may cause confusion.
        i.e. NA (nitrogenA) and Na (only sodium).
        '''
        badLabs = self.xdini.warnings       #No warnings are currently made upstream
        msg = QMessageBox()
        if len(badLabs) > 1:
            warnStr = 'The following atom labels may cause confusion as the first two characters of their labels could mean 2 different elements in the structure:'
            warnReq = 'Please change their labels to a number or a different letter before starting.'
        else:
            warnStr = 'This atom label may cause confusion as its first two characters could mean 2 different elements in the structure:'
            warnReq = 'Please the change the label to a number or a different letter before starting.'
        warningStr = '{}<br><br>{}<br><br>{}'.format(warnStr, ', '.join(badLabs).strip(', '), warnReq)
        msg.setText(warningStr)
        msg.exec_()


    def check4res(self):
        '''
        Check is xd.res exists.
        '''
        if os.path.isfile('xd.res'):
            self.toolbarRes2Inp.setEnabled(True)
        else:
            self.toolbarRes2Inp.setDisabled(True)

    def res2inpPress(self):
        '''
        Handle res2inp button press.
        '''
        if os.path.isfile('xd.res'):
            self.toolbarRes2Inp.setDisabled(True)
            res2inp()



    def refreshFolder(self):
        '''
        Reinitialise global variables, check xd.res and change user instructions.
        '''
        initialiseGlobVars()
        self.check4res()
        self.changeUserIns()

    def updateFolderTree(self):
        self.folderTree.setRootIndex(self.treeModel.index(QDir.currentPath()))

    def setFolder(self):
        '''
        Prompt user to choose a folder and change current working directory to the folder.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
        if folder:
            os.chdir(folder)
            self.updateFolderTree()
            self.changeWizBackup()
            self.setWindowTitle('{} - XD Toolkit {}'.format(os.path.basename(folder), self.versionNum))
            self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())
            self.resetLabels()
            self.check4res()

            hkl = True
            ins = True
            if not os.path.isfile('shelx.hkl'):
                hkl = False
            if not os.path.isfile('shelx.ins'):
                ins = False

            if not hkl or not ins:
                msg = QMessageBox()
                msg.setWindowTitle('Files not found')
                msgStr = ''
                if not ins and hkl:
                    msgStr = "Couldn't find <i>shelx.ins</i>."
                elif not hkl and ins:
                    msgStr = "Couldn't find <i>shelx.hkl</i>."

                elif not ins and not hkl:
                    msgStr = "Couldn't find <i>shelx.ins</i> or <i>shelx.hkl</i>."

                msg.setText(msgStr)
                msg.exec_()

            try:
                initialiseGlobVars()
                self.enableCheckNeebsButs()
                self.changeUserIns()
                self.resetWizInput()

            except Exception:
                self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())

    def openCwd(self):
        '''
        Open the current working directory in file explorer.
        '''
        if sys.platform == "win32":
            os.startfile(os.getcwd())
        else:
            opener ="open" if sys.platform == "darwin" else "xdg-open"
            subprocess.call([opener, os.getcwd()])

    def openMas(self):
        '''
        Open xd.mas in a text editor.
        '''
        global textEdAbsPath
        try:
            if sys.platform == 'win32':
                subprocess.call([textEdAbsPath,'xd.mas'])

            else:
                opener ="open" if sys.platform == "darwin" else "xdg-open"
                subprocess.call([opener, 'xd.mas'])
        except Exception as e:
            printExc('openMas',e)
            pass

    def loadIAM(self):
        '''
        Load shelx.ins and shelx.hkl from folder of structure solution.
        '''
        ins = False
        hkl = False

        folder = str(QFileDialog.getExistingDirectory(None, "Select Structure Solution Folder"))
        projectFolder = os.getcwd()

        for file in os.listdir(folder):
            if file[-4:] == '.res':
                copyfile((folder + '/' + file), (projectFolder  + '/shelx.ins'))
                ins = True
            elif file[-4:] == '.hkl':
                copyfile((folder + '/' + file), (projectFolder + '/shelx.hkl'))
                hkl = False

        if ins and hkl:
            self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())
            initialiseGlobVars()
            return True

        else:
            print('no files')
            msg = QMessageBox()
            msg.setWindowTitle('Missing  files')
            msgStr = ''
            if not ins:
                msgStr += 'res file not found.\n'
            if not hkl:
                msgStr += 'hkl file not found.'

            msg.setText(msgStr)
            msg.exec_()

    def checkNeebsPress(self):
        self.checkNeebs = checkNeebs()
        self.checkNeebs.show()

    def enableCheckNeebsButs(self):
        '''
        If all global variables exist, enable all check neighbour buttons, otherwise disable them.
        '''
        global globAtomLabs
        global globAtomTypes
        global globAtomPos
        global globAtomAngles

        if globAtomLabs and globAtomTypes and globAtomPos and globAtomAngles:
            for but in self.checkNeebsButs:
                but.setEnabled(True)

        else:
            for but in self.checkNeebsButs:
                but.setEnabled(False)



#########################################################################
#--------------------PREFERENCES-----------------------------------------
#########################################################################

    def initialiseSettings(self):
        '''
        Initialise all user settings.
        '''
        global molecoolQtAbsPath
        global globMercAbsPath
        global textEdAbsPath
        global xdlsmAbsPath
        global xdfourAbsPath
        global xdfftAbsPath
        global xdpropAbsPath
        global xdiniAbsPath
        global xdgeomAbsPath
        global xdgraphAbsPath
        global topxdAbsPath
        global xdpdfAbsPath
        global timeFileAbsPath
        global xdProgAbsPaths
        global manualAbsPath
        global cachePath
        global iconAbsPath

        projDir = os.getcwd()
        print('Installation folder: {}'.format(projDir))

        if sys.platform.startswith('win'):

            #Make folder in ProgramData on Windows if one doesn't already exist.
            XDTDataPath = 'C:/ProgramData/XDToolkit'
            if not os.path.isdir(XDTDataPath):
                os.makedirs(XDTDataPath)

        elif sys.platform.startswith('linux'):
            XDTDataPath = projDir

        #Make cache folder if one doesn't already exist.
        cachePath = XDTDataPath + '/cache'
        if not os.path.isdir(cachePath):
            os.makedirs(cachePath)

        iconAbsPath = projDir + '/res/flatearth.ico'
       # timeFileAbsPath = projDir + '/tools/lsmTimes.buckfast'
        manualAbsPath = projDir + '/res/XD Toolkit Manual.pdf'

#        if not os.path.isfile(timeFileAbsPath):
#            if not os.path.isdir('tools'):
#                os.mkdir('tools')
#            with open('{}/tools/lsmTimes.buckfast'.format(projDir),'w') as bucky:
#                pass

        if self.settings.value('xdpath'):

            xdExecAbsPath = self.settings.value('xdpath')

            if sys.platform.startswith('linux'):
                xdlsmAbsPath = '{}/xdlsm'.format(xdExecAbsPath)
                xdfourAbsPath = '{}/xdfour'.format(xdExecAbsPath)
                xdfftAbsPath = '{}/xdfft'.format(xdExecAbsPath)
                xdiniAbsPath = '{}/xdini'.format(xdExecAbsPath)
                xdgeomAbsPath = '{}/xdgeom'.format(xdExecAbsPath)
                xdpropAbsPath = '{}/xdprop'.format(xdExecAbsPath)
                xdpdfAbsPath = '{}/xdpdf'.format(xdExecAbsPath)
                topxdAbsPath = '{}/topxd'.format(xdExecAbsPath)

            elif sys.platform == 'win32':

                xdlsmAbsPath = '{}/xdlsm.exe'.format(xdExecAbsPath).replace('/','\\')
                xdfourAbsPath = '{}/xdfour.exe'.format(xdExecAbsPath).replace('/','\\')
                xdfftAbsPath = '{}/xdfft.exe'.format(xdExecAbsPath).replace('/','\\')
                xdiniAbsPath = '{}/xdini.exe'.format(xdExecAbsPath).replace('/','\\')
                xdgeomAbsPath = '{}/xdgeom.exe'.format(xdExecAbsPath).replace('/','\\')
                xdpropAbsPath = '{}/xdprop.exe'.format(xdExecAbsPath).replace('/','\\')
                xdpdfAbsPath = '{}/xdpdf.exe'.format(xdExecAbsPath).replace('/','\\')
                topxdAbsPath = '{}/topxd.exe'.format(xdExecAbsPath).replace('/','\\')
                xdgraphAbsPath = '{}/xdgraph.exe'.format(xdExecAbsPath).replace('/','\\')    #Can't find linux exec for xdgraph

        elif sys.platform == 'win32':

            if os.path.isdir('C:/WinXD/bin'):
                self.settings.setValue('xdpath', 'C:/WinXD/bin')
                xdExecAbsPath = 'C:/WinXD/bin'

                xdlsmAbsPath = '{}/xdlsm.exe'.format(xdExecAbsPath)
                xdfourAbsPath = '{}/xdfour.exe'.format(xdExecAbsPath)
                xdfftAbsPath = '{}/xdfft.exe'.format(xdExecAbsPath)
                xdiniAbsPath = '{}/xdini.exe'.format(xdExecAbsPath)
                xdgeomAbsPath = '{}/xdgeom.exe'.format(xdExecAbsPath)
                xdpropAbsPath = '{}/xdprop.exe'.format(xdExecAbsPath)
                xdpdfAbsPath = '{}/xdpdf.exe'.format(xdExecAbsPath)
                topxdAbsPath = '{}/topxd.exe'.format(xdExecAbsPath)
                xdgraphAbsPath = '{}/xdgraph.exe'.format(xdExecAbsPath)    #Can't find linux exec for xdgraph

        try:
            xdProgAbsPaths = {'xdlsm':xdlsmAbsPath, 'xdfour':xdfourAbsPath, 'xdfft':xdfftAbsPath,
                             'xdini':xdiniAbsPath, 'xdgeom':xdgeomAbsPath, 'xdprop':xdpropAbsPath,
                             'xdpdf':xdpdfAbsPath, 'topxd': topxdAbsPath}
        except Exception:
            pass

        if self.settings.value('molecoolpath'):
            molecoolQtAbsPath = self.settings.value('molecoolpath')
        elif sys.platform == 'win32':
            if os.path.isfile(os.path.expanduser('~') + '\\AppData\\Local\\MoleCoolQt\\molecoolqt.exe'):
                molecoolQtAbsPath = os.path.expanduser('~') + '\\AppData\\Local\\MoleCoolQt\\molecoolqt.exe'
            elif os.path.isfile(os.path.expanduser('~') + '\\AppData\\Local\\MoleCoolQt64\\molecoolqt.exe'):
                molecoolQtAbsPath = os.path.expanduser('~') + '\\AppData\\Local\\MoleCoolQt64\\molecoolqt.exe'

        if self.settings.value('mercurypath'):
            globMercAbsPath = self.settings.value('mercurypath')

        if self.settings.value('textedpath'):
            textEdAbsPath = self.settings.value('textedpath')
        else:
            if sys.platform.startswith('linux'):
                if os.path.isfile('/usr/bin/gedit'):
                    textEdAbsPath = '/usr/bin/gedit'
                    self.settings.setValue('textedpath','/usr/bin/gedit')
            elif sys.platform == 'win32':
                if os.path.isfile('C:/Windows/System32/notepad.exe'):
                    textEdAbsPath = 'C:/Windows/System32/notepad.exe'
                    self.settings.setValue('textedpath','C:/Windows/System32/notepad.exe')

        if not self.settings.value('senddata'):
            msg = '''Would you like to send anonymous usage data to help improve XD Toolkit?<br><br>
            You can always change this later in 'Preferences'.'''
            self.sendAnonData = QMessageBox.question(self, 'Send Anonymous Usage Data',
                            msg, QMessageBox.Yes, QMessageBox.No)

            if self.sendAnonData == QMessageBox.Yes:
                self.settings.setValue('senddata','yes')
            elif self.sendAnonData == QMessageBox.No:
                self.settings.setValue('senddata','no')

        try:
            os.chdir(self.settings.value('lastcwd'))
            self.setWindowTitle('{} - XD Toolkit {}'.format(os.path.basename(os.getcwd()), self.versionNum))
        except Exception:
            pass
        self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())



    def openPrefs(self):
        '''
        Open preferences window.
        '''
        self.prefs = prefGui()
        self.prefs.accepted.connect(self.updateNewPrefs)
        try:
            self.prefs.settingsXDLab.setText('Current path: ' + self.settings.value('xdpath', type = str))
        except Exception:
            pass
        try:
            self.prefs.settingsMoleCoolLab.setText('Current path: ' + self.settings.value('molecoolpath', type = str))
        except Exception:
            pass
        try:
            self.prefs.chooseMercPathLab.setText('Current path: ' + self.settings.value('mercurypath', type = str))
        except Exception:
            pass
        try:
            self.prefs.chooseTextPathLab.setText('Current path: ' + self.settings.value('textedpath', type = str))
        except Exception:
            pass
        try:
            if self.settings.value('senddata') == 'yes':
                self.prefs.anonDatBox.setChecked(True)
            elif self.settings.value('senddata') == 'no':
                self.prefs.anonDatBox.setChecked(False)
        except Exception:
            pass
        self.prefs.exec_()


    def updateNewPrefs(self):
        '''
        Update settings when user closes preferences window.
        '''
        try:
            self.settings.setValue('xdpath', self.prefs.xdpath)
        except Exception:
            pass

        try:
            self.settings.setValue('molecoolpath', self.prefs.molecoolpath)
        except Exception:
            pass

        try:
            self.settings.setValue('mercurypath', self.prefs.mercurypath)
        except Exception:
            pass

        try:
            self.settings.setValue('textedpath', self.prefs.textedpath)
        except Exception:
            pass

        if self.prefs.anonDatBox.isChecked():
            self.settings.setValue('senddata','yes')
        else:
            self.settings.setValue('senddata','no')

        self.initialiseSettings()


#########################################################################
#--------------------REFINEMENTS-----------------------------------------
#########################################################################

    def setupRef(self):
        '''
        Handle user clicking 'setup xd.mas' in 'Manual refinement' tab.
        Take input and setup xd.mas accordingly.
        '''
        refNum = self.manRefDropMenu.currentIndex() + 1
        successStr = 'xd.mas updated. Ready to run XDLSM.'
        failureStr = 'Error encountered. xd.mas not setup.'
        snlErrMsg = '<p>Please enter valid values for sin(&theta;/&lambda;) cutoffs.</p>'

        snlMin = str(self.manRefSnlMin.text()).strip()
        snlMax= str(self.manRefSnlMax.text()).strip()
        r = check4RB()
        rWarnMsg = 'No bond constraints detected.\n'
        c = check4CHEMCON()
        cWarnMsg = 'No chemical constraints detected.\n'
        procMsg = '\nProceed anyway?'

        validSnl = ((not snlMin and not snlMax) or (isfloat(snlMin) and isfloat(snlMax)))

        if not validSnl:
            self.manRefSetupLab.setText(snlErrMsg)

        else:

            if refNum == 1:                 #Scale factors
                try:
                    scaleFacRef()
                    self.manRefSetupLab.setText(successStr)
                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            elif refNum == 2:               #High angle non-H positions and ADPs

                try:
                    if snlMin and snlMax:               #Check snl isn't blank
                        if not r:
                            msg = rWarnMsg + procMsg
                            warningMsg = QMessageBox.question(self, 'Warning',
                            msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

                            if warningMsg == QMessageBox.Yes:
                                highAngleRef(float(snlMin), float(snlMax))
                                self.manRefSetupLab.setText(successStr)
                        else:
                            highAngleRef(float(snlMin), float(snlMax))
                            self.manRefSetupLab.setText(successStr)



                    else:
                        self.manRefSetupLab.setText(snlErrMsg)

                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            elif refNum == 3:               #Low angle H positions
                try:

                    if snlMin and snlMax:
                        if not r:
                            msg = rWarnMsg + procMsg
                            warningMsg = QMessageBox.question(self, 'Warning',
                            msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

                            if warningMsg == QMessageBox.Yes:
                                lowAngleRef(float(snlMin), float(snlMax))
                                self.manRefSetupLab.setText(successStr)

                        else:
                            lowAngleRef(float(snlMin), float(snlMax))
                            self.manRefSetupLab.setText(successStr)

                    else:
                        self.manRefSetupLab.setText(snlErrMsg)

                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            elif refNum == 4:               #Kappa and monopoles
                try:
                    if not c or not r:
                        msg = procMsg
                        if not c:
                            msg = c + msg
                        if not r:
                            msg = r + msg
                        warningMsg = QMessageBox.question(self, 'Warning', msg, QMessageBox.Yes, QMessageBox.No)

                        if warningMsg == QMessageBox.Yes:
                            kapMonRef()
                            self.manRefSetupLab.setText(successStr)

                    else:
                        kapMonRef()
                        self.manRefSetupLab.setText(successStr)

                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            elif refNum == 5:               #Dipoles
                self.manRefMultipoles(c,r,1)

            elif refNum == 6:               #Quadrupoles
                self.manRefMultipoles(c,r,2)

            elif refNum == 7:               #Octupoles
                self.manRefMultipoles(c,r,3)

            elif refNum == 8:               #Hexadecapoles
                self.manRefMultipoles(c,r,4)

            elif refNum == 9:               #Multipoles, non-H positions and ADPs
                try:
                    if not c or not r:
                        msg = procMsg
                        if not c:
                            msg = cWarnMsg + msg
                        if not r:
                            msg = rWarnMsg + msg
                        warningMsg = QMessageBox.question(self, 'Warning', msg, QMessageBox.Yes, QMessageBox.No)

                        if warningMsg == QMessageBox.Yes:
                            posADPMultRef()
                            missingSym = findUnaddedSym()

                            if missingSym:
                                for item in missingSym:
                                    if item in self.addedLocCoords.keys():
                                        x = self.addedLocCoords[item]
                                        addCustomLocCoords(x[0], x[1], x[2], x[3], x[4], x[5])
                                multipoleKeyTable()

                                missingSym = findUnaddedSym()
                                if missingSym:
                                    statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add local symmetry and local coordinate system manually in 'Tools' tab.'''
                                    self.manRefSetupLab.setText(statusStr)
                                else:
                                    self.manRefSetupLab.setText(successStr)
                            else:
                                self.manRefSetupLab.setText(successStr)

                    else:
                        posADPMultRef()
                        missingSym = findUnaddedSym()

                        if missingSym != []:
                            statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add local symmetry and local coordinate system manually in 'Tools' tab.'''
                            self.manRefSetupLab.setText(statusStr)
                        else:
                            self.manRefSetupLab.setText(successStr)

                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            elif refNum == 10:              #Lower symmetry
                try:

                    posADPMultRef()
                    lowerSym()
                    self.manRefSetupLab.setText(successStr)


                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            elif refNum == 11:              #Final refinement
                try:
                    checkMultRes()
                    self.manRefSetupLab.setText(successStr)

                except Exception as e:
                    printExc('setupRef',e)
                    self.manRefSetupLab.setText(failureStr)

            if snlMin:
                addSnlCutoff(float(snlMin), float(snlMax))


    def manRefMultipoles(self,c,r,l):
        '''
        Handle user setting up multipole refinements in 'Manual refinement' tab.
        Setup xd.mas according to user input.
        '''
        try:
            successStr = 'xd.mas updated. Ready to run XDLSM.'
            failureStr = 'Error encountered. xd.mas not setup.'
            rWarnMsg = 'No bond constraints detected.\n'
            cWarnMsg = 'No chemical constraints detected.\n'
            procMsg = '\nProceed anyway?'
            if not c or not r:
                msg = procMsg
                if not c:
                    msg = cWarnMsg + msg
                if not r:
                    msg = rWarnMsg + msg
                warningMsg = QMessageBox.question(self, 'Warning', msg, QMessageBox.Yes, QMessageBox.No)

                if warningMsg == QMessageBox.Yes:
                    multipoleMagician()
                    seqMultRef(l)
                    missingSym = findUnaddedSym()

                    if missingSym:
                        for item in missingSym:
                            if item in self.addedLocCoords.keys():
                                x = self.addedLocCoords[item]
                                addCustomLocCoords(x[0], x[1], x[2], x[3], x[4], x[5])
                        multipoleKeyTable()

                        missingSym = findUnaddedSym()
                        if missingSym:
                            statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add local symmetry and local coordinate system manually in 'Tools' tab.'''
                            self.manRefSetupLab.setText(statusStr)
                        else:
                            self.manRefStatusLab.setText(successStr)
                    else:
                        self.manRefSetupLab.setText(successStr)

            else:
                multipoleMagician()
                seqMultRef(l)
                missingSym = findUnaddedSym()

                if missingSym:
                    for item in missingSym:
                        if item in self.addedLocCoords.keys():
                            x = self.addedLocCoords[item]
                            addCustomLocCoords(x[0], x[1], x[2], x[3], x[4], x[5])
                    multipoleKeyTable()

                    missingSym = findUnaddedSym()
                    if missingSym:
                        statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add local symmetry and local coordinate system manually in 'Tools' tab.'''
                        self.manRefSetupLab.setText(statusStr)
                    else:
                        self.manRefStatusLab.setText(successStr)
                else:
                    self.manRefSetupLab.setText(successStr)

        except:
            self.manRefSetupLab.setText(failureStr)


    def refChosen(self):
        '''
        Update GUI labels based on what refinement is selected in drop down menu.
        '''
        refName = self.manRefDropMenu.currentText().strip()
        self.manRefBackupInput.setText(refName)
        refNum = self.manRefDropMenu.currentIndex() + 1

        #Set snl boxes accordingly for low/high angle refinement
        if refNum == 2:
            self.manRefSnlMin.setText('0.7')
            self.manRefSnlMax.setText('2.0')

        elif refNum == 3:
            self.manRefSnlMin.setText('0.0')
            self.manRefSnlMax.setText('0.7')

        else:
            self.manRefSnlMin.setText('')
            self.manRefSnlMax.setText('')

        self.changeUserIns()


    def changeUserIns(self):
        '''
        Change user instructions in 'Manual refinement' tab based on xd.mas and selected refinement.
        '''
        refNum = self.manRefDropMenu.currentIndex() + 1
        reqStr = ''

        if refNum in (4,5,6,7,8,9,10):

            c = check4CHEMCON()

            if not c:
                reqStr += '''- Add chemical constraints before setting up xd.mas.\n'''

        if refNum in (2,3,4,5,6,7,8,9,10,11):

            r = check4RB()

            if not r:
                reqStr += '''- Add bond constraints before setting up xd.mas.'''

        self.manRefSetupLab.setText(reqStr)


    def showResRef(self):
        '''
        Show results from xd_lsm.out in 'Manual refinement' tab.
        '''
        try:
            dmsda = getDMSDA('xd_lsm.out')
            if getConvergence('xd_lsm.out'):
                convStr = 'Yes'
            else:
                convStr = 'No'
            resStr = 'RF<sup>2</sup> = {0: .2f} %<br>Convergence - {1}<br>Average DMSDA = {2}<br>Max DMSDA = {3}'.format(getRF2(), convStr, dmsda[0], dmsda[1])

            i = 0
            if len(dmsda[2]) > 14:
                resStr += '<br><br>{0:23}{1:24}{0:23}{1}<br>'.format('INTERATOMIC VECTOR','DMSDA')
                while i+14 < len(dmsda[2]):

                    resStr += '<br>{0:7} ---> {1:7}  {2:< 25.0f}{3:7} ---> {4:7}  {5: .0f}'.format(dmsda[2][i][0], dmsda[2][i][1], dmsda[2][i][2], dmsda[2][i+14][0], dmsda[2][i+14][1], dmsda[2][i+14][2])
                    i += 1

                while i < 14:
                    resStr += '<br>{0:7} ---> {1:7}  {2:< 25.0f}'.format(dmsda[2][i][0], dmsda[2][i][1], dmsda[2][i][2])
                    i+=1
            else:
                resStr += '<br><br>{0:23}{1:24}<br>'.format('INTERATOMIC VECTOR','DMSDA')
                while i < len(dmsda[2]):
                    resStr += '<br>{0:7} ---> {1:7}  {2:< 25.0f}'.format(dmsda[2][i][0], dmsda[2][i][1], dmsda[2][i][2])
                    i += 1
            resStr = '<pre>' + resStr +'</pre>'

        #If full results can't be printed, try to print without DMSDA and if that doesn't work print error message
        except Exception:
            try:
                resStr = 'RF<sup>2</sup> = {0}<br>{1}<br>No DMSDA results found.'.format(getRF2(), getConvergence('xd_lsm.out'))
            except Exception:
                resStr = 'An error occurred.'

        self.manRefResLab.setText(resStr)


#########################################################################
#--------------------RESET BOND------------------------------------------
#########################################################################

    def armRBPress(self):
        '''
        Handle user clicking "Enable all reset bond instructions".
        '''
        try:
            armRBs()
            self.editRBLab.setText('All reset bond instructions enabled.')
        except Exception as e:
            print(e)
            self.editRBLab.setText('An error occurred.')

    def disarmRBPress(self):
        '''
        Handle user clicking "Disable all reset bond instructions".
        '''
        try:
            disarmRBs()
            self.editRBLab.setText('All reset bond instructions disabled.')
        except Exception as e:
            print(e)
            self.editRBLab.setText('An error occurred.')

    def delResetBondPress(self):
        '''
        Handle user clicking "Remove all reset bond instructions".
        '''
        try:
            delResetBond()
            self.editRBLab.setText('All reset bond instructions removed from xd.mas.')
        except Exception as e:
            print(e)
            self.editRBLab.setText('An error occurred.')

    def autoResetBondPress(self):
        '''
        Handle user clicking automatically add reset bond instructions.
        '''
        try:
            missedAtoms = autoResetBond(copy.copy(globAtomLabs), copy.copy(globAtomTypes))

            if len(missedAtoms) > 0:
                resStr = 'Appropriate bond constraints not found for atoms:' + '\n' + '\n'

                i = 0
                for atom in missedAtoms:
                    resStr += atom + ', '
                    i += 1

                    if i == 7:
                        resStr += '\n'
                        i = 0

                resStr = resStr[:-2]

                resStr += '\n' + '\n' + ' Please add these bond constraints manually before continuing.'
            else:
                resStr =  'Bond constraints added for all H atoms.'
            self.autoResetBondStatusLab.setText(resStr)
            self.resetBondStatusLab.setText('')
        except Exception as e:
            print(e)
            self.autoResetBondStatusLab.setText('An error occurred.')
        self.editRBLab.setText('')
        self.changeUserIns()
        self.wizTest()

    def resetBondInput(self):
        '''
        Handle user manually adding reset bond instructions.
        '''
        statusStr = ''
        try:
            #If 'All' is unchecked get atoms labels from the input and pass them to resetBond()
            if  self.allResetBox.isChecked() != True:

                resetAtomRaw = str(self.resetHInput.text().upper())
                if resetAtomRaw.strip():
                    if ',' in resetAtomRaw:
                        resetAtomRaw = resetAtomRaw.replace(' ','').strip(',')
                        resetAtomLabels = resetAtomRaw.split(',')
                    elif ' ' in resetAtomRaw:
                        resetAtomLabels = resetAtomRaw.split()
                    else:
                        resetAtomLabels = [resetAtomRaw.strip()]
                    resetAtomRawList = [label for label in resetAtomLabels]
                    resetAtomList = []
                    resetAtomStr = ''

                    for atomLab in resetAtomRawList:
                        if '(' not in atomLab:
                            if atomLab.startswith('H') and len(atomLab) > 1:
                                newAtomLab = 'H({})'.format(atomLab[1:])
                            else:
                                newAtomLab = 'H({})'.format(atomLab)
                            resetAtomList.append(newAtomLab)

                            resetAtomStr += (newAtomLab + ', ')
                        else:
                            resetAtomList.append(atomLab)
                    resetAtomList = sorted(list(set(resetAtomList)))
                    atomLabs = copy.copy(globAtomLabs)
                    atomLabs = atomLabs.keys()
                    wrongAtoms = [atom for atom in resetAtomList if atom not in atomLabs or atom[:2] != 'H(']
                    resetAtomList = [atom for atom in resetAtomList if atom not in wrongAtoms]

                    resetAtomStr = listjoin(resetAtomList, ', ')

                    resetBond(str(self.resetLengthInput.text()),resetAtomList,False)

                    if resetAtomList:
                        if len(resetAtomList)>1:
                            statusStr += 'Bond constraints added for {}<br>Bond length: {}<br><br>'.format(resetAtomStr, str(self.resetLengthInput.text()))
                        else:
                            statusStr += 'Bond constraint added for {}<br>Bond length: {}<br><br>'.format(resetAtomStr, str(self.resetLengthInput.text()))

                    if wrongAtoms:
                        statusStr += 'Following atom labels are incorrect: {}'.format(listjoin(wrongAtoms, ', '))

                else:
                    statusStr += 'No atoms selected'

            #If 'All' is checked run resetBond() with allornot True
            else:
                resetBond(str(self.resetLengthInput.text()),None,True)
                statusStr += 'Bond constraints added for all atoms<br>Bond length: {}'.format(str(self.resetLengthInput.text()))

            self.resetBondStatusLab.setText(statusStr)
            self.autoResetBondStatusLab.setText('')

        except Exception as e:
            print(e)
            self.resetBondStatusLab.setText('An error occurred.')

        self.changeUserIns()
        self.wizTest()

    def disableResetText(self):
        '''
        If 'All' box is checked, disable line edit for manually adding reset bonds.
        '''
        if self.allResetBox.isChecked():
            self.resetHInput.setEnabled(False)
        else:
            self.resetHInput.setEnabled(True)


#########################################################################
#--------------------CHEMICAL CONSTRAINTS--------------------------------
#########################################################################

    def autoCHEMCONPress(self):
        '''Handle user clicking automatic chemical constraint button.'''
        global globAtomEnv
        try:
            self.autoCHEMCONRunning = FindCHEMCONThread()
            self.autoCHEMCONRunning.finishedSignal.connect(
                    lambda: self.autoCHEMCONLab.setText(
                            'Chemical constraints added and xd.mas udpated.'))
            self.autoCHEMCONRunning.start()
            self.msg = QMessageBox()
            self.autoCHEMCONRunning.finishedSignal.connect(lambda: self.msg.accept())
            self.autoCHEMCONRunning.finishedSignal.connect(
                    lambda: writeCHEMCON(self.autoCHEMCONRunning.chemcon))
            self.msg.addButton(QMessageBox.Cancel)
            self.msg.setText('Finding chemical equivalent atoms...')
            self.msg.setFont(self.programFont)
            self.msg.setWindowTitle('Chemical constraints')
            self.msg.show()
    
            if self.msg == QMessageBox.Cancel:
                self.autoCHEMCONRunning.finishedSignal.disconnect(
                        lambda: writeCHEMCON(self.autoCHEMCONRunning.chemcon))
                self.autoCHEMCONRunning.finishedSignal.disconnect(
                        lambda: self.autoCHEMCONLab.setText(
                                'Chemical constraints added and xd.mas udpated.'))
        except PermissionError:
            self.autoCHEMCONLab.setText(self.permErrorMsg)
        except Exception as e:
            print(e)
            self.autoCHEMCONLab.setText('An error occurred.')
        self.CHEMCONStatusLab.setText('')

    def runCHEMCON(self):
        '''
        Handle user clicking 'Add chemical constraints'.
        '''
        global globAtomEnv
        chemcon = {}
            
        if self.chooseElementCHEMCON.isChecked() == True:

            inputText = str(self.inputElementCHEMCON.text().upper())
            inputElementList = labels2list(inputText)
            try:
                chemcon = findCHEMCONbyInputElement(inputElementList)
                writeCHEMCON(chemcon)
                self.CHEMCONStatusLab.setText('{} {}'.format('Chemical constraints added and xd.mas updated for elements', ', '.join(inputElementList).strip(', ')))
            except PermissionError:
                self.CHEMCONStatusLab.setText(self.permErrorMsg)
            except Exception as e:
                print(e)
                self.CHEMCONStatusLab.setText('An error occurred.')

            self.inputElementCHEMCON.setText('')

        elif self.chooseAtomCHEMCON.isChecked() == True:

            inputText = str(self.inputAtomCHEMCON.text().upper())

            if 'SHOW ME THE TRUTH' in inputText:
                webbrowser.open('https://www.youtube.com/watch?v=ZuToYSYdJS0')

            elif 'RESETMAS' in inputText:
                resetmas()
                self.changeUserIns()
                self.CHEMCONStatusLab.setText('xd.mas reset to test file')

            else:
                inputAtomList = rawInput2labels(inputText)

                atomLabs = copy.copy(globAtomLabs)
                atomLabs = atomLabs.keys()

                inputAtomList = sorted(list(set(inputAtomList)))
                wrongLabs = [atom for atom in inputAtomList if atom not in atomLabs]

                inputAtomList = [atom for atom in inputAtomList if atom not in wrongLabs]

                chemconAddedStr = ', '.join([lab for lab in inputAtomList]).strip(', ')

                try:
                    statusLabStr = ''
                    if inputAtomList:
                        chemcon = findCHEMCONbyInputAtoms(inputAtomList)
                        print(chemcon)
                        writeCHEMCON(chemcon)
                        statusLabStr = ('{} {}'.format('Chemical constraints added and xd.mas updated for', chemconAddedStr))
                    if wrongLabs:
                        statusLabStr += '<br><br>{}{}'.format('\nFollowing atom labels are incorrect: ', ', '.join(wrongLabs).strip(', '))
                    self.CHEMCONStatusLab.setText(statusLabStr)
                except PermissionError:
                    self.CHEMCONStatusLab.setText(self.permErrorMsg)
                except Exception as e:
                    print(e)
                    self.CHEMCONStatusLab.setText('An error occurred.')
            self.inputAtomCHEMCON.setText('')

        if chemcon:
            envSig = 0
            for atom, equis in chemcon.items():
                while str(envSig) in globAtomEnv.values():
                    envSig += 1
                globAtomEnv[atom] = str(envSig)
                for item in equis:
                    globAtomEnv[item] = str(envSig)
        self.autoCHEMCONLab.setText('')
        self.changeUserIns()
        self.wizTest()


#########################################################################
#--------------------TOOLS-----------------------------------------------
#########################################################################

    def addDUMPress(self):
        '''
        Add dummy atom from user input.
        '''
        dumInputBoxes = [self.DUMnum, self.Atom1, self.Atom2]

        inputDUMI = str(self.DUMnum.text())

        if inputDUMI.upper().startswith('DUM'):
            inputDUMI = inputDUMI[3:]

        if not inputDUMI.isdigit():
            self.addDUMLab.setText('Invalid input given.')
            return

        atom1 = rawInput2labels(str(self.Atom1.text()))[0]
        atom2 = rawInput2labels(str(self.Atom2.text()))[0]

        try:
            x = addDUM(atom1, atom2, inputDUMI) #Add dummy atom and store index actually used to x
        except Exception:
            self.addDUMLab.setText('An error occurred. Check input.')

        if x == inputDUMI:
            self.addDUMLab.setText('DUM{0} added to xd.mas.'.format(inputDUMI))
        else:
            self.addDUMLab.setText('DUM{0} already exists.\nDUM{1} added to xd.mas.'.format(inputDUMI,x))

        for item in dumInputBoxes:
            item.setText('')

    def alcs(self):
        '''
        Add local  coordinate system from user input.
        '''
        alcsInput = [self.alcsPAtomInput, self.alcsAtom1Input, self.alcsAtom2Input,
                     self.alcsLocSymInput]

        try:
            Patom = str(self.alcsPAtomInput.text())
            Patom = rawInput2labels(Patom)[0]
            axis1 = str(self.alcsAxis1Input.currentText()).upper()
            atom1 = str(self.alcsAtom1Input.text())
            atom1 = rawInput2labels(atom1)[0]
            axis2 = str(self.alcsAxis2Input.currentText()).upper()
            atom2 = str(self.alcsAtom2Input.text())
            atom2 = rawInput2labels(atom2)[0]
            sym = str(self.alcsLocSymInput.text()).upper()

            x = addCustomLocCoords(Patom, atom1, axis1, atom2, axis2, sym)
            if x == 'aw':
                self.addedLocCoords[Patom] = (Patom, atom1, axis1, atom2, axis2, sym)
                self.alcsStatusLab.setText('xd.mas updated with local coordinate system for ' + Patom)

                for item in alcsInput:
                    item.setText('')

            elif x == 'a':
                self.alcsStatusLab.setText('Atom table updated, but program unable to update key table.')
            elif x == 'w':
                self.alcsStatusLab.setText('Key table updated, but program unable to update atom table.')
            else:
                self.alcsStatusLab.setText('An error occurred. Check input.')
        except Exception:
            self.alcsStatusLab.setText('An error occurred. Check input.')


    def multKeyPress(self):
        '''
        Handle user pressing 'Add multipoles to key table button.
        '''
        try:
            multipoleKeyTable()
            self.multKeyStatusLab.setText('Key table updated in xd.mas.')

        except Exception:
            self.multKeyStatusLab.setText('An error occurred.')


#########################################################################
#--------------------RESULTS---------------------------------------------
#########################################################################

    def tabPropCheckInput(self):
        inputOk = True
        if not isfloat(str(self.tabPropCpRmin.text())) or not isfloat(str(self.tabPropCpRmin.text())):
            inputOk = False

        return inputOk

    def tabPropRun(self):

        if self.tabPropCheckInput():

            tabPropSetupMas(float(str(self.tabPropCpRmin.text())), float(str(self.tabPropCpRmax.text())))
            self.xdgeom.finishedSignal.connect(self.tabPropRunXDPROP)
            self.xdgeom.start()
            self.tabPropLab.setText('Running XDGEOM...')
        else:
            self.tabPropLab.setText('Invalid input given.')

    def tabPropRunXDPROP(self):
        self.xdprop.finishedSignal.connect(self.tabPropFinish)
        self.xdprop.start()
        self.tabPropLab.setText('Running XDPROP...')

    def tabPropFinish(self):

        self.tabPropLab.setText('Creating CSV file')

        if self.tabPropLength.isChecked():
            length = True
        if self.tabPropAngle.isChecked():
            angle = True
        if self.tabPropRho.isChecked():
            rho = True
        if self.tabPropD2rho.isChecked():
            d2rho = True
        if self.tabPropEllip.isChecked():
            ellip = True

        if angle:
            geom2Angles()

        self.tabPropLab.setText('Finished. Results saved to <i>"{}/bond_properties.csv"</i>'.format(os.getcwd()))

    def make3AtomLapmap(self):
        atoms = []
        try:
            atoms = rawInput2labels(str(self.taLapmapInput.text()))

            if len(atoms)>3:
                self.taLapmapLab.setText('More than 3 atoms given.')
                return
            elif len(atoms)<3:
                self.taLapmapLab.setText('Less than 3 atoms given.')
                return

        except Exception as e:
            print(e)
            self.taLapmapLab.setText("Can't read atom labels. Try space separated labels i.e. 'C(1) C(2) C(3)'")

        stepsize = str(self.taLapmapStepsize.text())
        npoints = str(self.taLapmapNpoints.text())

        npointsFine = npoints.isdigit()
        stepsizeFine = isfloat(stepsize)

        if not npointsFine:
            self.taLapmapLab.setText('Invalid value given for number of points.')
            return
        elif not stepsizeFine:
            self.taLapmapLab.setText('Invalid value given for stepsize.')
            return
        elif not stepsizeFine and not npointsFine:
            self.taLapmapLab.setText('Invalid values given for number of points and stepsize.')
            return

        setup3AtomLapmap(atoms, npoints, stepsize)

        self.xdprop.finishedSignal.connect(lambda: self.showLapmap(self.taLapmapLab))
        self.taLapmapLab.setText('Running XDPROP...')
        self.xdprop.start()

    def showLapmap(self, statusLabel):
        '''
        Make 2D laplacian map from xd_d2rho.grd file.
        '''
        if os.path.isfile('xd_d2rho.grd'):
            statusLabel.setText('Making 2D laplacian map...')
            statusLabel.repaint()
            QApplication.processEvents()
            try:
                self.xdprop.finishedSignal.disconnect(self.showLapmap)
            except Exception:
                pass
            if os.path.isfile('xd_d2rho.grd'):
                self.lapmap = resmap('xd_d2rho.grd')
                self.lapmap.show()
                statusLabel.setText('')
                self.lapmap.setWindowTitle('Quickplot laplacian map')
        else:
            statusLabel.setText('No file called "xd_d2rho.grd" found.')

    def make3DLapGrd(self):
        '''
        WORK IN PROGRESS: Make a 3D laplacian grd file around a given atom.
        '''
        pass


    def makeResStr(self, lsmOutFile):
        '''
        Create string of formatted results from xd_lsm.out file. Return string.
        '''
        try:
            dmsda = getDMSDA(lsmOutFile)
            conv=''
            c = getConvergence(lsmOutFile)
            if c:
                conv = 'Yes'
            else:
                conv = 'No'

            resStr = 'RF<sup>2</sup> = {0: .2f} %<br>Convergence – {1}<br>Average DMSDA = {2}<br>Max DMSDA = {3}'.format(getRF2(lsmOutFile), conv, dmsda[0], dmsda[1])
            print(resStr)
            SUs = readSUs(lsmOutFile)
            resStr+='<br><br>Largest standard uncertainties<br>'
            for item in reversed(SUs[-10:]):
                resStr+= '{0:8}{1:8}{2}<br>'.format(item[0], item[1], item[2])
            resStr=resStr[:-4]
            i = 0

            if len(dmsda[2]) > 14:
                resStr += '<br><br>{0:23}{1:24}{0:23}{1}<br>'.format('INTERATOMIC VECTOR','DMSDA')
                while i+14 < len(dmsda[2]):

                    resStr += '<br>{0:7} ---> {1:7}  {2:< 25.0f}{3:7} ---> {4:7}  {5: .0f}'.format(dmsda[2][i][0], dmsda[2][i][1], dmsda[2][i][2], dmsda[2][i+14][0], dmsda[2][i+14][1], dmsda[2][i+14][2])
                    i += 1

                while i < 14:
                    resStr += '<br>{0:7} ---> {1:7}  {2:< 25.0f}'.format(dmsda[2][i][0], dmsda[2][i][1], dmsda[2][i][2])
                    i+=1
            else:
                resStr += '<br><br>{0:23}{1:24}<br>'.format('INTERATOMIC VECTOR','DMSDA')
                while i < len(dmsda[2]):
                    resStr += '<br>{0:7} ---> {1:7}  {2:< 25.0f}'.format(dmsda[2][i][0], dmsda[2][i][1], dmsda[2][i][2])
                    i += 1
            resStr = '<pre>' + resStr +'</pre>'



        #If full results can't be printed, try to print without DMSDA and if that doesn't work print error message
        except Exception as e:
            print(e)
            try:
                resStr = 'RF<sup>2</sup> = {0}<br>{1}<br>No DMSDA results found.'.format(getRF2(lsmOutFile), getConvergence('xd_lsm.out'))
            except Exception:
                resStr = 'An error occurred.'

        return resStr


    def showResBackup(self):
        '''
        Show results from backup folder.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Choose backup folder to load results from"))

        resStr = ''

        if folder:
            resStr = self.makeResStr(folder + '/xd_lsm.out')

        self.showResLab.setText(resStr)
        self.lsmResLine.setVisible(True)

    def resBackupSum(self):
        '''
        Show wizard style summary of multiple refinements, from backup folder.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Choose backup folder to load results from"))
        results = []
        resStr = '<br>{0:47}{1:23}{2:14}{3:16}{4}<br>'.format(' Refinement', 'RF<sup>2</sup>', 'Convergence', 'Average DMSDA', 'Max DMSDA')

        for item in os.listdir(folder):
            c = ''
            rf2 = getRF2(folder + '/' + item + '/xd_lsm.out')
            print('got rf2')
            conv = getConvergence(folder + '/' + item + '/xd_lsm.out')
            print(folder + '/' + item + '/xd_lsm.out')
            if conv:
                c = 'Yes'
            else:
                c = 'No'
            try:
                dmsda = getDMSDA(folder + '/' + item + '/xd_lsm.out')
            except:
                dmsda = ['','']

            ref = item

            print(ref)
            refSplit = ref.split()
            refNum = ' '.join(refSplit[:2])
            refName = ' '.join(refSplit[2:])
            newRes = '<br>{0:>4} {1:40}{2:6.2f}{3:8}{4:14}{5:16}{6}'.format(refNum, refName, rf2, ' %', c, dmsda[0], dmsda[1])
            results.append(newRes)

        for res in sorted(results):
            resStr += res
        resStr = '<pre>' + resStr + '</pre>'
        self.showResLab.setText(resStr)
        self.lsmResLine.setVisible(True)


    def getResPress(self):
        '''
        Get results from xd_lsm.out and show them in results tab.
        '''
        resStr = self.makeResStr('xd_lsm.out')
        self.showResLab.setText(resStr)
        self.lsmResLine.setVisible(True)


    def FFT2FOUR(self):
        '''
        Run XDFOUR for atom nearest the largest peak in xd_fft.out.
        '''
        try:
            atom = FFTDetective()
            FOU3atoms(atom[0], copy.copy(globAtomLabs), copy.copy(globAtomPos))
            self.xdfour.finishedSignal.connect(self.showResDensMap)
            self.xdfour.start()
            self.xdfft.finishedSignal.disconnect(self.FFT2FOUR)
        except PermissionError:
            self.setupFOURStatusLab.setText(self.permErrorMsg)
        except FileNotFoundError:
            self.setupFOURStatusLab.setText(self.lstMissingErrorMsg)
        except Exception:
            self.setupFOURStatusLab.setText('An error occurred.')

    def showResDensMap(self):
        '''
        Show residual density map from xd_fou.grd.
        '''
        print('showresdensmap')
        try:
            self.xdfour.finishedSignal.disconnect(self.showResDensMap)
        except Exception:
            pass
        self.resmap = resmap('xd_fou.grd')
        self.resmap.show()

    def addFOURIns(self):
        '''
        Handle 'Run XDFOUR' button press in 'Results' tab.
        '''
        if self.setupFOURInputBox.isChecked() == True:
            rawInput = str(self.setupFOURInputText.text()).strip().upper()

            if rawInput:
                atom = rawInput2labels(rawInput)[0]

                try:
                    FOU3atoms(atom, copy.copy(globAtomLabs), copy.copy(globAtomPos))
                    self.xdfour.finishedSignal.connect(self.showResDensMap)
                    self.xdfour.start()
                except PermissionError:
                    self.setupFOURStatusLab.setText(self.permErrorMsg)
                except FileNotFoundError:
                    self.setupFOURStatusLab.setText(self.lstMissingErrorMsg)
                except KeyError as e:
                    self.setupFOURStatusLab.setText(self.atomNoExistErrorMsg)

            else:
                self.setupFOURStatusLab.setText('No atom selected.')

        elif self.setupFOURFFTBox.isChecked() == True:
            try:
                self.xdfft.start()
                self.xdfft.finishedSignal.connect(self.FFT2FOUR)

            except Exception as e:
                print(e)
                self.setupFOURStatusLab.setText('An error occurred.')

        elif self.quickplotGrdBox.isChecked():
            try:
                self.setupFOURStatusLab.setText('Making residual density map...')
                self.setupFOURStatusLab.repaint()
                QApplication.processEvents()
                self.showResDensMap()
                self.setupFOURStatusLab.setText('')
            except Exception as e:
                print(e)
                self.setupFOURStatusLab.setText('An error occurred.')


    def makeNPP(self):
        '''
        Handle user pressing 'Make normal probability plot' button.
        '''
        try:
            FOUcell()
            self.xdfour.finishedSignal.connect(self.showNPP)
            self.xdfour.start()
            self.resNPPLab.setText('Running XDFOUR...')
            self.resNPPBut.setText('Cancel')
            self.resNPPBut.clicked.connect(self.killXDFOUR)


        except Exception:
            self.resNPPLab.setText('An error occurred.')

    def showNPP(self):
        '''
        Show normal probability plot.
        '''
        try:
            self.xdfour.finishedSignal.disconnect(self.showNPP)
        except Exception:
            pass
        if not hasattr(self, 'npp'):
            self.npp = NPP()
            self.npp.show()
        elif not hasattr(self, 'npp2'):
            self.npp2 = NPP()
            self.npp2.show()
        else:
            self.npp3 = NPP()
            self.npp3.show()


    def getDpops(self):
        '''
        Run XDPROP with instructions to find d-orbital populations.
        '''
        try:
            print('getdpops')
            setupPROPDpops()
            print('setuppro')
            self.xdprop.finishedSignal.connect(self.showDpops)
            self.xdprop.start()

        except PermissionError:
            self.getDpopsStatusLab.setText(self.permErrorMsg)
        except Exception:
            self.getDpopsStatusLab.setText('An error occurred.')

    def showDpops(self):
        '''
        Show d-orbital populations from xd_pro.out.
        '''
        try:
            self.xdprop.finishedSignal.disconnect(self.showDpops)
        except Exception:
            pass
        try:
            dpopList = getDorbs()
            print(dpopList)
            dpopStr = 'd-orbital populations\n\n'
            print(dpopStr)
            for item in dpopList:
                dpopStr += item[0].rjust(7) + ':   ' + item[1] + ' %\n'

            if dpopStr:
                self.getDpopsStatusLab.setText(dpopStr)
            else:
                self.getDpopsStatusLab.setText('d-orbital populations not found in xd_pro.out.\nSetup xd.mas to get d-orbital populations with button above, run XDPROP and try again.')

        except:
            self.getDpopsStatusLab.setText('xd_pro.out missing or corrupted. Try running XDPROP again.')


    def autoTopxd(self):

        invalidLabels = []
        trueAtomList = getAtomList()
        if self.autoTopxdAll.isChecked():
            atomList = trueAtomList

        else:
            atomList = rawInput2labels(str(self.autoTopxdInput.text()))
            if atomList:
                for item in atomList:
                    if item not in trueAtomList:
                        invalidLabels.append(item)
                if invalidLabels:
                    self.autoTopxdLab.setText('Following atoms not in structure: {}'.format(listjoin(invalidLabels, ', ')))
                    return
            else:
                self.autoTopxdLab.setText('No atoms chosen.')
                return

        phiInput = str(self.autoTopxdPhi.text())
        thetaInput = str(self.autoTopxdTheta.text())

        if not phiInput.isdigit() or not thetaInput.isdigit():
            self.autoTopxdLab.setText('Invalid input given for phi/theta.')
            return

        estRunningTime = datetime.timedelta(minutes = (((int(phiInput)*int(thetaInput))/768)*30)*len(atomList))

        rawTimeStr = str(estRunningTime)

        days = 0
        hours = 0
        minutes = 0
        if '.' in rawTimeStr:
            rawTimeStr = rawTimeStr.split('.')[0]
        if ',' in rawTimeStr:
            rawTimeStrSplit = rawTimeStr.split(',')
            days = rawTimeStrSplit[0].split(' ')[0]
            time = rawTimeStrSplit[1].strip().split(':')

        else:
            time = rawTimeStr.strip().split(':')
        hours = time[0].lstrip('0')
        minutes = time[1].lstrip('0')
        print(time)

        estRunningTimeMsg = 'Estimated running time: '
        if days:
            if days=='1':
                estRunningTimeMsg += '{} day, '.format(days)
            else:
                estRunningTimeMsg += '{} days, '.format(days)
        if hours:
            if hours=='1':
                estRunningTimeMsg += '{} hour, '.format(hours)
            else:
                estRunningTimeMsg += '{} hours, '.format(hours)
        if minutes:
            if minutes=='1':
                estRunningTimeMsg += '{} minute, '.format(minutes)
            else:
                estRunningTimeMsg += '{} minutes, '.format(minutes)

        estRunningTimeMsg = estRunningTimeMsg.strip(', ')
        estRunningTimeMsg += '\n\nDo you want to start?'

        if not os.path.isdir('topxd'):
            os.makedirs('topxd')


        warningMsg = QMessageBox.question(self, 'Auto TOPXD', estRunningTimeMsg, QMessageBox.Yes, QMessageBox.No)

        if warningMsg == QMessageBox.Yes:
            self.autoTopxd = AutoTOPXD(atomList, phiInput, thetaInput)
            self.autoTopxd.finishedSignal.connect(self.autoTopxdFinished)
            self.autoTopxd.terminatedSignal.connect(self.autoTopxdTerminated)
            self.autoTopxd.show()
            self.autoTopxd.run()
        else:
            pass

    def autoTopxdInputHandler(self):
        '''
        Disables Auto TOPXD "All atoms" checkbox is there is input of specific atoms.
        '''
        rawInput = str(self.autoTopxdInput.text())
        if not rawInput:
            self.autoTopxdAll.setEnabled(True)
        else:
            self.autoTopxdAll.setEnabled(False)

    def autoTopxdFinished(self):
        '''
        When auto TOPXD finishes create csv file of Q and L from result files.
        '''
        with open('topxd/topxd_results.csv','w') as topxdcsv:
            topxdcsv = csv.writer(topxdcsv)
            topxdcsv.writerow(['Atom', 'Q', 'L'])

            for file in os.listdir('topxd'):
                splitFile = os.path.splitext(file)
                if splitFile[1] == '.out' and splitFile[0].startswith('topxd_'):
                    atom = splitFile[0].split('_')[1]
                    newRow = [atom]
                    with open('topxd/' + file, 'r', encoding='utf-8', errors = 'ignore') as topxdout:

                        for line in topxdout:
                            if line.startswith('              Q'):
                                newRow.append(float(line.split()[1]))

                            elif line.startswith('              L'):
                                newRow.append(float(line.split()[1]))
                                break
                        topxdcsv.writerow(newRow)
        self.autoTopxdLab.setText('csv file of atom Q and L values saved to "<i>{}</i>"'.format(os.path.abspath('topxd/topxd_results.csv')))


    def autoTopxdTerminated(self):
        pass

#########################################################################
#--------------------BACKUP----------------------------------------------
#########################################################################

    def genBackupPress(self,folder,label):
        '''
        Handle user backing up from different parts of GUI.
        '''
        nameAllowed = True
        for item in self.forbiddenChars:
            if item in folder:
                nameAllowed = False

        if nameAllowed == True and folder:
            try:
                backup(folder)
                label.setText(self.backupConfirmStr + '"' + os.getcwd().replace('\\','/') + '/Backup/' + folder + '"')
            except PermissionError:
                label.setText(self.permErrorMsg)
            except Exception:
                label.setText('An error occurred.')
        else:
            label.setText(self.invalidFolderNameStr)


    def loadBackupPress(self):
        '''
        Handle user pressing 'Load backup'.
        '''
        folderStr = os.getcwd() + '/Backup'

        if os.path.isdir(folderStr):
            folder = str(QFileDialog.getExistingDirectory(None, "Select Backup Folder",folderStr))
        else:
            folder = str(QFileDialog.getExistingDirectory(None, "Select Backup Folder",os.getcwd()))
        if folder:
            try:
                loadBackup(folder)
                self.loadBackupLab.setText('Backup loaded from: <i>"' + folder + '"</i>')
            except PermissionError:
                self.loadBackupLab.setText(self.permErrorMsg)
            except Exception:
                self.loadBackupLab.setText('An error occurred.')
        else:
            pass

        self.changeUserIns()


#########################################################################
#--------------------MISC------------------------------------------------
#########################################################################

    def openAbout(self):
        '''
        Open About XD Toolkit window.
        '''
        self.about = aboutBox()
        self.about.show()


    def openManual(self):

        if sys.platform=='win32':
            os.startfile(manualAbsPath)
        else:
            subprocess.call(['xdg-open', manualAbsPath])

    def updateHelpText(self):
        tab = self.tabWidget.currentIndex()
        helpText = ''
        if tab == 5:
            toolbox = self.toolsToolbox         
            try:
                helpText = self.helpTexts['{}{}'.format(tab, toolbox.currentIndex())]
            except KeyError:
                helpText = ''
        elif tab == 0:
            helpText = self.helpTexts['1']
        else:
            self.helpLabel.setText('')
        
        self.helpLabel.setText(helpText)
    
    def closeEvent(self, event):
        '''
        When program closes, kill XDLSM, save current working directory to settings and send lsmTimes.
        '''
        global timeFileAbsPath

        self.settings.setValue('lastcwd', os.getcwd())

#        if self.settings.value('senddata') == 'yes':
#
#            with open(timeFileAbsPath,'r') as lsmTimes:
#                x = len(lsmTimes.readlines())
#
#            if x > 100:
#                try:
#                    sendEmail(subject = 'XD Toolkit: XDLSM times', attachments = [timeFileAbsPath])
#                    with open(timeFileAbsPath, 'w') as lsmTimes:
#                        lsmTimes.write(' ')
#                except Exception:
#                    pass
        try:
            self.xdlsm.xdlsmRunning.terminate()
        except Exception:
            pass
        event.accept()


    def resetLabels(self):
        '''
        Clear all labels in GUI.
        '''
        for label in self.labList:
            label.setText('')
        self.lsmResLine.setVisible(False)
        try:
            self.nppCanvas.setParent(None)
        except Exception:
            pass

def customExceptHook(Type, value, traceback):
    print('Exception:')
    #sendEmail('<br><br>'.join(format_exception(Type, value, traceback)),toaddr = 'mcrav@chem.au.dk')
    print_exception(Type, value, traceback)
    pass

##Run GUI
if __name__ == '__main__':
    if sys.platform.startswith('win'):
        # On Windows calling this function is necessary.
        mp.freeze_support()
    print(os.getcwd())
    sys.excepthook = customExceptHook               #Accept any errors so GUI doesn't quit.
    app = QApplication(sys.argv)

    #Splash screen
    splash_pix = QPixmap('res/splash.png')
    splash = QSplashScreen(splash_pix)
    splash.setMask(splash_pix.mask())

    splash.setFont(getProgramFont())
    splash.showMessage('Initializing...',
                           Qt.AlignBottom | Qt.AlignLeft,
                           Qt.white)
    splash.show()
    app.processEvents()
    prog = XDToolGui()
    splash.finish(prog)
    app.aboutToQuit.connect(app.deleteLater)  #Fixes Anaconda bug where program only works on every second launch

    prog.show()

    sys.exit(app.exec_())
##
#import timeit
#os.chdir('/home/matt/dev/XDTstuff/test/carba')
#def wrapper(func, *args, **kwargs):
#     def wrapped():
#         return func(*args, **kwargs)
#     return wrapped
#wrapped = wrapper(ins2all)
#print(timeit.timeit(wrapped, number=1))
#print('\n\n~~~~~~~~~~~~~~\n\n')
#x = ins2all()
#from initfuncs import carbaTest
#if str(x) != carbaTest:
#    print('problems')
#else:
#    print('fine')

#x = ins2all()
#findSITESYM(trackAtom = 'C(01A)')
#x = ins2all()
#findCHEMCON()
#multipoleMagician()
#x = makeLCS()
#print(x)
