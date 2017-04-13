from quickplot import makeResMap
from asym2unit import applySymOps
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.stats import probplot
import matplotlib.pyplot as plt
from ast import literal_eval
import random
import copy
from traceback import print_exception
import subprocess
import hashlib
import time
import itertools
import numpy as np
import os
import smtplib
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import sys
import webbrowser
from shutil import copyfile, rmtree
from collections import Counter
from xdcool import Ui_MainWindow
from pref import Ui_pref
from about import Ui_aboutBox
from wizard import Ui_wizard
from resmap import Ui_resmap
from sendBug import Ui_sendBug
from sendSugg import Ui_sendSugg
from splash import Ui_splash
from PyQt5.QtWidgets import QWidget, QMessageBox, QLabel, QGridLayout, QDialogButtonBox, QSplashScreen, QPushButton, QApplication, QDialog, QLineEdit, QFileDialog, QMainWindow
from PyQt5.QtCore import QCoreApplication, QSettings, QThread, pyqtSignal, Qt
from PyQt5.QtGui import QPixmap, QFont

'''
#-------------------DEV TOOLS----------------------------------------  
'''
 
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
        tzero = time.time()
        rtn = f(*args, **kwargs)
        tfin = time.time()
        print('{0:40}{1:.7f} s'.format(f.__name__, tfin-tzero))
        return rtn
    return timeFunc
	
'''	
#-------------------UTILITIES----------------------------------------
'''

def lab2type(atomLabel):
    '''
    Convert atom label to atom type, e.g. 'C(2)' to 'C'. Return type.
    '''
    return atomLabel.split('(')[0].upper()

def spec2norm(atomLabel):
    '''
    Convert special atom label to normal atom label, e.g. 'C(2),asym' to 'C(2)'. Return label.
    '''
    return atomLabel.split(',')[0]

def rawInput2Labels(rawInput):
    '''
    Convert raw input of unspecified format atom labels to list of correctly formatted labels. Return list.
    '''
    return formatLabels(labels2list(rawInput))

def labels2list(inputText):
    '''
    Convert raw input of atom labels to list of unformatted atom labels. Return list.
    '''
    inputText = inputText.upper()
    if ',' in inputText:
        inputText = inputText.replace(' ','').strip(',')
        inputAtomList = inputText.split(',')
    elif ' ' in inputText:
        inputAtomList = inputText.split()
    else:
        inputAtomList = [inputText]
        
    return inputAtomList

def formatLabels(inputAtomList):
    '''
    Convert labels c1 to c(1).
    '''
    elements = findElements()
    inputNewAtomList = []
    
    for atomLab in inputAtomList:
        newAtomLab = ''
        if '(' not in atomLab:
            if atomLab[:2] in elements:
                newAtomLab = '{0}({1})'.format(atomLab[:2], atomLab[2:])
            elif atomLab[:1] in elements:
                newAtomLab = '{0}({1})'.format(atomLab[:1], atomLab[1:])
            elif atomLab[:3] == 'DUM':
                newAtomLab = atomLab
            
            if newAtomLab:
                inputNewAtomList.append(newAtomLab)

        else:
            inputNewAtomList.append(atomLab)
    
    return inputNewAtomList



def sendEmail(body = '', email = '', attachments = [], subject = ''):
    '''
    Send email to my email account from xdtoolkit@gmail.com
    '''
    fromaddr = "xdtoolkit@gmail.com"
    toaddr = "mcrav@chem.au.dk"
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = subject
    
    for f in attachments:
                
        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=os.path.basename(f)
            )
            part['Content-Disposition'] = 'attachment; filename="%s"' % os.path.basename(f)
            msg.attach(part)
     
    bodyText = '{0}<br><br>Email Address: {1}'.format(body,email)
    msg.attach(MIMEText(bodyText, 'html'))
    
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login("xdtoolkit@gmail.com", '***')
    text = msg.as_string()
    server.sendmail("xdtoolkit@gmail.com", "mcrav@chem.au.dk", text)
    
    server.quit()
    

def isfloat(x):
    '''
    Check if unknown value can be converted to a float. Return result as bool.
    '''
    floatBool = True
    
    try:
        float(x)
    except Exception:
        floatBool = False
        
    return floatBool
    

def fixLsmCif():
    '''
    Add unit cell paramaters to xd_lsm.cif.
    '''
    cellParams = getCellParams()
    
    with open('xd_lsm.cif','r') as cif, open('xd_lsmnew.cif','w') as newcif:
        
        for line in cif:
            
            if line.startswith('# Refinement details'):
                
                cellStr = '\n#Unit Cell Parameters\n\n'
                cellStr += '{0:35}{1}\n'.format('_cell_length_a', cellParams[0])
                cellStr += '{0:35}{1}\n'.format('_cell_length_b', cellParams[1])
                cellStr += '{0:35}{1}\n'.format('_cell_length_c', cellParams[2])
                cellStr += '{0:35}{1}\n'.format('_cell_angle_alpha', cellParams[3])
                cellStr += '{0:35}{1}\n'.format('_cell_angle_beta', cellParams[4])
                cellStr += '{0:35}{1}\n'.format('_cell_angle_gamma', cellParams[5])
                cellStr += '\n\n'
                newcif.write(cellStr)
                
            newcif.write(line)
            
    os.remove('xd_lsm.cif')
    os.rename('xd_lsmnew.cif','xd_lsm.cif')
                
def getCellParams():
    '''
    Get unit cell parameters from shelx.ins or xd.mas.
    '''
    if os.path.isfile('shelx.ins'):
        with open('shelx.ins','r') as ins:
            ins = open('shelx.ins','r')
            
            for line in ins:
                    
                if line.startswith('CELL'):
                    row = str.split(line)
                    a = float(row[2])
                    b = float(row[3])
                    c = float(row[4])
                    alpha = float(row[5])
                    beta = float(row[6])
                    gamma = float(row[7])
            
    elif os.path.isfile('xd.mas'):
        with open('xd.mas','r') as mas:
            
            for line in mas:
                
                if line.startswith('CELL '):
                
                    row = str.split(line)
                    
                    a = float(row[1])
                    b = float(row[2])
                    c = float(row[3])
                    alpha = float(row[4])
                    beta = float(row[5])
                    gamma = float(row[6])

    return [a, b, c, alpha, beta, gamma]        
    

def findElements():
    '''
    Find elements in compound from shelx.ins or xd.inp. Return list of elements.
    '''
    elements = []
    
    if os.path.isfile('shelx.ins'):
        with open('shelx.ins','r') as ins:
            for line in ins:
                if line.startswith('SFAC'):
                    row = str.split(line)
                    elements = [item.upper() for item in row[1:]]
                    break
                
    elif os.path.isfile('xd.inp'):
        with open('xd.inp','r') as inp:
            
            for line in inp:
                if line.startswith('END SCAT'):
                    scatTab = False
                    break
                    
                elif scatTab:
                    row = str.split(line)
                    if row[0].isalpha():
                        elements.append(row[0])
                        
                elif line.startswith('SCAT '):
                    scatTab = True
                
    return elements


def addSnlCutoff(snlmin = 0.0, snlmax = 2.0):
    '''
    Add sin(theta/lambda) cutoffs to xd.mas.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        for line in mas:
            
            if line.startswith('SKIP'):
                row = str.split(line)
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOn}  {snlMin:<5.3f} {snlMax:<5.3f}'.format(*row, snlOn = '*sinthl', snlMin = snlmin, snlMax = snlmax)
                newmas.write(rowStr + '\n')
            
            else:
                newmas.write(line)
                
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def removePhantomAtoms():
    '''
    Remove faulty atoms with format '(e' from xd.mas and xd.inp.
    '''
    phantoms = False
    
    with open('xd.inp','r') as inp, open('xdnew.inp','w') as newinp:
        i = 3
        phantomCount = 0
        
        for line in inp:
            if line.startswith('('):
                phantoms = True
                phantomCount += 1
                i = 0
            if i > 2:
                newinp.write(line)
            i+=1

   
    
    if phantoms:
        with open('xd.mas', 'r') as mas, open('xdnew.mas', 'w') as newmas:
            
            
            atomTab = False
            keyTab = False
            scatTab = False
            kappaPrinted = False
            k = 0
            
            for line in mas:
                if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                    atomTab = False
                elif line.startswith('KAPPA'):
                    keyTab = False
                    if kappaPrinted == False:
                        if k < phantomCount:
                            k += 1
                            continue
                            
                elif line.startswith('END SCAT'):
                    scatTab = False
                    
                if atomTab:
                    if line.startswith('('):
                        continue
                    else:
                        newmas.write(line)
                    
     
                elif scatTab:
                    if line.startswith('('):
                        continue
                    else:
                        newmas.write(line)
               
                elif keyTab:
                    if line.startswith('('):
                        continue
                    else:
                        newmas.write(line)
                    
                else:
                    newmas.write(line)
                    
                if line.startswith('ATOM     ATOM0'):
                    atomTab = True
                    
                elif line.startswith('KEY'):
                    keyTab = True
                
                elif line.startswith('SCAT'):
                    scatTab = True

        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
      
        numElements = 0

        with open('xdnew.inp','r') as inp, open('xd.inp','w') as newinp:
            i = 0
            j = 1
            for line in inp:
                if line.startswith('USAGE'):
                    row = str.split(line)
                    numAtoms = int(row[1])
                    numElements = int(row[4])
                    row[1] = str(numAtoms - phantomCount)
                    row[4] = str(numElements - phantomCount)
                    rowStr = '{0:10}{1:6}{2:3}{3:4}{4:4}{5:4}{6:4}{7:4}{8:4}{9:4}{10:4}{11:6}{12:3}{13:4}{14}\n'.format(*row)
                    newinp.write(rowStr)
                else:                                               #Code to remove extra line(s) of kappa table in inp
                    if j <= (numElements - phantomCount) and i > 2:
                        newinp.write(line)

                    elif i <= 2:
                        newinp.write(line)
                        
                    else:
                        i = -1
                    
                    
                if line[:1] == ' ':
                    i+= 1
                else:
                    i = 0
                    
                if i > 2:
                    j+=1

        
        os.remove('xdnew.inp')

    

def fixBrokenLabels():
    '''
    Change problematic atom labels in shelx.ins.
    '''
    try:  
        neebs = copy.copy(globAtomLabs)
        i = 0 
        atomBool = False
        elements = []
            
        with open('shelx.ins','r') as ins, open('shelxnew.ins','w') as newins:
            
            for line in ins:
                i = 0
                if line.startswith('HKLF'):
                    atomBool = False
                
                elif line.startswith('SFAC'):
                    row = str.split(line)
                    elements = [item.upper() for item in row[1:]]                
                
                if atomBool:
                    if line[:1].isalpha() and not line.startswith('AFIX'):
                        row = str.split(line)
                        atomLab = row[0]
                        newLab = atomLab.replace('(','').replace(')','')
                        #Fix labels that are just N or C or Co with no number
                        if newLab in elements:       #N Na error
                            newerLab = newLab + str(i)
                            while newerLab in neebs:
                                i += 1
                                newerLab = newLab + str(i)
                            newLab = newerLab
                        newLine = '{0:6}{1}'.format(newLab.upper(), line[6:])
                        newins.write(newLine)
                    else:
                        newins.write(line)
                
                else:
                    newins.write(line)
                    
                if line.startswith('FVAR'):
                    atomBool = True
        
        os.remove('shelx.ins')
        os.rename('shelxnew.ins', 'shelx.ins')
    
    finally:
        try:
            os.remove('shelxnew.ins')
        except Exception:
            pass
    

def getNumAtoms():
    '''
    Get total number of atoms in compound from shelx.ins or xd.inp.
    '''
    i = 0
    
    if os.path.isfile('xd.inp'):
        with open('xd.inp','r') as inp:
            
            for line in inp:
                
                if line.startswith('USAGE'):
                    row = str.split(line)
                    i = int(row[1])
                    break
                
    elif os.path.isfile('shelx.ins'):
        atomBool = False
        with open('shelx.ins','r') as ins:

            for line in ins:
                if line.startswith('HKLF') or line.startswith('REM'):
                    atomBool = False
        
                if atomBool and line[:1].isalpha() and not line.startswith('AFIX'):
                    i+=1
        
                if line.startswith('FVAR'):
                    atomBool = True
    return (int(i))

        
def totalEstTime():
    '''
    Estimate total estimated time for XD wizard to run. Return this value.
    '''
    numAtoms = getNumAtoms()
    s = 0.066*numAtoms +1.944
    ha = 0.342*numAtoms -3.861
    la = 0.050*numAtoms +1.348
    mk = 0.351*numAtoms -3.016
    m = 1.052*numAtoms -19.025
    nhpam = 0.882*numAtoms -6.180
 
    totalEstTime = s + ha + la + mk + m + nhpam + 5
 
    return totalEstTime


def check4errors():
    '''
    Check for errors in xd_lsm.out. Return errors.
    '''
    try:
        lsm = open('xd_lsm.out', 'r')
        complete = False
        error = ''
        
        for line in lsm.readlines()[-20:]:

            if line.startswith(' | Total program runtime'):
                complete = True
            elif line.strip().startswith('* * * '):
                error = line
                complete = False
            elif line.startswith('Noble gas configuration not recognized for element Cu'):
                error = line
                complete = False
    finally:
        lsm.close()
    return(complete,error)   


def ins2fracPos(insFile):
    '''
    Get fractonal coordinates of every atom in shelx.ins. Return these and unit cell parameters.
    '''
    ins = open(insFile,'r')
    atomBool = False
    atomPos = {}
    elements = []
    
    atomLab = ''
    
    for line in ins:
        
        if line.startswith('HKLF') or line.startswith('REM'):
            atomBool = False
            
        elif line.startswith('SFAC'):
            row = str.split(line)
            elements = [item.upper() for item in row[1:]]
            
        elif line.startswith('CELL'):
            row = str.split(line)
            a = float(row[2])
            b = float(row[3])
            c = float(row[4])
            alpha = float(row[5])
            beta = float(row[6])
            gamma = float(row[7])

        if atomBool and line[:1].isalpha() and not line.startswith('AFIX'):
            row = str.split(line)
            if row[0][:2] in elements:
                if len(row[0]) > 2:
                    atomLab = '{0}({1})'.format(row[0][:2],row[0][2:])
                else:
                    atomLab = row[0]
            elif row[0][:1] in elements:
                if len(row[0]) > 1:
                    atomLab = '{0}({1})'.format(row[0][:1],row[0][1:])
                else:
                    atomLab = row[0]
                
            atomPos[atomLab.upper()] = [np.array([float(row[2]), float(row[3]), float(row[4])])]

        if line.startswith('FVAR'):
            atomBool = True
            
    ins.close()

    return (atomPos, a, b, c, alpha, beta, gamma)


def inp2fracPos():
    '''
    Get fractional coordinates from xd.inp. Return these and unit cell parameters.
    '''
    with open('xd.inp','r') as inp, open('xd.mas','r') as mas:
        scatTab = False
        atomPos = {}
        elements = []
        
        for line in mas:
            
            if line.startswith('END SCAT'):
                scatTab = False
                break
            
            elif scatTab:
                row = str.split(line)
                if row[0].isalpha():
                    elements.append(row[0])
                    
            elif line.startswith('CELL '):
                
                row = str.split(line)
                
                a = float(row[1])
                b = float(row[2])
                c = float(row[3])
                alpha = float(row[4])
                beta = float(row[5])
                gamma = float(row[6])
        
            elif line.startswith('SCAT'):
                scatTab = True
            
        for line in inp:
            row = str.split(line)
            
            if '(' in row[0]:
                    
                atomPos[row[0].upper()] = [(np.array([float(row[12]), float(row[13]), float(row[14])]), row[0].upper())]   
    
    return (atomPos, a, b, c, alpha, beta, gamma)


def getBondDist(atom1c, atom2c, a, b, c, alpha, beta, gamma):
    '''
    Get distance between 2 atoms. Return distance.
    '''
    delta1 = a*(atom1c[0] - atom2c[0])
    delta2 = b*(atom1c[1] - atom2c[1])
    delta3 = c*(atom1c[2] - atom2c[2])

    d2 = delta1**2 + delta2**2 + delta3**2 + (2*delta1*delta2*np.cos(gamma)) + (2*delta1*delta3*np.cos(beta)) + (2*delta2*delta3*np.cos(alpha))
    
    bondDist = np.sqrt(d2)
    
    return bondDist
    

@timeDec
def ins2all():
    '''
    Find information about compound from shelx.ins input and calculations.
    Return nearest neighbours, atomic positions, bond distances, and bond angles.
    '''
    if os.path.isfile('shelx.ins'):
        atomPosData = applySymOps(ins2fracPos('shelx.ins')[0])
        atomPos = atomPosData[0]
        specAtomPos = atomPosData[1]
        insInfo = ins2fracPos('shelx.ins')

    else:
        atomPos = insInfo[0]

    a = insInfo[1]
    b = insInfo[2]
    c = insInfo[3]
    alpha = np.radians(insInfo[4])
    beta = np.radians(insInfo[5])
    gamma = np.radians(insInfo[6])
    atomPairs = []
    neebPairs = []
    neebTypes = {}
    neebLabs = {}
    asymNeebs = {}          #Dict to store atoms in asym. unit and their neighbours with special labels
    neebSpecialPairs = []
    neebSpecLabs = {}
    distances = {}
    specDistances = {}
    specAngles = {}
    oneCoordLst = ('H(','F(','CL','I(','BR')
    covradii = {'LU': 1.56, 'PO': 1.46, 'YB': 1.74, 'RU': 1.25, 'U': 1.42, 'TI': 1.32, 'GD': 1.61, 'SE': 1.16, 'SN': 1.41, 'P': 1.06, 'IN': 1.44, 'TC': 1.27, 'K': 2.03, 'I': 1.33, 'PR': 1.65, 'CL': 0.99, 'LA': 1.69, 'TB': 1.59, 'HG': 1.49, 'CO': 1.16, 'PT': 1.3, 'NI': 1.15, 'AL': 1.18, 'TL': 1.48, 'NB': 1.34, 'F': 0.72, 'HF': 1.44, 'SR': 1.91, 'MG': 1.36, 'Y': 1.62, 'PD': 1.28, 'TM': 1.56, 'O': 0.73, 'AT': 1.45, 'BR': 1.14, 'PB': 1.47, 'N': 0.75, 'SC': 1.44, 'SB': 1.41, 'NA': 1.54, 'BA': 1.98, 'FE': 1.17, 'MO': 1.3, 'GA': 1.26, 'V': 1.22, 'XE': 1.31, 'BE': 0.9, 'RH': 1.25, 'ZN': 1.25, 'CA': 1.74, 'SM': 1.62, 'RB': 2.16, 'NE': 0.71, 'CR': 1.18, 'OS': 1.26, 'BI': 1.46, 'ND': 1.64, 'LI': 1.23, 'C': 0.77, 'CE': 1.65, 'S': 1.02, 'KR': 1.12, 'ER': 1.57, 'CU': 1.17, 'B': 0.82, 'SI': 1.11, 'W': 1.3, 'TH': 1.65, 'CD': 1.48, 'IR': 1.27, 'H': 0.32, 'AG': 1.34, 'AS': 1.2, 'CS': 2.35, 'MN': 1.17, 'TE': 1.36, 'PM': 1.63, 'RE': 1.28, 'AR': 0.98, 'EU': 1.85, 'TA': 1.34, 'DY': 1.59, 'HE': 0.93, 'HO': 1.58, 'GE': 1.22, 'ZR': 1.45, 'AU': 1.34}
    
    #Make all possible pairs of atoms in atomPairs
    for atom in atomPos.keys():
        for atom2 in atomPos.keys():
            if atom != atom2 and {atom, atom2} not in atomPairs and not (atom[:2] in oneCoordLst and atom2[:2] in oneCoordLst):
                atomPairs.append({atom, atom2})
    i=0
    #Make list of pairs within bonding distances
    for pair in atomPairs:

        tupPair = tuple(pair)
        atom1 = tupPair[0]
        atom2 = tupPair[1]
        
        pos1SpecialLab = ''
        pos2SpecialLab = ''

        for pos in atomPos[atom1]:

            atom1c = np.array((pos[0][0], pos[0][1], pos[0][2]))
            pos1SpecialLab = pos[1]
            pos1InvSym = specAtomPos[pos1SpecialLab][1]
            
            for pos2 in atomPos[atom2]:

                atom2c = np.array((pos2[0][0], pos2[0][1], pos2[0][2]))
                pos2SpecialLab = pos2[1]
                
                bondDist = getBondDist(atom1c, atom2c, a, b, c, alpha, beta, gamma)
                    
                if bondDist < (covradii[atom1[:2].strip('(')] + covradii[atom2[:2].strip('(')] + 0.5):
                    distances[frozenset(tupPair)] = round(bondDist, 4)
                    specDistances[frozenset([pos1SpecialLab, pos2SpecialLab])] = round(bondDist, 4)
                    neebPairs.append(tupPair)   
                    neebSpecialPairs.append((pos1SpecialLab, pos2SpecialLab))
                    
                else: 
                    
                    #Get all x,y,z+-1 combos
                    for combo in getCombos(atom1c):
                        bondDist = getBondDist(combo[0], atom2c, a, b, c, alpha, beta, gamma)
                        
                        if bondDist < (covradii[atom1[:2].strip('(')] + covradii[atom2[:2].strip('(')] + 0.5):
                            tupPair = (atom1, atom2)
                            
                            combo1SpecialLab = atom1 + ',combo' + str(i)
                            i+=1
                            distances[frozenset(tupPair)] = round(bondDist, 4)
                            specDistances[frozenset([combo1SpecialLab, pos2SpecialLab])] = round(bondDist, 4)
                            
                            neebPairs.append(tupPair)   
                            neebSpecialPairs.append((combo1SpecialLab, pos2SpecialLab))
    
                            combo1InvSym = copy.copy(pos1InvSym)
                            combo1InvSym.insert(0, combo[1])
 
                            specAtomPos[combo1SpecialLab] = (combo[0], combo1InvSym)
                        
            
    #Get neighbours label dictionary       
    for pair in neebPairs:
        
        if pair[0] in neebLabs.keys():
            if not pair[1] in neebLabs[pair[0]]:
                neebLabs[pair[0]].append(pair[1])
        else:
            neebLabs[pair[0]] = [pair[1]]
            
        if pair[1] in neebLabs.keys():
            if not pair[0] in neebLabs[pair[1]]:
                neebLabs[pair[1]].append(pair[0])
        else:
            neebLabs[pair[1]] = [pair[0]]
    
    for atom, neighbours in neebLabs.items():
        neebLabs[atom] = sorted(neighbours)
        
    
    #Get neighbour label dictionary of actual neighbours with special labels
    for pair in neebSpecialPairs:

        if pair[0] in neebSpecLabs.keys():
            if not pair[1] in neebSpecLabs[pair[0]]:
                neebSpecLabs[pair[0]].append(pair[1])
        else:
            neebSpecLabs[pair[0]] = [pair[1]]
            
        if pair[1] in neebSpecLabs.keys():
            if not pair[0] in neebSpecLabs[pair[1]]:
                neebSpecLabs[pair[1]].append(pair[0])
        else:
            neebSpecLabs[pair[1]] = [pair[0]]
    
    for atom, neighbours in neebSpecLabs.items():
        neebSpecLabs[atom] = sorted(neighbours)
    
    asymNeebs = {}
    addedPosDict = {atom: [] for atom in neebLabs.keys()}
     
    #Generate all neighbours of parent atoms around that atom positions in asym unit
    for pair in neebSpecialPairs:
        splitLab = pair[0].split(',')
        tupPos = coords2tuple(specAtomPos[pair[1]][0])
        if splitLab[1] == 'asym':
            
            addedPos = addedPosDict[splitLab[0]]

            if tupPos not in addedPos:
                asymNeebs.setdefault(splitLab[0],[]).append((specAtomPos[pair[1]][0], pair[1]))
                addedPosDict[splitLab[0]].append(tupPos)

        else:          

            invSymOp = specAtomPos[pair[0]][1]              #Get inv symOp of asym atom
            pos = specAtomPos[pair[1]][0]                   #Position of atom not asym

            x = pos[0]
            y = pos[1]
            z = pos[2]
            for item in invSymOp:

                x = eval(item[0])                                  #Generate position in relation to asym atom
                y = eval(item[1])
                z = eval(item[2])
            newPos = np.array([x,y,z])
            tupPos = coords2tuple(newPos)
            
            addedPos = addedPosDict[splitLab[0]]
            if tupPos not in addedPos:

                newLab = pair[1] +'.dum' + str(i)
                asymNeebs.setdefault(splitLab[0], []).append((newPos, newLab))     #(pos, specLabel)
                specAtomPos[newLab] = (newPos, invSymOp)

                specDistances[frozenset([newLab, splitLab[0] + ',asym'])] = round(getBondDist(newPos, specAtomPos[splitLab[0] + ',asym'][0], a,b,c,alpha,beta,gamma), 4)
                addedPosDict[splitLab[0]].append(tupPos)
                i+=1
            
        splitLab = pair[1].split(',')
        tupPos = coords2tuple(specAtomPos[pair[0]][0])
        if splitLab[1] == 'asym':
            
            addedPos = addedPosDict[splitLab[0]]

            if tupPos not in addedPos:
                asymNeebs.setdefault(splitLab[0],[]).append((specAtomPos[pair[0]][0], pair[0]))
                addedPosDict[splitLab[0]].append(tupPos)
                
            
        else:                                                               #Atom not in asym so need to transform it
            invSymOp = specAtomPos[pair[1]][1]             #Get inv symOp of asym atom
            pos = specAtomPos[pair[0]][0]                   #Position of atom not asym

            x = pos[0]
            y = pos[1]
            z = pos[2]
            for item in invSymOp:
                x = eval(item[0])                                  #Generate position in relation to asym atom
                y = eval(item[1])
                z = eval(item[2])
            newPos = np.array([x,y,z])
            tupPos = coords2tuple(newPos)

            addedPos = addedPosDict[splitLab[0]]
            if tupPos not in addedPos:
                newLab = pair[0] + '.dum' + str(i)
                asymNeebs.setdefault(splitLab[0], []).append((newPos, newLab))     #(pos, specLabel)
                specAtomPos[newLab] = (newPos, invSymOp)

                specDistances[frozenset([newLab, splitLab[0] + ',asym'])] = getBondDist(newPos, specAtomPos[splitLab[0] + ',asym'][0],a,b,c,alpha,beta,gamma)
                addedPosDict[splitLab[0]].append(tupPos)
                i+=1

    #Get neighbour type dictionary
    for atom, neebs in asymNeebs.items():
        neebTypes[atom] = []
        for neeb in neebs:
            neebTypes[atom].append(neeb[1][:2].split(',')[0].strip('('))

    metricMatrix = np.array([[a**2, a*b*np.cos(gamma), a*c*np.cos(beta)],
                            [b*a*np.cos(gamma), b**2, b*c*np.cos(alpha)],
                            [c*a*np.cos(beta), c*b*np.cos(alpha), c**2]]) 
    

    #Work out angles based on positions generated by symmetry
    for atom, neighbours in asymNeebs.items():
        
        atom2 = atom
        neighbourLabs = [neighbour[1] for neighbour in neighbours]
        
        specAngles[atom2] = []
        pairs  = itertools.combinations(neighbourLabs, 2)
        
        
        for pair in pairs:

            atomsInAngle = (pair[0],atom2 + ',asym',pair[1])
            angle = atoms2angle(atomsInAngle, specAtomPos, specDistances, metricMatrix)

            specAngles[atom2].append((round(angle,2), atomsInAngle))
               

    #Remove atom positions from atomLab dict
    atomLabs = {atom : [neeb[1] for neeb in neighbours] for atom, neighbours in asymNeebs.items()}
    
    return (atomLabs, neebTypes, specAngles, specDistances, specAtomPos)


def getCombos(coord):
    '''
    Create all combinations of x, x+1, x-1 for x,y,z. Return combinations.
    '''
    combos = []
    
    xminus = coord[0] - 1
    yminus = coord[1] - 1
    zminus = coord[2] - 1
    xplus = coord[0] + 1
    yplus = coord[1] + 1
    zplus = coord[2] + 1
    x = coord[0]
    y = coord[1]
    z = coord[2]
    
    xs = [x, xminus, xplus]
    ys = [y, yminus, yplus]
    zs = [z, zminus, zplus]
    
    xSym = ['x','x+1','x-1']
    ySym = ['y','y+1','y-1']
    zSym = ['z','z+1','z-1']
    
    for i,xval in enumerate(xs):
        for j,yval in enumerate(ys):
            for k,zval in enumerate(zs):
                combos.append((np.array([xval, yval, zval]), (xSym[i], ySym[j], zSym[k])))

    return(combos)


def coords2tuple(coords):
    '''
    Convert numpy array of 3d coordinates to tuple with all numbers to 3 decimal places.
    Return tuple.
    '''
    tupPos = []

    for item in coords:
        stritem = '{0:.3f}'.format(item).strip()    #Round all numbers to 3 dp and make a string
        if item == 0:
            if stritem[0] == '-':           #Change -0. to 0.0
                stritem = stritem[1:]

        tupPos.append(stritem)

    return tuple(tupPos)
            

def atoms2angle(atoms, atomPos, distances, metricMatrix):
    '''
    Find angle between 3 atoms. Return angle.
    '''
    atom1c = atomPos[atoms[0]][0]
    atom2c = atomPos[atoms[1]][0]
    atom3c = atomPos[atoms[2]][0]
    
    r =  distances[frozenset((atoms[1],atoms[0]))]
    s =  distances[frozenset((atoms[1],atoms[2]))]

    X1 = np.array(atom2c-atom1c)
    X2 = np.array(atom2c-atom3c)
    
    cosphi = (np.dot(np.dot(np.transpose(X1),metricMatrix),X2)) / (r*s)

    angle = np.arccos(cosphi)
    
    degAngle = np.degrees(angle)
    roundAngle = round(degAngle, 2)
    return roundAngle       


def getmd5Hash(file):
    '''
    Create md5 hash of shelx.hkl file. Return hash as string.
    '''
    hklHash = ''
    
    if os.path.isfile(file):
        with open(file,'r') as hkl:
            hklTxt = hkl.read()
            hklHash = hashlib.md5(bytes(hklTxt, 'utf-8')).hexdigest()
            
    return hklHash


def inCache(file):
    '''
    Find out if project is in cache. Return result as bool and md5 hash of shelx.hkl.
    '''
    global cachePath
    inCache = False
    
    hklHash = getmd5Hash(file)
    
    if os.path.isdir('{}/{}'.format(cachePath, hklHash)):
        inCache = True
    
    return (inCache, hklHash)


def clearCache():
    '''
    Delete everything from cache.
    '''
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
                raise Exception
            else:
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
                

def res2inp():
    '''
    Rename xd.res to xd.inp.
    '''
    os.remove('xd.inp')
    os.rename('xd.res','xd.inp')
	
  
def backup(folderName):
    '''
    Copy refinement files to given folder.
    '''
    folder = 'Backup/' + folderName 
    if not os.path.isdir('{}{}'.format(os.getcwd(), '/Backup')):                   #Check to see if backup folder exists
        os.makedirs('Backup/')                          #If it doesn't exist, make it
        
                        #Create string of new folder path
    
    if not os.path.isdir('{}{}{}'.format(os.getcwd(), '/',folder)):                      #Check if new folder exists
        os.makedirs(folder)                             #If it doesn't exist, make it
    
    try:                                                #Copy all files to backup folder
        copyfile('xd.mas',folder + '/xd.mas')           #try and except used in case one of the files
    except:                                             #doesn't exist i.e. xd.res
        pass
    try:
        copyfile('xd.inp',folder + '/xd.inp')
    except:
        pass
    try:
        copyfile('xd.res',folder + '/xd.res')
    except:
        pass
    try:
        copyfile('xd_lsm.out',folder + '/xd_lsm.out')
    except:
        pass 
        

def initialiseMas():
    '''
    Fix 'noble gas configuration for CU' error.
    '''
    mas = open('xd.mas','r')
    newmas = open('xdnew.mas','w')
    
    scatTab = False
    
    for line in mas:
        
        if line.startswith('END SCAT'):
            scatTab = False
            
        if scatTab:
            if line[:2].upper() == 'CU':
                row = str.split(line)
                row[11] = '-10'

                i = 0
                for item in row:
                    if 3 < i < 25:
                        row[i] = float(item)
                    i += 1  
                rowStr = '{0:5}{1:5}{2:5}{3:7}{4:< 4.0f}{5:< 4.0f}{6:< 4.0f}{7:< 4.0f}{8:< 4.0f}{9:< 4.0f}{10:< 4.0f}{11:< 4.0f}{12:< 4.0f}{13:< 4.0f}{14:< 4.0f}{15:< 4.0f}{16:< 4.0f}{17:< 4.0f}{18:< 4.0f}{19:< 4.0f}{20:< 4.0f}{21:< 3.0f}{22:< 8.4f}{23:< 8.4f}{24: .3f}\n'.format(*row)
                newmas.write(rowStr)
            else:
                newmas.write(line)
        
        else:
            newmas.write(line)
        if line.startswith('SCAT'):
            scatTab = True
            
    mas.close()
    newmas.close()
    
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
    
    inp = open('xd.inp','r')
    newinp = open('xdnew.inp','w')
    
    i = 0
    cuFound = False
    
    for line in inp:
        if line[:2].upper() == 'CU':
            cuFound = True
        
        if cuFound:
            i+=1
        
        if i == 3:
            row = str.split(line)
            row[0] = 10.0000
            rowStr = '{0:< 9.4f}{1}\n'.format(*row)
            newinp.write(rowStr)
            cuFound = False
            i = 0
        
        else:
            newinp.write(line)
    
    inp.close()
    newinp.close()
    
    os.remove('xd.inp')
    os.rename('xdnew.inp','xd.inp')


def loadBackup(folderName):
    '''
    Copy files from given folder to project foler.
    '''
    folder = folderName                                 #folderName is the absolute path to the backup folder
    
    try:                                                #Copy files from backup folder to working directory
        copyfile(folder + '/xd_lsm.out','xd_lsm.out')   #try and except used in case one of the files doens't exist
    except:
        pass
    try:
        copyfile(folder + '/xd.mas','xd.mas')
    except:
        pass
    try:
        copyfile(folder + '/xd.inp','xd.inp')
    except:
        pass
    try:
        copyfile(folder + '/xd.res','xd.res')
    except:
        pass


def addNCST():
    '''
    Fix 'paramter/ncst size should be increased' error.
    '''
    mas = open('xd.mas', 'r')                       #Open xd.mas and xdnew.mas
    newmas = open('xdnew.mas','w') 
    
    i = 40
    
    for line in mas:                                #Go through xd.mas line by line

        #Add SIZE ncst 30 line at appropriate place
        if line.startswith('SAVE'):
            newmas.write(line)
            newmas.write('SIZE ncst {}\n'.format(i))    
            
        elif not line.startswith('SIZE ncst {}\n'.format(i)):
            newmas.write(line)
            
    #Close files
    newmas.close()
    mas.close()
        
    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')  


def addDUM(atom1, atom2, dumI = 1):
    '''
    Add dummy atom directly in between two atoms. Return number in name of new dummy atom.
    '''
    #Open "xd.inp" to read
    inp = open('xd.inp', 'r') 
    mas = open('xd.mas','r')
    
    addedDUMs = {}                          #Variable to know what dummy atoms are already added
 
    #Go through mas file line by line
    #Find any dummy atoms already there
    #[:-1] removes \n from end, save every dummy atom coords as key with index as value
    for line in mas:
        if line.startswith('DUM'):
            row = str.split(line)
            addedDUMs[tuple(row[1:])] = row[0][3:]
            if row[0] == atom1:
                atom1c = [float(row[1]), float(row[2]), float(row[3])]
            elif row[0] == atom2:
                atom2c = [float(row[1]), float(row[2]), float(row[3])]
            
        elif line.startswith('END ATOM'):
            break
          
    mas.close()                             #Close mas file
             
    atom1Line=[]                            #Lists to save inp lines for each atom
    atom2Line=[]
    #***Should add code to search file and index atom labels and dummy atoms to check if labels given are legit.
    
    #Find lines in xd.inp containing info on chosen atoms, and convert the lines to lists of numbers
    for line in inp:
        row = str.split(line)
        if row[0].upper() == atom1.upper(): 
            atom1Line = (str.split(line))
            atom1c = [float(atom1Line[12]), float(atom1Line[13]), float(atom1Line[14])]
        elif row[0].upper() == atom2.upper():
            atom2Line = (str.split(line))
            atom2c = [float(atom2Line[12]), float(atom2Line[13]), float(atom2Line[14])]
           
    inp.close()
    

    
    #Find the average of the x,y and z fractional coordinates of the 2 atoms
    dumX = (atom1c[0] + atom2c[0]) / 2
    dumY = (atom1c[1] + atom2c[1]) / 2   
    dumZ = (atom1c[2] + atom2c[2]) / 2

    #Open xd.mas and create new mas file xdnew.mas
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas', 'w')
    
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
                dumAtom = '!DUM{0} is halfway between {4} and {5}\nDUM{0:<5}{1:8.4f}{2:8.4f}{3:8.4f}\n'.format(dumI,dumX,dumY,dumZ,atom1,atom2)   #When a valid DUM index has been found write the dummy atom in the correct format
                newmas.write(dumAtom)
                    
        newmas.write(line)                          #Write every other line unchanged
        
    #Close both files, delete original xd.mas and rename xdnew.mas to xd.mas
    mas.close()
    newmas.close()
    
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
    
    return (str(dumI))                                   #Return dumI so that getLocalCoordSys knows what dummy atom label to use
    
def addCustomLocCoords(Patom, atom1, axis1, atom2, axis2, sym):
    '''
    Add custom local coordinate system. 
    Return letters 'a' and 'w' to confirm the atom and key tables were updated, respectively.
    '''
    try:
        mas = open('xd.mas', 'r')
        newmas = open('xdnew.mas', 'w')
        
        atomTab = False
        keyTab = False
        result = ''
        
        multipoleBank = {'NO': '00 000 00000 0000000 000000000', '1': '10 111 11111 1111111 111111111', 
                     'CYL': '10 001 00000 0000000 000000000', 'CYLX': '10 001 10000 1000000 000000000', 
                     'CYLXD': '10 001 10000 0000000 000000000', '2': '10 001 10010 1001000 100100010', 
                     'M': '10 110 10011 0110011 100110011', 'MM2': '10 001 10010 1001000 100100010', 
                     '4': '10 001 10000 1000000 100000010', '4MM': '10 001 10000 1000000 100000010', 
                     '3': '10 001 10000 1000010 100001000', '3M': '10 001 10000 1000010 100001000',  
                     '6': '10 001 10000 1000000 100000000', '6MM': '10 001 10000 1000000 100000000'}
        
        for line in mas:
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            elif line.startswith('KAPPA'):
                keyTab = False
                
            if atomTab:
                row = str.split(line)
                if row[0].upper() == Patom.upper():
                 
                    row[1] = atom1
                    row[2] = axis1
                    row[4] = atom2
                    row[5] = axis2
                    if len(row) > 12:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}\n'.format(*row)
                    else:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}\n'.format(*row)
                    newmas.write(rowStr)
                    result += 'a'
                else:
                    newmas.write(line)
            
                
            elif keyTab:
                row = str.split(line)
                if row[0].upper() == Patom.upper():
                    try:
                        rowStr = '{0:8}{1} {2}\n'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[sym])
                        newmas.write(rowStr)
                        result += 'w'
                    except Exception:
                        newmas.write(line)
                else:
                    newmas.write(line)
            else:
                newmas.write(line)
                
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                
            elif line.startswith('KEY'):
                keyTab = True
        
        mas.close()
        newmas.close()
        
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
    
    finally:
        mas.close()
        newmas.close()
    
    return result
    
'''	
#-------------------RESET BOND----------------------------------------
'''

def armRBs():
    '''
    Enable all reset bond instructions.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        for line in mas:
            if line.startswith('!RESET BOND'):
                newmas.write(line[1:])
            else:
                newmas.write(line)
    
    os.remove('xd.mas')
    os.rename('xdnew.mas', 'xd.mas')
    

def disarmRBs():
    '''
    Disable all reset bond instructions.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        for line in mas:
            if line.startswith('RESET BOND'):
                newmas.write('!' + line)
            else:
                newmas.write(line)
    
    os.remove('xd.mas')
    os.rename('xdnew.mas', 'xd.mas')
    
    
def check4RB():
    '''
    Check if reset bond instructions have been added. Return result as bool.
    '''
    try:
        mas = open('xd.mas','r')
        RB = False
        H = False
    
        for line in mas:
            #Find place in file to write RESET BOND instructions
            if line.startswith('RESET BOND'):
                RB = True
                break
            
            #If there are no Hs then function should return True as no RESET BOND instructions are required.
            elif line.startswith('H('):
                H = True
        mas.close()
        
        if not H:
            RB = True
            
        return RB
    
    except:
        return False


def resetBond(length,atoms,allornot = False):
    '''
    Add reset bond instruction to xd.mas for given atoms and length.
    '''
    mas = open('xd.mas', 'r')
    atomTab = False
    resetBonds = {}
    addedRBs = {}
    
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        #If 'All' is unchecked, make dictionary with H atoms inputted by user, their connected C atoms and inputted bond distance (1.09 default)
        if atomTab and line.startswith('H(') and allornot==False:
            row = str.split(line)
            row = [item.upper() for item in row]
              
            if row[0] in atoms:
                    resetBonds[row[0]] = [row[1],row[0],length]
        
        #If 'All' is checked make dictionary with all H atoms, their connected C atoms and inputted bond distance (1.09 default)
        elif atomTab and line.startswith('H('):
            row = str.split(line)
              
            resetBonds[row[0]] = [row[1],row[0],length]
        
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
        
        #Add RESET BOND instructions that aren't in the new list to the list so they are written
        #Any duplicate instructions aren't added and so only the new instruction is added to xd.mas
        if line.startswith('RESET BOND') or line.startswith('!RESET BOND'):
            row = str.split(line)
            i=0
            for item in row:
                
                if item[:1].isdigit():
                    if row[i-1] not in list(resetBonds.keys()):
                        resetBonds[row[i-1]] = [row[i-2],row[i-1],row[i]]
                i+=1                     
    
    mas.close()
    
    #Open test file to write
    newmas = open('xdnew.mas','w') 
    mas = open('xd.mas','r')
    j = 0
    startCount = False
    
    for line in mas:
        #Find place in file to write RESET BOND instructions
        if line.startswith('END ATOM'):
            
            startCount = True
        
        if j == 5:                      #Print RESET BOND instructions 5 lines after end of atom table
            newmas.write(line)
            newmas.write('RESET BOND  ')
            
            #variable to split RESET BOND onto several lines for easier reading
            i=1
            #Go through resetBonds and print the instructions moving to a new line every 7
            for key in resetBonds:
                if key not in list(addedRBs.keys()):
                    if (i%7==0):
                        resetIns = '\nRESET BOND  {0} {1} {2}  '.format(resetBonds[key][0], resetBonds[key][1], resetBonds[key][2])
                        newmas.write(resetIns)
                        
                    else:
                        resetIns = '{0} {1} {2}  '.format(resetBonds[key][0], resetBonds[key][1], resetBonds[key][2])
                        newmas.write(resetIns)
                    i += 1 
                
                else:
                    if resetBonds[key][2] != addedRBs[key][2]:
                        if (i%7==0):
                            resetIns = '\nRESET BOND  {0} {1} {2}  '.format(resetBonds[key][0], resetBonds[key][1], resetBonds[key][2])
                            newmas.write(resetIns)
                        
                        else:
                            resetIns = '{0} {1} {2}  '.format(resetBonds[key][0], resetBonds[key][1], resetBonds[key][2])
                            newmas.write(resetIns)
                        i += 1
        
            newmas.write('\n')
        
        elif line.startswith('RESET BOND') or line.startswith('!RESET BOND'):
            pass
        
        else:
            newmas.write(line)
        
        if startCount:
            j += 1
        
    newmas.close()
    mas.close()
    
    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def autoAddResetBond():
    '''
    Automatically add reset bond instructions based on detected H environments.
    '''
    #Get neighbours by type and label
    neighboursType = copy.copy(globAtomTypes)
    neighbourLabsRaw = copy.copy(globAtomLabs)
    CmethylList=[]
    XmethylList = []
    primaryHList = []
    secondaryHList = []
    aromaticHList = []
    acidHList = []
    alcoholHList = []
    nitro3HList = []
    nitro4HList = []
    alkeneHList = []
    XXCHList = []
    H2OHList = []
    addedRBs = []
     
    neighbourLabs = {}
    

    neighbourLabs = {atom:[neeb.split(',')[0] for neeb in neebs] for atom, neebs in neighbourLabsRaw.items()}
    
    #Go through neighbours by type
    for atom,neebs in neighboursType.items():

        sig = ''.join(sorted(neebs))
        #If a C is found if it has the CHHH neebor signature then find the labels of the Hs and add them to methyl list
        if atom[0:2] == 'C(':

            if sig == 'CHHH':
                atomNeebLabs = neighbourLabs[atom]
                for item in atomNeebLabs:
                    if item[0:2] == 'H(':
                        CmethylList.append(item)
                         
                #If a C is found if it has the CHHH neebor signature then find the labels of the Hs and add them to methyl list
            elif sorted(neebs).count('H') == 3 and 'C' not in neebs:
                atomNeebLabs = neighbourLabs[atom]
                for item in atomNeebLabs:
                    if item[0:2] == 'H(':
                        XmethylList.append(item)
                        
            elif 'HH' in sig and not 'HHH' in sig:
                atomNeebLabs = neighbourLabs[atom]
                for item in atomNeebLabs:
                    if item[0:2] == 'H(':
                        primaryHList.append(item)
                        
            elif 'H' in sig and not 'HH' in sig and not 'HHH' in sig and len(sig) == 4:
                atomNeebLabs = neighbourLabs[atom]
                for item in atomNeebLabs:
                    if item[0:2] == 'H(':
                        secondaryHList.append(item)
                        
            elif sig == 'CCH':
                
                aromaticConfirmed = confirmAromaticity(atom + ',asym')

                if aromaticConfirmed:   

                     atomNeebLabs = neighbourLabs[atom]
                     for item in atomNeebLabs:
                        if item[0:2] == 'H(':
                            aromaticHList.append(item)
                            
                elif len(neebs) == 3 and neebs.count('H') == 1 and not aromaticConfirmed:
                    for atomLab in neighbourLabs[atom]:
                        if atomLab.startswith('H('):
                           alkeneHList.append(atomLab)
                           
                           break
                       
            elif len(neebs) == 3 and neebs.count('H') == 1:
                for atomLab  in neighbourLabs[atom]:
                    if atomLab.startswith('H('):
                        XXCHList.append(atomLab)
                           
        elif atom[0:2] == 'N(':

            if 'H' in neebs:
                
                if len(neebs) == 4:
                    for atomLab in neighbourLabs[atom]:
                        if atomLab[:2] == 'H(':
                            nitro4HList.append(atomLab)
                elif len(neebs) == 3:
                    for atomLab in neighbourLabs[atom]:
                        
                        if atomLab[:2] == 'H(':
                           
                           nitro3HList.append(atomLab)
                           
        elif atom[:2] == 'O(':
           if 'H' in neebs and 'C' in neebs:
               for atomLab in neighbourLabs[atom]:
                   if atomLab[:2] == 'C(':
                       if neighboursType[atomLab].count('O') == 2:
                           for atomLab2 in neighbourLabs[atom]:
                               if atomLab2[:2] == 'H(':
                                   acidHList.append(atomLab2)
                                   break
                       
                       if neighboursType[atomLab].count('O') == 1:
                           for atomLab2 in neighbourLabs[atom]:
                               if atomLab2[:2] == 'H(':
                                   alcoholHList.append(atomLab2)
                                   break
       
           elif ''.join(neebs) == 'HH':

               for atomLab in neighbourLabs[atom]:
                   H2OHList.append(atomLab)
                   
    if XXCHList:
        resetBond('1.08', XXCHList)
        addedRBs.extend(XXCHList)               
    
    if H2OHList:
        resetBond('0.96',H2OHList)
        addedRBs.extend(H2OHList)
    
    if alcoholHList:        
        resetBond('0.967',alcoholHList)
        addedRBs.extend(alcoholHList)
        
    if acidHList:        
        resetBond('1.015',acidHList)
        addedRBs.extend(acidHList)

    if nitro4HList:        
        resetBond('1.033',nitro4HList)
        addedRBs.extend(nitro4HList)
        
    if nitro3HList:        
        resetBond('1.009',nitro3HList)
        addedRBs.extend(nitro3HList)
     
    if aromaticHList != []:        
        resetBond('1.083',aromaticHList)
        addedRBs.extend(aromaticHList)
                         
    if secondaryHList != []:          
        resetBond('1.099',secondaryHList)
        addedRBs.extend(secondaryHList)

    if primaryHList != []:            
        resetBond('1.092',primaryHList)
        addedRBs.extend(primaryHList)

    if XmethylList != []:            
        resetBond('1.066',XmethylList)
        addedRBs.extend(XmethylList)
                         
    if CmethylList != []:          
        resetBond('1.059',CmethylList)
        addedRBs.extend(CmethylList)
        
    if alkeneHList:
        resetBond('1.077',alkeneHList)
        addedRBs.extend(alkeneHList)

    addedRBs = [rb.upper() for rb in addedRBs]
    return addedRBs

def getPathAromatic(atom, atomNeebDict, usedBranches, lastPath=[], pathStr = ''):
    '''
    Find a path through the structure from a given atom, limitied by a list of branches already visited.
    Return a pipe separated string of path in the format 'C(1),asym|N(1),asym|H(1),asym'.
    Return the last branch in the path.
    Return the list of visited branches.
    '''
    currAtom = atom             #Current atom as program walks through structure.
    passedAtoms = []            #Atoms already visited in the path.
    
    newUsedBranches = []        #List of all branches explored in the structure.
                                #Branches are format 'ATOM~number of steps through path' i.e. 'C(1),asym~2'
                                #means C(1) was passed on the second step through the structure.
    
    aromaticConfirmed = False     
    steps = -1                   #Tracks number of steps taken down path.
    
    while currAtom not in passedAtoms:        #This condition will be satisfied until path can't go anywhere unvisited.
        
        passedAtoms.append(currAtom)          #Add current atoms to visited atoms list.
        steps += 1                            #Add step to total number of steps in path
        
        atomNeebs = atomNeebDict[currAtom.split(',')[0]]        #Get neighbours of current atom.
        
        if steps == 5:
            
            for neeb in atomNeebs:
                if neeb == atom:
                    aromaticConfirmed = True
                    return (pathStr, newUsedBranches[-1:], usedBranches, aromaticConfirmed)    #Return the path string, last new branch, and edited list of used branches.
        
        for neeb in atomNeebs:                              #Go through neighbours of current atom one by one.
            branchTag = '{0}~{1}'.format(neeb, str(steps))  #Make 'C(1),asym~2' format branch tag for current atom and number of steps.

            #If a neighbour of current atom is found that hasn't been visited,
            #and hasn't been branched too in other trips down the same path.
            if lab2type(neeb) == 'C' and neeb not in passedAtoms and branchTag not in usedBranches:

                try:
                    #If this path changes from the previous path at this atom,
                    #remove all branch tags downstream from current number of steps.
                    #i.e. if number of steps is 3 only keep branch tags with number of steps 0,1,2,3.
                    if lastPath[steps] != neeb:
                        usedBranches = [item for item in usedBranches if int(item.split('~')[1]) <= steps]
                        
                except IndexError:
                    #If the last path didn't go this far, or it was empty also remove downstream branch tags.
                    usedBranches = [item for item in usedBranches if int(item.split('~')[1]) <= steps]
                
                #Add branch tag of every atom in path to list.
                newUsedBranches.append(branchTag)

                pathStr += neeb + '|'   #Add new atom in path to pipe separated path string.
                currAtom = neeb         #Make the current atom the new atom, quit the for loop
                break                   #and look for the next atom along the path.
    
    pathStr = pathStr.strip('|')        #Remove the | from the end of the path string.
    
    return (pathStr, newUsedBranches[-1:], usedBranches, aromaticConfirmed)    #Return the path string, last new branch, and edited list of used branches.


def confirmAromaticity(atom):
    '''
    Check if atom is part of a phenyl ring. Return result as bool.
    '''
    lastPath = []
    paths = []
    usedBranches = []
    atomNeebs = copy.copy(globAtomLabs)
    aromaticConfirmed = False
    
    while True:
        
        pathRes = getPathAromatic(atom, atomNeebs, usedBranches, lastPath)  #Get a path.
        
        if pathRes[3] == True:
            aromaticConfirmed = True
            break
        
        #If no path is returned all paths have been found and the while loop is ended.
        if not pathRes[0]:
            break
        
        else:  

            #Store path as a list, to be the last path for the next time getPath is called.
            lastPath = pathRes[0].split('|')    

            #Update usedBranches to list returned from getPath.
            usedBranches = pathRes[2]
            
            #Format path to string of atom types e.g. 'CCCNCCCH'
            pathsFormatted = ''.join([item.split('(')[0] for item in pathRes[0].split('|')])

            paths.append(pathsFormatted)    #Add this path string to list of paths.

            usedBranches.extend(pathRes[1]) #Add the last branch from the current path to the list of used branches.
       
    #When no more paths can be found or aromatic path has been found return aromatic results.    
    return aromaticConfirmed


def autoResetBond():
    '''
    Automatically add reset bond instructions to xd.mas and find out what H atoms have been missed.
    Return missed atoms.
    '''
    resetBondsAdded = autoAddResetBond()

    with open('xd.mas', 'r') as mas:
        missedAtoms = []
        atomTab = False
        
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
                
            if atomTab and line.startswith('H('):
                row = str.split(line)
                
                if row[0].upper() not in resetBondsAdded:
                    missedAtoms.append(row[0])
                
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
    print(missedAtoms)
    return missedAtoms

  
def delResetBond():
    '''Delete all reset bond instructions from xd.mas.'''
    #Open xd.mas to read and xdnew.mas to write new mas file
    mas = open('xd.mas','r')
    newmas = open('xd.newmas','w')
    
    #Go through xd.mas line by line
    #If line has reset bond instructions don't write it, otherwise write line unchanged
    for line in mas:
        if line.startswith('RESET BOND') or line.startswith('!RESET BOND'):
            pass
        else:
            newmas.write(line)
    
    #Close files and rename xdnew.mas to xd.mas
    mas.close()
    newmas.close()
    
    os.remove('xd.mas')
    os.rename('xd.newmas','xd.mas')

'''	
#-------------------CHEMCON-------------------------------------------
'''

def findCHEMCON():
    '''
    Find chemical equivalency in structure. Return CHEMCON dictionary.
    '''    
    
    global globAtomEnv          #Global dictionary of atoms and their chemical environment hash values i.e. {'C(1)':'f11390f0a9cadcbb4f234c8e8ea8d236'}
    globAtomEnv = {}
    
    with open('xd.mas', 'r') as mas:
        
        envs = {}
        CHEMCON = {}
        atomTab = False
        
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            
            #Make dictionary of environment hash values and atoms with that environment
            #i.e. {'f11390f0a9cadcbb4f234c8e8ea8d236' : ['H(4)', 'H(3)']}
            
            if atomTab:
                row = str.split(line)
                atom = row[0].upper()
                print(atom)
                atomEnv = getEnvSig(atom + ',asym')
                globAtomEnv[atom] = atomEnv         
                envs.setdefault(atomEnv,[]).append(atom)
                
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                
        #Organise hash value dictionary into dictionary of parent atoms that appear first in ATOM table,
        #and children that appear further down.
        for env, atoms in envs.items():
            if len(atoms) > 1:
                CHEMCON[atoms[0]] = atoms[1:]
            else:
                CHEMCON[atoms[0]] = []
    
    return CHEMCON
            
'''
Short guide to path finding algorithm:
    
    findAllPaths will start with the first neighbour to starting atom, and exhaust all possibilities,
    before moving to the next neighbour. This applies no matter how many steps you are from the
    starting atom. For example, all branches 5+ steps down would be explored for one 5 step neighbour, before
    any other of the 5 step neighbours were explored.
    
    Only adding the last branch to usedBranches causes this behaviour. If the current path differs from 
    the last path at any stage then everything downstream is removed from usedBranches, as these branches
    are no longer used as something has been changed upstream.
'''

def getPath(atom, atomNeebDict, usedBranches, lastPath=[], pathStr = ''):
    '''
    Find a path through the structure from a given atom, limitied by a list of branches already visited.
    Return a pipe separated string of path in the format 'C(1),asym|N(1),asym|H(1),asym'.
    Return the last branch in the path.
    Return the list of visited branches.
    '''
    currAtom = atom             #Current atom as program walks through structure.
    passedAtoms = []            #Atoms already visited in the path.
    
    newUsedBranches = []        #List of all branches explored in the structure.
                                #Branches are format 'ATOM~number of steps through path' i.e. 'C(1),asym~2'
                                #means C(1) was passed on the second step through the structure.
                                
    steps = -1                   #Tracks number of steps taken down path.
    
    while currAtom not in passedAtoms:        #This condition will be satisfied until path can't go anywhere unvisited.
        
        passedAtoms.append(currAtom)          #Add current atoms to visited atoms list.
        steps += 1                            #Add step to total number of steps in path
        
        atomNeebs = atomNeebDict[currAtom.split(',')[0]]        #Get neighbours of current atom.

        for neeb in atomNeebs:                              #Go through neighbours of current atom one by one.
            branchTag = '{0}~{1}'.format(neeb, str(steps))  #Make 'C(1),asym~2' format branch tag for current atom and number of steps.

            #If a neighbour of current atom is found that hasn't been visited,
            #and hasn't been branched too in other trips down the same path.
            if neeb not in passedAtoms and branchTag not in usedBranches:

                try:
                    #If this path changes from the previous path at this atom,
                    #remove all branch tags downstream from current number of steps.
                    #i.e. if number of steps is 3 only keep branch tags with number of steps 0,1,2,3.
                    if lastPath[steps] != neeb:
                        usedBranches = [item for item in usedBranches if int(item.split('~')[1]) <= steps]
                        
                except IndexError:
                    #If the last path didn't go this far, or it was empty also remove downstream branch tags.
                    usedBranches = [item for item in usedBranches if int(item.split('~')[1]) <= steps]
                
                #Add branch tag of every atom in path to list.
                newUsedBranches.append(branchTag)

                pathStr += neeb + '|'   #Add new atom in path to pipe separated path string.
                currAtom = neeb         #Make the current atom the new atom, quit the for loop
                break                   #and look for the next atom along the path.
    
    pathStr = pathStr.strip('|')        #Remove the | from the end of the path string.
    
    return (pathStr, newUsedBranches[-1:], usedBranches)    #Return the path string, last new branch, and edited list of used branches.


def findAllPaths(atom):
    '''
    Find all possible paths through the molecule starting from a given atom. Return all paths in list.
    '''
    lastPath = []
    paths = []
    usedBranches = []
    atomNeebs = copy.copy(globAtomLabs)
    
    while True:
        
        pathRes = getPath(atom, atomNeebs, usedBranches, lastPath)  #Get a path.
        
        #If no path is returned all paths have been found and the while loop is ended.
        if not pathRes[0]:
            break
        
        else:  

            #Store path as a list, to be the last path for the next time getPath is called.
            lastPath = pathRes[0].split('|')    

            #Update usedBranches to list returned from getPath.
            usedBranches = pathRes[2]
            
            #Format path to string of atom types e.g. 'CCCNCCCH'
            pathsFormatted = ''.join([item.split('(')[0] for item in pathRes[0].split('|')])

            paths.append(pathsFormatted)    #Add this path string to list of paths.

            usedBranches.extend(pathRes[1]) #Add the last branch from the current path to the list of used branches.
       
    #When no more paths can be found return list of paths.    
    return paths

def getEnvSig(atom):
    '''
    Create a unique hash value for the chemical environment of a given atom. Return this hash value.
    '''
    #Sort list of paths alphabetically so it is the same no matter what order paths were found in
    #Make string with starting atom types followed by , joined sorted list of all paths.
    pathString = atom.split('(')[0].upper() + ','.join(sorted(findAllPaths(atom)))
    hashObj = hashlib.md5(bytes(pathString,'utf-8'))        #Generate unique hash value of paths.
    return hashObj.hexdigest()                              #Return digest of hash value.

def removeCHEMCON():
    '''
    Remove all CHEMCON from xd.mas.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        atomTab = False
    
            #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            #Detect end of ATOM table
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            
            if atomTab:
                row = str.split(line)
                if len(row) == 13:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}\n'.format(*row)
                    newmas.write(rowStr)
                else:
                    newmas.write(line)
            else:
                newmas.write(line)

            #Detect start of ATOM table    
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                            
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
            

def findCHEMCONbyElement():
    '''
    Find atoms of the same element and group them. Return dictionary of groupings.
    '''
    mas = open('xd.mas','r')
    
    atomTab = False
    prevElement = ''
    CHEMCON = {}
    currentAtom=''
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        #Detect end of ATOM table
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        #In atom table make CHEMCON dictionary
        if atomTab:
            row = str.split(line)
            
            #If first two characters are letters i.e. Co, check to see if it is the first time that element appears
            if line[0:2].isalpha():
                if line[0:2] != prevElement:
                    currentAtom = row[0].upper() 
                    CHEMCON[currentAtom] = []    #Make entry in dictionary for first instance of element in table
                      
                else:
                    CHEMCON[currentAtom].append(row[0].upper())    #Add atoms of same element to dictionary
                    
                prevElement = row[0][0:2]      #Update previous element to current element
                
            else:
                if line[0:1] != prevElement:
                    currentAtom = row[0].upper()
                    CHEMCON[currentAtom] = []        #Make entry in dictionary for first instance of element in table
                else:
                    CHEMCON[currentAtom].append(row[0].upper())     #Add atoms of same element to dictionary
                    
                prevElement = row[0][0:1]       #Update previous element to current element
        #Detect start of ATOM table    
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
                        
    mas.close()
    
    return CHEMCON


def findCHEMCONbyInputElement(inputElementList):
    '''
    Find atoms of given element and group them. Return grouping.
    '''
    mas = open('xd.mas','r')
    
    atomTab = False
    prevElement = ''
    CHEMCON = {}
    currentAtom=''
    eleList = inputElementList
    #Convert elements to upper case
    eleList = [item.upper() for item in eleList]

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        #Detect end of ATOM table
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        #In atom table make CHEMCON dictionary
        if atomTab:
            row = str.split(line)
            
            #If first two characters are letters i.e. Co, check to see if it is the first time that element appears
            if line[0:2].isalpha():
                if line[0:2].upper() in eleList:
                    if line[0:2] != prevElement:
                        currentAtom = row[0].upper()  
                        CHEMCON[currentAtom] = []    #Make entry in dictionary for first instance of element in table
                         
                    else:
                        CHEMCON[currentAtom].append(row[0].upper())    #Add atoms of same element to dictionary
                else: 
                    CHEMCON[row[0].upper()] = []  
                prevElement = row[0][0:2]      #Update previous element to current element
                
            else:
                if line[0:1].upper() in eleList:
                    if line[0:1] != prevElement:
                        currentAtom = row[0].upper()
                        CHEMCON[currentAtom] = []        #Make entry in dictionary for first instance of element in table
                    else:
                        CHEMCON[currentAtom].append(row[0].upper())     #Add atoms of same element to dictionary
                else: 
                    CHEMCON[row[0].upper()] = []    
                prevElement = row[0][0:1]       #Update previous element to current element
        #Detect start of ATOM table    
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
                        
    mas.close()

    return CHEMCON
        

def findCHEMCONbyNeebors():
    '''
    Group atoms of the same element with the same nearest neighbours. Return groupings.
    '''
    mas = open('xd.mas','r')

    atomTab = False                         #Bool to detect if the for loop is in the atom table when going through xd.mas
    CHEMCON = {}                            #Dict to return CHEMCON for writing to xd.mas with writeCHEMCON()
    neebsDict = copy.copy(globAtomTypes)    #Dict of nearest neighbours {'C(2):['C','C','H','H']}
    usedSigs = []                           #List of sigs of chemical environment i.e. CCHN = N bonded to CCH
    SigParentDict = {}                      #dict with sigs as keys and parent atom labels of that sig as values
    chemconSig=''                           

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        #Detect end of ATOM table
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        #In atom table make CHEMCON dictionary
        if atomTab:
            row = str.split(line)
            row = [item.upper() for item in row]
            if line[0:1] != 'H':
                
                
                if line[0:2].isalpha():
                    
                    neebList = sorted(neebsDict[row[0].upper()])
                    neebStr=''
                    
                    for item in neebList:
                        neebStr += item
                    chemconSig = neebStr + line[0:2].upper()
                    
                    #If first two characters are letters i.e. Co, check to see if it is the first time that chemconSig has appeared
                    if chemconSig not in usedSigs:
                        usedSigs.append(chemconSig)
                        SigParentDict[chemconSig] = row[0].upper()
                        CHEMCON[row[0].upper()] = []    #Make entry in dictionary for first instance of element in table    
                    else:
                        CHEMCON[SigParentDict[chemconSig]].append(row[0].upper())    #Add atoms of same element to dictionary
                    
                else:
                    neebList = sorted(neebsDict[row[0].upper()])

                    neebStr=''
                    
                    for item in neebList:
                        neebStr += item
                    chemconSig = neebStr + line[0:1]
                    
                    #If first two characters are letters i.e. Co, check to see if it is the first time that chemconSig has appeared
                    if chemconSig not in usedSigs:
                        usedSigs.append(chemconSig)
                        SigParentDict[chemconSig] = row[0]
                        CHEMCON[row[0]] = []    #Make entry in dictionary for first instance of element in table    
                    else:
                        CHEMCON[SigParentDict[chemconSig]].append(row[0])    #Add atoms of same element to dictionary
        #Code for Hs
            else:
                if row[1][0:2].isalpha():
                    neebStr = row[1][0:2].upper()
                else:
                    neebStr = row[1][0:1]
                chemconSig = neebStr + line[0:1]

                    #If first two characters are letters i.e. Co, check to see if it is the first time that chemconSig has appeared
                if chemconSig not in usedSigs:
                    usedSigs.append(chemconSig)
                    SigParentDict[chemconSig] = row[0]
                    CHEMCON[row[0]] = []    #Make entry in dictionary for first instance of element in table   
                else:
                    CHEMCON[SigParentDict[chemconSig]].append(row[0])    #Add atoms of same element to dictionary
            
        #Detect start of ATOM table    
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
                        
    mas.close()
    
    return CHEMCON


def findCHEMCONbyInputAtoms(atomList):
    '''
    Group given atoms. Retun grouping.
    '''
    mas = open('xd.mas','r')

    atomTab = False
    CHEMCON = {}
    atomList = atomList
    parent = ''
    parent2 = ''          #Parent2 is a new parent if an old parent is overwritten
    parentFound = False
    parent2Found = False
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        #Detect end of ATOM table
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        #In atom table make CHEMCON dictionary
        if atomTab:
            
            row = str.split(line)
            
            if row[0].upper() in atomList: 
                
                if parentFound == False:
                
                    CHEMCON[row[0].upper()] = []
                    parentFound = True
                    parent = row[0].upper()
                
                else:
                    CHEMCON[parent].append(row[0].upper())
                    
            elif len(row) == 13:
                print(row[0])
                if not parent2Found:
                    if row[12] == parent:
                        parent2 = row[0].upper()
                        parent2Found = True
                        CHEMCON[parent2] = []
                    
                elif row[12] == parent:
                    CHEMCON[parent2].append(row[0].upper())
                    
                
        #Detect start of ATOM table    
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
                        
    mas.close()
    
    return CHEMCON


def writeCHEMCON(CHEMCONdict):
    '''
    Write CHEMCON column in xd.mas with given groupings.
    '''
    atomTab = False
    CHEMCON = CHEMCONdict

    written = False
    
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas','w') 
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        if atomTab:
            written = False
            row = str.split(line)
            
            #If atom is in dictionary it is added with local coordinate system
            if row[0].upper() in CHEMCON.keys():
                
                #If clause is so that if CHEMCON has already been added for all Cs, you can add CHEMCON for a few of the Cs and the parent will lose its previous CHEMCON 
                if len(row) != 13:
                    newmas.write(line)
                    written = True
                    
                else:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                    newmas.write(rowStr + '\n')
                    written = True
                    
            else:
                for atom,equi in CHEMCON.items():
                    if row[0].upper() in equi:
                        rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                        newmas.write(rowStr + atom + '\n')  
                        written = True
            
            if written == False:
                newmas.write(line)
            
        else:     
            newmas.write(line)
        
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
        
    newmas.close()
    mas.close()
        
    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')  


def check4CHEMCON():
    '''
    Check if CHEMCON has been added to xd.mas. Return result as bool.
    '''
    atomTab = False
    chemcon = False

    try:
        mas = open('xd.mas', 'r')

    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            
            if atomTab:
                row = str.split(line)
                #If atom is in dictionary it is added with local coordinate system
                if len(row) == 13:
                    
                    chemcon = True
                    break
                
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
        mas.close()
    except FileNotFoundError:
        pass
    
    return chemcon
	
'''		
#-------------------RESULTS--------------------------------------------
'''

def FFTDetective():
    '''
    Find closest atom to highest peak in xd_fft.out. Return atom label.
    '''
    with open('xd_fft.out','r') as fft:    #Open output file from XDFFT
    
        table = False                   #Bool to find start of table with highest peaks
        i=-1                            #Index starts on -1 so that when the for loop reaches the first line of the table i == 1
        j = 0
        suspectAtom = ('','')           #Variable initialised for atom at top of table
        peakDict = {}
        
        for line in fft:                #For loop goes through the table line by line
            
            if table and line.startswith(' -----'):
                j+= 1
                if j == 2:
                    table = False
                    break
        
            if table == True:           #If line is in the table i += 1 to move to the next line
                row = str.split(line)
                if len(row) == 10:
                    peakDict[abs(float(row[9]))] = row[5].upper()
                
            if line.startswith('  no  peak'):
                table = True
    
    biggestPeak = max(peakDict.keys())
    suspectAtom = (peakDict[biggestPeak], biggestPeak)
    return suspectAtom


def FOUcell():
    '''
    Setup XDFOUR instructions in xd.mas to run over the entire unit cell.
    '''
    mas = open('xd.mas','r')                    #Open xd.mas
    newmas = open('xdnew.mas','w')              #Open new file to write then copy to xd.mas
    
    fouBool = False                     #Bool to detect start of XDFOUR instructions
        
    for line in mas:                            #Go through xd.mas line by line
        
        if line.startswith('CELL '):            #Find line in xd.mas with dimensions of unit cell
            
            row = str.split(line)               #Split line into list
            
            cellA = float(row[1])               #Get the a,b,c dimensions of unit cell and store them as floats
            cellB = float(row[2])
            cellC = float(row[3])
            
            maxCell = max([cellA,cellB,cellC])  #Find the largest of a,b,c
            
            x = 100/maxCell                     #Get a factor for a,b,c that will make the largest == 100
            
            fouA = int(cellA*x)                 #Multiply a,b,c by this factor to get nx,y,z values for XDFOUR
            fouB = int(cellB*x)
            fouC = int(cellC*x)
            
            newmas.write(line)
            
        elif line.startswith('   MODULE *XDFOUR') or line.startswith('   MODULE XDFOUR'):           #Find start of XDFOUR section
            
            fouBool = True
            
            newmas.write('   MODULE *XDFOUR' + '\n')                                                #Write header and following lines as whole unit cell instructions
            
            newmas.write('SELECT *fobs *fmod1  fmod2  print  snlmin  0.000  snlmax  2.000' + '\n')
            newmas.write('GRID   3-points  perp *cryst' +'\n')
                           
            newmas.write('LIMITS xmin  0.0 xmax  1.0 nx  ' + str(fouA) +'\n')                       #Write lines with n values worked out previously.
            newmas.write('LIMITS ymin  0.0 ymax  1.0 ny  ' + str(fouB) +'\n')                       #These values determine how many slices are taken in each direction
            newmas.write('LIMITS zmin  0.0 zmax  1.0 nz  ' + str(fouC) +'\n')                       #so they depend on the lengths in each direction of the unit cell

        elif not fouBool:
            newmas.write(line)                  #Write every other line in xd.mas unchanged.
        
        if line.startswith('   END XDFOUR'):    #Change fouBool at end of XDFOUR instructions to resume writing lines unchanged.
            fouBool = False
            newmas.write(line)
            
    mas.close()                     #Close both files
    newmas.close()
    
    os.remove('xd.mas')             #Delete original xd.mas file
    os.rename('xdnew.mas','xd.mas') #Rename xdnew.mas to xd.mas

 
def getDumNeebs():
    '''
    Find neighbouring dummy atoms. Return all neighbours as asym or dummy atoms.
    '''
    global globAtomLabs
    
    atomLabsRaw = copy.copy(globAtomLabs)
    dumNeebs = {}
    
    for atom, neebs in atomLabsRaw.items():
        dumNeebs[atom] = [spec2masLab(neeb) for neeb in neebs]
        
    return dumNeebs
        
    
def FOU3atoms(atom):
    '''
    Setup XDFOUR instructions in xd.mas to run on the plane of a given atom and 2 of its neighbours.
    '''
    atom = atom.upper()                 #Convert atom given in argument to uppercase to avoid problems with upper/lower
    fouBool = False                     #Bool to detect start of XDFOUR instructions
    neebsRaw = copy.copy(globAtomLabs)   #Get dict of nearest neighbours with labels for each atom
    neebs = {}
    
    for lab, neebors in neebsRaw.items():
        neebs[lab] = [item.split(',')[0] for item in neebors if item.split(',')[1] == 'asym'] 
    
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        for line in mas:                    #Go through mas file line by line
            
            if line.startswith('   MODULE *XDFOUR') or line.startswith('   MODULE XDFOUR'):                 #Find the start of XDFOUR section.
                fouBool = True                                                                              #Change fouBool to True so that program knows not to write XDFOUR instructions twice
                newmas.write('   MODULE *XDFOUR\n')                                                    #Write XDFOUR instructions to make residual map around 3 atoms
                newmas.write('SELECT *fobs *fmod1  fmod2  print  snlmin  0.000  snlmax  2.000\n')
                newmas.write('GRID  *3-points  perp  cryst\n')
                
                rowStr = '{0}{1:7}{2}'.format('ATOM  label ', atom, 'symm  1 trans 0 0 0 *mark on plot')               #Add atom inputted to XDFOUR instructions
                newmas.write(rowStr + '\n') 
                
                if len(neebs[atom]) > 1:
                    print(atom)
                    print(neebs[atom][0])
                    print(neebs[atom][1])
                    rowStr = '{0}{1:7}{2}\n'.format('ATOM  label ', neebs[atom][0], 'symm  1 trans 0 0 0 *mark on plot') #If H isn't the target atom add first 2 nearest neighbours in neighbour dict as they will be most likely non-H
                    newmas.write(rowStr) 
                                                                                                  
                    rowStr = '{0}{1:7}{2}\n'.format('ATOM  label ', neebs[atom][1], 'symm  1 trans 0 0 0 *mark on plot')
                    newmas.write(rowStr) 
    
                elif len(neebsRaw[atom]) > 1:
                    print(atom)
                    print('elif')
                    
                    i = 0
                    
                    for neeb in neebsRaw[atom][:2]:
                        splitLab = neeb.split(',')
        
                        pos = globAtomPos[neeb][0]
                        print(pos)
                        x = pos[0]
                        y = pos[1]
                        z = pos[2]
                        rowStr = 'XYZ label {0} {1:.4f} {2:.4f} {3:.4f} symm  1 trans 0 0 0 *mark on plot\n'.format(splitLab[0],x,y,z)
                        newmas.write(rowStr) 
                        if i == 1:
                            break
                        i+=1
                
                #1 neighbour atoms
                else:
                    print('else')
                    neighbour = neebsRaw[atom][0]
                    pos = globAtomPos[neighbour][0]
                    x,y,z = pos[0], pos[1], pos[2]
                    splitLab = neighbour.split(',')
                    rowStr = 'XYZ label {0} {1:.4f} {2:.4f} {3:.4f} symm  1 trans 0 0 0 *mark on plot\n'.format(splitLab[0],x,y,z)
                    newmas.write(rowStr) 
                    for neeb in neebsRaw[neighbour.split(',')[0]]:
                        if neeb.split(',')[0] != atom:
                            nextNeighbour = neeb
                            break
                    pos = globAtomPos[nextNeighbour][0]
                    x,y,z = pos[0],pos[1],pos[2]
                    splitLab = nextNeighbour.split(',')
                    rowStr = 'XYZ label {0} {1:.4f} {2:.4f} {3:.4f} symm  1 trans 0 0 0 *mark on plot\n'.format(splitLab[0],x,y,z)
                    newmas.write(rowStr) 
                        
                newmas.write('LIMITS xmin -2.0 xmax  2.0 nx  50\n')     #Finish writing appropriate XDFOUR instructions
                newmas.write('LIMITS ymin -2.0 ymax  2.0 ny  50\n')
                newmas.write('LIMITS zmin  0.0 zmax  0.0 nz  1\n')
                newmas.write('   END XDFOUR\n')
        
            elif not fouBool:
                newmas.write(line)                  #Write every other line in xd.mas unchanged.
            
            if line.startswith('   END XDFOUR'):    #Change fouBool at end of XDFOUR instructions to resume writing lines unchanged.
                fouBool = False
            
    os.remove('xd.mas')                 #Remove xd.mas
    os.rename('xdnew.mas','xd.mas')     #Rename xdnew.mas to xd.mas


def grd2values():
    '''
    Get a list of values from xd_fou.grd. Return list.
    '''
    values = []
    with open('xd_fou.grd','r') as grd:
        x = ''
        while not x.startswith('! Values'):
            x = grd.readline()
        
        for line in grd.readlines():

            values.extend([float(item) for item in str.split(line)])

    return values


def setupPROPDpops():
    '''
    Setup XDPROP instructions in xd.mas to find d-orbital populations.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
    
        for line in mas:                    #Go through mas file line by line
            if line.startswith('!D-POP'):   #If dpop line is turned off write it without the ! to turn it on
                newmas.write(line[1:])
            else:
                newmas.write(line)          #Write every other line unchanged

    os.remove('xd.mas')                 #Remove xd.mas
    os.rename('xdnew.mas','xd.mas')     #Rename xdnew.mas to xd.mas

def getDorbs():
    '''
    Find d-orbital populations in xd_pro.out. Return populations.
    '''
    with open('xd_pro.out','r', encoding = 'utf-8', errors = 'ignore') as pro:      #errors ignore stops an error when the profile out is unreadable because some part isn't utf-8 decodable.
        
        orbTable = False                                    #Bool to detect d-orbital population section
        orbPops = []                                        #Initialise list to add orbital populations
        i = -1                                              #i == -1 so that i == 1 at the first line of the table
        
        for line in pro:                                    #Go through pro file line by line
            
            if line.startswith('   z2/xz'):                 #Change orbTable to False once orbital populations have been passed
                orbTable = False
                
            elif orbTable:                                  #If in the orbital population table i += 1
                i+=1
                
                if i >= 1:                                  #If in the right part of the orbital population table
                    row = str.split(line)                   #split the line into a list
                    orbPops.append((row[0],row[3]))         #and append the orbital and the percentage population as a tuple to the orbPops list
                    
                
            elif line.startswith(' Orbital populations'):   #If at the start of the orbital population table set orbTable to True
                orbTable = True
        
    

    return orbPops  #Return list of tuples with orbitals and their populations


def readSUs(file):
    '''
    Get all SUs from xd_lsm.out and sort them. Return sorted list of tuples of form (atom, parameter, SU).
    '''
    i = 0
    paramsList = []
    paramTab = False
    
    with open(file,'r') as lsm:

        for line in lsm:
            
            if line.startswith(' ATOM #   1'):
                paramTab = True
            
            if paramTab:
                if line.startswith(' ATOM #'):
                    row = line.split()
                    atom = row[4]
                    i=0
                    
                elif i > 5 and not line.startswith('------'):
                   
                    row = line.split()
                    if len(row) == 7:
                         SU = line[48:56]
                         param = row[0]
  
                         paramsList.append((atom, param, float(SU)))
                    
            if line.startswith(' SCALE') and len(line.split()) == 7:
                break    
            
            i+=1

    paramsList.sort(key=lambda item: item[2])
        
    return paramsList


def getRF2(folder=None):
    '''
    Get final RF2 value from xd_lsm.out. Return value.
    '''
    if not folder:
        lsmout = open('xd_lsm.out','r')
    else:
        filePath = folder + '/xd_lsm.out'
        lsmout = open(filePath,'r')
    finalCycle = False
    rf2 = ''
    rFound= False
    
    for line in lsmout:
        #Find final cycle and line with R-value
        if finalCycle:
            if line.startswith('  R{F^2} ='):
                row = str.split(line)
                rf2 = float(row[2])
                rFound = True
                finalCycle = False
        elif line.startswith('                      Residuals after final cycle'):
            finalCycle = True
        
    lsmout.close()
    
    #Only return anything if R-value has been found, for error handling in the button function
    if rFound == True:
        return (rf2*100)

def getNumMultipoles():
    '''
    Get number of multipoles with significant populations from xd.res.
    '''
    i = 30
    j = 0
    
    with open('xd.res','r') as res:
             
        for line in res:

            row = str.split(line)
            if '(' in row[0]:
                i = 0
                
            if 1 < i < 5:
                for item in row:
                    if float(item) > 0.02:
                        j+=1
                
            i += 1
                
    return j

def getNumLowAngRefl():
    '''
    Get number of reflections with sin(theta/lambda) from xd_lsm.out.'
    '''
    i = 0
    
    with open('xd_lsm.out','r') as lsm:
        obs = False
        
        for line in lsm:
            if line.startswith('   NO.   H   K   L SINTHL'):
                obs = True
                
            if obs:
                row = str.split(line)
                if row[0].isdigit():
                    if float(row[4]) < 0.5:
                        i+=1
                
            if line.startswith(' Condition(s) met:'):
                obs = False
                
    return i

def getKrauseParam():
    '''
    Return Krause parameter. If there are no low angle reflections, return None.
    '''
    refl = getNumLowAngRefl()
    mult = getNumMultipoles()
    
    if refl > 0:
        return float(mult/refl) 	
    
def getConvergence(lsmFile):
    '''
    Find out if refinement in xd_lsm.out converged. Return result as bool.
    '''
    convCritMet = False
    
    with open(lsmFile,'r') as lsmout:
    
        for line in lsmout:
                                
            if 'convergence criterion was met' in line:
                convCritMet = True
                
    return convCritMet


def getDMSDA(lsmOut):
    '''
    Get DMSDA results from xd_lsm.out. Return DMSDA results.
    '''
    lsm = open(lsmOut,'r')          #Open xd_lsm.out to read

    dmsdaBool = False               #Bool to detect start of DMSDA section
    dmsdaList = []                  #Initialise list to store DMSDA results
    dmsdaFull = []                  #Dictionary to store 2 atoms as keys and DMSDA as values
    parentAtom = ''                 #Initialise string to store starting atom of interatomic vector
                
    for line in lsm:                #Go through xd_lsm.out line by line
            
        if line.startswith('-') and dmsdaBool == True:  #Detect end of DMSDA section
            dmsdaBool = False                           
                
        if dmsdaBool == True:                           #If in DMSDA section split line into list
            row = str.split(line)                
            
            if not line.startswith('      '):   #If starting atom of interatomic vectors is at the start of the current line
                    
                parentAtom = row[0] #parentAtom = starting atom of interatomic vectors
                i = 1               #i = 1 so as not to include the parent atom
                                                         
                for item in row[1:]:            #Go through row and find atom labels
                    if item[0:1].isalpha(): 
                        if row[i+1] == '*':            #If there is a star DMSDA value is row[i+3] otherwise it is row[i+2]
                            dmsdaFull.append((parentAtom, row[i], int(row[i+3])))         #Append DMSDA value to list
                        else:
                            dmsdaFull.append((parentAtom, row[i], int(row[i+2])))
                    i+=1                                                    #Increment i to track current position in line
                    
            else:       #If starting atom of interatomic vectors is at the start of previous line
                    
                i = 0   #i = 0 this time as first item in row isn't starting atom
                    
                for item in row:            #Go through row and find atom labels
                    if item[0:1].isalpha():
                        if row[i+1] == '*':         #If there is a star DMSDA value is row[i+3] otherwise it is row[i+2]
                            dmsdaFull.append((parentAtom, row[i], int(row[i+3])))     #Append DMSDA value to list
                        else:
                            dmsdaFull.append((parentAtom, row[i], int(row[i+2])))
                    i += 1                                                  #Increment i to track current position in line
                    
        if line.startswith(' ATOM-->  ATOM    /  DIST  DMSDA ATOM'):        #Detect start of DMSDA section
                dmsdaBool = True
    
    dmsdaList = [item[2] for item in dmsdaFull]
    #Get list of dmsda values without interatomic vectors
    averageDMSDA = sum([abs(int(x)) for x in dmsdaList])/len(dmsdaList)    #Calculate average dmsda value
        
    lsm.close()                 #close xd_lsm.out

    return (str('%.1f' % averageDMSDA),str(max(dmsdaList)),dmsdaFull) #return tuple of average dmsda, max dmsda and the dmsda dict)

	
'''	   
#--------------------REFINEMENTS----------------------------------------
'''

def scaleFacRef():
    '''
    Setup xd.mas to refine scale factors.
    '''
    setupmas()
    resetKeyTable()
        
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas','w') 
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
                   
        if line.startswith('SELECT cycle'):
            row = str.split(line)
            row[2] = '-10'
            rowStr = ' '.join(row)
            newmas.write(rowStr + '\n')
                        
        else:     
            newmas.write(line)   
    
    newmas.close()
    mas.close()
        
    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')  
    
    
def highAngleRef(sinthlMin,sinthlMax):
    '''
    Setup xd.mas to refine high angle non-H positions and ADPs.
    '''
    setupmas()
    resetKeyTable()
    
    keyTab = False
        
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas','w') 
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('KAPPA'):
            keyTab = False
            
        if line.startswith('SKIP'):
            row = str.split(line)
            rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOn}  {snlMin:<5.3f} {snlMax:<5.3f}'.format(*row, snlOn = '*sinthl', snlMin = sinthlMin, snlMax = sinthlMax)
            newmas.write(rowStr + '\n')
                        
        elif keyTab:
            row = str.split(line)
            
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
    
    newmas.close()
    mas.close()
        
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
                row = str.split(line)
                
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
                row = str.split(line)
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
    try:
        mas = open('xd.mas', 'r')
        newmas = open('xdnew.mas','w') 
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            
            if atomTab:
                
                row = str.split(line)
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
                
                row = str.split(line)
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
                
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
        

        kapInpRes('xd.res', inpTable, i, Hpresent)
        kapInpRes('xd.inp', inpTable, i, Hpresent)
            

        
        #Create new xd.mas file
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')  
        
        return i
    finally:
        mas.close()
        newmas.close()


def kapInpRes(fileName, inpTable, i):
    '''
    Setup xd.inp or xd.res with correctly formatted kappa parameters.
    '''
    inpTableBool = False
    Hpresent = False
    
    #If user hasn't pressed res2inp also update res file so it is copied over correctly to inp
    if os.path.exists(fileName):
        try:    
            res = open(fileName,'r')
            newres = open('xdnew.buckfast','w')
            
            for line in res:

                if line.startswith('H('):
                    Hpresent = True 

                    
                if line.startswith('USAGE'):
                    row = str.split(line)
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
              
            res.close()
            newres.close()
            
            os.remove(fileName)
            os.rename('xdnew.buckfast',fileName)
        finally:
            res.close()
            newres.close()
    
    
def getEleNum():
    '''
    Get the number of different elements in the compound. Return number.
    '''
    try:
        mas = open('xd.mas','r')
        
        scat = False
        i = 0
        
        for line in mas:
            if line.startswith('END SCAT'):
                break
            
            if scat:
                i += 1
            
            if line.startswith('SCAT'):
                scat = True
            
    finally:
        mas.close()
        
    return i
            

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
    try:    
        mas = open('xd.mas', 'r')
        newmas = open('xdnew.mas','w') 
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            
            if line.startswith('EXTCN'):
                keyTab = False
                            
            if keyTab:
                row = str.split(line)
                
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
                row = str.split(line)
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOff}  {8} {9}'.format(*row, snlOff = 'sinthl')
                newmas.write(rowStr + '\n')
                
            elif line.startswith('!RESET BOND'):
                newLine = line[1:]
                newmas.write(newLine)
            
            else:     
                newmas.write(line)   
                
            if line.startswith('KEY     XYZ'):
                keyTab = True        
        
        newmas.close()
        mas.close()
            
        #Create new xd.mas file
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
    finally:
        mas.close()
        newmas.close()

 
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
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas','w') 
    keyTab = False
    
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('KAPPA'):
                keyTab = False
        
        if keyTab:
            row = str.split(line)
            
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
    
    newmas.close()
    mas.close()
        
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


def mm2Tom():
    '''
    Convert mm2 and cyl to m in atom table.
    '''
    XAtoms = ('F(','CL','BR','I(', 'O(', 'N(')
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        atomTab = False
        
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('!DUM') or line.startswith('DUM'):
                atomTab = False
     
            if atomTab:
                row = str.split(line)

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
            
            if line.startswith('ATOM     ATOM0'):
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
                row = str.split(line)
                
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
            
            elif line.startswith('ATOM     ATOM0'):
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

            
            row = str.split(line)
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
                row = str.split(line)
                
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
#--------------------MISC-----------------------------------------------	
'''

def setupmas():
    '''
    Setup details in xd.mas for refinement.
    '''
    try:
        mas = open('xd.mas', 'r')
        newmas = open('xdnew.mas','w') 
        
        atomTab = False
        
        for line in mas:
            
            row = str.split(line)
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            
            if atomTab:
                row = str.split(line)
                
                row[10] = '4'
                
                if line.startswith('H('):
                    row[7] = '1'
                
                if len(row)==13:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)
                else:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                newmas.write(rowStr + '\n')
                
            #Set model to multipoles max l = 4
            elif row[0:2] == ['SELECT','*model']:
                row = str.split(line)
                rowStr = row[0] + ' ' + row[1] + '  4  ' + row[3] + '  ' + row[4] + '  ' + row[5] + ' ' + row[6] + ' ' + 'F^2' + '  ' + row[8] + ' ' + row[9] + ' ' + row[10]
                newmas.write(rowStr + '\n')
                
            elif line.startswith('!RESET    bond C(1) H(1) 1.09 ...'):
                pass
            
            #Add appropriate fmod1
            elif line.startswith('FOUR') or line.startswith('!FOUR'):
                row = str.split(line)
                rowStr = 'FOUR  fmod1  4  2  0  0   fmod2 -1  2  0  0'
                newmas.write(rowStr + '\n')
                
            #Add convcrit    
            elif row[0:2] == ['SELECT','cycle']:
                
                row[2] = '25'
                row[11] = '*convcrit' 
                row[12] = '0.1E-3'
                rowStr = ' '.join(row)
                newmas.write(rowStr + '\n')
            
            #Turn on rigid bond test
            elif line.startswith('!DMSDA'):
                row = str.split(line)
                rowStr = 'DMSDA    ' + row[1] + '   ' + row[2] 
                newmas.write(rowStr + '\n')
            
            elif line.startswith('SKIP'):
                row = str.split(line)
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOff}  {8} {9}'.format(*row, snlOff = 'sinthl')
                newmas.write(rowStr + '\n')
            
            else:     
                newmas.write(line)
            
                #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                
        #Close files
        newmas.close()
        mas.close()
            
        #Create new xd.mas file
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')  
    
    finally:
        mas.close()
        newmas.close()


def findNeeborType(lstFile):
    '''
    Get nearest neighbour types for every atom from lst file. Return nearest neighbour types.
    '''
    #neighboursSym uses format {'O1':['C','C']}
    neighboursSym = {}
    #Opens lst file from SHELXL refinement. Need to update so user chooses lst file from file explorer.
    lstFile = open(lstFile,'r')
    bondTab = False 
    i=0
    j=0
    atom=''
    
    for line in lstFile:
        
        #Find start of bond table
        if line.startswith(' Bond lengths and angles'):
            bondTab=True
        
        #End of bond table
        if bondTab==True and i==2:
            bondTab=False
        
        #If in right part of bond table and line not empty split line into list
        if bondTab and line.strip():
            row = str.split(line)
            
            #If not on the first atom or last line of a block add connected atoms to dictionary
            if j > 0 and not line.startswith('      '):
                
                
                if row[0][:2].isalpha():
                    neighboursSym[atom].append((row[0][:2]).upper())
                else:
                    neighboursSym[atom].append((row[0][:1]).upper())
            
            #If on the first line of a block add a new key to dictionary for that atom
            if 'Distance' in row:
                if row[0][:2].isalpha():
                    atom = (row[0][:2] + '(' + row[0][2:] + ')').upper()
                else:
                    atom = (row[0][:1] + '(' + row[0][1:] + ')').upper()
                neighboursSym[atom] = []
            j += 1
            
        #Find double empty line at bottom of bond table.
        if not line.strip():
            i+=1
            j=0
        elif line.strip():
            i=0
                
    lstFile.close()
 
    return neighboursSym


def findNeeborLabels(lstFile):
    '''
    Get nearest neighbour labels from lst file. Return nearest neighbour labels.
    '''
    #neighbours uses format {'O1':['C1','C2']}
    neighbours = {}
    
    #Opens lst file from SHELXL refinement. Need to update so user chooses lst file from file explorer.
    lstFile = open(lstFile,'r')
    bondTab = False 
    i=0
    j=0
    atom=''
    
    for line in lstFile:
        
        #Find start of bond table
        if line.startswith(' Bond lengths and angles'):
            bondTab=True
        
        #End of bond table
        if bondTab==True and i==2:
            bondTab=False
        
        #If in right part of bond table and line not empty split line into list
        if bondTab and line.strip():
            row = str.split(line)
            
            #If not on the first atom or last line of a block add connected atoms to dictionary
            if j > 0 and not line.startswith('      '):
                
                
                if row[0][:2].isalpha():
                    neighbours[atom].append((row[0][:2].upper() + '(' + row[0][2:] + ')').upper())
                else:
                    neighbours[atom].append((row[0][:1].upper() + '(' + row[0][1:] + ')').upper())
            
            #If on the first line of a block add a new key to dictionary for that atom
            if 'Distance' in row:
                if row[0][:2].isalpha():
                    atom = (row[0][:2] + '(' + row[0][2:] + ')').upper()
                else:
                    atom = (row[0][:1] + '(' + row[0][1:] + ')').upper()
                neighbours[atom] = []
            j += 1
            
        #Find double empty line at bottom of bond table.
        if not line.strip():
            i+=1
            j=0
        elif line.strip():
            i=0
                
    lstFile.close()

    return neighbours


def convert2XDLabel(atomLabel):
    '''
    Convert C1 label to C(1) label
    '''
    if atomLabel[:2].isalpha():
        newLabel = (atomLabel[:2] + '(' + atomLabel[2:] + ')').upper()
    else:
        newLabel = (atomLabel[:1] + '(' + atomLabel[1:] + ')').upper()
        
    return newLabel


def getBondAngles(lstFile):
    '''
    Get bond angles from lst file. Return bond angles.
    '''
    lstFile = open(lstFile, 'r')
    bondAngleDict = {}
    secondAtomList = []  #List to store list of other nearest neighbours to get the atoms involved in angles
    i = 0
    j = 0
    bondTab = False
    
    for line in lstFile:
        
        if j == 2 and bondTab == True:
            bondTab = False
    
        if bondTab:
            if line.strip():
                j = 0
                row = str.split(line)
                
                if 'Distance' in row:
                    parentAtom = convert2XDLabel(row[0])
                    bondAngleDict[parentAtom] = []
                    secondAtomList.clear()
                            
                else:
                    
                    distanceBeen = False            #Bool to detect if distance number has already been passed in line
                    i = 0
                    #Convert atom to XD label format
                    atomLabStr = convert2XDLabel(row[0])
                    secondAtomList.append(atomLabStr.upper())
                    for item in row:
                        
                        if item[0:1].isdigit():
                            if distanceBeen:
                                bondAngleDict[parentAtom].append((item,[atomLabStr, parentAtom, secondAtomList[i]]))
                                i+=1
                            else:
                                distanceBeen = True
            else:
                j+=1

        
                #Find start of bond table
        if line.startswith(' Bond lengths and angles'):
            bondTab=True
    
    return bondAngleDict  
        

def findSITESYM(trackAtom = None):
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
                        
 
#--------------------------------4 NEIGHBOURS-------------------------------------------------------------
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
        
                        
def roundNum(x, base=5):
    '''
    Rounds number to nearest multiple of 5, or other given base.
    '''
    return int(base * round(float(x)/base))


def writeSITESYM(atomSymDict):
    '''
    Adds SITESYM column to atom table in xd.mas. Returns any atoms for which SITESYM wasn't added.
    '''
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas','w') 
    atomTab = False
    SITESYMs = atomSymDict
    unaddedSym = []
    
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        if atomTab:
            row = str.split(line)
            
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
            
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
            
    newmas.close()
    mas.close()
    
    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
    
    return unaddedSym
        

def findUnaddedSym():
    '''
    Find any atoms with 'NO' in the SITESYM column of the atom table. Return atoms.
    '''
    mas = open('xd.mas', 'r')

    atomTab = False
    noSymAtoms = []
    
    #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
            atomTab = False
        
        if atomTab:
            row = str.split(line)
            
            #If sym label is NO add to list of atoms with unadded SITESYM
            if row[11] == 'NO':
                noSymAtoms.append(row[0])                            
            
        if line.startswith('ATOM     ATOM0'):
            atomTab = True
            
    mas.close()

    return(noSymAtoms)

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
                
                row = str.split(line)
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

def findMasCHEMCON():
    '''
    Get CHEMCON from mas file. Return as dictionary of children and their parents.
    '''
    with open('xd.mas','r') as mas:
        
        atomTab = False
        chemcon = {}
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
                break
            
            if atomTab:
                row = line.upper().split()
                if len(row) == 13:
                    chemcon[row[0]] = row[12]
            
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                
    return (chemcon)

def makeLCS():
    '''
    Write local coordinate systems to mas file for all atoms. Return atoms for which local coordinate system could not be found as list.
    '''
    global globAtomEnv
    
    x = getLocalCoordSys()
    
    if globAtomEnv:
        y = writeLocalCoordSys(getCHEMCONLocalCoordSys(x[0], x[1]))
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
    try:  
        mas = open('xd.mas', 'r')
        newmas = open('xdnew.mas','w') 
        
        atomTab = False
        coordSystems = atomLocCoordsDict

        unaddedLocCoords = []
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
            
            if atomTab:
                row = str.split(line)
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
            
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                    
        newmas.close()
        mas.close()
        
        #Create new xd.mas file
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
    
    finally:
        mas.close()
        newmas.close()
        
    return unaddedLocCoords                 #List of atoms for which local coordinates haven't been added
            

def multipoleKeyTable():
    '''
    Write key table with multipoles based on SITESYM column of atom table.
    '''
    #Initialise dictionary containing multipole KEY table settings for common local symmetries.
    multipoleBank = {'NO': '00 000 00000 0000000 000000000', '1': '10 111 11111 1111111 111111111', 
                     'CYL': '10 001 00000 0000000 000000000', 'CYLX': '10 001 10000 1000000 000000000', 
                     'CYLXD': '10 001 10000 1000000 000000000', '2': '10 001 10010 1001000 100100010', 
                     'M': '10 110 10011 0110011 100110011', 'MM2': '10 001 10010 1001000 100100010', 
                     '4': '10 001 10000 1000000 100000010', '4MM': '10 001 10000 1000000 100000010', 
                     '3': '10 001 10000 1000010 100001000', '3M': '10 001 10000 1000010 100001000',  
                     '6': '10 001 10000 1000000 100000000', '6MM': '10 001 10000 1000000 100000000'}
    
    XAtoms = ('CL', 'F(', 'BR','I(', 'O(', 'N(')
    noHexaAtoms = ['C(', 'O(', 'N(']

    with open('xd.mas','r') as mas:
   
        atomTab = False             #Initialise atomTab and keyTab bools to detect atom table and key table
        keyTab = False
        newMultipoles = {}          #Initialise dict to store atom labels and their corresponding multipoles
        atomSyms = {}
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
        #In atom table, if a row has no chemcon (12 items in row) add to dictionary for printing to new file. 
        #Atoms with CHEMCON will inherit multipole configuration if everything is set to 0.
            if atomTab:
                row = str.split(line.upper())
                
                if len(row)==12:
                    atomSyms[row[0]] = row[11]
                    if line[:2] not in XAtoms:
                        if line[:2] not in noHexaAtoms:
                            #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                            newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                        else:
                            #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                            newMultStr = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                            #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                            newMultipoles[row[0]] = newMultStr[:-9] + '000000000'
                    #If cyl is sym label for a halogen include z2 quadrupole
                    else:
                        if row[11] == 'CYL':
                            if line[:2] != 'F(':
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLX'])
                            else:
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLXD'])
                        else:
                            newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                
                elif len(row) == 13:
                    if row[11] != atomSyms[row[12]]:
                        if line[:2] not in XAtoms:
                            if line[:2] not in noHexaAtoms:
                                #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
    
                            else:
                                #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                                newMultStr = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                                #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                                newMultipoles[row[0]] = newMultStr[:-9] + '000000000'
                        
                        #If cyl is sym label for a halogen include z2 quadrupole
                        else:
                            if row[11] == 'CYL':
                                if line[:2] != 'F(':
                                    newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLX'])
                                else:
                                    newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLXD'])
                            else:
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                            print(newMultipoles[row[0]])
                    
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                
        
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        #Find key table the same as with atomTab
        for line in mas:
            
            if line.startswith('KAPPA'):
                keyTab = False
            
            #In the key table, if the atom needs multipoles they are written
            #otherwise write all 0s and atom will inherit multipoles from CHEMCON atom.
            if keyTab:
                row = str.split(line)
                if row[0].upper() in newMultipoles:
                    newmas.write(newMultipoles[row[0].upper()] + '\n')
                else:
                    newmas.write('{0:8}{1}\n'.format(row[0].upper(), '000 000000 0000000000 000000000000000 00 000 00000 0000000 000000000'))
            else:
                newmas.write(line)
            
            if line.startswith('KEY     XYZ'):
                keyTab = True
    
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
 
    
def resetKeyTable():
    '''
    Reset key table to all 0s.
    '''
    keyTab = False
        
    mas = open('xd.mas', 'r')
    newmas = open('xdnew.mas','w') 
    #Go through xd.mas and flip keyTab to true when you reach the start of the atom table and false when you reach the end of the atom table
    for line in mas:
        
        if line.startswith('EXTCN'):
            keyTab = False
                        
        if keyTab:
            if line[:5] != 'KAPPA':
                row = str.split(line)
            
                rowStr = '{0:8}{1}'.format(row[0], '000 000000 0000000000 000000000000000 00 000 00000 0000000 000000000')
                newmas.write(rowStr + '\n')
            else:
                newmas.write('KAPPA   000000\n')
    
        else:     
            newmas.write(line)   
        
        if line.startswith('KEY     XYZ'):
            keyTab = True        
    
    newmas.close()
    mas.close()
        
    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')  
	
#---------------------------------------------WIZARD------------------------------------------------------


def seqMultRef(l):
    '''
    Add multipoles up to a given l value to key table.
    '''
    l+=5       #Make l correspond to row index in key table

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        keyTab = False
        
        for line in mas:
        
            if line.startswith('KAPPA'):
                keyTab = False
    
            if keyTab:
                row = str.split(line)
                rowStr = line[:48]
                for i in range(6,10):
                    if i <= l:
                        rowStr += (' ' + row[i])
                    else:
                        rowStr += (' ' + ''.join(['0' for char in row[i]]))
                rowStr += '\n'
                newmas.write(rowStr)
                
            else:
                newmas.write(line)
            
            if line.startswith('KEY     XYZ'):
                keyTab = True
                
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def wizAddResetBond():
    '''
    Add reset bond instructions in xdwiz.mas to xd.mas.
    '''
    try:
        rbstr = ''
        with open('xdwiz.mas', 'r') as rb:
            
            
            for line in rb:
                if line.startswith('RESET BOND'):
                    rbstr += line
                
                elif line.startswith('KEY'):
                    break
                
    
            
        with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
            for line in mas:
                if line.startswith('WEIGHT'):
                    newmas.write(line)
                    newmas.write(rbstr)
                else:
                    newmas.write(line)
        
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
    except Exception:
        pass


def wizAddLocCoords():
    '''
    Add local coordinate systems in xdwiz.mas to xd.mas
    '''
    try:
        rb = open('xdwiz.mas', 'r')
        ccstr = ''
        atomTab = False
        
        for line in rb:
            if line.startswith('END ATOM'):
                atomTab = False
                
            if atomTab:
                ccstr += line
                    
            
            if line.startswith('ATOM     ATOM0'):
                atomTab = True 
        rb.close()
        
        mas = open('xd.mas','r')
        newmas = open('xdnew.mas','w')
        
        for line in mas:
    
            if line.startswith('END ATOM'):
                atomTab = False
                
            if atomTab:
                pass
            
            else:
                newmas.write(line)
                    
            if line.startswith('ATOM     ATOM0'):
                atomTab = True
                newmas.write(ccstr)
                
            
        mas.close()
        newmas.close()
        
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
        
    finally:
        mas.close()
        newmas.close()


def wizAddCHEMCON():
    '''
    Add CHEMCON from xdwiz.mas to xd.mas.
    '''
    with open('xdwiz.mas','r') as rb:
       
        cc = []
        atomTab = False
        
        for line in rb:
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
                
            if atomTab:
                row = str.split(line)
                if len(row) == 13:
                    cc.append(row[12])
                else:
                    cc.append(' ')
                    
            
            if line.startswith('ATOM     ATOM0'):
                atomTab = True 

        
        with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

            i=0
            for line in mas:
        
                if line.startswith('END ATOM') or line.startswith('!DUM') or line.startswith('DUM'):
                    atomTab = False
                    
                if atomTab:
                    row = str.split(line)
                    if len(row) == 13:
                        row[12] = cc(i)
                    else:
                        row.append(' ')
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}\n'.format(*row)
                    newmas.write(rowStr)
                
                else:
                    newmas.write(line)
                        
                if line.startswith('ATOM     ATOM0'):
                    atomTab = True
                    
                i+=1
        
        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')
        


'''
#--------------------GUI-------------------------------------------------
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
            self.xdlsmRunning = subprocess.Popen(xdlsmAbsPath, shell = False, cwd = os.getcwd())
            self.startSignal.emit()
            self.xdlsmRunning.wait()
            try:
                fixLsmCif()
            except Exception:
                pass
            self.finishedSignal.emit()
        
        except Exception:
            self.warningSignal.emit()
       
class XDProg(QThread):
    '''
    Run XD Program in QThread.
    '''        
    startSignal = pyqtSignal()
    finishedSignal = pyqtSignal()
    warningSignal = pyqtSignal()
    
    def __init__(self, prog):
        QThread.__init__(self)
        self.xdProgName = prog
        
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
                self.xdProgRunning = subprocess.Popen(self.xdProg, shell = False, cwd = os.getcwd())
                self.startSignal.emit()
                self.xdProgRunning.wait()
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
        fixBrokenLabels()
        
        self.xdiniRunning = subprocess.Popen([xdiniAbsPath, ''.join(compoundID4XDINI.split()), 'shelx'], shell = False, cwd = os.getcwd())
        self.startSignal.emit()
        
        try:
            initialiseGlobVars()
            
        except Exception:
            pass
        
        self.xdiniRunning.wait() 
        removePhantomAtoms()
        self.finishedSignal.emit()
        
        
class checkLSMOUT(QThread):
    '''
    *DOESN'T WORK: Checks xd_lsm.out while XDLSM is running.
    '''
    updateSignal = pyqtSignal()
    
    def __init__(self):
        QThread.__init__(self)
        
        
    def __del__(self):
        self.wait()
    
    def stop(self):
        self.xdlsmIsRunning = False
    
    def run(self):

        i = 0
        self.xdlsmIsRunning = True

        while self.xdlsmIsRunning:
            time.sleep(1)
            k = 4
            lsmout = open('xd_lsm.out','r')
            
            for line in lsmout:

                if line.strip().startswith('Residuals after cycle'):
 
                    row = str.split(line)
                    
                    if int(row[3]) > i:
                        i += 1
                        cycle = row[3]

                        k = 0
                        
                elif k == 3:
                    row = str.split(line)
                    RF2 = float(row[2])*100
                    self.statusMsg = 'Cycle {0} - RF<sup>2</sup> = {1:.2f} %'.format(cycle, RF2)
                    print(self.statusMsg)
                    self.updateSignal.emit()
                    
                k+=1
            lsmout.close()
            
            
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
            

class resmap(QWidget, Ui_resmap):
    '''
    Residual density map window.
    '''
    def __init__(self, parent=None):
        super(resmap, self).__init__(parent)
        self.setupUi(self)
        self.setup()
            
    def setup(self):
        '''
        Create residual density map and display it with save button.
        '''
        fig = makeResMap()
        canvas = FigureCanvas(fig)
        canvas.setParent(self)       
        saveBut = QPushButton('Save PNG file')
        saveBut.clicked.connect(self.savePng)
        saveBut.setFixedWidth(150)
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
        filename = QFileDialog.getSaveFileName(self,"PNG of residual density map",os.getcwd(), "PNG Files (*.png)")
        if filename[0]:
            makeResMap(filename[0])
            self.saveLab.setText('Residual map saved to <i>"{}"</i>'.format(filename[0]))

class NPP(QWidget, Ui_resmap):
    '''
    Show normal probability plot in new window.
    '''
    def __init__(self, parent=None):

        super(NPP, self).__init__(parent)
        self.setupUi(self)
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

        
class aboutBox(QWidget, Ui_aboutBox):
    '''
    About XD Toolkit window.
    '''
    def __init__(self, parent=None):
        super(aboutBox, self).__init__(parent)
        self.setupUi(self)
        
        
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
        os.makedirs('Backup/' + self.folder)
        print('backup made')
        if os.path.isfile('shelx.ins') and os.path.isfile('shelx.hkl'):
            print('is shelx.ins')
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

        if self.collectDat:
            global timeFileAbsPath
            if self.i > 0 and self.i < len(self.refList):
                try:
                    with open(timeFileAbsPath, 'a') as timeFile:
                        
                        self.finishingTime = time.time()
                       
                        runningTime = self.finishingTime - self.startingTime
                        numAtoms = getNumAtoms()
                        ref = self.refList[self.i]
                        refSplit = ref.split()
                        refName = ''.join(refSplit[2:])
    
                        try:
                            refName = refName[:13]
                        except Exception:
                            pass
                        
                        timeStr = '{0:21}{1:10}{2:<13.2f}{3}\n'.format(refName, str(numAtoms), runningTime, sys.platform)
                        timeFile.write(timeStr)
                    
                except Exception:
                    pass
        
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
                    initialiseMas()
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
                multipoleMagician()
                
            elif refName.upper() == ('DIPOLES'):
                multipoleMagician()
                seqMultRef(1)
                
            elif refName.upper().startswith('QUADRUPOLES'):
                multipoleMagician()
                seqMultRef(2)
            
            elif refName.upper().startswith('OCTUPOLES'):
                multipoleMagician()
                seqMultRef(3)
                
            elif refName.upper().startswith('HEXADECAPOLES'):
                multipoleMagician()
                seqMultRef(4)
                
            elif refName.upper().startswith('MULTIPOLES,'):
                posADPMultRef()
            
            elif refName.upper().startswith('LOWER'):
                lowerSym()
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
            
        except Exception:
            self.wizStatusLab.setText('''Couldn't setup xd.mas for ''' + self.refList[self.i][4:].lower().replace('-h','-H'))
    
    
class sendBug(QDialog, Ui_sendBug):
    '''
    Send bug report window.
    '''
    def __init__(self, parent=None):
        super(sendBug, self).__init__(parent)
        self.setupUi(self)
        self.buttonBox.button(QDialogButtonBox.Ok).setText("Send")
        
        
class sendSugg(QDialog, Ui_sendSugg):
    '''
    Send suggestion window.
    '''
    def __init__(self, parent=None):
        super(sendSugg, self).__init__(parent)
        self.setupUi(self)
        
        self.buttonBox.button(QDialogButtonBox.Ok).setText("Send")
    
    
class loadingScreen(QThread):
    '''
    Loading popup while compound initializes.
    '''
    def __init__(self):
        QThread.__init__(self)
        
    def __del__(self):
        self.wait()
        
    def run(self):
        pass
        
    
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

        self.cwdStatusLab = QLabel()
        self.settings = QSettings('prefs')
        self.initialiseSettings()	#Load user preferences
        
        initialiseGlobVars()
    
        self.changeUserIns()
        
        
        #Check for XD files and if they are not there prompt user to find directory.
        if not self.settings.value('xdpath'):
            self.findXD()
                
        self.addedLocCoords = {}
        self.ins = ''
        self.forbiddenChars = ['*', '?', '"', '/', '\\', '<', '>', ':', '|']
        self.invalidFolderNameStr = ('Invalid folder name given.\n'
                                     'Remove the following characters: * ? " / \ < > : |')
        self.permErrorMsg = 'Close MoleCoolQT and try again.'
        self.lstMissingErrorMsg = 'Select add shelx.ins to project folder and try again.'
        self.backupConfirmStr = 'Current files backed up to: '
        self.labList = [self.resNPPLab, self.loadBackupLab, self.customBackupLab, self.autoResetBondStatusLab, self.resetBondStatusLab, self.CHEMCONStatusLab, self.resBackupLab, self.getResLab, self.setupFOURStatusLab, self.getDpopsStatusLab]
        self.tabWidget.setCurrentIndex(0)
        toolboxes = [self.rbToolbox, self.backupToolbox, self.resToolbox, self.toolsToolbox]
        for item in toolboxes:
            item.setCurrentIndex(0)
        #Display current working directory on startup.
        self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())
        self.statusbar.addWidget(self.cwdStatusLab)
        self.toolbarRefreshCwd.triggered.connect(self.refreshFolder)
        self.toolbarSetFolder.triggered.connect(self.setFolder)
        self.toolbarOpenmas.triggered.connect(self.openMas)
        self.menuSendBug.triggered.connect(self.bugRepDialog)
        self.menuSendSugg.triggered.connect(self.sendSuggDialog)
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
        self.menuManual.triggered.connect(lambda: os.startfile(manualAbsPath))
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
        #self.lsmProgUpdate = checkLSMOUT()
        #self.lsmProgUpdate.updateSignal.connect(lambda: self.XDLSMLab.setText(self.lsmProgUpdate.statusMsg))
        self.xdWizardTestBut.clicked.connect(self.wizTest)
        self.wiz2ndStageObjects = [self.wizBackupInput, self.xdWizardBut, self.wizReuseMasBut, self.wizUniSnlMax, self.xdWizardTestBut, self.wizHighSnlMin, self.wizHighSnlMax, self.wizLowSnlMin, self.wizLowSnlMax]
        self.wizSnlInput = [self.wizHighSnlMin, self.wizHighSnlMax, self.wizLowSnlMin, self.wizLowSnlMax]
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
        self.XDLSMLabs = [self.XDLSMLab, self.pkgXDLab]
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
        
#        self.xdgraph = XDGRAPH()
#        self.xdgraph.startSignal.connect(self.startXDGRAPH)
#        self.xdgraph.finishedSignal.connect(self.finishedXDGRAPH)
#        self.XDGRAPHLabs = [self.pkgXDLab]
#        self.XDGRAPHButs = [self.pkgXDGRAPHBut]
#        self.pkgXDGRAPHBut.clicked.connect(lambda: self.xdgraph.start())
        
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
        
        for prog in self.xdProgsRunning:
            prog.warningSignal.connect(self.findXD)
            
#--------------------EMAIL-----------------------------------------------

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

#--------------------WIZARD----------------------------------------------

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
                
        elif os.path.isfile('xd.inp'):
            self.wizCheckIni()
    
    #Check that XDINI has created xd.mas, xd.hkl and xd.inp.
    def wizCheckIni(self):
        '''
        Check that XDINI has created all files and test them.
        Unlock second stage of XD Wizard if all files are present, regardless of test result.
        '''        
        if os.path.isfile('xd.mas') and os.path.isfile('xd.inp') and os.path.isfile('xd.hkl'):
            self.xdWizINILab.setText('Compound initialized successfully. Follow instructions below and click "Test".')
            try:
                self.wizTest()
            except Exception:
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
        r = check4RB()
        s = False
        f = os.path.isfile('xd.mas') and os.path.isfile('xd.inp') and os.path.isfile('xd.hkl')

        multipoleMagician()
        missingSym = findUnaddedSym()

        if missingSym:
            m = False
            
        else:
            m = True
        
        try:
            for item in self.wizSnlInput:
                a = float(str(item.text()))     #checks if input can be converted to float
            
            if str(self.wizUniSnlMax.text()).strip():
                a = float(str(item.text()))
                
            s = True
            
        except ValueError:
            s = False
        print(c)
        print(r)
        print(s)
        print(f)
        print(m)
        if c and r and s and f and m:
            
            self.xdWizardBut.setEnabled(True)
            estTime = totalEstTime()
            minutes, seconds = divmod(estTime, 60)
            hours, minutes = divmod(minutes, 60)
            timeStr = '{0:.0f}hrs {1:.0f}mins {2:.0f}secs'.format(hours, minutes, seconds)
            #wizStr = '{0} Estimated running time is {1}'.format('Ready to run XD Wizard.', timeStr)
            wizStr = 'Ready to run XD Wizard.'
            self.wizTestStatusLab.setText(wizStr)
        
        else: 

            wizStr = ''
            
            if not c:
                wizStr += '-Add chemical equivalency in CHEMCON tab.<br>'

            if not r:
                wizStr += '-Add reset bond instructions in RESET BOND tab.<br>'

            if not m:
                wizStr += '{0}{1}'.format('-Add local coordinate system and SITESYM manually for ', ', '.join(missingSym).strip(', '))

            if not s:
                wizStr += '-Invalid input for sin(&theta;/&lambda;) cutoffs.<br>'

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
        rWarnMsg = 'No reset bond instructions detected.\n'
        cWarnMsg = 'No chemical constraints detected.\n'
        mWarnMsg = 'No local coordinate system added for {}'.format(', '.join(testRes[2][1]).strip(', '))
        testPassed = True
        
        for item in testRes:
            if not item:
                testPassed = False
        
        if testPassed:
            return True
                
        else:
            if testRes[3] and testRes[4]:
                
                msg = procMsg
                if not testRes[0]:
                    msg = cWarnMsg + msg
                if not testRes[1]:
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
        print(self.settings.value('senddata'))
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
    
            if self.wizSeqMultBox.isChecked():
                self.xdWizRunning.refList = ['0 - XDINI', '1 - Scale factors', '2 - High angle non-H positions and ADPs',
                            '3 - Low angle H positions and isotropic ADPs', '4 - Kappa and monopoles', '5 - Dipoles', 
                            '6 - Quadrupoles', '7 - Octupoles', '8 - Hexadecapoles',
                            '9 - Multipoles, and non-H positions and ADPs', '10 - Lower symmetry', 
                            '11 - Final refinement']
    
            if not self.lowerSymBox.isChecked():
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

            copyfile('xd.mas','xdwiz.mas')
            
            self.xdWizRunning.show()
            self.xdWizRunning.xdWiz(str(self.wizBackupInput.text()))        
            
            
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
    #        saveBut = QPushButton('Save PNG file')
    #        saveBut.clicked.connect(self.savePng)
    #        saveBut.setFixedWidth(150)
    #        self.saveLab = QLabel()
            self.nppWizLayout.addWidget(self.nppCanvas, 0)
    #        self.wizLayout.addWidget(saveBut,1)
    #        self.wizLayout.addWidget(self.saveLab,2)
            self.setLayout(self.nppWizLayout)
            # prevent the canvas to shrink beyond a point
            # original size looks like a good minimum size
            self.nppCanvas.setMinimumSize(self.nppCanvas.size())
            self.nppCanvas.setMaximumSize(700,200000)
            self.xdWizardStatusLab.setText(results)
            
            if self.settings.value('senddata') == 'yes':
                self.xdWizRunning.tfin = time.time()
                runningTime = self.xdWizRunning.tfin - self.xdWizRunning.tzero
                with open(timeFileAbsPath,'a') as lsmTimes:
                    lsmTimes.write('{0:10}{1:<15}{2:<13.2f}{3:<13}{4}\n'.format('WIZ', getNumAtoms(), runningTime,  len(self.xdWizRunning.refList), sys.platform))
        
        except Exception:   #Handle user cancelling halfway through
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

#--------------------UTILITIES-------------------------------------------
    
    def multKeyPress(self):
        '''
        Handle user pressing 'Add multipoles to key table button.
        '''
        try:
            multipoleKeyTable()
            self.multKeyStatusLab.setText('Key table updated in xd.mas.')
        
        except Exception:
            self.multKeyStatusLab.setText('An error occurred.')

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
                
                if not os.environ.get['XD_DATADIR']:
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


    def disableXDButs(self, survivors = []):
        '''
        Disable all buttons that run XD programs, except given survivors.
        '''
        for but in self.xdProgButs:
            if but not in survivors:
                but.setEnabled(False)
        self.xdWizCmpID.returnPressed.disconnect()
                
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
           
                        initialiseMas()
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
        
        
#    def startXDGRAPH(self):
#        '''
#        Handle XDGRAPH starting.
#        '''
#        for lab in self.XDGRAPHLabs:
#            lab.setText('XDGRAPH running')
#        
#        for but in self.XDGRAPHButs:
#            but.setText('Cancel')
#            but.disconnect()
#            but.clicked.connect(self.killXDGRAPH)
#            
#        self.xdgraph.finishedSignal.disconnect()
#        self.xdgraph.finishedSignal.connect(self.finishedXDGRAPH)
#        
#    def finishedXDGRAPH(self):
#        '''
#        Handle XDGRAPH finishing.
#        '''        
#        for but in self.XDGRAPHButs:
#            but.setText('Run XDGRAPH')                          
#            but.disconnect()
#            but.clicked.connect(lambda: self.xdgraph.start())
#            
#        for lab in self.XDGRAPHLabs:
#            lab.setText('XDGRAPH finished')
#    
#    def killXDGRAPH(self):
#        '''
#        Kill XDGRAPH.
#        '''
#        try:  
#            self.xdgraph.finishedSignal.disconnect()
#            self.xdgraph.xdProgRunning.terminate()
#            for but in self.XDGRAPHButs:
#                but.setText('Run XDGRAPH')
#                but.disconnect()
#                but.clicked.connect(lambda: self.xdgeom.start())
#            for lab in self.XDGRAPHLabs:
#                lab.setText('XDGRAPH terminated')            
#        except Exception:
#            pass
        
        
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
            self.xdini.finishedSignal.connect(self.enableXDButs)
            self.xdini.start()
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
        self.XDINILab.setText('Compound initialized successfully. Ready to begin refinement.')                         #Sets status label to 'XDINI finished'
        self.refChosen()
        self.runXDINIBut.setText('Run XDINI')  
        self.check4res()
        self.enableXDButs()
        self.changeUserIns()
        removePhantomAtoms()
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
            res2inp()
            self.toolbarRes2Inp.setDisabled(True)
            
    def refreshFolder(self):
        '''
        Reinitialise global variables, check xd.res and change user instructions.
        '''
        initialiseGlobVars()
        self.check4res()
        self.changeUserIns()
      

    def setFolder(self):
        '''
        Prompt user to choose a folder and change current working directory to the folder.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Select Directory"))
        
        os.chdir(folder)
        self.changeWizBackup()
        self.cwdStatusLab.setText('Current project folder: ' + os.getcwd())
        self.resetLabels()
        self.check4res()
        
        try:
            initialiseGlobVars()
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
        except Exception:
            pass

    def loadIAM(self):
        '''
        Load shelx.ins and shelx.hkl from folder of structure solution.
        '''
        ins = False
        hkl = False
        
        folder = str(QFileDialog.getExistingDirectory(None, "Select Structure Solution Folder"))
        projectFolder = str(QFileDialog.getExistingDirectory(None, "Select Project Folder"))

        for file in os.listdir(folder):
            if file[-4:] == '.res':
                copyfile((folder + '/' + file), (projectFolder  + '/shelx.ins'))
                ins = True
            elif file[-4:] == '.hkl':
                copyfile((folder + '/' + file), (projectFolder + '/shelx.hkl'))
                hkl = False
        
        if ins and hkl:
            os.chdir(projectFolder)
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
        
        atom1 = rawInput2Labels(str(self.Atom1.text()))[0]
        atom2 = rawInput2Labels(str(self.Atom2.text()))[0]
        
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
            Patom = rawInput2Labels(Patom)[0]
            axis1 = str(self.alcsAxis1Input.currentText())
            atom1 = str(self.alcsAtom1Input.text())
            atom1 = rawInput2Labels(atom1)[0]
            axis2 = str(self.alcsAxis2Input.currentText())
            atom2 = str(self.alcsAtom2Input.text())
            atom2 = rawInput2Labels(atom2)[0]
            sym = str(self.alcsLocSymInput.text())

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
        
            
#--------------------PREFERENCES-----------------------------------------
        #Set default directory on startup to last working directory or to C:/Users/User
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
        
        
        cachePath = os.getcwd() + '/cache'
        
        if not os.path.isdir(cachePath):
            os.makedirs(cachePath)
            
        timeFileAbsPath = os.getcwd() + '/lsmTimes.buckfast'
        manualAbsPath = os.getcwd() + 'res/XD Toolkit Manual.pdf'
                                   
        if not os.path.isfile(timeFileAbsPath):
            with open('lsmTimes.buckfast','w') as bucky:
                pass

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
        print('senddata')
        print(self.settings.value('senddata'))
        try:
            os.chdir(self.settings.value('lastcwd'))
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
        

#--------------------REFINEMENTS-----------------------------------------
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
        rWarnMsg = 'No reset bond instructions detected.\n'
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
                except:
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
                    
                except:
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
                    
                except:
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
                    
                except:
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
                                    statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add SITESYM and local coordinate system manually in 'Tools' tab.'''
                                    self.manRefSetupLab.setText(statusStr)
                                else:
                                    self.manRefStatusLab.setText(successStr)
                            else:
                                self.manRefSetupLab.setText(successStr)
                        
                    else:
                        posADPMultRef()
                        missingSym = findUnaddedSym()
        
                        if missingSym != []:
                            statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add SITESYM and local coordinate system manually in 'Tools' tab.'''
                            self.manRefSetupLab.setText(statusStr)
                        else:
                            self.manRefSetupLab.setText(successStr)
                        
                except:
                    self.manRefSetupLab.setText(failureStr) 
                    
            elif refNum == 10:              #Lower symmetry
                try:
                
                    posADPMultRef()
                    lowerSym()
                    self.manRefSetupLab.setText(successStr)
                    
                  
                except:
                    self.manRefSetupLab.setText(failureStr)
                
            elif refNum == 11:              #Final refinement
                try:
                    checkMultRes()
                    self.manRefSetupLab.setText(successStr)
                    
                except:
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
            rWarnMsg = 'No reset bond instructions detected.\n'
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
                            statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add SITESYM and local coordinate system manually in 'Tools' tab.'''
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
                        statusStr = 'Unable to find local coordinate system for atoms: ' + ', '.join(missingSym) + '''\nPlease add SITESYM and local coordinate system manually in 'Tools' tab.'''
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
                reqStr += '''- Add chemical equivalency in 'CHEMCON' tab before setting up xd.mas.\n'''
        
        if refNum in (2,3,4,5,6,7,8,9,10,11):

            r = check4RB()
            
            if not r:
                reqStr += '''- Add reset bond instructions in 'RESET BOND' tab before setting up xd.mas.'''
        
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
        
        

#--------------------RESET BOND------------------------------------------

    def armRBPress(self):
        '''
        Handle user clicking enable all reset bond instructions.
        '''
        try:
            armRBs()
            self.armRBLab.setText('All reset bond instructions enabled.')
            self.disarmRBLab.setText('')
        except Exception:
            self.armRBLab.setText('An error occurred.')
            
    def disarmRBPress(self):
        '''
        Handle user clicking disable all reset bond instructions.
        '''
        try:
            disarmRBs()
            self.disarmRBLab.setText('All reset bond instructions disabled.')
            self.armRBLab.setText('')
        except Exception:
            self.disarmRBLab.setText('An error occurred.')
    
    def autoResetBondPress(self):
        '''
        Handle user clicking automatically add reset bond instructions.
        '''
        try:
            missedAtoms = autoResetBond()
            resStr = 'Auto RESET BOND finished.' + '\n' 
            if len(missedAtoms) > 0:
                resStr = 'Auto RESET BOND finished.' + '\n' + 'Appropriate bond lengths not found for atoms:' + '\n' + '\n'    
            
                i = 0
                for atom in missedAtoms:
                    resStr += atom + ', '
                    i += 1
                    
                    if i == 7:
                        resStr += '\n'
                        i = 0
                        
                resStr = resStr[:-2]
                
                resStr += '\n' + '\n' + ' Please add these bond lengths manually before continuing.'
            else:
                resStr = 'Auto RESET BOND finished.' + '\n' + 'RESET BOND instructions added for all H atoms.'
            self.autoResetBondStatusLab.setText(resStr)        
        except:
            self.autoResetBondStatusLab.setText('An error occurred. Please check ins file is in project folder.')        
        self.changeUserIns()
        self.wizTest()
        
    def resetBondInput(self):
        '''
        Handle user manually adding reset bond instructions.
        '''
        try:
            #If 'All' is unchecked get atoms labels from the input and pass them to resetBond()
            if  self.allResetBox.isChecked() != True:
               
                resetAtomRaw = str(self.resetHInput.text().upper())
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
                            newAtomLab = 'H({0})'.format(atomLab[1:])
                        else:
                            newAtomLab = 'H({0})'.format(atomLab)
                        resetAtomList.append(newAtomLab)
                
                        resetAtomStr += (newAtomLab + ', ')
                    else:
                        resetAtomList.append(atomLab)
                resetAtomList = sorted(list(set(resetAtomList)))
                atomLabs = copy.copy(globAtomLabs)
                atomLabs = atomLabs.keys()
                wrongAtoms = [atom for atom in resetAtomList if atom not in atomLabs or atom[:2] != 'H(']
                resetAtomList = [atom for atom in resetAtomList if atom not in wrongAtoms]
                resetAtomStr = ', '.join(lab for lab in resetAtomList).strip(', ')

                resetBond(str(self.resetLengthInput.text()),resetAtomList,False)
            #If 'All' is checked run resetBond() with allornot True
            else:
                resetBond(str(self.resetLengthInput.text()),None,True)
            statusStr = ''
            if resetAtomList:
                statusStr += '{}{}<br>Bond length: {}'.format('RESET BOND instructions added for atoms ',resetAtomStr, str(self.resetLengthInput.text()))
            if wrongAtoms:
                statusStr += '<br><br>{}{}'.format('Following atom labels are incorrect: ', ', '.join(wrongAtoms).strip(', '))
            self.resetBondStatusLab.setText(statusStr)
        
        except:
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
     
    def delResetBondPress(self):
        '''
        Handle user clicking delete all reset bond instructions.
        '''
        try:
            delResetBond()
            self.delResetBondLab.setText('All reset bond instructions removed from xd.mas.')
        except:
            self.delResetBondLab.setText('An error occurred.')
            
            
#--------------------CHEMCON---------------------------------------------
		
    def runCHEMCON(self):
        '''
        Handle user clicking 'Add CHEMCON'.
        '''
        if self.autoCHEMCON.isChecked():
            try:
                writeCHEMCON(findCHEMCON())
                self.CHEMCONStatusLab.setText('CHEMCON added and xd.mas updated')
            except PermissionError:
                self.CHEMCONStatusLab.setText(self.permErrorMsg)
            except Exception:
                self.CHEMCONStatusLab.setText('An error occurred.')
            
        elif self.chooseElementCHEMCON.isChecked() == True:
            
            inputText = str(self.inputElementCHEMCON.text().upper())
            inputElementList = labels2list(inputText)
            try:
                writeCHEMCON(findCHEMCONbyInputElement(inputElementList))
                self.CHEMCONStatusLab.setText('{} {}'.format('CHEMCON added and xd.mas updated for elements', ', '.join(inputElementList).strip(', ')))
            except PermissionError:
                self.CHEMCONStatusLab.setText(self.permErrorMsg)
            except Exception:
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
                inputAtomList = rawInput2Labels(inputText)
                
                atomLabs = copy.copy(globAtomLabs)  
                atomLabs = atomLabs.keys()

                inputAtomList = sorted(list(set(inputAtomList)))
                wrongLabs = [atom for atom in inputAtomList if atom not in atomLabs]

                inputAtomList = [atom for atom in inputAtomList if atom not in wrongLabs]

                chemconAddedStr = ', '.join([lab for lab in inputAtomList]).strip(', ')
                        
                try:
                    statusLabStr = ''
                    if inputAtomList:
                        writeCHEMCON(findCHEMCONbyInputAtoms(inputAtomList))
                        statusLabStr = ('{} {}'.format('CHEMCON added and xd.mas updated for', chemconAddedStr))
                    if wrongLabs:
                        statusLabStr += '<br><br>{}{}'.format('\nFollowing atom labels are incorrect: ', ', '.join(wrongLabs).strip(', '))
                    self.CHEMCONStatusLab.setText(statusLabStr)
                except PermissionError:
                    self.CHEMCONStatusLab.setText(self.permErrorMsg)
                except Exception:
                    self.CHEMCONStatusLab.setText('An error occurred.')
            self.inputAtomCHEMCON.setText('')
        self.changeUserIns()
        self.wizTest()

		

#--------------------RESULTS---------------------------------------------
        
    def makeResStr(self, lsmOutFile):
        '''
        Create string of formatted results from xd_lsm.out file. Return string.
        '''
        try:
            dmsda = getDMSDA(lsmOutFile)
            c = getConvergence(lsmOutFile)
            if c:
                conv = 'Yes'
            else:
                conv = 'No'
                
            resStr = 'RF<sup>2</sup> = {0: .2f} %<br>Convergence - {1}<br>Average DMSDA = {2}<br>Max DMSDA = {3}'.format(getRF2(), conv, dmsda[0], dmsda[1]) 
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
        except Exception:
            try:
                resStr = 'RF<sup>2</sup> = {0}<br>{1}<br>No DMSDA results found.'.format(getRF2(), getConvergence('xd_lsm.out'))
            except Exception:
                resStr = 'An error occurred.'
        
        return resStr
    
        
    def showResBackup(self):
        '''
        Show results from backup folder.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Choose backup folder to load results from"))
        
        if folder:
            resStr = self.makeResStr(folder + '/xd_lsm.out')

        self.resBackupLab.setText(resStr)
        
        
    def resBackupSum(self):
        '''
        Show wizard style summary of multiple refinements, from backup folder.
        '''
        folder = str(QFileDialog.getExistingDirectory(None, "Choose backup folder to load results from"))
        results = []
        resStr = '<br>{0:47}{1:23}{2:14}{3:16}{4}<br>'.format(' Refinement', 'RF<sup>2</sup>', 'Convergence', 'Average DMSDA', 'Max DMSDA')
        
        for item in os.listdir(folder):
            c = ''
            rf2 = getRF2(folder + '/' + item)
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
        self.resBackupLab.setText(resStr)
        
        
    def getResPress(self):
        '''
        Get results from xd_lsm.out and show them in results tab.
        '''
        resStr = self.makeResStr('xd_lsm.out')
        self.getResLab.setText(resStr)
    
    
    def FFT2FOUR(self):
        '''
        Run XDFOUR for atom nearest the largest peak in xd_fft.out.
        '''
        try:
            atom = FFTDetective()
            FOU3atoms(atom[0])
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
        self.resmap = resmap()
        self.resmap.show()
        
    def addFOURIns(self):
        '''
        Handle 'Run XDFOUR' button press in 'Results' tab.
        '''
        if self.setupFOURInputBox.isChecked() == True:
            rawInput = str(self.setupFOURInputText.text()).strip().upper()
            atom = rawInput2Labels(rawInput)[0]
            try:
                FOU3atoms(atom)
                self.xdfour.finishedSignal.connect(self.showResDensMap)
                self.xdfour.start()
            except PermissionError:
                self.setupFOURStatusLab.setText(self.permErrorMsg)
            except FileNotFoundError:
                self.setupFOURStatusLab.setText(self.lstMissingErrorMsg)
            except Exception:
                self.setupFOURStatusLab.setText('An error occurred.')
                    
        elif self.setupFOURFFTBox.isChecked() == True:
            try:
                self.xdfft.start()
                self.xdfft.finishedSignal.connect(self.FFT2FOUR)  
                
            except Exception:
                self.setupFOURStatusLab.setText('An error occurred.')

        elif self.quickplotGrdBox.isChecked():
            try:
                self.setupFOURStatusLab.setText('')
                self.showResDensMap()
            except Exception:
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


#--------------------BACKUP------------------------------------------------
    
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
                self.loadBackupLab.setText('Backup loaded from: "' + folder + '"') 
            except PermissionError:
                self.loadBackupLab.setText(self.permErrorMsg)
            except Exception:
                self.loadBackupLab.setText('An error occurred.')
        else:
            pass
        
        self.changeUserIns()
     

#--------------------MISC---------------------------------------------

    def openAbout(self):
        '''
        Open About XD Toolkit window.
        '''
        self.about = aboutBox()
        self.about.show()
        
        
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
       
def customExceptHook(Type, value, traceback):
    print_exception(Type, value, traceback)
    pass

##Run GUI
if __name__ == '__main__':
    
    sys.excepthook = customExceptHook               #Accept any errors so GUI doesn't quit.
    app = QApplication(sys.argv)
    
    #Splash screen
    splash_pix = QPixmap('res/flatearth.png')
    splash = QSplashScreen(splash_pix)
    splash.setMask(splash_pix.mask())
    font = QFont()
    font.setFamily("Bitstream Vera Sans Mono")
    splash.setFont(font)
    splash.showMessage('Initializing...',
                           Qt.AlignBottom | Qt.AlignLeft,
                           Qt.white)
    splash.show()
    
    prog = XDToolGui()
    app.aboutToQuit.connect(app.deleteLater)  #Fixes Anaconda bug where program only works on every second launch
    splash.finish(prog)
    prog.show()
    
    sys.exit(app.exec_())

#os.chdir('/home/matt/dev/XDTstuff/test/data/urea')
#x = ins2all()
#multipoleMagician()
#x = ins2all()
#findCHEMCON()
#multipoleMagician()
#x = makeLCS()
#print(x)


  



