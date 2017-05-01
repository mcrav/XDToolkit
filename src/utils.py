'''
#####################################################################
#-------------------UTILITIES----------------------------------------
#####################################################################
'''

import os

def atomTableEnds(line):
    return (line.startswith('END ATOM') or line.startswith('!') or line.startswith('DUM'))

def atomTableBegins(line):
    return (line.startswith('ATOM     ATOM0'))

def listjoin(listItem, divider):
    return divider.join(listItem).strip(divider)

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

def convert2XDLabel(atomLabel):
    '''
    Convert C1 label to C(1) label
    '''
    if atomLabel[:2].isalpha():
        newLabel = (atomLabel[:2] + '(' + atomLabel[2:] + ')').upper()
    else:
        newLabel = (atomLabel[:1] + '(' + atomLabel[1:] + ')').upper()

    return newLabel

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

def rawInput2labels(rawInput):
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

def getCellParams():
    '''
    Get unit cell parameters from shelx.ins or xd.mas.
    '''
    a = 0
    b = 0
    c = 0
    alpha = 0
    beta = 0
    gamma = 0
    
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
    return int(i)

def res2inp():
    '''
    Rename xd.res to xd.inp.
    '''
    os.remove('xd.inp')
    os.rename('xd.res','xd.inp')

def getEleNum():
    '''
    Get the number of different elements in the compound. Return number.
    '''
    with open('xd.mas','r') as mas:

        scat = False
        i = 0

        for line in mas:
            if line.startswith('END SCAT'):
                break

            if scat:
                i += 1

            if line.startswith('SCAT'):
                scat = True

    return i

def findMasCHEMCON():
    '''
    Get CHEMCON from mas file. Return as dictionary of children and their parents.
    '''
    with open('xd.mas','r') as mas:

        atomTab = False
        chemcon = {}

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.upper().split()
                if len(row) == 13:
                    chemcon[row[0]] = row[12]

            if atomTableBegins(line):
                atomTab = True

    return chemcon

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

def getAtomList():
    '''
    Get list of atoms from xd.mas. Return list.
    '''
    with open('xd.mas','r') as mas:  
        
        atomTab = False
        atoms = []
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.split()
                atoms.append(row[0].upper())

            if atomTableBegins(line):
                atomTab = True
                
    return atoms


