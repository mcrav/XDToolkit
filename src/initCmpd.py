'''
#####################################################################
#-------------------COMPOUND INTIIALIZATION--------------------------
#####################################################################
'''

import numpy as np
import itertools
import os
import hashlib
from ast import literal_eval
from copy import copy
from shutil import rmtree
from devTools import timeDec
from asym2unit import applySymOps

def ins2fracPos(insFile):
    '''
    Get fractonal coordinates of every atom in shelx.ins. Return these and unit cell parameters.
    '''
    with open(insFile,'r') as ins:
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

                            combo1InvSym = copy(pos1InvSym)
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
