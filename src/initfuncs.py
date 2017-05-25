from utils import coords2tuple
from asym2unit import applySymOps
from devtools import timeDec, atomsInPair
from databank import covradii
import os
import numpy as np
import itertools
import copy
import multiprocessing as mp
from functools import partial
#from time import time

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

def getBondDist(atom1c, atom2c, a, b, c, alpha, beta, gamma):
    '''
    Get distance between 2 atoms. Return distance.
    '''
    delta1 = a*(atom1c[0] - atom2c[0])
    delta2 = b*(atom1c[1] - atom2c[1])
    delta3 = c*(atom1c[2] - atom2c[2])

    return np.sqrt(delta1**2 + delta2**2 + delta3**2 + (2*delta1*delta2*np.cos(gamma)) + (2*delta1*delta3*np.cos(beta)) + (2*delta2*delta3*np.cos(alpha)))

def getCombos(coord, xSymOps, ySymOps, zSymOps):
    '''
    Create all combinations of x, x+1, x-1 for x,y,z. Return combinations.
    '''
    combos = []

    xs = [coord[0], coord[0]-1, coord[0]+1]
    ys = [coord[1], coord[1]-1, coord[1]+1]
    zs = [coord[2], coord[2]-1, coord[2]+1]

    for i,xval in enumerate(xs):
        for j,yval in enumerate(ys):
            for k,zval in enumerate(zs):
                combos.append((np.array([xval, yval, zval]), (xSymOps[i], ySymOps[j], zSymOps[k])))

    return(combos)

def getCombosMp(xSymOps, ySymOps, zSymOps, coord):
    '''
    Create all combinations of x, x+1, x-1 for x,y,z. Return combinations.
    '''
    combos = []

    xs = [coord[0], coord[0]-1, coord[0]+1]
    ys = [coord[1], coord[1]-1, coord[1]+1]
    zs = [coord[2], coord[2]-1, coord[2]+1]

    for i,xval in enumerate(xs):
        for j,yval in enumerate(ys):
            for k,zval in enumerate(zs):
                combos.append((np.array([xval, yval, zval]), (xSymOps[i], ySymOps[j], zSymOps[k])))

    return(tuple(coord), combos)

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
    roundAngle = np.round(degAngle, 2)
    return roundAngle

def getAtomPairs(atomList):
    '''
    Make all pair combinations of atoms, given a list of atom labels. Return list.
    '''
    atomPairs = []

    for atom in atomList:
        for atom2 in atomList:
            if (atom, atom2) not in atomPairs and (atom2, atom) not in atomPairs:
                atomPairs.append((atom, atom2))

    return atomPairs

def trackBond(trackBond, atom1, atom2, bondDist, cutoffDist):
    if atomsInPair(trackBond, (atom1, atom2)) and bondDist < 5:
        print(atom1)
        print(atom2)
        print(bondDist)
        print(cutoffDist)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~')

def getCutoffDist(atom1, atom2):
    return covradii[atom1[:2].strip('(')] + covradii[atom2[:2].strip('(')] + 0.5

def makeNeebLabDict(neebPairs):
    '''
    Create dictionary of atoms and their nearest neighbours given a list of bonded pairs. Return dictionary.
    '''
    neebLabs = {}
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

    return neebLabs

def makeNeebSpecLabDict(neebSpecialPairs):
    '''
    Create dictionary of atoms and their nearest neighbours all with 'C(1),asym.combo3' type labels.
    Return dictionary.
    '''
    neebSpecLabs = {}

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

    return neebSpecLabs

def makeNeebTypeDict(neebLabDict):
    '''
    Create dictionary of atoms and their nearest neighbour types, from dictionary of nearest neighbour labels.
    Return dictionary.
    '''
    neebTypes = {}

    for atom, neebs in neebLabDict.items():
        neebTypes[atom] = []
        for neeb in neebs:
            neebTypes[atom].append(neeb[1][:2].split(',')[0].strip('('))

    return neebTypes

def makeMetricMatrix(a, b, c, alpha, beta, gamma):
    return np.array([[a**2, a*b*np.cos(gamma), a*c*np.cos(beta)],
                    [b*a*np.cos(gamma), b**2, b*c*np.cos(alpha)],
                    [c*a*np.cos(beta), c*b*np.cos(alpha), c**2]])
    
@timeDec
def findBondingAtoms(unitCellParams, atomPairs, atomPos, specAtomPos, trackBondAtoms = None):
    '''
    Find atoms within bonding distance.
    Return...
    '''
    xSymOps = ['x','x+1','x-1']
    ySymOps = ['y','y+1','y-1']
    zSymOps = ['z','z+1','z-1']
    specDistances = {}
    distances = {}
    atomList = []
    neebPairs = []
    neebSpecialPairs = []
    i=0                 #i is counter for creating .combo1 type labels
    comboDict = {}
    
    for pair in atomPairs:
        if pair[0] not in atomList:
            atomList.append(pair[0])
        if pair[1] not in atomList:
            atomList.append(pair[1])    
    
    for atom in atomList:
        for pos in atomPos[atom]:
            comboDict[tuple(pos[0])] = getCombos(pos[0], xSymOps, ySymOps, zSymOps)

    for pair in atomPairs:

        atom1 = pair[0]
        atom2 = pair[1]

        for pos in atomPos[atom1]:

            atom1c = np.array((pos[0][0], pos[0][1], pos[0][2]))
            pos1SpecialLab = pos[1]
            pos1InvSym = specAtomPos[pos1SpecialLab][1]

            for pos2 in atomPos[atom2]:

                atom2c = np.array((pos2[0][0], pos2[0][1], pos2[0][2]))
                pos2SpecialLab = pos2[1]

                bondDist = getBondDist(atom1c, atom2c, *unitCellParams)
                cutoffDist = getCutoffDist(atom1, atom2)
                #Testing
                if trackBondAtoms:
                    trackBond(trackBondAtoms, atom1, atom2, bondDist, cutoffDist)

                if bondDist < cutoffDist and bondDist != 0:
                    bondDist = np.round(bondDist,4)
                    distances[frozenset(pair)] = bondDist
                    specDistances[frozenset([pos1SpecialLab, pos2SpecialLab])] = bondDist
                    neebPairs.append(pair)
                    neebSpecialPairs.append((pos1SpecialLab, pos2SpecialLab))
                    break

                else:
                    #Get all x,y,z+-1 combos
                    #for combo in getCombos(atom1c, xSymOps, ySymOps, zSymOps):
                    for combo in comboDict[tuple(pos[0])]:
                        
                        bondDist = getBondDist(combo[0], atom2c, *unitCellParams)
                        #Testing
                        if trackBondAtoms:
                            trackBond(trackBondAtoms, atom1, atom2, bondDist, cutoffDist)

                        if bondDist < cutoffDist and bondDist !=0:
                            bondDist = np.round(bondDist, 4)

                            combo1SpecialLab = atom1 + ',combo' + str(i)
                            i+=1
                            #OPTIMIZATION COMMENT distances[frozenset(atom1, atom2)] = bondDist
                            specDistances[frozenset([combo1SpecialLab, pos2SpecialLab])] = bondDist

                            neebPairs.append((atom1, atom2))
                            neebSpecialPairs.append((combo1SpecialLab, pos2SpecialLab))

                            combo1InvSym = copy.copy(pos1InvSym)
                            combo1InvSym.insert(0, combo[1])

                            specAtomPos[combo1SpecialLab] = (combo[0], combo1InvSym)
                            break
                
    return (specDistances, specAtomPos, neebPairs, neebSpecialPairs)

def makeAsymNeebDict(neebLabs, neebSpecialPairs, specAtomPos, specDistances, unitCellParams):
    '''
    Make dictionary of atoms in asymmetric unit and their nearest neighbours with 'C(1),asym.combo3' type labels.
    Update specAtomPos and specDistances with new positions and bond distances.
    Return asymNeebs, specAtomPos and specDistances.
    '''
    asymNeebs = {}
    i=0
    addedPosDict = {atom: [] for atom in neebLabs.keys()}

    for pair in neebSpecialPairs:
        j=0
        while j < 2:
            splitLab = pair[j].split(',')
            tupPos = coords2tuple(specAtomPos[pair[j-1]][0])
            if splitLab[1] == 'asym':

                addedPos = addedPosDict[splitLab[0]]

                if tupPos not in addedPos:
                    asymNeebs.setdefault(splitLab[0],[]).append((specAtomPos[pair[j-1]][0], pair[j-1]))
                    addedPosDict[splitLab[0]].append(tupPos)

            else:
                invSymOp = specAtomPos[pair[j]][1]              #Get inv symOp of asym atom
                pos = specAtomPos[pair[j-1]][0]                   #Position of atom not asym

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

                    newLab = pair[j-1] +'.dum' + str(i)
                    asymNeebs.setdefault(splitLab[0], []).append((newPos, newLab))     #(pos, specLabel)
                    specAtomPos[newLab] = (newPos, invSymOp)

                    specDistances[frozenset([newLab, splitLab[0] + ',asym'])] = np.round(getBondDist(newPos, specAtomPos[splitLab[0] + ',asym'][0], *unitCellParams), 4)
                    addedPosDict[splitLab[0]].append(tupPos)
                    i+=1
            j+=1

    return (asymNeebs, specAtomPos, specDistances)

def getBondAngles(asymNeebs, specAtomPos, specDistances, metricMatrix):
    '''
    Get all bond angles in the structure. Return bond angles in dictionary.
    '''
    specAngles = {}

    for atom, neighbours in asymNeebs.items():
        atom2 = atom
        neighbourLabs = [neighbour[1] for neighbour in neighbours]

        specAngles[atom2] = []
        pairs  = itertools.combinations(neighbourLabs, 2)

        for pair in pairs:
            atomsInAngle = (pair[0],atom2 + ',asym',pair[1])
            angle = atoms2angle(atomsInAngle, specAtomPos, specDistances, metricMatrix)
            specAngles[atom2].append((np.round(angle,2), atomsInAngle))

    return specAngles

@timeDec
def ins2all(trackBondAtoms = None):
    '''
    Find information about compound from shelx.ins input and calculations.
    Return nearest neighbours, atomic positions, bond distances, and bond angles.
    '''
    #print('ins2all\n~~~~~~~~~~~~~~~~~~~~~~~~')
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
    unitCellParams = (a, b, c, alpha, beta, gamma)
    
    #Make all possible pairs of atoms in atomPairs
    atomPairs = getAtomPairs(atomPos.keys())

    #Find bonded atom pairs
    findBondsRes = findBondingAtoms(unitCellParams, atomPairs, atomPos, specAtomPos)
    specDistances = findBondsRes[0]
    specAtomPos = findBondsRes[1]
    neebPairs = findBondsRes[2]
    neebSpecialPairs = findBondsRes[3]
    
    #Get neighbours label dictionary
    neebLabs = makeNeebLabDict(neebPairs)

    #Get neighbour label dictionary of actual neighbours with special labels
    #neebSpecLabs = makeNeebSpecLabDict(neebSpecialPairs)

    # #Generate all neighbours of parent atoms around that atom positions in asym unit
    asymNeebsRes = makeAsymNeebDict(neebLabs, neebSpecialPairs, specAtomPos, specDistances, unitCellParams)
    asymNeebs = asymNeebsRes[0]
    specAtomPos = asymNeebsRes[1]
    specDistances = asymNeebsRes[2]

    #Get neighbour type dictionary
    neebTypes = makeNeebTypeDict(asymNeebs)

    metricMatrix = makeMetricMatrix(*unitCellParams)

    #Work out angles based on positions generated by symmetry
    specAngles = getBondAngles(asymNeebs, specAtomPos, specDistances, metricMatrix)

    #Remove atom positions from atomLab dict
    atomLabs = {atom : [neeb[1] for neeb in neighbours] for atom, neighbours in asymNeebs.items()}
    
    return (atomLabs, neebTypes, specAngles, specDistances, specAtomPos)

#print(getBondDist(np.array([1,0,0]), np.array([0.9,-0.5,0]), 1, 1, 1, 90, 90, 90))
            
            
#carbaTest = '''({'O(0)': ['C(A),asym'], 'C(A)': ['O(0),asym', 'N(0),asym', 'N(A),asym'], 'N(0)': ['C(A),asym', 'C(1),asym', 'C(14),asym'], 'C(1)': ['N(0),asym', 'C(6),asym', 'C(2),asym'], 'C(14)': ['N(0),asym', 'C(9),asym', 'C(13),asym'], 'N(A)': ['H(AA),asym', 'C(A),asym', 'H(AB),asym'], 'H(AA)': ['N(A),asym'], 'H(AB)': ['N(A),asym'], 'C(6)': ['C(1),asym', 'C(7),asym', 'C(5),asym'], 'C(2)': ['C(1),asym', 'H(2),asym', 'C(3),asym'], 'C(9)': ['C(14),asym', 'C(10),asym', 'C(8),asym'], 'C(13)': ['C(14),asym', 'H(13),asym', 'C(12),asym'], 'C(7)': ['C(6),asym', 'H(7),asym', 'C(8),asym'], 'C(5)': ['C(6),asym', 'H(5),asym', 'C(4),asym'], 'H(2)': ['C(2),asym'], 'C(3)': ['C(2),asym', 'H(3),asym', 'C(4),asym'], 'H(7)': ['C(7),asym'], 'C(8)': ['C(7),asym', 'C(9),asym', 'H(8),asym'], 'C(10)': ['C(9),asym', 'H(10),asym', 'C(11),asym'], 'H(5)': ['C(5),asym'], 'C(4)': ['C(5),asym', 'C(3),asym', 'H(4),asym'], 'H(13)': ['C(13),asym'], 'C(12)': ['C(13),asym', 'C(11),asym', 'H(12),asym'], 'H(3)': ['C(3),asym'], 'H(10)': ['C(10),asym'], 'C(11)': ['C(10),asym', 'H(11),asym', 'C(12),asym'], 'H(8)': ['C(8),asym'], 'H(11)': ['C(11),asym'], 'H(4)': ['C(4),asym'], 'H(12)': ['C(12),asym']}, {'O(0)': ['C'], 'C(A)': ['O', 'N', 'N'], 'N(0)': ['C', 'C', 'C'], 'C(1)': ['N', 'C', 'C'], 'C(14)': ['N', 'C', 'C'], 'N(A)': ['H', 'C', 'H'], 'H(AA)': ['N'], 'H(AB)': ['N'], 'C(6)': ['C', 'C', 'C'], 'C(2)': ['C', 'H', 'C'], 'C(9)': ['C', 'C', 'C'], 'C(13)': ['C', 'H', 'C'], 'C(7)': ['C', 'H', 'C'], 'C(5)': ['C', 'H', 'C'], 'H(2)': ['C'], 'C(3)': ['C', 'H', 'C'], 'H(7)': ['C'], 'C(8)': ['C', 'C', 'H'], 'C(10)': ['C', 'H', 'C'], 'H(5)': ['C'], 'C(4)': ['C', 'C', 'H'], 'H(13)': ['C'], 'C(12)': ['C', 'C', 'H'], 'H(3)': ['C'], 'H(10)': ['C'], 'C(11)': ['C', 'H', 'C'], 'H(8)': ['C'], 'H(11)': ['C'], 'H(4)': ['C'], 'H(12)': ['C']}, {'O(0)': [], 'C(A)': [(121.61, ('O(0),asym', 'C(A),asym', 'N(0),asym')), (122.28, ('O(0),asym', 'C(A),asym', 'N(A),asym')), (116.11, ('N(0),asym', 'C(A),asym', 'N(A),asym'))], 'N(0)': [(120.93000000000001, ('C(A),asym', 'N(0),asym', 'C(1),asym')), (121.55, ('C(A),asym', 'N(0),asym', 'C(14),asym')), (117.19, ('C(1),asym', 'N(0),asym', 'C(14),asym'))], 'C(1)': [(119.12, ('N(0),asym', 'C(1),asym', 'C(6),asym')), (120.11, ('N(0),asym', 'C(1),asym', 'C(2),asym')), (120.76000000000001, ('C(6),asym', 'C(1),asym', 'C(2),asym'))], 'C(14)': [(119.76000000000001, ('N(0),asym', 'C(14),asym', 'C(9),asym')), (119.2, ('N(0),asym', 'C(14),asym', 'C(13),asym')), (120.97, ('C(9),asym', 'C(14),asym', 'C(13),asym'))], 'N(A)': [(109.92, ('H(AA),asym', 'N(A),asym', 'C(A),asym')), (130.18000000000001, ('H(AA),asym', 'N(A),asym', 'H(AB),asym')), (119.44, ('C(A),asym', 'N(A),asym', 'H(AB),asym'))], 'H(AA)': [], 'H(AB)': [], 'C(6)': [(123.59999999999999, ('C(1),asym', 'C(6),asym', 'C(7),asym')), (118.04000000000001, ('C(1),asym', 'C(6),asym', 'C(5),asym')), (118.36, ('C(7),asym', 'C(6),asym', 'C(5),asym'))], 'C(2)': [(119.89, ('C(1),asym', 'C(2),asym', 'H(2),asym')), (120.22, ('C(1),asym', 'C(2),asym', 'C(3),asym')), (119.89, ('H(2),asym', 'C(2),asym', 'C(3),asym'))], 'C(9)': [(117.68000000000001, ('C(14),asym', 'C(9),asym', 'C(10),asym')), (123.20999999999999, ('C(14),asym', 'C(9),asym', 'C(8),asym')), (119.09, ('C(10),asym', 'C(9),asym', 'C(8),asym'))], 'C(13)': [(119.86, ('C(14),asym', 'C(13),asym', 'H(13),asym')), (120.28, ('C(14),asym', 'C(13),asym', 'C(12),asym')), (119.86, ('H(13),asym', 'C(13),asym', 'C(12),asym'))], 'C(7)': [(116.06999999999999, ('C(6),asym', 'C(7),asym', 'H(7),asym')), (127.84999999999999, ('C(6),asym', 'C(7),asym', 'C(8),asym')), (116.06999999999999, ('H(7),asym', 'C(7),asym', 'C(8),asym'))], 'C(5)': [(119.33, ('C(6),asym', 'C(5),asym', 'H(5),asym')), (121.34999999999999, ('C(6),asym', 'C(5),asym', 'C(4),asym')), (119.33, ('H(5),asym', 'C(5),asym', 'C(4),asym'))], 'H(2)': [], 'C(3)': [(120.0, ('C(2),asym', 'C(3),asym', 'H(3),asym')), (120.0, ('C(2),asym', 'C(3),asym', 'C(4),asym')), (120.0, ('H(3),asym', 'C(3),asym', 'C(4),asym'))], 'H(7)': [], 'C(8)': [(126.17, ('C(7),asym', 'C(8),asym', 'C(9),asym')), (116.92, ('C(7),asym', 'C(8),asym', 'H(8),asym')), (116.92, ('C(9),asym', 'C(8),asym', 'H(8),asym'))], 'C(10)': [(119.25, ('C(9),asym', 'C(10),asym', 'H(10),asym')), (121.51000000000001, ('C(9),asym', 'C(10),asym', 'C(11),asym')), (119.23999999999999, ('H(10),asym', 'C(10),asym', 'C(11),asym'))], 'H(5)': [], 'C(4)': [(119.62, ('C(5),asym', 'C(4),asym', 'C(3),asym')), (120.19, ('C(5),asym', 'C(4),asym', 'H(4),asym')), (120.19, ('C(3),asym', 'C(4),asym', 'H(4),asym'))], 'H(13)': [], 'C(12)': [(119.44, ('C(13),asym', 'C(12),asym', 'C(11),asym')), (120.28, ('C(13),asym', 'C(12),asym', 'H(12),asym')), (120.28, ('C(11),asym', 'C(12),asym', 'H(12),asym'))], 'H(3)': [], 'H(10)': [], 'C(11)': [(119.98, ('C(10),asym', 'C(11),asym', 'H(11),asym')), (120.05, ('C(10),asym', 'C(11),asym', 'C(12),asym')), (119.97, ('H(11),asym', 'C(11),asym', 'C(12),asym'))], 'H(8)': [], 'H(11)': [], 'H(4)': [], 'H(12)': []}, {frozenset({'C(A),asym', 'O(0),asym'}): 1.24, frozenset({'C(A),1', 'O(0),1'}): 1.24, frozenset({'C(A),3', 'O(0),3'}): 1.24, frozenset({'O(0),4', 'C(A),4'}): 1.24, frozenset({'C(A),asym', 'N(0),asym'}): 1.3817999999999999, frozenset({'C(A),1', 'N(0),1'}): 1.3817999999999999, frozenset({'C(A),3', 'N(0),3'}): 1.3817999999999999, frozenset({'N(0),4', 'C(A),4'}): 1.3817999999999999, frozenset({'C(1),asym', 'N(0),asym'}): 1.4322999999999999, frozenset({'N(0),1', 'C(1),1'}): 1.4322999999999999, frozenset({'C(1),3', 'N(0),3'}): 1.4322999999999999, frozenset({'C(1),4', 'N(0),4'}): 1.4322999999999999, frozenset({'C(14),asym', 'N(0),asym'}): 1.4331, frozenset({'C(14),1', 'N(0),1'}): 1.4331, frozenset({'C(14),3', 'N(0),3'}): 1.4331, frozenset({'N(0),4', 'C(14),4'}): 1.4331, frozenset({'N(A),asym', 'H(AA),asym'}): 0.88360000000000005, frozenset({'N(A),1', 'H(AA),1'}): 0.88360000000000005, frozenset({'H(AA),3', 'N(A),3'}): 0.88360000000000005, frozenset({'H(AA),4', 'N(A),4'}): 0.88360000000000005, frozenset({'N(A),asym', 'C(A),asym'}): 1.3574999999999999, frozenset({'C(A),1', 'N(A),1'}): 1.3574999999999999, frozenset({'C(A),3', 'N(A),3'}): 1.3574999999999999, frozenset({'C(A),4', 'N(A),4'}): 1.3574999999999999, frozenset({'N(A),asym', 'H(AB),asym'}): 0.79790000000000005, frozenset({'N(A),1', 'H(AB),1'}): 0.79790000000000005, frozenset({'H(AB),3', 'N(A),3'}): 0.79790000000000005, frozenset({'H(AB),4', 'N(A),4'}): 0.79790000000000005, frozenset({'C(1),asym', 'C(6),asym'}): 1.4048, frozenset({'C(1),1', 'C(6),1'}): 1.4048, frozenset({'C(6),3', 'C(1),3'}): 1.4048, frozenset({'C(1),4', 'C(6),4'}): 1.4048, frozenset({'C(1),asym', 'C(2),asym'}): 1.3934, frozenset({'C(2),1', 'C(1),1'}): 1.3934, frozenset({'C(2),3', 'C(1),3'}): 1.3934, frozenset({'C(1),4', 'C(2),4'}): 1.3934, frozenset({'C(9),asym', 'C(14),asym'}): 1.4027000000000001, frozenset({'C(14),1', 'C(9),1'}): 1.4027000000000001, frozenset({'C(14),3', 'C(9),3'}): 1.4027000000000001, frozenset({'C(14),4', 'C(9),4'}): 1.4027000000000001, frozenset({'C(14),asym', 'C(13),asym'}): 1.3971, frozenset({'C(14),1', 'C(13),1'}): 1.3971, frozenset({'C(14),3', 'C(13),3'}): 1.3971, frozenset({'C(14),4', 'C(13),4'}): 1.3971, frozenset({'C(7),asym', 'C(6),asym'}): 1.4635, frozenset({'C(6),1', 'C(7),1'}): 1.4635, frozenset({'C(6),3', 'C(7),3'}): 1.4635, frozenset({'C(6),4', 'C(7),4'}): 1.4635, frozenset({'C(5),asym', 'C(6),asym'}): 1.4095, frozenset({'C(6),1', 'C(5),1'}): 1.4095, frozenset({'C(5),3', 'C(6),3'}): 1.4095, frozenset({'C(5),4', 'C(6),4'}): 1.4095, frozenset({'C(2),asym', 'H(2),asym'}): 0.94999999999999996, frozenset({'C(2),1', 'H(2),1'}): 0.94999999999999996, frozenset({'C(2),3', 'H(2),3'}): 0.94999999999999996, frozenset({'H(2),4', 'C(2),4'}): 0.94999999999999996, frozenset({'C(3),asym', 'C(2),asym'}): 1.3924000000000001, frozenset({'C(2),1', 'C(3),1'}): 1.3924000000000001, frozenset({'C(2),3', 'C(3),3'}): 1.3924000000000001, frozenset({'C(2),4', 'C(3),4'}): 1.3924000000000001, frozenset({'C(7),asym', 'H(7),asym'}): 0.94999999999999996, frozenset({'H(7),1', 'C(7),1'}): 0.94999999999999996, frozenset({'H(7),3', 'C(7),3'}): 0.94999999999999996, frozenset({'C(7),4', 'H(7),4'}): 0.94999999999999996, frozenset({'C(8),asym', 'C(7),asym'}): 1.3492, frozenset({'C(7),1', 'C(8),1'}): 1.3492, frozenset({'C(8),3', 'C(7),3'}): 1.3492, frozenset({'C(7),4', 'C(8),4'}): 1.3492, frozenset({'C(9),asym', 'C(10),asym'}): 1.4068000000000001, frozenset({'C(10),1', 'C(9),1'}): 1.4068000000000001, frozenset({'C(9),3', 'C(10),3'}): 1.4068000000000001, frozenset({'C(10),4', 'C(9),4'}): 1.4068000000000001, frozenset({'C(9),asym', 'C(8),asym'}): 1.4633, frozenset({'C(8),1', 'C(9),1'}): 1.4633, frozenset({'C(9),3', 'C(8),3'}): 1.4633, frozenset({'C(8),4', 'C(9),4'}): 1.4633, frozenset({'H(5),asym', 'C(5),asym'}): 0.94999999999999996, frozenset({'H(5),1', 'C(5),1'}): 0.94999999999999996, frozenset({'H(5),3', 'C(5),3'}): 0.94999999999999996, frozenset({'C(5),4', 'H(5),4'}): 0.94999999999999996, frozenset({'C(4),asym', 'C(5),asym'}): 1.3885000000000001, frozenset({'C(4),1', 'C(5),1'}): 1.3885000000000001, frozenset({'C(4),3', 'C(5),3'}): 1.3885000000000001, frozenset({'C(4),4', 'C(5),4'}): 1.3885000000000001, frozenset({'H(13),asym', 'C(13),asym'}): 0.94999999999999996, frozenset({'H(13),1', 'C(13),1'}): 0.94999999999999996, frozenset({'C(13),3', 'H(13),3'}): 0.94999999999999996, frozenset({'H(13),4', 'C(13),4'}): 0.94999999999999996, frozenset({'C(12),asym', 'C(13),asym'}): 1.3919999999999999, frozenset({'C(12),1', 'C(13),1'}): 1.3919999999999999, frozenset({'C(13),3', 'C(12),3'}): 1.3919999999999999, frozenset({'C(12),4', 'C(13),4'}): 1.3919999999999999, frozenset({'H(3),asym', 'C(3),asym'}): 0.94999999999999996, frozenset({'C(3),1', 'H(3),1'}): 0.94999999999999996, frozenset({'H(3),3', 'C(3),3'}): 0.94999999999999996, frozenset({'C(3),4', 'H(3),4'}): 0.94999999999999996, frozenset({'C(4),asym', 'C(3),asym'}): 1.3983000000000001, frozenset({'C(3),1', 'C(4),1'}): 1.3983000000000001, frozenset({'C(4),3', 'C(3),3'}): 1.3983000000000001, frozenset({'C(4),4', 'C(3),4'}): 1.3983000000000001, frozenset({'C(10),asym', 'H(10),asym'}): 0.94999999999999996, frozenset({'H(10),1', 'C(10),1'}): 0.94999999999999996, frozenset({'H(10),3', 'C(10),3'}): 0.94999999999999996, frozenset({'C(10),4', 'H(10),4'}): 0.94999999999999996, frozenset({'C(10),asym', 'C(11),asym'}): 1.3862000000000001, frozenset({'C(11),1', 'C(10),1'}): 1.3862000000000001, frozenset({'C(11),3', 'C(10),3'}): 1.3862000000000001, frozenset({'C(10),4', 'C(11),4'}): 1.3862000000000001, frozenset({'H(8),asym', 'C(8),asym'}): 0.94999999999999996, frozenset({'H(8),1', 'C(8),1'}): 0.94999999999999996, frozenset({'H(8),3', 'C(8),3'}): 0.94999999999999996, frozenset({'H(8),4', 'C(8),4'}): 0.94999999999999996, frozenset({'H(11),asym', 'C(11),asym'}): 0.94999999999999996, frozenset({'C(11),1', 'H(11),1'}): 0.94999999999999996, frozenset({'C(11),3', 'H(11),3'}): 0.94999999999999996, frozenset({'H(11),4', 'C(11),4'}): 0.94999999999999996, frozenset({'C(12),asym', 'C(11),asym'}): 1.3977999999999999, frozenset({'C(12),1', 'C(11),1'}): 1.3977999999999999, frozenset({'C(12),3', 'C(11),3'}): 1.3977999999999999, frozenset({'C(12),4', 'C(11),4'}): 1.3977999999999999, frozenset({'C(4),asym', 'H(4),asym'}): 0.94999999999999996, frozenset({'H(4),1', 'C(4),1'}): 0.94999999999999996, frozenset({'H(4),3', 'C(4),3'}): 0.94999999999999996, frozenset({'C(4),4', 'H(4),4'}): 0.94999999999999996, frozenset({'H(12),asym', 'C(12),asym'}): 0.94999999999999996, frozenset({'C(12),1', 'H(12),1'}): 0.94999999999999996, frozenset({'C(12),3', 'H(12),3'}): 0.94999999999999996, frozenset({'H(12),4', 'C(12),4'}): 0.94999999999999996}, {'O(0),asym': (array([-0.010687,  0.125029,  0.406055]), [('x', 'y', 'z')]), 'O(0),1': (array([ 0.510687,  0.625029,  0.093945]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'O(0),3': (array([ 0.489313,  0.374971,  0.906055]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(0),4': (array([ 0.010687, -0.125029, -0.406055]), [('-x', '-y', '-z')]), 'N(0),asym': (array([ 0.205975,  0.268633,  0.437844]), [('x', 'y', 'z')]), 'N(0),1': (array([ 0.294025,  0.768633,  0.062156]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'N(0),3': (array([ 0.705975,  0.231367,  0.937844]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'N(0),4': (array([-0.205975, -0.268633, -0.437844]), [('-x', '-y', '-z')]), 'N(A),asym': (array([ 0.180961,  0.100682,  0.538785]), [('x', 'y', 'z')]), 'N(A),1': (array([ 0.319039,  0.600682, -0.038785]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'N(A),3': (array([ 0.680961,  0.399318,  1.038785]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'N(A),4': (array([-0.180961, -0.100682, -0.538785]), [('-x', '-y', '-z')]), 'H(AA),asym': (array([ 0.111039,  0.037887,  0.549421]), [('x', 'y', 'z')]), 'H(AA),1': (array([ 0.388961,  0.537887, -0.049421]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(AA),3': (array([ 0.611039,  0.462113,  1.049421]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(AA),4': (array([-0.111039, -0.037887, -0.549421]), [('-x', '-y', '-z')]), 'C(A),asym': (array([ 0.118291,  0.161877,  0.458263]), [('x', 'y', 'z')]), 'C(A),1': (array([ 0.381709,  0.661877,  0.041737]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(A),3': (array([ 0.618291,  0.338123,  0.958263]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(A),4': (array([-0.118291, -0.161877, -0.458263]), [('-x', '-y', '-z')]), 'C(1),asym': (array([ 0.168697,  0.331654,  0.347829]), [('x', 'y', 'z')]), 'C(1),1': (array([ 0.331303,  0.831654,  0.152171]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(1),3': (array([ 0.668697,  0.168346,  0.847829]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(1),4': (array([-0.168697, -0.331654, -0.347829]), [('-x', '-y', '-z')]), 'C(14),asym': (array([ 0.350449,  0.313577,  0.500372]), [('x', 'y', 'z')]), 'C(14),1': (array([  1.49551000e-01,   8.13577000e-01,  -3.72000000e-04]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(14),3': (array([ 0.850449,  0.186423,  1.000372]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(14),4': (array([-0.350449, -0.313577, -0.500372]), [('-x', '-y', '-z')]), 'C(6),asym': (array([ 0.298932,  0.333164,  0.277616]), [('x', 'y', 'z')]), 'C(6),1': (array([ 0.201068,  0.833164,  0.222384]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(6),3': (array([ 0.798932,  0.166836,  0.777616]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(6),4': (array([-0.298932, -0.333164, -0.277616]), [('-x', '-y', '-z')]), 'C(2),asym': (array([ 0.005834,  0.391513,  0.331057]), [('x', 'y', 'z')]), 'C(2),1': (array([ 0.494166,  0.891513,  0.168943]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(2),3': (array([ 0.505834,  0.108487,  0.831057]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(2),4': (array([-0.005834, -0.391513, -0.331057]), [('-x', '-y', '-z')]), 'H(2),asym': (array([-0.080096,  0.391125,  0.379696]), [('x', 'y', 'z')]), 'H(2),1': (array([ 0.580096,  0.891125,  0.120304]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(2),3': (array([ 0.419904,  0.108875,  0.879696]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(2),4': (array([ 0.080096, -0.391125, -0.379696]), [('-x', '-y', '-z')]), 'C(7),asym': (array([ 0.474285,  0.274849,  0.291454]), [('x', 'y', 'z')]), 'C(7),1': (array([ 0.025715,  0.774849,  0.208546]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(7),3': (array([ 0.974285,  0.225151,  0.791454]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(7),4': (array([-0.474285, -0.274849, -0.291454]), [('-x', '-y', '-z')]), 'H(7),asym': (array([ 0.524443,  0.243814,  0.234333]), [('x', 'y', 'z')]), 'H(7),1': (array([-0.024443,  0.743814,  0.265667]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(7),3': (array([ 1.024443,  0.256186,  0.734333]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(7),4': (array([-0.524443, -0.243814, -0.234333]), [('-x', '-y', '-z')]), 'C(9),asym': (array([ 0.527304,  0.302774,  0.471944]), [('x', 'y', 'z')]), 'C(9),1': (array([-0.027304,  0.802774,  0.028056]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(9),3': (array([ 1.027304,  0.197226,  0.971944]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(9),4': (array([-0.527304, -0.302774, -0.471944]), [('-x', '-y', '-z')]), 'C(5),asym': (array([ 0.258142,  0.39445 ,  0.189123]), [('x', 'y', 'z')]), 'C(5),1': (array([ 0.241858,  0.89445 ,  0.310877]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(5),3': (array([ 0.758142,  0.10555 ,  0.689123]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(5),4': (array([-0.258142, -0.39445 , -0.189123]), [('-x', '-y', '-z')]), 'H(5),asym': (array([ 0.343596,  0.395267,  0.140237]), [('x', 'y', 'z')]), 'H(5),1': (array([ 0.156404,  0.895267,  0.359763]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(5),3': (array([ 0.843596,  0.104733,  0.640237]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(5),4': (array([-0.343596, -0.395267, -0.140237]), [('-x', '-y', '-z')]), 'C(13),asym': (array([ 0.313113,  0.362341,  0.591056]), [('x', 'y', 'z')]), 'C(13),1': (array([ 0.186887,  0.862341, -0.091056]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(13),3': (array([ 0.813113,  0.137659,  1.091056]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(13),4': (array([-0.313113, -0.362341, -0.591056]), [('-x', '-y', '-z')]), 'H(13),asym': (array([ 0.192352,  0.372359,  0.607952]), [('x', 'y', 'z')]), 'H(13),1': (array([ 0.307648,  0.872359, -0.107952]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(13),3': (array([ 0.692352,  0.127641,  1.107952]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(13),4': (array([-0.192352, -0.372359, -0.607952]), [('-x', '-y', '-z')]), 'C(3),asym': (array([-0.031392,  0.451897,  0.243245]), [('x', 'y', 'z')]), 'C(3),1': (array([ 0.531392,  0.951897,  0.256755]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(3),3': (array([ 0.468608,  0.048103,  0.743245]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(3),4': (array([ 0.031392, -0.451897, -0.243245]), [('-x', '-y', '-z')]), 'H(3),asym': (array([-0.143066,  0.4919  ,  0.231779]), [('x', 'y', 'z')]), 'H(3),1': (array([ 0.643066,  0.9919  ,  0.268221]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(3),3': (array([ 0.356934,  0.0081  ,  0.731779]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(3),4': (array([ 0.143066, -0.4919  , -0.231779]), [('-x', '-y', '-z')]), 'C(10),asym': (array([ 0.665549,  0.337877,  0.53947 ]), [('x', 'y', 'z')]), 'C(10),1': (array([-0.165549,  0.837877, -0.03947 ]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(10),3': (array([ 1.165549,  0.162123,  1.03947 ]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(10),4': (array([-0.665549, -0.337877, -0.53947 ]), [('-x', '-y', '-z')]), 'H(10),asym': (array([ 0.786532,  0.329826,  0.522409]), [('x', 'y', 'z')]), 'H(10),1': (array([-0.286532,  0.829826, -0.022409]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(10),3': (array([ 1.286532,  0.170174,  1.022409]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(10),4': (array([-0.786532, -0.329826, -0.522409]), [('-x', '-y', '-z')]), 'C(8),asym': (array([ 0.573326,  0.259828,  0.375515]), [('x', 'y', 'z')]), 'C(8),1': (array([-0.073326,  0.759828,  0.124485]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(8),3': (array([ 1.073326,  0.240172,  0.875515]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(8),4': (array([-0.573326, -0.259828, -0.375515]), [('-x', '-y', '-z')]), 'H(8),asym': (array([ 0.683353,  0.217352,  0.371748]), [('x', 'y', 'z')]), 'H(8),1': (array([-0.183353,  0.717352,  0.128252]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(8),3': (array([ 1.183353,  0.282648,  0.871748]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(8),4': (array([-0.683353, -0.217352, -0.371748]), [('-x', '-y', '-z')]), 'C(11),asym': (array([ 0.62928 ,  0.383945,  0.630331]), [('x', 'y', 'z')]), 'C(11),1': (array([-0.12928 ,  0.883945, -0.130331]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(11),3': (array([ 1.12928 ,  0.116055,  1.130331]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(11),4': (array([-0.62928 , -0.383945, -0.630331]), [('-x', '-y', '-z')]), 'H(11),asym': (array([ 0.724911,  0.40713 ,  0.674673]), [('x', 'y', 'z')]), 'H(11),1': (array([-0.224911,  0.90713 , -0.174673]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(11),3': (array([ 1.224911,  0.09287 ,  1.174673]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(11),4': (array([-0.724911, -0.40713 , -0.674673]), [('-x', '-y', '-z')]), 'C(4),asym': (array([ 0.095614,  0.453652,  0.17194 ]), [('x', 'y', 'z')]), 'C(4),1': (array([ 0.404386,  0.953652,  0.32806 ]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(4),3': (array([ 0.595614,  0.046348,  0.67194 ]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(4),4': (array([-0.095614, -0.453652, -0.17194 ]), [('-x', '-y', '-z')]), 'H(4),asym': (array([ 0.070842,  0.495111,  0.112086]), [('x', 'y', 'z')]), 'H(4),1': (array([ 0.429158,  0.995111,  0.387914]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(4),3': (array([ 0.570842,  0.004889,  0.612086]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(4),4': (array([-0.070842, -0.495111, -0.112086]), [('-x', '-y', '-z')]), 'C(12),asym': (array([ 0.452107,  0.396308,  0.656755]), [('x', 'y', 'z')]), 'C(12),1': (array([ 0.047893,  0.896308, -0.156755]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'C(12),3': (array([ 0.952107,  0.103692,  1.156755]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(12),4': (array([-0.452107, -0.396308, -0.656755]), [('-x', '-y', '-z')]), 'H(12),asym': (array([ 0.426765,  0.427632,  0.719032]), [('x', 'y', 'z')]), 'H(12),1': (array([ 0.073235,  0.927632, -0.219032]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(12),3': (array([ 0.926765,  0.072368,  1.219032]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(12),4': (array([-0.426765, -0.427632, -0.719032]), [('-x', '-y', '-z')]), 'H(AB),asym': (array([ 0.262283,  0.128288,  0.571565]), [('x', 'y', 'z')]), 'H(AB),1': (array([ 0.237717,  0.628288, -0.071565]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z')]), 'H(AB),3': (array([ 0.762283,  0.371712,  1.071565]), [['0.5-x', '-0.5+y', '0.5-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(AB),4': (array([-0.262283, -0.128288, -0.571565]), [('-x', '-y', '-z')])})'''