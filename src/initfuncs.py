from utils import coords2tuple
from asym2unit import applySymOps
from devtools import timeDec, atomsInPair
from databank import covradii
import os
import numpy as np
import itertools
import copy

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

    d2 = delta1**2 + delta2**2 + delta3**2 + (2*delta1*delta2*np.cos(gamma)) + (2*delta1*delta3*np.cos(beta)) + (2*delta2*delta3*np.cos(alpha))

    bondDist = np.sqrt(d2)

    return bondDist

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

def getAtomPairs(atomList):
    '''
    Make all pair combinations of atoms, given a list of atom labels. Return list.
    '''
    atomPairs = []

    for atom in atomList:
        for atom2 in atomList:
            if (atom, atom2) not in atomPairs:
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

def findBondingAtoms(unitCellParams, atomPairs, atomPos, specAtomPos, trackBondAtoms = None):
    '''
    Find atoms within bonding distance.
    Return...
    '''
    specDistances = {}
    distances = {}
    neebPairs = []
    neebSpecialPairs = []
    i=0                 #i is counter for creating .combo1 type labels

    for pair in atomPairs:

        atom1 = pair[0]
        atom2 = pair[1]

        pos1SpecialLab = ''
        pos2SpecialLab = ''

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
                    distances[frozenset(pair)] = round(bondDist, 4)
                    specDistances[frozenset([pos1SpecialLab, pos2SpecialLab])] = round(bondDist, 4)
                    neebPairs.append(pair)
                    neebSpecialPairs.append((pos1SpecialLab, pos2SpecialLab))

                else:
                    #Get all x,y,z+-1 combos
                    for combo in getCombos(atom1c):
                        bondDist = getBondDist(combo[0], atom2c, *unitCellParams)
                        #Testing
                        if trackBondAtoms:
                            trackBond(trackBondAtoms, atom1, atom2, bondDist, cutoffDist)

                        if bondDist < cutoffDist and bondDist !=0:
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

    return (specDistances, specAtomPos, neebPairs, neebSpecialPairs, distances)

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

                    specDistances[frozenset([newLab, splitLab[0] + ',asym'])] = round(getBondDist(newPos, specAtomPos[splitLab[0] + ',asym'][0], *unitCellParams), 4)
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
            specAngles[atom2].append((round(angle,2), atomsInAngle))

    return specAngles

@timeDec
def ins2all(trackBondAtoms = None):
    '''
    Find information about compound from shelx.ins input and calculations.
    Return nearest neighbours, atomic positions, bond distances, and bond angles.
    '''
    print('ins2all\n~~~~~~~~~~~~~~~~~~~~~~~~')
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
