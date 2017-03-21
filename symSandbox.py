import os
import numpy as np
import copy
from XDToolkit import ins2all, ins2fracPos
    
def readINS():
    
    centeringList = {1:'P', 2:'I', 3:'R', 4:'F', 5:'A', 6:'B', 7:'C'}
    centrosym = True
    symOps = []
    
    with open('shelx.ins','r') as ins:
        
        for line in ins:
            if line.startswith('LATT'):
                row = str.split(line)
                if row[1].isdigit():
                    centering = [centeringList[int(item)] for item in row[1:]]
                else:
                    centering = [item for item in row[1:]]
                if row[1].startswith('-'):
                    centrosym = False
            
            elif line.startswith('SYMM'):

                symOp = line[4:].lower().split(',')
                symOp = [item.strip() for item in symOp]
                symOps.append(symOp)
            
            elif symOps and centering and not line.startswith('SYMM'):
                break
                
    return (centering, centrosym, symOps)

def genAtomsbySym(insInfo, x,y,z):
    centSymOps = []
    symOps = []
    
    LATT = insInfo[0]
    SYMM = insInfo[2]
    
    if 'P' in LATT:
        centSymOps.append(np.array([x, y, z]))

    if 'I' in LATT:
        centSymOps.append(np.array([0.5 + x, 0.5 + y, 0.5 + z]))
            
    if 'R' in LATT:
        pass
    
    if 'F' in LATT:
        centSymOps.append(np.array([x, 0.5 + y, 0.5 + z]))
        centSymOps.append(np.array([0.5 + x, y, 0.5 + z]))
        centSymOps.append(np.array([0.5 + x, 0.5 + y, z]))
        
    if 'A' in LATT:
        centSymOps.append(np.array([x, 0.5 + y, 0.5 + z]))
        
    if 'B' in LATT:
        centSymOps.append(np.array([0.5 + x, y, 0.5 + z]))
        
    if 'C' in LATT:
        centSymOps.append(np.array([0.5 + x, 0.5 + y, z]))   
    
    for insSymOp in SYMM:
        for symOp in centSymOps:
            x = symOp[0]
            y = symOp[1]
            z = symOp[2]
            
            symOps.append(np.array([eval(insSymOp[0]), eval(insSymOp[1]), eval(insSymOp[2])]))
        
    return symOps
        

def getSymOps(roundNum = 4):
    asymUnit = ins2fracPos('shelx.ins')[0]
    atoms = copy.copy(asymUnit)
    
    atoms2 = {}
    atoms3 = {}
    specAtoms = {}
    
    insInfo = readINS()
    print(insInfo)
    centering = insInfo[0]
    centrosym = insInfo[1]
    SYMM = insInfo[2]
    
    for lab, atomPos in atoms.items():
        
        i = 0

        x = atomPos[0][0]
        y = atomPos[0][1]
        z = atomPos[0][2]
        
        if lab not in atoms2.keys():
            atoms2[lab] = [(atomPos[0], lab + ',' + str(i))]
            i+=1
        
        for item in genAtomsbySym(insInfo, x,y,z):
            atoms2[lab].append(  (item, lab + ',' + str(i))  ) 
            i+=1
        
        for item in SYMM:
            coords = np.array([eval(item[0]), eval(item[1]), eval(item[2])])
            atoms2[lab].append((coords, lab + ',' + str(i)))
            i+=1
            
        if centrosym:
            
            for item in genAtomsbySym(insInfo, -x,-y,-z):

                atoms2[lab].append((item, lab + ',' + str(i)))   
                i+=1
            
            atoms2[lab].append((np.array([-x,-y,-z]), lab + ',' + str(i)))
            i+=1
            x = -x
            y = -y
            z = -z
            for item in SYMM:
                coords = np.array([eval(item[0]), eval(item[1]), eval(item[2])])
                atoms2[lab].append(  (coords, lab + ',' + str(i))  )

            
    usedPos = []    

    for atom, posList in atoms2.items():         #remove duplicates

        if atom not in atoms3.keys():
            atoms3[atom] = []
        for pos in posList:
            
            tupPos = (pos[0][0].round(2), pos[0][1].round(2), pos[0][2].round(2))

            if tupPos not in usedPos:
                usedPos.append(tupPos)
                atoms3[atom].append(pos)
            else:
                pass
            
    
        
        
    print(len(atoms.values()))
    print(getNumPos(atoms2))
    del atoms2
    
    atoms4 = catchEdgeAtoms(atoms3)
    
    for atom, posList in atoms4.items():
        for pos in posList:
            specAtoms[pos[1]] = pos[0]
    
    print(getNumPos(atoms3))
    print(getNumPos(atoms4))

    return (atoms4, specAtoms)

def getNumPos(atomPos):
    i = 0
    for lab, posList in atomPos.items():
        for pos in posList:
            i += 1
    return i
    
def catchEdgeAtoms(atomPos):
    xEdge = False
    yEdge = False
    zEdge = False
    
    atoms2 = {}

    
    for atom, posList in atomPos.items():
        i = 0
        for pos in posList:
            xEdge = False
            yEdge = False
            zEdge = False

            x = pos[0][0]
            y = pos[0][1]
            z = pos[0][2]
    
            if atom not in atoms2.keys():
                atoms2[atom] = []
            atoms2[atom].append(pos)
            for i,coord in enumerate(pos[0]):
                
                if coord > 0.8:
                    
                    if i == 0:
                        xEdge = True
                        atoms2[atom].append((np.array([x-1, y, z]), atom + ',edge' + str(i)))
                        i+=1
                    elif i == 1:
                        yEdge = True
                        atoms2[atom].append((np.array([x, y-1, z]), atom + ',edge' + str(i)))
                        i+=1
                    elif i ==2:
                        zEdge = True
                        atoms2[atom].append((np.array([x, y, z-1]), atom + ',edge'+ str(i)))
                        i+=1
            if xEdge and yEdge:
                atoms2[atom].append((np.array([x-1, y-1, z]), atom + ',edge'+ str(i)))
                i+=1
            elif xEdge and zEdge:
                atoms2[atom].append((np.array([x-1, y, z-1]), atom + ',edge'+ str(i)))
                i+=1
            elif yEdge and zEdge:
                atoms2[atom].append((np.array([x, y-1, z-1]), atom + ',edge'+ str(i)))
                i+=1
            elif xEdge and yEdge and zEdge:
                atoms2[atom].append((np.array([x-1, y-1, z-1]), atom + ',edge'+ str(i)))
                i+=1
    return atoms2

os.chdir('/home/matt/Dev/XDToolkit/test/data/Mn-dimer')
#y = getSymOps()
xxx = ins2all(file = 'shelx.ins', atomPosInput = getSymOps())

