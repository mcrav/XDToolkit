import numpy as np
import copy
    
def rearrangeSYMM(symOp):
    newSymOp = ['x','y','z']
    for item in symOp:
        newItem = item.replace(' ','')
        if 'x' in newItem:
            newSymOp[0] = newItem
        elif 'y' in newItem:
            newSymOp[1] = newItem
        else:
            newSymOp[2] = newItem
    
    return newSymOp

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

                symOp = line[4:].lower().strip('\n').split(',')
                symOp = rearrangeSYMM(symOp)
                symOps.append(symOp)
            
            elif symOps and centering and not line.startswith('SYMM'):
                break
                
    return (centering, centrosym, symOps)

def genAtomsbySym(insInfo, x,y,z):
    centSymOps = []
    symOps = []
    
    LATT = insInfo[0]
    SYMM = insInfo[2]
    
    #For each centering append (position, inverseSymOp)
    if 'P' in LATT:
        centSymOps.append((np.array([x, y, z]), ('x', 'y', 'z')))

    if 'I' in LATT:
        centSymOps.append((np.array([0.5 + x, 0.5 + y, 0.5 + z]), ('x - 0.5', 'y - 0.5', 'z - 0.5')))
            
    if 'R' in LATT:
        centSymOps.append((np.array([(2/3) + x, (1/3) + y, (1/3) + z]), ('x - (2/3)', 'y - (1/3)', 'z - (1/3)')))  
        centSymOps.append((np.array([(1/3) + x, (2/3) + y, (2/3) + z]), ('x - (1/3)', 'y - (2/3)', 'z - (2/3)')))
        
    if 'F' in LATT:
        centSymOps.append((np.array([x, 0.5 + y, 0.5 + z]), ('x', 'y -0.5', 'z-0.5')))
        centSymOps.append((np.array([0.5 + x, y, 0.5 + z]), ('x-0.5', 'y', 'z-0.5')))
        centSymOps.append((np.array([0.5 + x, 0.5 + y, z]), ('x-0.5', 'y-0.5', 'z')))
        
    if 'A' in LATT:
        centSymOps.append((np.array([x, 0.5 + y, 0.5 + z]), ('x', 'y -0.5', 'z-0.5')))
        
    if 'B' in LATT:
        centSymOps.append((np.array([0.5 + x, y, 0.5 + z]), ('x-0.5', 'y', 'z-0.5')))
        
    if 'C' in LATT:
        centSymOps.append((np.array([0.5 + x, 0.5 + y, z]), ('x-0.5', 'y-0.5', 'z')))
    
    for insSymOp in SYMM:
        for symOp in centSymOps:
            x = symOp[0][0]
            y = symOp[0][1]
            z = symOp[0][2]
            
            invSymOp = invertSYMMOp(insSymOp)

            #Coords given as well as tuple of SYMM op and centering op
            symOps.append((np.array([eval(insSymOp[0]), eval(insSymOp[1]), eval(insSymOp[2])]), [invSymOp, symOp[1]]))
            
    return symOps
 
def check4minus(symOp):
    if '-x' in symOp or '-y' in symOp or '-z' in symOp:
        return True
    else:
        return False
       
def invertSYMMOp(insSymOp):
    invSymOp = []
    newItem = ''

    for i,item in enumerate(insSymOp):
        itemPacked = item.replace(' ','')

        if not check4minus(itemPacked):

            if itemPacked[0].isdigit():
                newItem = '-' + itemPacked
            elif itemPacked[0] == '-' and itemPacked[1].isdigit():
                newItem = itemPacked[1:]
            else:
                newItem = itemPacked
        else:
            newItem = itemPacked

        invSymOp.append(newItem)

    return invSymOp
    

def applySymOps(fracPos):
    asymUnit = fracPos
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
            atoms2[lab] = [(atomPos[0], lab + ',asym', [('x','y','z')])]
            i+=1
        
        for item in genAtomsbySym(insInfo, x,y,z):                              #CENT, CENT+SYMM
            atoms2[lab].append(  (item[0], lab + ',' + str(i), item[1])  )      #item[0] is position, item[1] is inverse symOp 
            i+=1
        
        for item in SYMM:
            coords = np.array([eval(item[0]), eval(item[1]), eval(item[2])])
            invSymOp = invertSYMMOp(item)
            atoms2[lab].append((coords, lab + ',' + str(i), [invSymOp]))   #SYMM
            i+=1
            
        if centrosym:

            for item in genAtomsbySym(insInfo, -x,-y,-z):
                centSym = item[1][1]
                insSym = item[1][0]
                atoms2[lab].append((item[0], lab + ',' + str(i), [insSym, centSym, ('-x','-y','-z'  )]))               #INV, CENT, CENT+SYM
                i+=1
            
            atoms2[lab].append((np.array([-x,-y,-z]), lab + ',' + str(i), [('-x','-y','-z')]))    #INV
            i+=1
            
            x = -x
            y = -y
            z = -z
            
            for item in SYMM:
                coords = np.array([eval(item[0]), eval(item[1]), eval(item[2])])
                
                invSymOp = invertSYMMOp(item)
                atoms2[lab].append(  (coords, lab + ',' + str(i), [invSymOp, ('-x','-y','-z')])   )     #SYMM + inv
                
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

#    atoms4 = catchEdgeAtoms(atoms3)
#    for atom, posList in atoms3.items():
#        
#        for pos in posList:
#            for pos2 in getCombos(pos[0]):
#                atoms4.setdefault(atom,[]).append((pos2, 'ub' + str(i)))
#                i+=1
    
    for atom, posList in atoms3.items():

        for pos in posList:

            specAtoms[pos[1]] = (pos[0], pos[2])

    
    print(getNumPos(atoms3))
#    print(getNumPos(atoms4))
#    print(atoms4['SI(0),3'])
#    print('atoms4------------------------------------------------------')
    return (atoms3, specAtoms)

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

    j = 0
    for atom, posList in atomPos.items():
        
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
                    invSym = copy.copy(pos[2])
                    if i == 0:
                        xEdge = True
                        
                        invSym.append(('x+1','y','z'))
                        atoms2[atom].append((np.array([x-1, y, z]), atom + ',edge' + str(j), invSym))
                        j+=1
                    elif i == 1:
                        yEdge = True

                        invSym.append(('x','y+1','z'))
                        atoms2[atom].append((np.array([x, y-1, z]), atom + ',edge' + str(j), invSym))
                        j+=1
                    elif i ==2:
                        zEdge = True

                        invSym.append(('x','y','z+1'))
                        atoms2[atom].append((np.array([x, y, z-1]), atom + ',edge'+ str(j), invSym))
                        j+=1
                        
            if xEdge and yEdge:

                invSym.append(('x+1','y+1','z'))
                atoms2[atom].append((np.array([x-1, y-1, z]), atom + ',edge'+ str(j), invSym))
                j+=1
            elif xEdge and zEdge:

                invSym.append(('x+1','y','z+1'))
                atoms2[atom].append((np.array([x-1, y, z-1]), atom + ',edge'+ str(j), invSym))
                j+=1
            elif yEdge and zEdge:

                invSym.append(('x','y+1','z+1'))
                atoms2[atom].append((np.array([x, y-1, z-1]), atom + ',edge'+ str(j), invSym))
                j+=1
            elif xEdge and yEdge and zEdge:
                invSym.append(('x+1','y+1','z+1'))
                atoms2[atom].append((np.array([x-1, y-1, z-1]), atom + ',edge'+ str(j), invSym))
                j+=1

    return atoms2


def getCombos(coord):
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
    
    for xval in xs:
        for yval in ys:
            for zval in zs:
                combos.append((xval, yval, zval))
    
    print(len(combos))
    return(combos)

