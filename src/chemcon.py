'''
######################################################################
#-------------------CHEMCON-------------------------------------------
######################################################################
'''

import os
import hashlib
from subprocess import PIPE, Popen
from utils import listjoin, atomTableBegins, atomTableEnds
from devtools import timeDec
from numba import jit
'''
######################################
Short guide to path finding algorithm:
######################################

findAllPaths will start with the first neighbour to starting atom, and exhaust all possibilities,
before moving to the next neighbour. This applies no matter how many steps you are from the
starting atom. For example, all branches 5+ steps down would be explored for one 5 step neighbour, before
any other of the 5 step neighbours were explored.

Only adding the last branch to usedBranches causes this behaviour. If the current path differs from
the last path at any stage then everything downstream is removed from usedBranches, as these branches
are no longer used as something has been changed upstream.
'''

testDict = {'O(1)': ['C(1),asym'], 'C(1)': ['O(1),asym', 'O(2),asym', 'C(2),asym'], 'O(2)': ['C(1),asym'], 'O(3)': ['C(3),asym', 'H(4),asym'], 'C(3)': ['O(3),asym', 'C(2),asym', 'H(31),asym', 'H(32),asym'], 'H(4)': ['O(3),asym'], 'N(1)': ['C(2),asym', 'H(11),asym', 'H(12),asym', 'H(13),asym'], 'C(2)': ['N(1),asym', 'C(1),asym', 'C(3),asym', 'H(2),asym'], 'H(11)': ['N(1),asym'], 'H(12)': ['N(1),asym'], 'H(13)': ['N(1),asym'], 'H(2)': ['C(2),asym'], 'H(31)': ['C(3),asym'], 'H(32)': ['C(3),asym']}
testAtom = 'H(31),asym'

usedBranches = []
lastPath = []
pathStr = ''

#Mn-dimer with CPPgetEnvSig -                      ~730s
#Mn-dimer unoptimized getEnvSig -                   580.7068872s (9mins 40secs)
#Mn-dimer optimized getEnvSig -                     312.5087068s (5mins 10secs)
#Mn-dimer optimized w/ multiprocessing map          149.7411444s (2mins 30secs)
#Mn-dimer optimized w/ multiprocessing map_async    146.1650229s (2mins 25secs)

def CPPgetEnvSig(atomNeebDict, atom):
    '''
    Create an md5 hash value for the chemical environment of a given atom. Return this hash value.
    This function calls a compiled C++ script.
    '''
    inputStr = atom
    for item, neebs in atomNeebDict.items():
        inputStr+= '&' + item + '&' + listjoin(neebs, '*')
    os.chdir('/home/matt/dev/XDToolkit/src')
    cppFC = Popen(['./a.out'], shell=False, stdout=PIPE, stdin=PIPE)
    res = cppFC.communicate(bytes(inputStr, 'utf-8'))[0]
    cppFC.terminate()
    x = res.decode('utf-8')
    x = atom.split('(')[0].upper() +','.join(sorted(x.split('^')))
    y = hashlib.sha256(bytes(x, 'utf-8'))
    envHash = y.hexdigest()
    return (atom, envHash)

def findAllPaths(atom, atomLabsDict):
    '''
    Find all possible paths through the molecule starting from a given atom. Return all paths in list.
    '''
    paths = []
    atomNormLabs = {}
    pathList =[]
    passedAtoms=[]
    usedBranches = []
    lastPath = []

    for neebs in atomLabsDict.values():
        for neebLab in neebs:
            if neebLab not in atomNormLabs:
                atomNormLabs[neebLab] = neebLab.split(',')[0]
    
    while True:
        currAtom = atom             #Current atom as program walks through structure.
        passedAtoms=[]
        pathList=[]
        pathAppend = pathList.append
        passAppend = passedAtoms.append

        steps=-1    
                                    #Branches are format 'ATOM~number of steps through path' i.e. 'C(1),asym~2'
                                    #means C(1) was passed on the second step through the structure.
    
        while currAtom not in passedAtoms:        #This condition will be satisfied until path can't go anywhere unvisited.
    
            passAppend(currAtom)          #Add current atoms to visited atoms list.
            steps += 1                            #Add step to total number of steps in path
            
            atomNeebs = atomLabsDict[atomNormLabs[currAtom]]        #Get neighbours of current atom.
    
            for neeb in atomNeebs:                              #Go through neighbours of current atom one by one.
                branchTag = (steps, neeb)
                #If a neighbour of current atom is found that hasn't been visited,
                #and hasn't been branched too in other trips down the same path.
                if  branchTag not in usedBranches: 
                    
                    if neeb not in passedAtoms:
                        newUsedBranch = branchTag
                        currAtom = neeb         #Make the current atom the new atom, quit the for loop
                        pathAppend(neeb)   #Add new atom in path to pipe separated path string.
                        
                        try:
                            #If this path changes from the previous path at this atom,
                            #remove all branch tags downstream from current number of steps.
                            #i.e. if number of steps is 3 only keep branch tags with number of steps 0,1,2,3.
                            if lastPath[steps] != neeb:
                                usedBranches = [item for item in usedBranches if item[0] <= steps]
        
                        except IndexError:
                            #If the last path didn't go this far, or it was empty also remove downstream branch tags.
                            usedBranches = [item for item in usedBranches if item[0] <= steps]
        
                        #Add branch tag of every atom in path to list.
                        break                   #and look for the next atom along the path.
        
        usedBranches.append(newUsedBranch)
        #If no path is returned all paths have been found and the while loop is ended.
        if not pathList:
            break

        else:
            lastPath = pathList
            #Format path to string of atom types e.g. 'CCCNCCCH'
            pathStr=''
            for item in pathList:
                pathStr+=item[:2]
            paths.append(pathStr)  
            
    #When no more paths can be found return list of paths.
    return paths

def getEnvSig(atomLabsDict, atom):
    '''
    Create a sha256 hash value for the chemical environment of a given atom. 
    Return this hash value.
    '''
    #Sort list of paths alphabetically so it is the same no matter what order paths were found in
    #Make string with starting atom types followed by , joined sorted list of all paths.
    pathString = atom.split('(')[0].upper() + ','.join(sorted(findAllPaths(atom, atomLabsDict)))
    hashObj = hashlib.sha256(bytes(pathString,'utf-8'))        #Generate unique hash value of paths.
    envHash = (atom, hashObj.hexdigest())
    return envHash    

def removeCHEMCON():
    '''
    Remove all CHEMCON from xd.mas.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        atomTab = False

            #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            #Detect end of ATOM table
            if atomTableEnds(line):
                atomTab = False

            if atomTab:
                row = line.split()
                if len(row) == 13:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}\n'.format(*row)
                    newmas.write(rowStr)
                else:
                    newmas.write(line)
            else:
                newmas.write(line)

            #Detect start of ATOM table
            if atomTableBegins(line):
                atomTab = True
    
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def findCHEMCONbyElement():
    '''
    Find atoms of the same element and group them. Return dictionary of groupings.
    '''
    with open('xd.mas','r') as mas:

        atomTab = False
        prevElement = ''
        CHEMCON = {}
        currentAtom=''

            #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            #Detect end of ATOM table
            if atomTableEnds(line):
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
            if atomTableBegins(line):
                atomTab = True

    return CHEMCON


def findCHEMCONbyInputElement(inputElementList):
    '''
    Find atoms of given element and group them. Return grouping.
    '''
    with open('xd.mas','r') as mas:

        atomTab = False
        elementParents = {}
        CHEMCON = {}
        eleList = inputElementList
        #Convert elements to upper case
        eleList = [item.upper() for item in eleList]
        eleList = [item + '(' for item in eleList if len(item)==1]

            #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            #Detect end of ATOM table
            if atomTableEnds(line):
                atomTab = False

            #In atom table make CHEMCON dictionary
            if atomTab:
                row = line.upper().split()

                if line[:2] in eleList:
                    if line[:2] in elementParents:
                        CHEMCON[elementParents[line[:2]]].append(row[0])
                    else:
                        CHEMCON[row[0]] = []
                        elementParents[line[:2]] = row[0]
                    
            #Detect start of ATOM table
            if atomTableBegins(line):
                atomTab = True
        
    return CHEMCON

def findCHEMCONbyNeebors():
    '''
    Group atoms of the same element with the same nearest neighbours. Return groupings.
    '''
    with open('xd.mas','r') as mas:

        atomTab = False                         #Bool to detect if the for loop is in the atom table when going through xd.mas
        CHEMCON = {}                            #Dict to return CHEMCON for writing to xd.mas with writeCHEMCON()
        neebsDict = copy.copy(globAtomTypes)    #Dict of nearest neighbours {'C(2):['C','C','H','H']}
        usedSigs = []                           #List of sigs of chemical environment i.e. CCHN = N bonded to CCH
        SigParentDict = {}                      #dict with sigs as keys and parent atom labels of that sig as values
        chemconSig=''

            #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            #Detect end of ATOM table
            if atomTableEnds(line):
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
            if atomTableBegins(line):
                atomTab = True

    return CHEMCON

def findCHEMCONbyInputAtoms(atomList):
    '''
    Group given atoms. Retun grouping.
    '''
    with open('xd.mas','r') as mas:

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
            if atomTableEnds(line):
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
                    if not parent2Found:
                        if row[12] == parent:
                            parent2 = row[0].upper()
                            parent2Found = True
                            CHEMCON[parent2] = []

                    elif row[12] == parent:
                        CHEMCON[parent2].append(row[0].upper())


            #Detect start of ATOM table
            if atomTableBegins(line):
                atomTab = True

    return CHEMCON


def writeCHEMCON(CHEMCONdict):
    '''
    Write CHEMCON column in xd.mas with given groupings.
    '''
    atomTab = False
    CHEMCON = CHEMCONdict
    written = False

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if atomTableEnds(line):
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

            if atomTableBegins(line):
                atomTab = True

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
        with open('xd.mas','r') as mas:

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
            for line in mas:

                if atomTableEnds(line):
                    atomTab = False

                if atomTab:
                    row = str.split(line)
                    #If atom is in dictionary it is added with local coordinate system
                    if len(row) == 13:
                        chemcon = True
                        break

                if atomTableBegins(line):
                    atomTab = True

    except FileNotFoundError:
        pass

    return chemcon



import time

#def CPPwrapper(func, *args, **kwargs):
#     def CPPwrapped():
#         return func(*args, **kwargs)
#     return CPPwrapped
#CPPwrapped = CPPwrapper(CPPgetEnvSig, testAtom, testDict)
#print(timeit.timeit(CPPwrapped, number=1000))
#x = CPPgetEnvSig(testAtom, testDict)
#print(x)
#print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
#def wrapper(func, *args, **kwargs):
#     def wrapped():
#         return func(*args, **kwargs)
#     return wrapped
#wrapped = wrapper(getEnvSig, testDict, testAtom)
#print(timeit.timeit(wrapped, number=10000))
#s = time.time()
#x = getEnvSig(testDict, testAtom)
#f = time.time()
#print(x)
#print(f-s)
#print('~'*50)

#import multiprocessing as mp
#from functools import partial
#import time
#os.chdir('/home/matt/dev/XDToolkit/src')
#y=time.time()
#
#pool = mp.Pool()
#z = pool.map(partial(getEnvSig, atomLabsDict=testDict), testAtom)
#
#x=time.time()
#print(x-y)
