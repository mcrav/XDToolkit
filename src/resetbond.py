'''
######################################################################
#-------------------RESET BOND----------------------------------------
######################################################################
'''
from utils import lab2type
import os
import copy


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

def getHList():
    '''
    Get list of all hydrogen atom labels from xd.mas. Return list.
    '''
    with open('xd.mas','r') as mas:
        Hs = []
        for line in mas:
            atomTab = False
    
            #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
            for line in mas:
    
                if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                    atomTab = False
    
                #If 'All' is unchecked, make dictionary with H atoms inputted by user, their connected C atoms and inputted bond distance (1.09 default)
                if atomTab and line.startswith('H('):
                    row = line.split()
                    Hs.append(row[0].upper())
              
    
                if line.startswith('ATOM     ATOM0'):
                    atomTab = True
                    
    return Hs

def getRBHList():
    '''
    Get list of all H atoms with reset bond instructions in xd.mas. Return list.
    '''
    with open('xd.mas','r') as mas:
        RBHs = []
        for line in mas:
            if line.startswith('RESET BOND') or line.startswith('!RESETBOND'):
                row = line.split()
                for item in row:
                    if item.startswith('H('):
                        RBHs.append(item)
                        
    return RBHs

def check4RBHs():
    '''
    Check if there are any hydrogen atoms without reset bond instructions.
    Return list of these atoms.
    '''
    Hs = getHList()
    RBHs = getRBHList()
    noRBHs = [item for item in Hs if item not in RBHs]
    return noRBHs

def check4RB():
    '''
    Check if reset bond instructions have been added. Return result as bool.
    '''
    try:
        with open('xd.mas','r') as mas:
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

        if not H:
            RB = True

        return RB

    except:
        return False


def resetBond(length,atoms,allornot = False):
    '''
    Add reset bond instruction to xd.mas for given atoms and length.
    '''
    with open('xd.mas', 'r') as mas:
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
    #Open test file to write
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

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

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def autoAddResetBond(atomLabDict, atomTypeDict):
    '''
    Automatically add reset bond instructions based on detected H environments.
    '''
    #Get neighbours by type and label
    neighboursType = atomTypeDict
    neighbourLabsRaw = copy.copy(atomLabDict)
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

                aromaticConfirmed = confirmAromaticity(atom + ',asym', atomLabDict)

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

    if aromaticHList:
        resetBond('1.083', aromaticHList)
        addedRBs.extend(aromaticHList)

    if secondaryHList:
        resetBond('1.099',secondaryHList)
        addedRBs.extend(secondaryHList)

    if primaryHList:
        resetBond('1.092',primaryHList)
        addedRBs.extend(primaryHList)

    if XmethylList:
        resetBond('1.066',XmethylList)
        addedRBs.extend(XmethylList)

    if CmethylList:
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


def confirmAromaticity(atom, atomLabsDict):
    '''
    Check if atom is part of a phenyl ring. Return result as bool.
    '''
    lastPath = []
    paths = []
    usedBranches = []
    aromaticConfirmed = False

    while True:

        pathRes = getPathAromatic(atom, atomLabsDict, usedBranches, lastPath)  #Get a path.

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

def autoResetBond(atomLabsDict, atomTypeDict):
    '''
    Automatically add reset bond instructions to xd.mas and find out what H atoms have been missed.
    Return missed atoms.
    '''
    resetBondsAdded = autoAddResetBond(atomLabsDict, atomTypeDict)

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

    return missedAtoms


def delResetBond():
    '''
    Delete all reset bond instructions from xd.mas.
    '''
    #Open xd.mas to read and xdnew.mas to write new mas file
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        #Go through xd.mas line by line
        #If line has reset bond instructions don't write it, otherwise write line unchanged
        for line in mas:
            if line.startswith('RESET BOND') or line.startswith('!RESET BOND'):
                pass
            else:
                newmas.write(line)

    os.remove('xd.mas')
    os.rename('xd.newmas','xd.mas')
    



