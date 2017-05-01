'''
#######################################################################
#-------------------RESULTS--------------------------------------------
#######################################################################
'''

import os
import subprocess as sub
from utils import getAtomList

def FFTDetective():
    '''
    Find closest atom to highest peak in xd_fft.out. Return atom label.
    '''
    with open('xd_fft.out','r') as fft:    #Open output file from XDFFT

        table = False                   #Bool to find start of table with highest peaks
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
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

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

    os.remove('xd.mas')             #Delete original xd.mas file
    os.rename('xdnew.mas','xd.mas') #Rename xdnew.mas to xd.mas

def FOU3atoms(atom, atomLabs, atomPos):
    '''
    Setup XDFOUR instructions in xd.mas to run on the plane of a given atom and 2 of its neighbours.
    '''
    atom = atom.upper()                 #Convert atom given in argument to uppercase to avoid problems with upper/lower
    fouBool = False                     #Bool to detect start of XDFOUR instructions
    neebsRaw = atomLabs   #Get dict of nearest neighbours with labels for each atom
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

                        pos = atomPos[neeb][0]
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
                    pos = atomPos[neighbour][0]
                    x,y,z = pos[0], pos[1], pos[2]
                    splitLab = neighbour.split(',')
                    rowStr = 'XYZ label {0} {1:.4f} {2:.4f} {3:.4f} symm  1 trans 0 0 0 *mark on plot\n'.format(splitLab[0],x,y,z)
                    newmas.write(rowStr)
                    for neeb in neebsRaw[neighbour.split(',')[0]]:
                        if neeb.split(',')[0] != atom:
                            nextNeighbour = neeb
                            break
                    pos = atomPos[nextNeighbour][0]
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


def getRF2(file = None):
    '''
    Get final RF2 value from xd_lsm.out. Return value.
    '''
    if not file:
        lsm = open('xd_lsm.out','r')

    else:
        filePath = file
        lsm = open(filePath,'r')

    try:
        finalCycle = False
        rf2 = ''
        rFound= False

        for line in lsm:
            #Find final cycle and line with R-value
            if finalCycle:
                
                if line.startswith('  R{F^2} ='):
                    row = str.split(line)
                    rf2 = float(row[2])
                    rFound = True
                    finalCycle = False

            elif line.startswith('                      Residuals after final cycle'):
                finalCycle = True

    finally:
        lsm.close()

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


def getDMSDA(lsmFile):
    '''
    Get DMSDA results from xd_lsm.out. Return DMSDA results.
    '''
    with open(lsmFile,'r')  as lsm:         #Open xd_lsm.out to read

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

    return (str('{0:.1f}'.format(averageDMSDA)),str(max(dmsdaList)),dmsdaFull) #return tuple of average dmsda, max dmsda and the dmsda dict)

def setupmasTOPXD(atom):
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        
        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:
            #Detect end of ATOM table
            if line.startswith('ATBP *atoms'):
                row = line.split()
                row[2] = atom
                newLine = ' '.join(row) + '\n'
                newmas.write(newLine)

            else:
                newmas.write(line)

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def runTOPXD(atomList):
    for atom in atomList:
        setupmasTOPXD(atom)
        file = open('topxd_' + atom, 'w')
        topxd = sub.Popen('/home/matt/dev/XDTstuff/xd-distr/bin/topxd', shell = False, cwd = os.getcwd(), stdout = file)
        topxd.wait()
        file.close()
        
        
os.chdir('/home/matt/dev/XDTstuff/topxdtest')
runTOPXD(getAtomList())