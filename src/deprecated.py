def wizAddLocCoords():
    '''
    Deprecated. Add local coordinate systems in xdwiz.mas to xd.mas
    '''
    with open('xdwiz.mas','r') as rb:
        ccstr = ''
        atomTab = False

        for line in rb:
            if line.startswith('END ATOM'):
                atomTab = False

            if atomTab:
                ccstr += line

            if line.startswith('ATOM     ATOM0'):
                atomTab = True

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
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

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def findNeeborType(lstFile):
    '''
    Get nearest neighbour types for every atom from lst file. Return nearest neighbour types.
    '''
    #neighboursSym uses format {'O1':['C','C']}
    neighboursSym = {}
    #Opens lst file from SHELXL refinement. Need to update so user chooses lst file from file explorer.
    with open(lstFile,'r') as lstFile:
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

    return neighboursSym

def findNeeborLabels(lstFile):
    '''
    Get nearest neighbour labels from lst file. Return nearest neighbour labels.
    '''
    #neighbours uses format {'O1':['C1','C2']}
    neighbours = {}

    #Opens lst file from SHELXL refinement. Need to update so user chooses lst file from file explorer.
    with open(lstFile,'r') as lstFile:

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

    return neighbours

def getBondAngles(lstFile):
    '''
    Get bond angles from lst file. Return bond angles.
    '''
    with open(lstFile, 'r')as lstFile:

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

            with open('xd_lsm.out','r') as lsm:

                for line in lsm:
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

#Works but unused
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
