import os
import csv
import re
import numpy as np

def getCPSearch():
    with open('xd_pro.out','r', encoding = 'utf-8', errors = 'ignore') as pro, open('xd_pro.csv','w') as proCsv:
        proCsv = csv.writer(proCsv, delimiter = ',')
        proCsv.writerow(['Atom 1', 'Atom 2', 'rho', 'd2rho', 'ellip'])
        cpTable = False
        i=0
        for line in pro:
            
            if not line.strip():
                cpTable = False  
                
            if cpTable:
                row = line.split()
                
                if i==0:
                    csvRow = [row[0], row[1][1:], row[2], row[3]]
                    
                else:
                    csvRow.append(row[3])
                    print(csvRow)
                    proCsv.writerow(csvRow)
                    
                i+=1
                if i==2:
                    i=0
            
            if line.startswith('                                                 Hessian Eigenvalues   ellip'):
                cpTable = True
            
            
def geom2Lengths():
    with open('xd_geo.out','r') as geo, open('bond_lengths.csv','w') as lenCsv:
        lenCsv = csv.writer(lenCsv, delimiter = ',')
        lenCsv.writerow(['Atom 1', 'Atom 2', 'Bond length'])
        bondLen = False
        for line in geo:
            
            if bondLen:
                if line.strip():
                    row = line.split()
                    csvRow = [line[:7].replace('-','').strip(), line[9:16].strip(), row[-1]]
                    print(csvRow)
                    lenCsv.writerow(csvRow)
                    
                else:
                    bondLen = False
                    break
                
            elif line.startswith(' Bond lengths'):
                bondLen = True   
                
def geom2Angles():
    with open('xd_geo.out','r') as geo, open('bond_angles_2.5.csv','w') as angCsv:
        angCsv = csv.writer(angCsv, delimiter = ',')
        angCsv.writerow(['Atom 1', 'Atom 2', 'Atom 3', 'Bond angle'])
        bondAng = False
        for line in geo:
            
            if bondAng:
                if line.strip():
                    row = line.split()
                    csvRow = [line[1:7].replace('-','').strip(), line[9:15].replace('-','').strip(), line[17:24].replace('-','').strip(), row[-1]]
                    print(csvRow)
                    angCsv.writerow(csvRow)
                    
                else:
                    bondAng = False
                    break
                
            elif line.startswith(' Bond angles'):
                bondAng = True  
    
def getBondLengths():
    
    lengthDict = {}
    
    with open('bond_lengths.csv','r') as lengths:
        lengths = csv.reader(lengths)
        
        for row in lengths:
            lengthDict[frozenset([row[0], row[1]])] = row[2]
            
    return lengthDict

def addLengthsToProOut():
    lengths = getBondLengths()
    
    with open('sorted_xd_pro.csv','r') as pro, open('xd_pro_lengths.csv','w') as proLengths:
        pro = csv.reader(pro)
        proLengths = csv.writer(proLengths)
        proLengths.writerow(['Atom 1', 'Atom 2', 'rho', 'd2rho', 'ellip', 'Bond Length'])
        i=0
        for row in pro:
            if i >0:
                newRow = row
                newRow.append(lengths[frozenset([row[0], row[1]])])
                proLengths.writerow(newRow)
                print(newRow)
            i+=1
            
def sortTables(atomOrder, file):
    rows = {}
    atomOrder = [item + '(' for item in atomOrder if len(item)==1]
    newfilename = 'sorted_' + file
    with open(file,'r') as file, open(newfilename,'w') as newFile:
        file = csv.reader(file)
        newFile = csv.writer(newFile)
        i=0
        for row in file:
            if i>0:
                rows.setdefault(row[0][:2],[]).append(row)
            else:
                newFile.writerow(row)
            i+=1
        
        for element in atomOrder:
            for row in sorted(rows[element], key = lambda row: int(re.sub("[^0-9]", "", row[0].split('(')[1].strip(')')))):
                newFile.writerow(row)

def readTOPXDFolder(folder):
    os.chdir(folder)
    with open('topxdRes.csv','w') as topxdcsv:
        topxdcsv = csv.writer(topxdcsv)
        topxdcsv.writerow(['Atom', 'Q', 'L'])
        
        for file in os.listdir(folder):
            splitFile = os.path.splitext(file)
            if splitFile[1] == '.out':
                atom = splitFile[0].split('_')[1]
                newRow = [atom]
                with open(file, 'r', encoding='utf-8', errors = 'ignore') as topxdout:
                    
                    for line in topxdout:
                        if line.startswith('              Q'):
                            newRow.append(float(line.split()[1]))
                            
                        elif line.startswith('              L'):
                            newRow.append(float(line.split()[1]))
                            break
                    topxdcsv.writerow(newRow)
                    print(newRow)
 


os.chdir('/home/matt/Work/cucf3_final')
#
readTOPXDFolder(os.path.abspath('cucf3_highangle_topxd'))
atomOrder = ['CU','F','N','C','H']
sortTables(atomOrder, 'topxdRes.csv')

#sortTables(atomOrder, 'xd_pro.csv')
#geom2Angles()
#geom2Lengths()
#addLengthsToProOut()
#getCPSearch()



def getBondDist(atom1c, atom2c, a, b, c, alpha, beta, gamma):
    '''
    Get distance between 2 atoms. Return distance.
    '''
    delta1 = a*(atom1c[0] - atom2c[0])
    delta2 = b*(atom1c[1] - atom2c[1])
    delta3 = c*(atom1c[2] - atom2c[2])

    return np.sqrt(delta1**2 + delta2**2 + delta3**2 + (2*delta1*delta2*np.cos(gamma)) + (2*delta1*delta3*np.cos(beta)) + (2*delta2*delta3*np.cos(alpha)))

def findClosestAtom2ResidualPeak(residualPeakCoords, xd_fft):  
    '''
    Find closest atom to residual peak in xd_fft.out. 
    Return atom and distance from atom to residual peak.
    '''           
    with open(xd_fft,'r') as fft:
        atomTab = False
        unitCellFound = False
        atomPos = {}
        i=0
        for line in fft:
            
            if atomTab:
                if i==2:
                    break
                elif line.startswith(' ------'):
                    i+=1
                    continue
                else:
                    row = line.split()
                    atomPos[row[1]] = np.array([float(row[3]), float(row[4]), float(row[5])])
                    
            elif not unitCellFound and line.startswith(' Unit cell parameters'):
                row = line.split()
                a = float(row[3])
                b = float(row[4])
                c = float(row[5])
                alpha = np.radians(float(row[6]))
                beta = np.radians(float(row[7]))
                gamma = np.radians(float(row[8]))
                unitCellParams = (a, b, c, alpha, beta, gamma)
                unitCellFound = True
            
            elif line.startswith('  no  atom'):
                atomTab = True
                
    minDist = 1000
    minAtom = ''
    resPeak = residualPeakCoords
    for atom, pos in atomPos.items():
        bondDist = getBondDist(resPeak, pos, *unitCellParams)
        if bondDist < minDist:
            minDist = bondDist
            minAtom = atom
    
    print(minAtom)
    print(minDist)
    return(minAtom, minDist)
