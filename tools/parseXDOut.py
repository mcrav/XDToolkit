import os
import csv
import re

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
 


os.chdir('/home/matt/Work/cosph_cov/topxd_files')
atomOrder = ['CO','S(','P(','C(','H(']
sortTables(atomOrder, 'topxdRes.csv')
#sortTables(atomOrder, 'xd_pro.csv')
#geom2Angles()
#geom2Lengths()
#addLengthsToProOut()
#getCPSearch()
                

    