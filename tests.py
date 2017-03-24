
import os
from XDToolkit import ins2all

os.chdir('/home/matt/dev/XDTstuff')

def getTestDict(bondFile):
    '''
    Get test dictionary from bond file.
    '''
    neebs = {}
    print(bondFile)
    with open(bondFile,'r') as bonds:
        
        for line in bonds:
            if line[0].isdigit():
                row = str.split(line)
                
                neebs.setdefault(row[1],[]).append(row[2])
                neebs.setdefault(row[2],[]).append(row[1])
                
    return neebs
                

def checkNeebs(testDictRaw):
    '''
    Test neebs against given test dictionary.
    '''
    x = ins2all()
    atomLabsRaw = x[0]
    atomLabs = {}
    testDict = {}
    
    problems = {}
    
    for atom,neebs in atomLabsRaw.items():
        atomLabs[atom.replace('(','').replace(')','')] = tuple(sorted([neeb.split(',')[0].replace('(','').replace(')','') for neeb in neebs]))
        
    for atom,neebs in testDictRaw.items():
        testDict[atom] = tuple(sorted(neebs))
    
    for atom,neebs in atomLabs.items():
        if not neebs == testDict[atom]:
            problems[atom] = neebs
            
    return problems


def checkAllNeebs():
    '''
    Test all neighbours in testData folder.
    '''
    problems = {}
    
    for folder in os.listdir('testData'):
        if folder != 'alanyl-methionine':
        
            os.chdir('/home/matt/dev/XDTstuff/testData/{}'.format(folder))
            print(os.getcwd())
            testDict = getTestDict('bonds.tsv')
            
            x = checkNeebs(testDict)
            
            if x:
                problems[folder] = x
                print(x)
    
    for item in problems:
        print('\n\n')
        print('------------------------------------------')
        print(item.upper())
        print('')
        print(problems[item])
    return problems


def checkFolder(folder):
    '''
    Test specific folder.
    '''
    os.chdir('/home/matt/dev/XDTstuff/testData/{}'.format(folder))
    print(checkNeebs(getTestDict('bonds.tsv')))
            

#checkAllNeebs()
checkFolder('carba') 
        
    