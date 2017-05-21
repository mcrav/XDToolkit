import os

def findLabels():
    
   file = open('XDToolkit.py','r')
   labList = []
   
   for line in file:
       row = str.split(line,'.')
       i=0
       for item in row:
           if len(item) > 4 and (i+2) < len(row):
               if item[-4:] == 'self' and row[i+2].startswith('setText('):
                   label = 'self.' + row[i+1]
                   if label not in labList and label != 'self.cwdLab' and 'UserIns' not in label:
                       labList.append(label)
       i+=1
       
   file.close()
   return(labList)
      
def writeLabList():

    labList = findLabels()

    file = open('XDToolkit.py','r')
    newfile = open('XDToolkitnew.py','w')
    
    for line in file:
        if line.startswith('        self.labList'):
            newfile.write('        self.labList = ' + str(labList) + '\n')
        else:
            newfile.write(line)
            
    file.close()
    newfile.close()

    file = open('XDToolkitnew.py','r')
    newfile = open('XDToolkitnew2.py','w')
    
    for line in file:
        if line.startswith('        self.labList'):
            for char in line:
                if char.isalpha() or char.isdigit() or char == '[' or char == ']' or char == '.' or char == ' ' or char == ',' or char == '=':
                    newfile.write(char)
            newfile.write('\n')
        else:
            newfile.write(line)
            
    file.close()
    newfile.close()

    os.remove('XDToolkitnew.py')
    os.remove('XDToolkit.py')
    os.rename('XDToolkitnew2.py','XDToolkit.py')

def getNumLinesCode():
    i = 0
    targetFiles = ['XDToolkit.py',
                   'xdfiletools.py',
                   'xderrfix.py',
                   'wizardfuncs',
                   'utils.py',
                   'results.py',
                   'resetbond.py',
                   'initfuncs.py',
                   'emailfuncs.py',
                   'chemcon.py',
                   'databank.py',
                   'backup.py',
                   'asym2unit.py']
    for file in os.listdir(os.getcwd()):
        if file in targetFiles:
            i+=len(open(file,'r').readlines())
            
    return i

#os.chdir('/home/matt/dev/XDToolkit/src')
#print(getNumLinesCode())
    
