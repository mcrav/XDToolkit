import os

def findLabels():
    
   file = open('XD Toolkit.py','r')
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

    file = open('XD Toolkit.py','r')
    newfile = open('XD Toolkitnew.py','w')
    
    for line in file:
        if line.startswith('        self.labList'):
            newfile.write('        self.labList = ' + str(labList) + '\n')
        else:
            newfile.write(line)
            
    file.close()
    newfile.close()

    file = open('XD Toolkitnew.py','r')
    newfile = open('XD Toolkitnew2.py','w')
    
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

    os.remove('XD Toolkitnew.py')
    os.remove('XD Toolkit.py')
    os.rename('XD Toolkitnew2.py','XD Toolkit.py')
    
