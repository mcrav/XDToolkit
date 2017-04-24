'''
########################################################################
#--------------------WIZARD---------------------------------------------
########################################################################
'''

import os
import shutil

def seqMultRef(l):
    '''
    Add multipoles up to a given l value to key table.
    '''
    l+=5       #Make l correspond to row index in key table

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        keyTab = False

        for line in mas:

            if line.startswith('KAPPA'):
                keyTab = False

            if keyTab:
                row = str.split(line)
                rowStr = line[:48]
                for i in range(6,10):
                    if i <= l:
                        rowStr += (' ' + row[i])
                    else:
                        rowStr += (' ' + ''.join(['0' for char in row[i]]))
                rowStr += '\n'
                newmas.write(rowStr)

            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def wizAddResetBond():
    '''
    Add reset bond instructions in xdwiz.mas to xd.mas.
    '''
    rbstr = ''
    with open('xdwiz.mas', 'r') as rb:
        for line in rb:
            if line.startswith('RESET BOND'):
                rbstr += line
            elif line.startswith('KEY'):
                break
            
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        for line in mas:
            if line.startswith('WEIGHT'):
                newmas.write(line)
                newmas.write(rbstr)
            else:
                newmas.write(line)

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')


def wizAddLocCoords():
    '''
    Add local coordinate systems in xdwiz.mas to xd.mas
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


def wizAddCHEMCON():
    '''
    Add CHEMCON from xdwiz.mas to xd.mas.
    '''
    #Read CHEMCON from xdwiz.mas
    with open('xdwiz.mas','r') as rb:

        cc = []
        atomTab = False

        for line in rb:
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False

            if atomTab:
                row = str.split(line)
                if len(row) == 13:
                    cc.append(row[12])
                else:
                    cc.append(' ')

            if line.startswith('ATOM     ATOM0'):
                atomTab = True
            
    print(cc)
    #Write CHEMCON to xd.mas
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        i=0
        for line in mas:

            if line.startswith('END ATOM') or line.startswith('!DUM') or line.startswith('DUM'):
                atomTab = False

            if atomTab:
                row = line.split()
                row.append(cc[i])
                rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}\n'.format(*row)
                newmas.write(rowStr)
                i+=1

            else:
                newmas.write(line)

            if line.startswith('ATOM     ATOM0'):
                atomTab = True

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
    
def wizAddMultipoles():
    '''
    Copies entire xdwiz.mas file to xd.mas.
    '''
    os.remove('xd.mas')
    shutil.copyfile('xdwiz.mas', 'xd.mas')
    
def wizAddCustomLCS(customLCS):
    '''
    Add dictionary of custom xd.mas atom table lines to current xd.mas file.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        atomTab = False
        
        for line in mas:

            if line.startswith('END ATOM') or line.startswith('!DUM') or line.startswith('DUM'):
                atomTab = False

            if atomTab:
                row = line.upper().split()
                
                if row[0] in customLCS:
                    newmas.write(customLCS[row[0]])
                else:
                    newmas.write(line)
                    
            else:
                newmas.write(line)

            if line.startswith('ATOM     ATOM0'):
                atomTab = True

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

