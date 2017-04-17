'''
######################################################################
#-------------------XD FILE TOOLS-------------------------------------
######################################################################
'''

import os

def addSnlCutoff(snlmin = 0.0, snlmax = 2.0):
    '''
    Add sin(theta/lambda) cutoffs to xd.mas.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        for line in mas:

            if line.startswith('SKIP'):
                row = str.split(line)
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOn}  {snlMin:<5.3f} {snlMax:<5.3f}'.format(*row, snlOn = '*sinthl', snlMin = snlmin, snlMax = snlmax)
                newmas.write(rowStr + '\n')

            else:
                newmas.write(line)

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def setupmas():
    '''
    Setup details in xd.mas for refinement.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        atomTab = False

        for line in mas:

            row = str.split(line)
            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False

            if atomTab:
                row = str.split(line)

                row[10] = '4'

                if line.startswith('H('):
                    row[7] = '1'

                if len(row)==13:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}{12}'.format(*row)
                else:
                    rowStr = '{0:9}{1:10}{2:3}{3:9}{4:9}{5:4}{6:4}{7:3}{8:4}{9:4}{10:3}{11:10}'.format(*row)
                newmas.write(rowStr + '\n')

            #Set model to multipoles max l = 4
            elif row[0:2] == ['SELECT','*model']:
                row = str.split(line)
                rowStr = row[0] + ' ' + row[1] + '  4  ' + row[3] + '  ' + row[4] + '  ' + row[5] + ' ' + row[6] + ' ' + 'F^2' + '  ' + row[8] + ' ' + row[9] + ' ' + row[10]
                newmas.write(rowStr + '\n')

            elif line.startswith('!RESET    bond C(1) H(1) 1.09 ...'):
                pass

            #Add appropriate fmod1
            elif line.startswith('FOUR') or line.startswith('!FOUR'):
                row = str.split(line)
                rowStr = 'FOUR  fmod1  4  2  0  0   fmod2 -1  2  0  0'
                newmas.write(rowStr + '\n')

            #Add convcrit
            elif row[0:2] == ['SELECT','cycle']:

                row[2] = '25'
                row[11] = '*convcrit'
                row[12] = '0.1E-3'
                rowStr = ' '.join(row)
                newmas.write(rowStr + '\n')

            #Turn on rigid bond test
            elif line.startswith('!DMSDA'):
                row = str.split(line)
                rowStr = 'DMSDA    ' + row[1] + '   ' + row[2]
                newmas.write(rowStr + '\n')

            elif line.startswith('SKIP'):
                row = str.split(line)
                rowStr = '{0:7}{1:5}{2} {3} {4:9}{5} {6} {snlOff}  {8} {9}'.format(*row, snlOff = 'sinthl')
                newmas.write(rowStr + '\n')

            else:
                newmas.write(line)

                #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
            if line.startswith('ATOM     ATOM0'):
                atomTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def resetKeyTable():
    '''
    Reset key table to all 0s.
    '''
    keyTab = False

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        #Go through xd.mas and flip keyTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('EXTCN'):
                keyTab = False

            if keyTab:
                if line[:5] != 'KAPPA':
                    row = str.split(line)

                    rowStr = '{0:8}{1}'.format(row[0], '000 000000 0000000000 000000000000000 00 000 00000 0000000 000000000')
                    newmas.write(rowStr + '\n')
                else:
                    newmas.write('KAPPA   000000\n')

            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def multipoleKeyTable():
    '''
    Write key table with multipoles based on SITESYM column of atom table.
    '''
    #Initialise dictionary containing multipole KEY table settings for common local symmetries.
    multipoleBank = {'NO': '00 000 00000 0000000 000000000', '1': '10 111 11111 1111111 111111111',
                     'CYL': '10 001 00000 0000000 000000000', 'CYLX': '10 001 10000 1000000 000000000',
                     'CYLXD': '10 001 10000 1000000 000000000', '2': '10 001 10010 1001000 100100010',
                     'M': '10 110 10011 0110011 100110011', 'MM2': '10 001 10010 1001000 100100010',
                     '4': '10 001 10000 1000000 100000010', '4MM': '10 001 10000 1000000 100000010',
                     '3': '10 001 10000 1000010 100001000', '3M': '10 001 10000 1000010 100001000',
                     '6': '10 001 10000 1000000 100000000', '6MM': '10 001 10000 1000000 100000000'}

    XAtoms = ('CL', 'F(', 'BR','I(', 'O(', 'N(')
    noHexaAtoms = ['C(', 'O(', 'N(']

    with open('xd.mas','r') as mas:

        atomTab = False             #Initialise atomTab and keyTab bools to detect atom table and key table
        keyTab = False
        newMultipoles = {}          #Initialise dict to store atom labels and their corresponding multipoles
        atomSyms = {}

        #Go through xd.mas and flip atomTab to true when you reach the start of the atom table and false when you reach the end of the atom table
        for line in mas:

            if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
                atomTab = False
        #In atom table, if a row has no chemcon (12 items in row) add to dictionary for printing to new file.
        #Atoms with CHEMCON will inherit multipole configuration if everything is set to 0.
            if atomTab:
                row = str.split(line.upper())

                if len(row)==12:
                    atomSyms[row[0]] = row[11]
                    if line[:2] not in XAtoms:
                        if line[:2] not in noHexaAtoms:
                            #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                            newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                        else:
                            #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                            newMultStr = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                            #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                            newMultipoles[row[0]] = newMultStr[:-9] + '000000000'
                    #If cyl is sym label for a halogen include z2 quadrupole
                    else:
                        if row[11] == 'CYL':
                            if line[:2] != 'F(':
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLX'])
                            else:
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLXD'])
                        else:
                            newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])

                elif len(row) == 13:
                    if row[11] != atomSyms[row[12]]:
                        if line[:2] not in XAtoms:
                            if line[:2] not in noHexaAtoms:
                                #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])

                            else:
                                #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                                newMultStr = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                                #Add new key table line to dict. Multipoles selected from multipoleBank based on SITESYM
                                newMultipoles[row[0]] = newMultStr[:-9] + '000000000'

                        #If cyl is sym label for a halogen include z2 quadrupole
                        else:
                            if row[11] == 'CYL':
                                if line[:2] != 'F(':
                                    newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLX'])
                                else:
                                    newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank['CYLXD'])
                            else:
                                newMultipoles[row[0]] = '{0:8}{1} {2}'.format(row[0], '000 000000 0000000000 000000000000000', multipoleBank[row[11]])
                            print(newMultipoles[row[0]])

            if line.startswith('ATOM     ATOM0'):
                atomTab = True

    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:
        #Find key table the same as with atomTab
        for line in mas:

            if line.startswith('KAPPA'):
                keyTab = False

            #In the key table, if the atom needs multipoles they are written
            #otherwise write all 0s and atom will inherit multipoles from CHEMCON atom.
            if keyTab:
                row = str.split(line)
                if row[0].upper() in newMultipoles:
                    newmas.write(newMultipoles[row[0].upper()] + '\n')
                else:
                    newmas.write('{0:8}{1}\n'.format(row[0].upper(), '000 000000 0000000000 000000000000000 00 000 00000 0000000 000000000'))
            else:
                newmas.write(line)

            if line.startswith('KEY     XYZ'):
                keyTab = True

    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')
