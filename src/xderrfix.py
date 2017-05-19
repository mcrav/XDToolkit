'''
#####################################################################
#-------------------XD ERROR FIXES-----------------------------------
#####################################################################
'''

import os
from utils import getCellParams, atomTableBegins, atomTableEnds

def check4errors():
    '''
    Check for errors in xd_lsm.out. Return errors.
    '''
    with open('xd_lsm.out','r') as lsm:
        complete = False
        error = ''

        for line in lsm.readlines()[-20:]:

            if line.startswith(' | Total program runtime'):
                complete = True
            elif line.strip().startswith('* * * '):
                error = line
                complete = False
            elif line.startswith('Noble gas configuration not recognized for element Cu'):
                error = line
                complete = False

    return(complete,error)

def fixLsmCif():
    '''
    Add unit cell paramaters to xd_lsm.cif.
    '''
    cellParams = getCellParams()

    with open('xd_lsm.cif','r') as cif, open('xd_lsmnew.cif','w') as newcif:

        for line in cif:

            if line.startswith('# Refinement details'):

                cellStr = '\n#Unit Cell Parameters\n\n'
                cellStr += '{0:35}{1}\n'.format('_cell_length_a', cellParams[0])
                cellStr += '{0:35}{1}\n'.format('_cell_length_b', cellParams[1])
                cellStr += '{0:35}{1}\n'.format('_cell_length_c', cellParams[2])
                cellStr += '{0:35}{1}\n'.format('_cell_angle_alpha', cellParams[3])
                cellStr += '{0:35}{1}\n'.format('_cell_angle_beta', cellParams[4])
                cellStr += '{0:35}{1}\n'.format('_cell_angle_gamma', cellParams[5])
                cellStr += '\n\n'
                newcif.write(cellStr)

            newcif.write(line)

    os.remove('xd_lsm.cif')
    os.rename('xd_lsmnew.cif','xd_lsm.cif')

def removePhantomAtoms():
    '''
    Remove faulty atoms with format '(e' from xd.mas and xd.inp.
    '''
    phantoms = False

    with open('xd.inp','r') as inp, open('xdnew.inp','w') as newinp:
        i = 3
        phantomCount = 0

        for line in inp:
            if line.startswith('('):
                phantoms = True
                phantomCount += 1
                i = 0
            if i > 2:
                newinp.write(line)
            i+=1

    if phantoms:
        with open('xd.mas', 'r') as mas, open('xdnew.mas', 'w') as newmas:

            atomTab = False
            keyTab = False
            scatTab = False
            kappaPrinted = False
            k = 0

            for line in mas:
                if atomTableEnds(line):
                    atomTab = False
                elif line.startswith('KAPPA'):
                    keyTab = False
                    if kappaPrinted == False:
                        if k < phantomCount:
                            k += 1
                            continue

                elif line.startswith('END SCAT'):
                    scatTab = False

                if atomTab:
                    if line.startswith('('):
                        continue
                    else:
                        newmas.write(line)

                elif scatTab:
                    if line.startswith('('):
                        continue
                    else:
                        newmas.write(line)

                elif keyTab:
                    if line.startswith('('):
                        continue
                    else:
                        newmas.write(line)

                else:
                    newmas.write(line)

                if atomTableBegins(line):
                    atomTab = True

                elif line.startswith('KEY'):
                    keyTab = True

                elif line.startswith('SCAT'):
                    scatTab = True

        os.remove('xd.mas')
        os.rename('xdnew.mas','xd.mas')

        numElements = 0

        with open('xdnew.inp','r') as inp, open('xd.inp','w') as newinp:
            i = 0
            j = 1
            for line in inp:
                if line.startswith('USAGE'):
                    row = str.split(line)
                    numAtoms = int(row[1])
                    numElements = int(row[4])
                    row[1] = str(numAtoms - phantomCount)
                    row[4] = str(numElements - phantomCount)
                    rowStr = '{0:10}{1:6}{2:3}{3:4}{4:4}{5:4}{6:4}{7:4}{8:4}{9:4}{10:4}{11:6}{12:3}{13:4}{14}\n'.format(*row)
                    newinp.write(rowStr)
                else:                                               #Code to remove extra line(s) of kappa table in inp
                    if j <= (numElements - phantomCount) and i > 2:
                        newinp.write(line)

                    elif i <= 2:
                        newinp.write(line)

                    else:
                        i = -1

                if line[:1] == ' ':
                    i+= 1
                else:
                    i = 0

                if i > 2:
                    j+=1

        os.remove('xdnew.inp')

def fixBrokenLabels(atomLabsDict):
    '''
    Change problematic atom labels in shelx.ins.
    '''
    try:
        neebs = atomLabsDict.keys()
        i = 0
        atomBool = False
        elements = []

        with open('shelx.ins','r') as ins, open('shelxnew.ins','w') as newins:

            for line in ins:
                i = 0
                if line.startswith('HKLF'):
                    atomBool = False

                elif line.startswith('SFAC'):
                    row = str.split(line)
                    elements = [item.upper() for item in row[1:]]

                if atomBool:
                    if line[:1].isalpha() and not line.startswith('AFIX'):
                        row = str.split(line)
                        atomLab = row[0]
                        newLab = atomLab.replace('(','').replace(')','')
                        #Fix labels that are just N or C or Co with no number
                        if newLab in elements:       #N Na error
                            newerLab = newLab + str(i)
                            while newerLab in neebs:
                                i += 1
                                newerLab = newLab + str(i)
                            newLab = newerLab
                        newLine = '{0:6}{1}'.format(newLab.upper(), line[6:])
                        newins.write(newLine)
                    else:
                        newins.write(line)

                else:
                    newins.write(line)

                if line.startswith('FVAR'):
                    atomBool = True

        os.remove('shelx.ins')
        os.rename('shelxnew.ins', 'shelx.ins')

    finally:
        try:
            os.remove('shelxnew.ins')
        except Exception:
            pass

def addNCST():
    '''
    Fix 'paramter/ncst size should be increased' error.
    '''
    with open('xd.mas','r') as mas, open('xdnew.mas','w') as newmas:

        i = 40

        for line in mas:                                #Go through xd.mas line by line

            #Add SIZE ncst 30 line at appropriate place
            if line.startswith('SAVE'):
                newmas.write(line)
                newmas.write('SIZE ncst {}\n'.format(i))

            elif not line.startswith('SIZE ncst {}\n'.format(i)):
                newmas.write(line)

    #Create new xd.mas file
    os.remove('xd.mas')
    os.rename('xdnew.mas','xd.mas')

def fixCuNobleGasMas(masFile):
    '''
    Give Cu 4s1 3d10 configuration in scattering table in given mas file.
    3d10 written as valence electrons.
    '''
    with open(masFile,'r') as mas, open('xdnew.mas','w') as newmas:

        scatTab = False

        for line in mas:

            if line.startswith('END SCAT'):
                scatTab = False

            if scatTab:
                if line[:2].upper() == 'CU':
                    row = str.split(line)
                    row[11] = '-10'

                    i = 0
                    for item in row:
                        if 3 < i < 25:
                            row[i] = float(item)
                        i += 1
                    rowStr = '{0:5}{1:5}{2:5}{3:7}{4:< 4.0f}{5:< 4.0f}{6:< 4.0f}{7:< 4.0f}{8:< 4.0f}{9:< 4.0f}{10:< 4.0f}{11:< 4.0f}{12:< 4.0f}{13:< 4.0f}{14:< 4.0f}{15:< 4.0f}{16:< 4.0f}{17:< 4.0f}{18:< 4.0f}{19:< 4.0f}{20:< 4.0f}{21:< 3.0f}{22:< 8.4f}{23:< 8.4f}{24: .3f}\n'.format(*row)
                    newmas.write(rowStr)
                else:
                    newmas.write(line)

            else:
                newmas.write(line)
            if line.startswith('SCAT'):
                scatTab = True

    os.remove(masFile)
    os.rename('xdnew.mas',masFile)
    
def fixCuNobleGasError():
    '''
    Fix 'noble gas configuration for CU' error.
    '''
    print('~~~~~~~~~~~~~~~~~\n\nFIXING CU NOBLE GAS CONFIGURATION ERROR\n\n~~~~~~~~~~~~~~~~~~~')
    fixCuNobleGasMas('xd.mas')
    fixCuNobleGasMas('xdwiz.mas')
    
    #Add 10 to Cu valence monopole population in inp file.
    with open('xd.inp','r') as inp, open('xdnew.inp','w') as newinp:

        i = 0
        cuFound = False

        for line in inp:
            if line[:2].upper() == 'CU':
                cuFound = True

            if cuFound:
                i+=1

            if i==3:
                row = str.split(line)
                row[0] = 10.0000
                rowStr = '{0:< 9.4f}{1}\n'.format(*row)
                newinp.write(rowStr)
                cuFound = False
                i = 0

            else:
                newinp.write(line)

    os.remove('xd.inp')
    os.rename('xdnew.inp','xd.inp')

