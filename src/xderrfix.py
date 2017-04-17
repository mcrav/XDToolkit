'''
#####################################################################
#-------------------XD ERROR FIXES-----------------------------------
#####################################################################
'''

import os

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
                if line.startswith('END ATOM') or line.startswith('DUM') or line.startswith('!'):
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

                if line.startswith('ATOM     ATOM0'):
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
