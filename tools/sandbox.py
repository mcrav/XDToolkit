

def moveLigandsCloserGrd(grdFile, fractionCloser, atoms):
    
    with open(grdFile,'r') as grd, open('new_' + grdFile,'w') as newGrd:
        i=-1000
        atomBool = False
        for line in grd:
            
            if line.startswith('! Connections'):
                atomBool = False
            
            if atomBool and i>1:

                row = line.split()
                if row[0] in atoms:
                    row[1] = float(row[1])*fractionCloser
                    row[2] = float(row[2])*fractionCloser
                    row[3] = float(row[3])*fractionCloser
                    print('{}\n{}\n{}\n'.format(row[1], row[2], row[3]))
                    newLine = '{0:<10}{1:< 10.5f}{2:< 10.5f}{3:< 9.5f}{4}\n'.format(row[0], row[1], row[2], row[3], row[4])
                    newGrd.write(newLine)
                else:
                    newGrd.write(line)
                continue
            else:
                newGrd.write(line)
            
            if line.startswith('! Objects'):
                atomBool = True
                i=0
                
            i+=1