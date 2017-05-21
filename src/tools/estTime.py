

def timeForRef():
    
    file = open('lsmTimes.buckfast','r')
    timeRef = {}
    
    for line in file:
        if not line.startswith('  '):
            row = str.split(line)
            time = float(row[2])
            atoms = int(row[1])
            if row[0] in timeRef.keys():
                timeRef[row[0]][0].append(time)
                timeRef[row[0]][1].append(atoms)
            else:
                timeRef[row[0]] = ([time],[atoms])
            
    print (timeRef)
    return(timeRef)
    file.close()
    
import numpy as np
import matplotlib.pyplot as plt

reflabs = ('s','ha','la','mk','m','nhpam')
data = timeForRef()
timeEst = []
for ref in reflabs:
    
    N = 50
    y = np.array(data[ref][0])
    x = np.array(data[ref][1])
    colors = ['000000']
    area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radii
    
    plt.scatter(x,y, s=area, c=colors, alpha=0.5)
    plt.ylim([0,100])
    plt.xlim([0,100])
    plt.show()
    
    m, b = np.polyfit(x, y, 1)
    
    plt.plot(x, y, '.')
    plt.plot(x, m*x + b, '-')
    estStr = '{0} = {1:.3f}*numAtoms {2:+.3f}'.format(ref,m,b)
    timeEst.append(estStr)

print(timeEst)

#Scale time = 0.837 * numAtoms - 7.174
#High angle time = 1.812 * numAtoms -18.0969
#Low angle time = 0.114 * numAtoms 1.626