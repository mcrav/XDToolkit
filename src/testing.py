import os
from test import Ui_MainWindow
from XDToolkit import rawInput2labels, customExceptHook
from PyQt5.QtWidgets import QWidget, QMessageBox, QLabel, QGridLayout, QDialogButtonBox, QSplashScreen, QPushButton, QApplication, QDialog, QLineEdit, QFileDialog, QMainWindow
from PyQt5.QtCore import QCoreApplication, QSettings, QThread, pyqtSignal, Qt
from PyQt5.QtGui import QPixmap, QFont
import PyQt5.QtCore
import sys
import shutil
from ast import literal_eval
from utils import spec2norm
from initfuncs import ins2all

def listjoin(listItem, splitter = ' '):
    return splitter.join(listItem).strip(splitter)

def addHTMLPre(string):
    return '<pre>' + string + '</pre>'

def removeHTMLPre(string):
    return string.strip('<pre>').strip('</pre>')

def addHTMLBold(string):
    return '<b>' + string + '</b>'
#
#os.chdir('/home/matt/dev/XDTstuff')
#
#def getTestDict(bondFile):
#    '''
#    Get test dictionary from bond file.
#    '''
#    neebs = {}
#    print(bondFile)
#    with open(bondFile,'r') as bonds:
#
#        for line in bonds:
#            if line[0].isdigit():
#                row = str.split(line)
#
#                neebs.setdefault(row[1],[]).append(row[2])
#                neebs.setdefault(row[2],[]).append(row[1])
#
#    return neebs
#
#

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
        atomLabs[atom] = tuple(sorted([item.split(',')[0] for item in neebs]))

    for atom,neebs in testDictRaw.items():
        testDict[atom] = tuple(sorted(neebs))

    for atom,neebs in atomLabs.items():
        if not neebs == testDict[atom]:
            problems[atom] = neebs

    print(problems)
    return problems

def getDistances(atom, dists):
    
    for pair, dist in dists.items():
        if atom in [spec2norm(item) for item in pair]:
            print(pair)
            print(dist)
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
#
#
#def checkAllNeebs():
#    '''
#    Test all neighbours in testData folder.
#    '''
#    problems = {}
#
#    for folder in os.listdir('testData'):
#        if folder != 'alanyl-methionine':
#
#            os.chdir('/home/matt/dev/XDTstuff/testData/{}'.format(folder))
#            print(os.getcwd())
#            testDict = getTestDict('bonds.tsv')
#
#            x = checkNeebs(testDict)
#
#            if x:
#                problems[folder] = x
#                print(x)
#
#    for item in problems:
#        print('\n\n')
#        print('------------------------------------------')
#        print(item.upper())
#        print('')
#        print(problems[item])
#    return problems
#
#
#def checkFolder(folder):
#    '''
#    Test specific folder.
#    '''
#    os.chdir('/home/matt/dev/XDTstuff/testData/{}'.format(folder))
#    print(checkNeebs(getTestDict('bonds.tsv')))

class TestMainWindow(QMainWindow, Ui_MainWindow):
    '''
    Main window of testing program.
    '''
    def __init__(self, parent=None):
        super(TestMainWindow, self).__init__(parent)
        self.setupUi(self)

        self.setFolderBut.clicked.connect(self.initialiseFolder)
        self.addNeebsBut.clicked.connect(self.addNeebs)
        self.neebInput.returnPressed.connect(self.addNeebs)
        self.testAllBut.clicked.connect(self.atomLabsTestAll)
        self.atomLabs = {}

        self.testBut.clicked.connect(self.runAtomLabsTest)

        self.atomsDrop.currentIndexChanged.connect(self.updateNeebInput)

    def updateNeebInput(self):
        atom = str(self.atomsDrop.currentText())
        if atom in self.atomLabs:
            self.neebInput.setText(' '.join(self.atomLabs[atom]).strip())
        else:
            self.neebInput.setText('')

    def initialiseFolder(self):
        folder = str(QFileDialog.getExistingDirectory(None, "Select Directory"))

        folderName = folder.split('/')[-1]
        testFolderAbsPath = '/home/matt/dev/XDTstuff/testData/{}'.format(folderName)

        if not os.path.isdir(testFolderAbsPath):
            os.makedirs(testFolderAbsPath)

            os.chdir(testFolderAbsPath)
            shutil.copyfile('{}/shelx.ins'.format(folder), 'shelx.ins'.format(folderName))
            shutil.copyfile('{}/shelx.hkl'.format(folder), 'shelx.hkl'.format(folderName))

            with open('atomLabs.test','w') as atomLabs:
                atomLabs.write('{}')
                self.atomLabs = {}

        else:
            os.chdir(testFolderAbsPath)
            if os.path.isfile('atomLabs.test'):
                with open('atomLabs.test','r') as atomLabs:
                    self.atomLabs = literal_eval(atomLabs.read())

        self.atomList = self.getAtomList()
        self.atomsDrop.clear()
        self.atomsDrop.addItems(self.atomList)

        self.addedNeebsLab.setText('')
        self.testLab.setText('')

    def getAtomList(self):
        '''
        Get list of atoms in asymmetric unit.
        '''
        atoms = []

        if os.path.isfile('shelx.ins'):
            atomBool = False

            with open('shelx.ins','r') as ins:

                for line in ins:
                    if line.startswith('HKLF') or line.startswith('REM'):
                        atomBool = False

                    if atomBool and line[:1].isalpha() and not line.startswith('AFIX'):
                        row = line.split()
                        atoms.extend(rawInput2labels(row[0]))

                    if line.startswith('FVAR'):
                        atomBool = True

        return atoms

    def addNeebs(self):
        '''
        Add user inputted neighbours to test dictionary file.
        '''
        with open('atomLabs.test','w') as atomLabs:

            addedNeebs = rawInput2labels(str(self.neebInput.text()))
            currAtom = str(self.atomsDrop.currentText())
            self.atomLabs[currAtom] = addedNeebs
            print(self.atomLabs)
            atomLabs.write(str(self.atomLabs))
            self.addedNeebsLab.setText('Atom: {}       Neighbours: {}'.format(currAtom, addedNeebs))
            self.atomsDrop.setCurrentIndex(self.atomsDrop.currentIndex() + 1)

    def runAtomLabsTest(self):
        '''
        Test user created test dict against ins2all generated atomLabs dict.
        '''
        tildaBreak = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~<br><br>'
        self.testLab.setText('Running test')
        results = checkNeebs(self.atomLabs)
        if results:
            resStr = 'Differences: {}<br><br>{}'.format(listjoin(results, ', '), tildaBreak)
            for atom, neebs in results.items():
                resStr += 'Atom: {}<br><br>'.format(atom)
                resStr += '{0:<30}{1}<br>'.format('ins2all neighbours:', listjoin(neebs, ', '))
                resStr += '{0:<30}{1}<br><br>'.format('User inputted neighbours:', listjoin(self.atomLabs[atom], ', '))
                resStr += tildaBreak

            resStr = '<pre>' + resStr + '</pre>'
            self.testLab.setText(resStr)

        else:
            self.testLab.setText('No problems')

    def atomLabsTestAll(self):
        '''
        Run atomLab test on everything in testData folder.
        '''
        problems = []
        resStr = ''
        testDataDir = '/home/matt/dev/XDTstuff/testData'
        for item in os.listdir(testDataDir):

            self.atomLabs = {}
            os.chdir(testDataDir + '/' + item)

            if os.path.isfile('atomLabs.test'):
                with open('atomLabs.test','r') as atomLabs:
                    self.atomLabs = literal_eval(atomLabs.read())
            result = checkNeebs(self.atomLabs)
            if result:
                problems.append(item)
                comment = 'Fail'
            else:
                comment = 'Pass'

            resStr = removeHTMLPre(resStr)
            newRes = '{0:20}  {1}'.format(item, comment)

            if comment == 'Fail':
                resStr += ' <br><b>{}</b> '.format(newRes)
            else:
                resStr += ' <br>{}'.format(newRes)

            resStr = addHTMLPre(resStr)
            self.testLab.setText(resStr)
            self.testLab.repaint()
            QApplication.processEvents()

        resStr += '<br>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~<br><i>TESTS FINISHED</i>'
        self.testLab.setText(resStr)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    sys.excepthook = customExceptHook
    prog = TestMainWindow()
    app.aboutToQuit.connect(app.deleteLater)
    prog.show()
    sys.exit(app.exec_())

#checkAllNeebs()
#checkFolder('kmno4')
