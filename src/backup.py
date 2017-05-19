import os
from shutil import copyfile

def backup(folderName):
    '''
    Copy refinement files to given folder.
    '''
    folder = 'Backup/' + folderName
    if not os.path.isdir('{}{}'.format(os.getcwd(), '/Backup')):                   #Check to see if backup folder exists
        os.makedirs('Backup/')                          #If it doesn't exist, make it

    #Create string of new folder path
    if not os.path.isdir('{}{}{}'.format(os.getcwd(), '/',folder)):                      #Check if new folder exists
        os.makedirs(folder)                             #If it doesn't exist, make it

    backupFiles = ('xd.mas','xd.inp','xd.res','xd_lsm.out','xd.cov', 'xd.fou')

    for file in backupFiles:
        try:
            copyfile(file, folder + '/' + file)
        except Exception:
            pass

def loadBackup(folderName):
    '''
    Copy files from given folder to project foler.
    '''
    folder = folderName                                 #folderName is the absolute path to the backup folder

    try:                                                #Copy files from backup folder to working directory
        copyfile(folder + '/xd_lsm.out','xd_lsm.out')   #try and except used in case one of the files doens't exist
    except:
        pass
    try:
        copyfile(folder + '/xd.mas','xd.mas')
    except:
        pass
    try:
        copyfile(folder + '/xd.inp','xd.inp')
    except:
        pass
    try:
        copyfile(folder + '/xd.res','xd.res')
    except:
        pass
