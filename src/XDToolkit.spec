
#Builds on windows on my account, with my python34 environment. 
#To build on another computer the data paths will need to be changed and the python environment must be <=3.4
#And have all necessary packages installed i.e. numpy, scipy, matplotlib, PyQt etc.
#On my computer it is necessary to include the tk and tcl folders from the python environment folder
#This may not always be necessary but you can try building without them and then drag them into the dist folder
#if the exe doesn't work.

block_cipher = None

a = Analysis(['XDToolkit.py'],
             pathex=[],
             binaries=[],
             datas=[('res','res'),
		                ('C:/Users/Matthew_2/Desktop/XDToolkit-master/LICENSE.txt','.'),
                        ('C:/Users/Matthew_2/Anaconda3/envs/python34/Lib/tkinter','tk'),
                        ('C:/Users/Matthew_2/Anaconda3/envs/python34/tcl','tcl')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='XDToolkit',
          debug=False,
          strip=False,
          upx=True,
          console=False )

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='XDToolkit')
