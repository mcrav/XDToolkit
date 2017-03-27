# -*- mode: python -*-

block_cipher = None


a = Analysis(['XDToolkit.py'],
             pathex=['C:\\Users\\Matthew_2\\Python\\XDToolkit-master','C:\\Users\\Matthew_2\\Anaconda3\\envs\\python34\\Lib\\site-packages','C:\\Users\\Matthew_2\\Anaconda3\\envs\\python34\\Lib','C:\\Users\\Matthew_2\\Anaconda3\\envs\\python34'],
             binaries=[],
             datas=[('tcl','tcl'),
					('tk','tk'),
					('res','res'),
					('LICENSE.txt','.')],
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
          console=False
		)
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='XDToolkit')
