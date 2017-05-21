# -*- mode: python -*-

block_cipher = None

a = Analysis(['XDToolkit.py'],
             pathex=[],
             binaries=[],
             datas=[('/home/matt/dev/XDToolkit/res','res'),
		                ('/home/matt/dev/XDToolkit/LICENSE.txt','.')],
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
