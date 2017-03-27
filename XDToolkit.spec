# -*- mode: python -*-

block_cipher = None


a = Analysis(['XDToolkit.py'],
             pathex=['/home/matt/dev/XDToolkit'],
             binaries=[],
             datas=[('res/flatearth.ico','res'),
		    ('res/folder.ico','res'),
                    ('res/mercury_48x48.png','res'),
                    ('res/molecoolico.png','res'),
                    ('res/notepadico.png','res'),
                    ('res/play.ico','res'),
                    ('res/refreshico.png','res'),
		    ('res/settings.png','res'),
                    ('res/XD Toolkit Manual.pdf','res')],
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
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='XDToolkit')
