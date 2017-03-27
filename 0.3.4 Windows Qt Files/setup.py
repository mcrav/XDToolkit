import sys
from cx_Freeze import setup, Executable

includefiles = [
        ('qwindows.dll','platforms\qwindows.dll'),
        ('flatearth.ico','flatearth.ico'),
        ('imageformats/','imageformats/'),
        ('XD Toolkit Manual.pdf','XD Toolkit Manual.pdf')
]

setup(
    name = "XD Toolkit",
    version = "1.0",
    description = "A toolkit for working with WinXD",
    options = {'build_exe': {'include_files':includefiles}}, 
    executables = [Executable("XD Toolkit.py", base = "Win32GUI")])