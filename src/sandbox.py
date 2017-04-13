import os
import codecs

os.chdir('/home/matt/dev/XDTstuff/test')

with open('xd_pro.out', "r",encoding='utf-8', errors='ignore') as fdata:

    z = fdata.read()
    
    print(z)

