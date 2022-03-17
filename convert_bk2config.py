#! /bin/bash
# convert a bk72s.k file into config.xyz and config.dat

import os
import numpy as np
import csv

file2open='bk7s2.k'
f=open(file2open,'r')
rawdata=f.readlines()
f.close()

file2open='config.xyz'
if os.path.isfile(file2open):
    os.remove(file2open)

fc=open('config.dat', 'w')

f=open(file2open, 'w')
noa=rawdata[1].strip()
comment='sample comment line'

f.write(noa+'\n')
f.write(comment+'\n')

comment='sample comment line'
label='C '
#print(noa)
#print(comment)
for i in rawdata[2:]:
    i=i.strip()
    line=label + i.split()[0]+' '+i.split()[1]+' '+i.split()[2]
    f.write(line+'\n')
    
    #linec=i.split()[0]+' '+i.split()[1]+' '+i.split()[2]+' '+i.split()[3]
    linec=i
    fc.write(linec+'\n')
    
    #print(line)
f.close()    
fc.close()


print('')
print('bk7s2.k has been converted into config.xyz and config.dat')
print('The content of config.xyz is as follows')
print('')
f=open(file2open,'r')
rawdata=f.readlines()
for i in rawdata:
    print(i.strip())
f.close()


print('')
print('The content of config.dat is as follows')
print('')
fc=open('config.dat','r')
rawdata=fc.readlines()
for i in rawdata:
    print(i.strip())
fc.close()






