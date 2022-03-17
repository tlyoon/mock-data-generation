## Mandatory requirement: the presence of correct versino config.dat

## This script creates a new file experimental.dat from experimental.template.dat
## and replace the real, imaginary index in it with that from config.dat

import numpy as np, os

### config.dat ####
f=open('config.dat','r')
rawdata=f.readlines()
f.close()
cdata={}
count=0
for i in rawdata:
    cdata[count]=i.strip()
    count=count+1
nr=cdata[0].split()[-2]    
ni=cdata[0].split()[-1]    
### config.dat ####

#### experimental.dat
f=open('experimental.template.dat','r')
rawdata=f.readlines()
f.close()

f=open('experimental.dat','w')
edata={}
count=0
for i in rawdata:
    edata[count]=i.strip()
    line=edata[count].split()
    if count>=2:
        wl=line[0]; Csca=line[4];Cabs=line[5]
        
        newline=wl + ' ' + str(nr) + ' ' + str(ni) + ' Cext    Csca    Cabs'
        f.write(newline+'\n')
    else:
        f.write(edata[count] +'\n')
    count=count+1
ledata=len(edata)-1
#print(edata[ledata])
f.close()
#### experimental.dat
