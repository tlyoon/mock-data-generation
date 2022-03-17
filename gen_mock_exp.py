#! /bin/bash
import os
import numpy as np
import csv
#import pickle


fn2='gmt_mock.dat'
f=open(fn2,'r')
cextdata=f.readlines()
lcextdata=np.shape(cextdata)[0]
cext={};count=0;
for i in range(lcextdata):
        if i>=2:
            cext[count]=cextdata[i].split()[1]
            count=count+1
f.close()            

## normlization constant
Norm=np.sum([ float(i) for i in cext.values() ])
##

#print(os.listdir())

fn1='experimental.dat'
f=open(fn1,'r')
expdata=f.readlines()
lexpdata=np.shape(expdata)[0]
nf=len(expdata[lexpdata-1].split())
f.close()

try:
    if os.path.isfile('experimental_mock.dat'):
        os.remove('experimental_mock.dat')
except:
    zero=0    

count=0    
try:
    if os.path.isfile('gmt_mock.norm.dat'):
        os.remove('gmt_mock.norm.dat')
except:
    zero=0    
f2=open('gmt_mock.norm.dat','w')
f2.write(cextdata[0])
f2.write(cextdata[1])


f3=open('experimental_mock.dat','w')
f3.write(expdata[0])
f3.write(expdata[1])
f3.close()

with open('experimental_mock.dat', "a",newline='') as file: 
    writer = csv.writer(file, delimiter=' ')
    
    for i in range(lexpdata):
        if i <=1:
            line=expdata[i]
            #writer.writerow(line)
            #print('i:',i, line)
        if i >=2:
            expdatasplit=[ i for i in expdata[i].split() ]
            l1=[ expdatasplit[i] for i in range(0,3) ]
            l4=[ cext[count] ]
            
            l2=[ expdatasplit[i] for i in range(4,nf) ]
            
            joined2= l1 + [ str(float(l4[0])/Norm) ] +  l2
            joined2=joined2[0] + ' ' + joined2[1] + ' ' + joined2[2] + ' ' + joined2[3] \
                + ' ' + joined2[4] + ' ' + joined2[5]
                
            count=count+1
            #joined=l1+l4+l2
            joined=l1+[ str(float(l4[0])/Norm) ]+l2
                
            #print(i,'***',joined2)
            #joined=[ float(i) for i in joined if isinstance(i, float) ]
            
            writer.writerow(joined)
            f2.write(joined2.split()[0] + ' '  + joined2.split()[3] + '\n')
    f.close()

f2.close()
print('The Cext column in', fn1,'has been replaced by that from gmt_mock.dat. The resultant file is saved as experimental_mock.dat.')

print('Copy experimental_mock.dat to overwrite experimental.dat')
import shutil
shutil.copyfile( 'experimental_mock.dat' , 'experimental.dat' )