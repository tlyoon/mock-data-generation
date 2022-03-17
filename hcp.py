## Mandatory requirement: presence of experimental.template.dat

## this script generates random hcp configuration. 
## Input required: number of atoms, noa
## Input required: lattice constatn, lc. 
## Input required: real part of the sphere's refractive index, n
## Input required: imaginary part of the sphere's refractive index, k
## 
## Refer to Piri, N., Shams-Nateri, A., & Mokhtari, J. (2015) for choosing default values for the size and reflective index of the sphere. 
## All sizes are to be expressed in unit of nm. 
## Default size (inferred from Piri): 100 nm
## Default reflective index (inferred from Piri): n=1.5; k=0.001



import itertools
import numpy as np
import random

noa=14  ## default 14
lc=100  ## default 200. Common diameter of the spheres
n=1.15   ## Real part of reflective index
nk=9.5 ## Imaginary  part of reflective index



### don't touch anything below 
f=open('experimental.template.dat','r')
nwl=len(f.readlines())-2
f.close

import time
strings = time.strftime("%m,%d,%H,%M,%S")
strings=strings.replace(',','')
state='nwl.'+str(nwl)+'.noa.'+str(noa)+'.lc.'+str(lc)+'.n.'+str(n)+'.nk.'+str(nk)+'.'+strings
#print(state)
f=open(state,'w')
f.close()



rc=lc/2

def hcp(i,j,k,r):
    x=3*i + (j+k)%2
    y=np.sqrt(3)*( j + (1/3)*(k%2) )
    z=(2/3)*np.sqrt(6)*k
    return x*r,y*r,z*r
    
def coor(i,j,k):
    return hcp(i,j,k,rc)

def dist(i,j,k,ip,jp,kp):
    dist=np.sqrt(np.sum([ (coor(i,j,k)[q]-coor(ip,jp,kp)[q])**2 for q in [0,1,2] ]))
    return dist

def dist1(i,j,k):
    return np.sqrt(np.sum([ coor(i,j,k)[q]**2 for q in [0,1,2] ]))


fill=[]
dstop=10000
#i=random.randint(-1,1)
ip=0;i=0
jp=0;j=0
kp=0;k=0


while len(fill) < noa:
    while round(dstop/lc,1) > 1.4:
        rchoice=random.choice([0,1,2])        
        a = [i,j,k]
        a[rchoice] = a[rchoice] + random.randint(-1,1)
        ip=a[0]
        jp=a[1]
        kp=a[2]
        #dstop=round(dist(i,j,k,ip,jp,kp),2)
        dummy=round(dist(i,j,k,ip,jp,kp),1)
        if dummy==0:
            dstop=10000
        else:
            dstop=dummy
            #print(i,j,k,ip,jp,kp,round(dist(i,j,k,ip,jp,kp)/(lc),1))
            
    #print(round(dist(i,j,k,ip,jp,kp)/lc,1)<=1.4)
#    print(i,j,k,ip,jp,kp,round(dist(i,j,k,ip,jp,kp)/lc,1))
    #if round(dist(i,j,k,ip,jp,kp)/lc,1)<=1.86:
    #i,j,k = a        #print('')
    #print(i,j,k,ip,jp,kp)
    #print('*** ',i,ip,j,' ', jp,k,kp)
    #print('exit: dstop',dstop)
    fill.append([i,j,k])
    fill.append([ip,jp,kp])
        #print(i,j,k,ip,jp,kp,round(dist(i,j,k,ip,jp,kp)/(lc),1))
    #print('fill, before',fill)
    fill.sort()
    fill=list(k for k,_ in itertools.groupby(fill))
    #print('fill, after',fill)
    dstop=10000
    
    rchoice=random.choice([0,1,2])        
    b = [i,j,k]
    b[rchoice] = b[rchoice] + random.randint(-1,1)
    i=b[0]
    j=b[1]
    k=b[2]
    #i,j,k = a

    #print('last:',i,j,k)
    #print('')
'''
print('len(fill)',len(fill))
print('fill:',fill)
'''
if len(fill)!=noa:
    fill=fill[:noa]
'''
print(' ')
print('fill:',fill)
print('len(fill):',len(fill))
print('')
'''
f=open('config.xyz','w')
f.write(str(noa)+'\n')
f.write('spheres in hcp'+'\n')

fc=open('config.dat', 'w')
for i in fill:
    x=coor(i[0],i[1],i[2])[0]
    y=coor(i[0],i[1],i[2])[1]
    z=coor(i[0],i[1],i[2])[2]
    #print(i, coor(i[0],i[1],i[2]))
    line= 'C ' + str(x) + ' ' + str(y) + ' ' +  str(z)
    #print(line)
    f.write(line +'\n')
    
    linec=str(x) + ' ' + str(y) + ' ' +  str(z) + ' ' + str(rc) + ' ' + str(n) + ' ' + str(nk)
    fc.write(linec +'\n')
f.close()
fc.close()

#print('')
f=open('config.xyz','r')
rawdata=f.readlines()
for i in rawdata:
    zero=0
#    print(i.strip())
f.close()



#print('')
f=open('config.dat','r')
rawdata=f.readlines()
for i in rawdata:
    zero=0
#    print(i.strip())
f.close()



