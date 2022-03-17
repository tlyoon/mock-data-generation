### this script convert a config.dat file to produce output files 'mycal.in' and 'radius.dat'.

import numpy as np
import pandas as pd

df=pd.read_csv('config.dat',header=None,skiprows=0)

radius=df.iloc[0].tolist()[0].split()[-3]

f=open('radius.dat','w')
f.write(str(radius+'\n'))
f.close()

f=open('mycal.in','w')
f.write('** total number of varaibles'+'\n')
f.write('** data begins from line 27 and onwards'+'\n')
f.write('** the rest are all dummies'+'\n')

for i in range(4,27):
    f.write('**'+'\n')
    
count=0    
for i in range(len(df)):
    line=df.iloc[i].tolist()[0]
    x=line.split()[0]
    y=line.split()[1]
    z=line.split()[2]
    count=count+1;
    state='variable '+ str(count) + ' ' + x
    #print(state)
    f.write(state+'\n')
    
    count=count+1;
    #print(count,y);
    state='variable '+ str(count) + ' ' + y
    #print(state)
    f.write(state+'\n')
    
    count=count+1
    #print(count,z);
    state='variable '+ str(count) + ' ' + z
    #print(state)
    f.write(state+'\n')
f.write('end'+'\n')	
f.close()
    

