import os 

fnout='experimental_trimmed.'
fnin='realistic_experimental.dat'

f=open(fnin,'r')
rawdata=f.readlines()
f.close()
lr=len(rawdata)
#print(lr)

f=open(fnout,'w')
f.write(rawdata[0].strip()+'\n')
f.write(rawdata[1].strip()+'\n')

count=0
nskip=10  ### nskip=10 gives 57 wavelengths. Increase nskip if you need less number of wavelengths
for i in range(2,len(rawdata),nskip):
    line=rawdata[i].strip()
    f.write(line+'\n')
    count=count+1
f.close()

fnout2=fnout+str(count)+'.dat'
os.rename(fnout, fnout2)

f=open(fnout2,'r')
rawdata=f.readlines()
f.close()

print('trimmed version of',fnin,'(',fnin,'contains',lr-2,'wavelengths) has been saved as',fnout2)

import shutil
shutil.copy(fnout2, 'experimental.template.dat')
print(fnout2,'has been saved as','experimental.template.dat')
