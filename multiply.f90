PROGRAM multiply

      IMPLICIT NONE

Character(Len=20) :: Input,Output
Real(Kind=8) :: wl,re,im
Real(Kind=8) :: nwl,nre,nim

Input="xdata"
Output="xdata.dat"

Open(Unit=1,File=Input,Status="Old")
Read(Unit=1,Fmt=*) wl,re,im
nwl=wl/1
nre=re/1
nim=im/1
Open(Unit=2,File=Output,Status="Replace")
Write(Unit=2,Fmt=*) nwl,nre,nim

      END PROGRAM multiply


