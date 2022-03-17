C
C  multisphere-scattering Fortran code gmm01s.f 
C  the scattering formulations and numerical techniques used as well as
C  the usage are exactly the same as the gmm01f.f code 
C  the only difference between this code and gmm01f.f is that this code 
C  is less computer-memory-demanding 
C  This code can calculate aggregates of a larger number of spheres 
C  than the gmm01f.f code for the same available amount of computer 
C  memory. Consequently, this gmm01s.f code is slower in computing 
C  speed than the gmm01f.f code. 
C  Yu-lin Xu   
C  released to public  8/2001
C  available online at http://www.astro.ufl.edu/~xu
C  last updated  5/2003
C
C----------------------------------------------------------------------  
C
C  For questions/comments/suggestions/bugs/problems please contact 
C  Yu-lin Xu at xu@astro.ufl.edu or yxu@erc.ufl.edu
C
C-----------------------------------------------------------------------
      PROGRAM gmm01s
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
C----------------------------------------------------------------------
C  Two parameters in "gmm01f.par": nLp, np
C  nLp --- the maximum number of spheres, must be equal to or greater 
C          than the sphere-number in an aggregate actually calculated  
C  np  --- the maximum scattering order in the incident and scattered
C          field expansions, must be equal to or greater than the 
C          highest scattering order required in actual calculations;
C          usually, this can be estimated by Wiscombe's criterion  
C          regarding the largest component sphere in an aggregate
C  an example of "gmm01f.par":   parameter (nLp=10000,np=20)
C
C  
C  ******** Computer memory required by this code is at the level of 
C                           np^4 + nLp*np^2
C  which is significantly less than (nLp**2*np**3 + np**4) for the code
C  gmm01f.f, especially for a large number of component particles. This 
C  code is written in double precision arithmetic. An individual sphere 
C  may reach a size parameter way beyond ~200. There is no limit to the 
C  overall dimension of an aggregate. 
C  Note that there may be a convergence problem for aggregates of very  
C  large sizes or a very large number of component particles. ********* 
C
C----------------------------------------------------------------------
C  
C  Revision history:
C
C  September 2001: A bug in subroutine transs.f was found and fixed 
C                    by Zhifang Lin. 
C  October 2001:   One major and two minor bugs in the calculation of
C                    internal field distributions were found and fixed.
C                    These involve subroutine internd.f.
C  November 2001:  In solving boundary condition linear equations for 
C                    partial scattering coefficients, van der Vorst's 
C                    BI-CGSTAB scheme is replaced by its version of 
C                    Martin H. Gutknecht [SIAM J. Sci. comput. 14, 
C                    pp.1020-1033 (1993)], as suggested by Zhifang Lin.
C  May 2003:  A bug in the calculation of amplitude scattering matrix 
C             elements of s1y(theta,phi) in backward scattering 
C             directions (i.e., theta>90 degrees) was found and fixed. 
C             The format of output for the amplitude and Mueller 
C             matrices has been changed.
C
C----------------------------------------------------------------------
      parameter (nmp=np*(np+2),nmp0=(np+1)*(np+4)/2)
      parameter (NXMAX=3000,nangmax=181,MOR=181,ncmax=360)
C---------------------------------------------------------------------- 
C  NXMAX - the maximum dimension of an auxiliary array in calculating 
C          Mie scattering coefficients, must not be less than 
C                  1.1*np*|m|+10 
C          where |m| is the largest refractive index of all spheres) 
C----------------------------------------------------------------------
      parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
      parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)	
      parameter (nrc=4*np*(np+1)*(np+2)/3+np)
      integer u,v,u0,nmax(nLp),uvmax(nLp),ind(nLp),ind2(nLp)
      double precision k,lnfacd,r0(6,nLp),x(nLp),dang(nangmax),
     +   r00(3,nLp),rsr0(NXMAX),rsi0(NXMAX),
     +   rsx0(NXMAX),px0(NXMAX),w1(np),w2(np),w3(np),w4(np),
     +   rsr(np,nLp),rsi(np,nLp),rsx(np,nLp),px(np,nLp),betar(MOR),
     +   thetr(MOR),phair(MOR),smue(4,4),mue(4,4,ncmax,nangmax),
     +   besj(0:2*np+1),besy(0:2*np+1),i11(nangmax),
     +   i21(nangmax),i22(nangmax),i12(nangmax),inat(nangmax),
     +   pol(nangmax),cscaxi(nLp),cscayi(nLp),cextxi(nLp),
     +   cextyi(nLp),cabsxi(nLp),cabsyi(nLp),cexti(nLp),cabsi(nLp),
     +   cscai(nLp),assymi(nLp),assymxi(nLp),assymyi(nLp),
     +   cprxi(nLp),cpryi(nLp),cpri(nLp),drot(nrc),
     +   c0i(nLp),c1i(nLp)
      complex*16 A,B,cmz,Aj,Bj,A2,B2,Aj2,Bj2,A0,B0,ephi,ci,cin,
     +   atr0(ni0),btr0(ni0),atr(2,np,nmp),at(nmp),bt(nmp),
     +   ek(np),ref(nLp),ref0(nLp),p0(nLp,nmp),
     +   q0(nLp,nmp),an(np),bn(np),aMie(nLp,np),bMie(nLp,np),
     +	 as(nLp,nmp),bs(nLp,nmp),as0(nLp,nmp),bs0(nLp,nmp),
     +   asc(nLp,nmp),bsc(nLp,nmp),as1(nLp,nmp),bs1(nLp,nmp),
     +   ast(nLp,nmp),bst(nLp,nmp),asp(nLp,nmp),bsp(nLp,nmp),
     +   asv(nLp,nmp),bsv(nLp,nmp),B2i(nLp),
     +   s2x(ncmax,nangmax),s4x(ncmax,nangmax),
     +	 s3y(ncmax,nangmax),s1y(ncmax,nangmax),
     +	 atj(nmp),btj(nmp),py0(NXMAX),py(NXMAX),dpy(NXMAX)
      CHARACTER FLNAME*20,fileout*20,fileout1*19,fileout2*21,
     +   tailn*3,cnr2*2,cnr3*3,cnr1*1,flout*22
      COMMON/MIESUB/ twopi,pih
      common/rot/bcof(0:np+2),dc(-np:np,0:nmp)
      common/fnr/fnr(0:2*(np+2))
      common/pitau/pi(nmp0),tau(nmp0)
      common/tran/atr
      common/ig0/iga0(ni0)
      common/g0/ga0(ng0)
      common/cofmnv0/cof0(ni0)
      common/crot/cofsr(nmp)
      pih = dacos(0.d0)
      twopi = 4.d0*pih
      pione = 2.d0*pih
      ci=dcmplx(0.d0,1.d0)
      cin=dcmplx(0.d0,-1.d0)
      gcs=0.d0
      gcv=0.d0
      idpq=0	
      OPEN(UNIT=1,FILE='gmm01f.in',STATUS='OLD')
      READ(1,'(a20)') FLNAME
C-----------------------------------------------------------------------
C  FLNAME - the input file name for the aggregate configuration 
C           the first line in this file is the incident wavelength, 
C           the second line is the number of spheres, 
C           rest lines provide coordinates of sphere-centers, radii, and 
C           refractive indexes of the spheres; each line contains six 
C           real numbers of x, y, z, r, Re(m), Im(m) for a sphere
C-----------------------------------------------------------------------
      write(6,'(a12,a20)') 'input file: ',FLNAME
      READ(1,*) nbeta,nthet,nphai
C-----------------------------------------------------------------------
C  The product of nbeta, nthet, and nphai is the total number of 
C  orientations to be averaged for the input sphere-aggregate.
C  In the Cartesian coordinate system, the direction of propagation of 
C  the incident plane wave defines the positive z-axis, the x-z plane is 
C  the scattering plane. 
C  nbeta - # of rotation around the z-axis,  
C          this rotation angle corresponds to the Euler angle of "alpha" 
C  nthet - # of rotation around the y-axis, 
C          corresponding to the Euler angle of "beta"
C  nphai - # of rotation by the Euler angle of "gamma"
C  The definitions of the three Euler angles follow the convention used 
C  by Edmonds ["Angular Momentum in Quantum Mechanics," pp.6-8 (1957)].
C  examples:
C  When only a single orientation needs to be calculated:
C    nbeta=nthet=nphai=1 (the input line will be 1,1,1 or 1 1 1)
C  When an aggregate rotates around the y-axis (i.e., rotates in the  
C  scattering plane), the number of orientations of the aggregate to be 
C  calculated is 19, then nbeta=1, nthet=19, nphai=1 (1,19,1 or 1 19 1) 
C
C  *** In addition to the output file of mueller.out for the scattering
C      Mueller matrix, there is also an output file of amp.out for the
C      amplitude scattering matrix when nbeta=nthet=nphai=1 (i.e., when
C      a single fixed-orientation is calculated).
C-----------------------------------------------------------------------
      write(6,'(a19,3i5)') 'nbeta,nthet,nphai: ',nbeta,nthet,nphai
      if(nbeta.gt.MOR.or.nthet.gt.MOR.or.nphai.gt.MOR) then
         write(6,*) '***  parameter MOR too small  ***'
         write(6,*) ' MOR must be >nbeta,nthet,nphai given above'
         write(6,*) ' Please change MOR in the parameter line of the' 
         write(6,*) '   main code, recompile, then try again'
         stop
      endif 
      if(nbeta*nphai.gt.1) then 
         idran=1
      else
         idran=0
      endif
      nram=nbeta*nthet*nphai
      if(nram.lt.1) then 
         write(6,*) 'please check (nbeta,nthet,nphai) in gmm01f.in'
         stop
      endif
      READ(1,*) betami,betamx,thetmi,thetmx,phaimi,phaimx
C-----------------------------------------------------------------------
C  betami,betamx,thetmi,thetmx,phaimi,phaimx are in degrees 
C  betami,betamx -- the range of the Euler angle of rotation (alpha)  
C                   to be calculated
C  thetmi,thetmx -- the range of the Euler angle of rotation (beta)  
C                   to be calculated
C  phaimi,phaimx -- the range of the Euler angle of rotation (gamma)  
C                   to be calculated
C  an example:
C  When only a fixed orientation needs to be calculated, betami=betamx,
C  thetmi=thetmx, and phaimi=phaimx, this line could be, e.g.,  
C  (0. 0. 0. 0. 0. 0.), (0. 0. 90. 90. 0. 0.), (30. 30. 40. 40. 45. 45.) 
C------------------------------------------------------------------------
      if(nbeta.eq.1) betamx=betami
      if(nthet.eq.1) thetmx=thetmi
      if(nphai.eq.1) phaimx=phaimi
      write(6,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +         betami,betamx,thetmi,thetmx,phaimi,phaimx
      READ(1,*) idMie
C------------------------------------------------------------------------
C  idMie=1: calculating only coherent Mie-scattering, no interaction 
C------------------------------------------------------------------------
      write(6,'(a7,i3)') 'idMie: ',idMie
      READ(1,*) idd
C------------------------------------------------------------------------
C  When calculating a single orientation, put idd=0
C  When idd=1, the number of orientations to be calculated is doubled. 
C  Each orientation is coupled with the orientation that the aggregate  
C  is rotated by 90 degrees around the z axis. This insures that the 
C  averaged polarizations are zero when the scattering angle is either 0 
C  or 180 degrees, as is suppossed to be for an average over random 
C  orientations.
C------------------------------------------------------------------------
      if(nram.eq.1) idd=0
      READ(1,*) idc,iseed
C------------------------------------------------------------------------
C  idc for choosing in the three schemes for an orientation average 
C  idc=1 to devide the range of [thetmi,thetmx] by the cos(theta) scheme
C        (theta is the scattering angle, corresponding to the Euler angle 
C        of ritation "beta") 
C  idc=0 to devide the range of [thetmi,thetmx] by "degrees"
C  idc=-1 to pick up all the three Euler angles of rotation randomly using 
C         a random number generator 
C  iseed is the seed number for the random number generator, used only 
C        when idc=-1 (i.e., when idc=1 or 0, iseed has no function)
C        iseed can be an arbitrary positive integer
C  when calculating only one single orientation, i.e., 
C       nram=nbeta*nthet*nphai=1, idc and iseed have no function
C------------------------------------------------------------------------
      if(idpq.eq.1) then
         nram=nphai+nthet
         idc=0
         idd=0
         nbeta=1
         betami=0.d0
         betamx=betami
         thetr(nthet+1)=0.d0
         phair(nphai+1)=0.d0
      endif
      write(6,'(a5,i3)') 'idd: ',idd
      if(idd.eq.1) nram=2*nram
      write(6,'(a36,i5)') '# of orientations to be averaged: ',nram	
      if(idc.lt.0) then 
         write(6,'(a11,i4,i12)') 'idc,iseed: ',idc,iseed
      else
         write(6,'(a5,i3)') 'idc: ',idc
      endif	
      READ(1,*) factor1,factor2,MXINT
C------------------------------------------------------------------------
C  factor1 and factor2 are numerical factors used for improving the 
C  convergence behavior of the iterative solution process in solving 
C  interacting equations, which are in the range of [0,1]
C  factor1 is for x-polarized incident plane wave and factor2 for
C  y-polarized incident wave 
C  MXINT is the maximum number of iterations allowed in the iterative 
C  solution process
C  This code uses two alternative methods in the iterative solution of 
C  interacting equations: iteration scheme [see Fuller and Kattawar, Opt. 
C  Lett. 13, 90 (1988); Xu, Appl. Opt. 34, 4573 (1995)] and BI-CGSTAB,  
C  the stabilized Bi-Conjugate Gradient [see H.A. van der Vorst, SIAM J. 
C  Sci. Stat. Comput. 13, 631, (1992); M.H. Gutknecht, SIAM J. Sci. 
C  comput. 14, pp.1020-1033 (1993)].  
C  When factor1=0 or factor2=0, the code directly goes to BI-CGSTAB 
C  without using the other iteration scheme.  
C  When factor1=factor2=1 it is equivalent to  Fuller and Kattawar's 
C  order-of-scattering method. 
C  In many cases, a divergence will occur when factor1=factor2=1. Then,  
C  setting factor1,factor2<1 may help to converge to a numerical 
C  solution. When MXINT is exceeded, it will automatically switch to 
C  BI-CGSTAB.
C  
C  ******** In most circumstances the BI-CGSTAB is much more efficient   
C  than the other iterative schemes so that factor1=factor2=0 is the best 
C  choice for the iterative solution of the boundary condition equations.
C  When a user would like to try other iterative schemes, just simply 
C  comment out the following two lines factor1=0 and factor2=0 and 
C  the user-specified input values of factor1 and factor2 will be used .
C  10/20/02
C-----------------------------------------------------------------------
      factor1=0
      factor2=0
      write(6,'(a,2f6.2)') 'Numerical factors for convergence:',
     +                        factor1,factor2
      write(6,'(a37,i5)') 'Maximum iterations allowed:',MXINT     
      READ(1,*) NADD
C-----------------------------------------------------------------------
C  NADD is the number of terms to add to the scattering orders required  
C  by the Wiscombe's criterion, which can be negative or positive in  
C  the range of [-9,99] 
C  Normally, set NADD=0
C-----------------------------------------------------------------------
      write(6,'(a41,i3)') 
     +   'Scat. orders added to Wiscombe criterion:',NADD
      READ(1,*) eps,small
      if(eps.gt.1.d-20) eps=1.d-20
C-----------------------------------------------------------------------
C  eps: error tolerance for determining single-sphere field-expansion 
C       truncation, default: 1.d-20 
C       (the default value of 1.d-20 allows to use Wiscombi's criterion) 
C  small: error tolerance for the iterative solution process for solving 
C         the interacting scattering coefficients (1.d-6)
C------------------------------------------------------------------------
      write(6,'(a35,e10.2)') 'error tolerance for Mie-expansions:',eps
      write(6,'(a22,e10.2)') 'Convergence criterion:',small
      READ(1,*) fint
C------------------------------------------------------------------------
C  fint is the interaction index in the range of [0,1] (default: 0.02)
C  In the scattering calculations, a quantity "f" for a pair of component 
C  spheres is defined by f=(r_i+r_j)/d_{ij}, where (r_i,r_j) are the 
C  radii of the two spheres and d_{ij} is the separation distance of the 
C  two sphere-centers. When f<fint, interaction between the pair is 
C  considered to be negligible and no interaction will be calculated 
C  between the two spheres. When fint=0, no sphere is excluded in 
C  interaction calculations.
C------------------------------------------------------------------------
      if(fint.lt.0.d0.or.fint.gt.1.d0) then
         fint=0.02d0
         write(6,'(a37)') 'Interaction index: using default 0.02'
      else
         write(6,'(a18,f7.3)') 'Interaction index:',fint
      endif
      READ(1,*) sang,pang
C------------------------------------------------------------------------
C  sang -- the scattering angle interval for output  
C  example: when sang=1, the results in output will be written for every 
C  degree of scattering angle from 0 to 180 
C  pang -- the azimuth angle interval for the two-dimensional Mueller 
C          matrix output
C      (1) For each azimuth angle, the number of scattering angles that 
C          will be calculated is the same as that determined by "sang", 
C          i.e., the number of scattering angles that will calculated is 
C          the same as calculated for the scattering plane of the azimuth 
C          angle=0 (180/sang +1).
C      (2) When pang = 0., no additional calculations for the scattering 
C          matrix map, i.e., calculating only the scattering plane of the 
C          azimuth angle = 0. 
C      (3) When pang>0, the number of azimuth angles in the range of 
C          (0,360) degrees that will be calculated in addition to 0 (360)  
C          degrees is npng-1 with npng=360/pang. For example, when 
C          pang=180, npng=2, the azimuth angles 0 (360) and 180 degrees 
C          will be calculated. In the output for the Mueller matrix at 
C          each azimuth angle, there are nang2 (=2*nang-1) sets of the 
C          16 elements. 
C------------------------------------------------------------------------
      write(6,'(a32,f7.3)') 'scat.-angle-interval in output: ',sang
      if(sang.le.0.d0) sang=1.d0
      nang=90.d0/sang+1.d0
      nang2=2*nang-1
      if(nang2.gt.nangmax) then
         write(6,*) 'sang too small'
         write(6,*) 'please increase sang in the input file gmm01f.in'
         write(6,*) '  and try again, or'
         write(6,*) '  increase nangmax in the parameter line of the'
         write(6,*) '  main code, recompile, then try again'
         stop
      endif
      write(6,'(a,f9.3)') 
     +   'azimuth-angle-interval in Mueller matrix output: ',pang
      if(pang.lt.0.0001d0) then
         npng=1
      else
         npng=360.0d0/pang
      endif
      if(npng.gt.ncmax) then
         write(6,*) 'pang too small'
         write(6,*) 
     +      'please increase pang in the input file gmm01f.in'
         write(6,*) 
     +      'and try again, or increase ncmax in the parameter'
         write(6,*) 
     +      ' line of the main code, recompile, then try again'
         stop
      endif
      READ(1,*) idphoto,nphoto,dphi,istart,iend,istep
C------------------------------------------------------------------------
C  idphoto=0 or nphoto=0 for not calculating internal field distributions
C  idphoto=1 for calculating internal field distributions in spherical 
C                coordinates r/a [0,1] and theta [0,360] for each cross 
C                section at angle phi
C  idphoto=2 for calculating internal field distributions in Cartesian 
C                coordinates in x'/a [-1,1] and z/a [-1,1] for each cross 
C                section at angle phi
C  nphoto+1: number of mesh points in (r/a,theta) or in (x'/a,z/a)
C            nphoto must be an even number
C  dphi:   dphi=0 for calculating phi=0 only
C          otherwise, the calculated number of phi's would be 180/dphi
C          (dphi >= 0.1)
C  The internal field distributions can be calculated for a selected 
C  component sphere or a selected range of spheres by specifying 
C  (istart,iend,istep). For example, only the internal field distribution 
C  of the 5th sphere will be calculated when (istart,iend,istep)=(5,5,1).
C  When (istart,iend,istep)=(1,5,2), the internal field distributions of  
C  the 1st, 3rd, and 5th spheres will be calculated. When istart=1,
C  iend=nL, and istep=1, the internal distributions will be calculated 
C  for all component spheres. 
C  Note that this code calculates the internal distributions for only a 
C  fixed orientation and does not average over orientations. Thus, it 
C  calculates the internal fields for only the first orientation.  
C  Note also that the output files for the internal field distributions  
C  have a large size, especially when nphoto is large.
C------------------------------------------------------------------------
      if(idphoto.gt.0) then 
         write(6,'(a16,i3,i6,f5.1,3i6)') 'idphoto,nphoto,dphi: ',
     +                     idphoto,nphoto,dphi,istart,iend,istep
      else 
         write(6,'(a9,i3)') 'idphoto: ',idphoto
      endif
      close(1)
      write(6,'(/)')
	
      OPEN(UNIT=2,FILE=FLNAME,STATUS='OLD')
      READ(2,*) w
C------------------------------------------------------------------------
C  w -- incident wavelength
C------------------------------------------------------------------------
      READ(2,*) nL
C------------------------------------------------------------------------
C  nL -- number of spheres in the aggregate
C------------------------------------------------------------------------
      if(nL.gt.nLp) then
         write(6,*) 'Parameter nLp too small, must be >', nL 
         write(6,*) 
     +      'Change nLp in gmm01f.par, recompile, then try again'
         stop
      endif
      if(nL.eq.1) then 
         idMie=1
         betami=0.d0
         betamx=0.d0
         thetmi=0.d0
         thetmx=0.d0
         phaimi=0.d0
         phaimx=0.d0
         nbeta=1
         nthet=1
         nphai=1
         nram=1
         idc=0
         idd=0
      endif
C------------------------------------------------------------------------
C  input the configuration and particle parameters for the aggregate
C  each line includes 6 numbers:
C  x-, y-, z-coordinates of the sphere-center, the radius of the sphere 
C  in the same unit of the incident wavelength, the real and imaginary 
C  parts of the refractive index 
C------------------------------------------------------------------------
      do 1 i=1,nL
         read(2,*,err=10) (r0(j,i),j=1,6)
         x0=r0(1,i)
         y0=r0(2,i)
         z0=r0(3,i)
         r00(1,i)=x0
         r00(2,i)=y0
         r00(3,i)=z0
         if(r0(6,i).gt.0.d0) r0(6,i)=-r0(6,i)
         if(r0(5,i).eq.1.d0.and.r0(6,i).eq.0.d0) goto 1 
         gcs=gcs+r0(4,i)*r0(4,i)
         gcv=gcv+r0(4,i)*r0(4,i)*r0(4,i)
 1    continue
      close(2)
      gcsr=dsqrt(gcs)
      gcvr=gcv**(1.d0/3.d0)
      goto 11
 10   write(6,*) 'fatal error in the input file'
      stop
 11   k=twopi/w
      xv=k*gcvr
      xs=k*gcsr
      write(6,'(a,f7.3,a,f7.3)') ' volume-equiv. xv: ',xv,
     +      '  surface-equiv. xs: ',xs
      write(6,'(/)')
        
      betami=betami*pih/90.d0
      betamx=betamx*pih/90.d0
      thetmi=thetmi*pih/90.d0
      thetmx=thetmx*pih/90.d0
      phaimi=phaimi*pih/90.d0
      phaimx=phaimx*pih/90.d0
      if(idc.gt.0) then
         call orientcd(betami,betamx,thetmi,thetmx,phaimi,phaimx,
     +           MOR,MOR,MOR,nbeta,nthet,nphai,betar,thetr,phair)
      else
         call orientud(betami,betamx,thetmi,thetmx,phaimi,phaimx,
     +           MOR,MOR,MOR,nbeta,nthet,nphai,betar,thetr,phair)
      endif
      fileout='gmm01s.out'
      fileout1=fileout
      fileout2='pq'//fileout1
      if(idMie.eq.1) then
         write(6,*) 
     +      '*** Calculating only coherent Mie-scattering ***'
         write(6,*) 
     +      '*** No interaction included ********************'
      endif

      do i=1,nL
         x(i)=k*r0(4,i)
         ref(i)=dcmplx(r0(5,i),r0(6,i))
         ref0(i)=dcmplx(r0(5,i),-r0(6,i))
      enddo
      do j=1,np
         do i=1,nL
            aMie(i,j)=0.d0
            bMie(i,j)=0.d0
         enddo
      enddo
      nmax0=1
      do i=1,nL
         if(i.eq.1) goto  12
         if(x(i).eq.x(i-1).and.ref(i).eq.ref(i-1)) then
            nmax(i)=nmax(i-1)
            uvmax(i)=uvmax(i-1)
            goto 15
         endif
 12      write(6,'(a,i3,a,f10.4)') 'sphere #',i, 
     +         '   individual size parameter: ',x(i)
c
c  calculating Mie-scattering coefficients for each spheres
c  the ratio method of Wang and van der Hulst is used in calculating  
c  Riccati-Bessel functions [see Wang and van der Hulst, Appl. Opt. 
c  30, 106 (1991), Xu, Gustafson, Giovane, Blum, and Tehranian, 
c  Phys. Rev. E 60, 2347 (1999)]
c 
         call abMiexud(x(i),ref(i),np,NXMAX,nmax(i),an,bn,NADD,
     +                 rsr0,rsi0,rsx0,px0,w1,w2,w3,w4,eps)
         if(nmax(i).gt.np) then
            write(6,*) ' Parameter np too small, must be >',nmax(i)
            write(6,*) ' Please change np in gmm01f.par, recompile,' 
            write(6,*) '   then try again'
            stop
         endif
         uvmax(i)=nmax(i)*(nmax(i)+2)
         write(6,'(a,1x,i4)') 
     +      ' Actual single-sphere expansion truncation:',nmax(i)
         do j=1,nmax(i)
            rsr(j,i)=rsr0(j)
            rsi(j,i)=rsi0(j)
            rsx(j,i)=rsx0(j)
            px(j,i)=px0(j)
            temp1=an(j)
            temp2=bn(j)
            if(j.eq.1.or.j.eq.nmax(i)) write(6,'(i10,4e15.7)')
     +         j,temp1,dimag(an(j)),temp2,dimag(bn(j))
         enddo
 15      do j=1,nmax(i)
            aMie(i,j)=an(j)
            bMie(i,j)=bn(j)
            rsr(j,i)=rsr0(j)
            rsi(j,i)=rsi0(j)
            rsx(j,i)=rsx0(j)
            px(j,i)=px0(j)
         enddo
         if(nmax(i).gt.nmax0) nmax0=nmax(i)
      enddo
      cextx=0.d0
      cexty=0.d0
      cabsx=0.d0
      cabsy=0.d0
      cscax=0.d0
      cscay=0.d0
      cprx=0.d0
      cpry=0.d0
      cbakx=0.d0
      cbaky=0.d0
      do i=1,nL
         cextxi(i)=0.d0
         cextyi(i)=0.d0
         cabsxi(i)=0.d0
         cabsyi(i)=0.d0
         cscaxi(i)=0.d0
         cscayi(i)=0.d0
         cprxi(i)=0.d0
         cpryi(i)=0.d0
      enddo
      do i=1,nang2
         i11(i)=0.d0
         i21(i)=0.d0
         i22(i)=0.d0
         i12(i)=0.d0
         do jc=1,npng
            do j=1,4
               do m=1,4
                  mue(j,m,jc,i)=0.d0
               enddo
            enddo
         enddo
      enddo
      iram=0
      if(idpq.eq.1) then 
         open(12,file=fileout2,status='unknown')
         write(12,*) 'input file: ',FLNAME 
      endif
      write(6,'(/)')
      write(6,*) 'original input sphere-positions: '
      i=1
      write(6,'(i5,3f14.5)') i,r0(1,1),r0(2,1),r0(3,i)
      i=nL
      write(6,'(i5,3f14.5)') i,r0(1,i),r0(2,i),r0(3,i)
      if(idpq.eq.1) then
         nphaic=nphai+1
         phair(nphaic)=0.d0
      else
         nphaic=nphai
      endif

C----------------------------------------------
C  calculating constants and Gaunt coefficients
C----------------------------------------------
      n0=nmax0+2
      fnr(0)=0.d0
      do n=1,2*n0
         fnr(n)=dsqrt(dble(n))
      enddo
      bcof(0)=1.d0
      do n=0,n0-1
         bcof(n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n)/fnr(n+1)/fnr(n+1)
      enddo
C
C  the formulation used here for the calculation of Gaunt coefficients 
C  can be found in Bruning and Lo, IEEE Trans. Anten. Prop. Ap-19, 378 
C  (1971) and Xu, J. Comput. Appl. Math. 85, 53 (1997), J. Comput. Phys. 
C  139, 137 (1998)
C
      call cofsrd(nmax0)	
      call cofd0(nmax0)
      call cofnv0(nmax0)
      call gau0(nmax0)
	
      do ibeta=1,nbeta
         do iphai=1,nphaic
            if(idpq.eq.1.and.iphai.lt.nphaic) then
               nthetc=1
            else
               nthetc=nthet
            endif
            do ithet=1,nthetc
               if(idc.lt.0) then
                  betar(ibeta)=(betamx-betami)*ran1d(iseed)
                  phair(iphai)=(phaimx-phaimi)*ran1d(iseed)
                  thetr(ithet)=(thetmx-thetmi)*ran1d(iseed)
               endif 
               do 19 irot=1,2
                  if(idd.ne.1.and.irot.eq.2) goto 19
                  iram=iram+1
                  if(irot.eq.1) then
                     alph=0.d0
                  else
                     alph=pih
                  endif
                  ca=dcos(alph)
                  sa=dsin(alph)
                  beta=betar(ibeta)
                  cb=dcos(beta)
                  sb=dsin(beta)	   
                  do i=1,nL
                     x0=r00(1,i)
                     y0=r00(2,i)
                     r0(1,i)=cb*x0-sb*y0
                     r0(2,i)=sb*x0+cb*y0
                  enddo
                  phai=phair(iphai)
                  thet=thetr(ithet)
                  if(idpq.eq.1.and.nthetc.eq.1) thet=0.d0
                  cb=dcos(phai)
                  sb=dsin(phai)
                  cz=dcos(thet)
                  sz=dsin(thet)
                  if(iram.eq.1.or.iram/50*50.eq.iram) then 
                     write(6,'(a,2i5)') 'iram & nram: ', 
     +                                   iram,nram
                  endif     
                  do i=1,nL
                     x0=r0(1,i)
                     y0=r0(2,i)
                     z0=r00(3,i)
                     r0(1,i)=ca*cz*x0-(ca*sz*sb+sa*cb)*y0
     +                               +(ca*sz*cb-sa*sb)*z0
                     r0(2,i)=sa*cz*x0-(sa*sz*sb-ca*cb)*y0
     +                               +(sa*sz*cb+ca*sb)*z0
                     r0(3,i)=-sz*x0-cz*sb*y0+cz*cb*z0	
                  enddo
                  if(iram.eq.1.or.iram/50*50.eq.iram) then 
                     i=1	         
                     write(6,'(i5,3f14.5)') 
     +                     i,r0(1,i),r0(2,i),r0(3,i)
                     i=nL
                     write(6,'(i5,3f14.5)') 
     +                     i,r0(1,i),r0(2,i),r0(3,i)
                  endif

                  indpol=0
                  factor=factor1
                  write(6,'(a8,i4,a32)') ' orien.#',iram,
     +               '  Solving for x-pol. inci. state'
 18               do imn=1,nmp
                     do i=1,nL
                        p0(i,imn)=0.d0
                        q0(i,imn)=0.d0
                        as(i,imn)=0.d0
                        bs(i,imn)=0.d0
                     enddo
                  enddo
C--------------------------------------------------
C  calculating incident wave expansion coefficients
C--------------------------------------------------
                  do 20 i=1,nL
                     cz=dcos(k*r0(3,i))
                     sz=dsin(k*r0(3,i))
                     cmz=0.5d0*dcmplx(cz,sz)
                     do 21 n=1,nmax(i)
                        imn=n*n+n+1
                        A=fnr(2*n+1)*cmz
                        p0(i,imn)=aMie(i,n)*A
                        q0(i,imn)=bMie(i,n)*A
                        p0(i,imn-2)=-p0(i,imn)
                        q0(i,imn-2)=q0(i,imn)
                        if(indpol.gt.1) then
                           p0(i,imn)=p0(i,imn)*cin
                           q0(i,imn)=q0(i,imn)*cin
                           p0(i,imn-2)=p0(i,imn)
                           q0(i,imn-2)=-q0(i,imn)
                        endif
                        as(i,imn)=p0(i,imn)
                        bs(i,imn)=q0(i,imn)
                        as(i,imn-2)=p0(i,imn-2)
                        bs(i,imn-2)=q0(i,imn-2)
 21                  continue
 20               continue
C--------------------------------------------------------------
c  begins iteration process to solve the interaction equations
c  for partial interacting scatteing coefficients
c  factor1,factor2=0: BI-CGSTAB [see H.A. van der Vorst, SIAM,  
c  J. Sci. Stat. Comput. 13, 631 (1992)]
c  factor1=factor2=1: order-of-scattering [Fuller and Kattawar, 
c  Opt. Lett. 13, 90 & 1063 (1988)]
c  0<factor1,factor2<1: an iteration scheme [Xu, Appl. Opt. 34,
c  4573 (1995)] 
c  When the number of iterations exceeds the maximum MXINT, it 
c  switches to use BI-CGSTAB.   
C---------------------------------------------------------------
                  if(idMie.eq.1) goto 490
                  if(nL.eq.1) goto 490
                  do i=1,nL
                     ind(i)=0
                     ind2(i)=0
                     c0i(i)=0.d0
                     do n=1,nmax(i)
                        imn=n*n+n+1
                        c0i(i)=c0i(i)
     +                         +p0(i,imn)*dconjg(p0(i,imn)) 
                        c0i(i)=c0i(i)
     +                         +q0(i,imn)*dconjg(q0(i,imn))
                        c0i(i)=c0i(i)
     +                     +p0(i,imn-2)*dconjg(p0(i,imn-2)) 
                        c0i(i)=c0i(i)
     +                     +q0(i,imn-2)*dconjg(q0(i,imn-2))
                     enddo
                  enddo
                  icase=1
                  niter=1
                  if(factor1.lt.0.001d0.or.factor2.lt.0.001d0) 
     +               goto 61
                  if(iram.eq.1) write(6,*) 
     +               'Starting iteration solution process'	
                  do i=1,nL
                     do imn=1,uvmax(i)
                        as0(i,imn)=p0(i,imn)
                        bs0(i,imn)=q0(i,imn)
                     enddo
                  enddo
 60               call transs(k,nL,r0,nmax0,nmax,uvmax,fint,ek,
     +                besj,besy,drot,as0,bs0,as1,bs1,ind,icase)
                  do 601 i=1,nL
                     if(ind(i).gt.0) goto 601
                     if(ind2(i).gt.0) goto 601
                     c1i(i)=0.d0
                     do imn=1,uvmax(i)
                        n=dsqrt(dble(imn))
                        as0(i,imn)=p0(i,imn)-aMie(i,n)*as1(i,imn)
                        bs0(i,imn)=q0(i,imn)-bMie(i,n)*bs1(i,imn)
                        A=as0(i,imn)-as(i,imn)
                        B=bs0(i,imn)-bs(i,imn)
                        c1i(i)=c1i(i)+A*dconjg(A)
                        c1i(i)=c1i(i)+B*dconjg(B)
                        as0(i,imn)=as(i,imn)+factor*A
                        bs0(i,imn)=bs(i,imn)+factor*B
                        as(i,imn)=as0(i,imn)
                        bs(i,imn)=bs0(i,imn) 
                     enddo
 601              continue
                  cext0=0.d0
                  cext1=0.d0
                  do 602 i=1,nL
                     if(ind(i).gt.0) goto 602
                     if(ind2(i).gt.0) goto 602
                     cext0=cext0+c0i(i)
                     cext1=cext1+c1i(i)
                     temp=c1i(i)/c0i(i)
                     if(temp.lt.small) ind2(i)=1
 602              continue
                  temp=cext1/cext0
                  if(temp.lt.small) goto 490
                  if(iram.eq.1.or.iram.eq.nram) then
                     if(niter.eq.1.or.niter/20*20.eq.niter) then 
                        write(6,'(a11,i4,2x,e15.7)') 
     +                        'iteration #',niter,temp
                     endif
                  endif
                  if(niter.gt.MXINT) then
                     write(6,*) 
     +               '*** Maximum iterations exceeded  ***'
                     write(6,*) 
     +               '*** Switched to Bi-CGSTAB method ***'
                     do i=1,nL
                        ind(i)=0
                        ind2(i)=0
                        do imn=1,uvmax(i)
                           as(i,imn)=p0(i,imn)
                           bs(i,imn)=q0(i,imn)
                        enddo
                     enddo
                     niter=1
                     goto 61
                  endif
                  niter=niter+1
                  goto 60
 61               if(iram.eq.1) write(6,*) 
     +               'Starting Bi-CGSTAB solution process'
                  call transs(k,nL,r0,nmax0,nmax,uvmax,fint,
     +               ek,besj,besy,drot,as,bs,as1,bs1,ind,icase)
                  do 611 i=1,nL
                     c1i(i)=0.d0
                     if(ind(i).gt.0) goto 611
                     do imn=1,uvmax(i)
                        n=dsqrt(dble(imn))
                        as1(i,imn)=-aMie(i,n)*as1(i,imn)
                        bs1(i,imn)=-bMie(i,n)*bs1(i,imn)
                        c1i(i)=c1i(i)
     +                         +as1(i,imn)*dconjg(as1(i,imn))
                        c1i(i)=c1i(i)
     +                         +bs1(i,imn)*dconjg(bs1(i,imn))
                     enddo
 611              continue
                  temp=0.d0
                  B0=0.d0
                  do 612 i=1,nL
                     cext0=c1i(i)/c0i(i)
                     if(cext0.lt.small) ind(i)=1
                     if(ind(i).gt.0) goto 612
                     if(cext0.gt.temp) temp=cext0
                     B0=B0+c1i(i)
 612              continue
                  if(temp.lt.small) goto 490
                  do 613 i=1,nL
                     if(ind(i).gt.0) goto 613
                     do imn=1,uvmax(i)
                        asp(i,imn)=as1(i,imn)
                        bsp(i,imn)=bs1(i,imn)
                        as0(i,imn)=as1(i,imn)
                        bs0(i,imn)=bs1(i,imn)
                     enddo
 613              continue
                  call transs(k,nL,r0,nmax0,nmax,uvmax,fint,ek,
     +               besj,besy,drot,asp,bsp,ast,bst,ind,icase)
                  A0=0.d0
                  do 614 i=1,nL
                     if(ind(i).gt.0) goto 614
                     do imn=1,uvmax(i)
                        n=dsqrt(dble(imn))
                        ast(i,imn)=aMie(i,n)*ast(i,imn)
     +                             +asp(i,imn)
                        bst(i,imn)=bMie(i,n)*bst(i,imn)
     +                             +bsp(i,imn)
                        A0=A0+dconjg(as0(i,imn))*ast(i,imn)
                        A0=A0+dconjg(bs0(i,imn))*bst(i,imn)
                     enddo
 614              continue
                  if(cdabs(A0).lt.1.d-100) goto 490
                  Aj=B0/A0	
  62              do 621 i=1,nL
                     if(ind(i).gt.0) goto 621
                     do imn=1,uvmax(i)
                        asv(i,imn)=asp(i,imn)-Aj*ast(i,imn)
                        bsv(i,imn)=bsp(i,imn)-Aj*bst(i,imn)
                     enddo
 621              continue
                  call transs(k,nL,r0,nmax0,nmax,uvmax,fint,ek,
     +               besj,besy,drot,asv,bsv,asc,bsc,ind,icase)
 	
                  A2=0.d0
                  B2=0.d0
                  do 622 i=1,nL
                     if(ind(i).gt.0) goto 622
                     do imn=1,uvmax(i)
                        n=dsqrt(dble(imn))
                        asc(i,imn)=aMie(i,n)*asc(i,imn)
     +                             +asv(i,imn)
                        bsc(i,imn)=bMie(i,n)*bsc(i,imn)
     +                             +bsv(i,imn)
                        A2=A2+dconjg(asc(i,imn))*asv(i,imn)
                        A2=A2+dconjg(bsc(i,imn))*bsv(i,imn)
                        B2=B2+dconjg(asc(i,imn))*asc(i,imn)
                        B2=B2+dconjg(bsc(i,imn))*bsc(i,imn)
                     enddo
 622              continue
                  if(cdabs(B2).lt.1.d-100) goto 490
                  Bj=A2/B2
                  do 623 i=1,nL
                     if(ind(i).gt.0) goto 623
                     do imn=1,uvmax(i)
                        asp(i,imn)=asv(i,imn)-Bj*asc(i,imn)
                        bsp(i,imn)=bsv(i,imn)-Bj*bsc(i,imn)
                     enddo
 623              continue
                  do 624 i=1,nL
                     c1i(i)=0.d0
                     if(ind(i).gt.0) goto 624
                     do imn=1,uvmax(i)
                        Aj2=Aj*as1(i,imn)+Bj*asv(i,imn)
                        Bj2=Aj*bs1(i,imn)+Bj*bsv(i,imn)
                        as(i,imn)=as(i,imn)+Aj2
                        bs(i,imn)=bs(i,imn)+Bj2
                        c1i(i)=c1i(i)+Aj2*dconjg(Aj2)
                        c1i(i)=c1i(i)+Bj2*dconjg(Bj2)
                     enddo
 624              continue
                  cext0=0.d0
                  cext1=0.d0
                  do 625 i=1,nL
                     if(ind(i).gt.0) goto 625
                     cext0=cext0+c0i(i)
                     cext1=cext1+c1i(i)
 625              continue
                  temp=cext1/cext0
                  if(temp.lt.small) goto 490
                  if(niter.gt.MXINT) then
                     write(6,*) 'Caution:'
                     write(6,*) 
     +               '*** Maximum iterations exceeded ***'
                     write(6,*) 
     +               '***  Results may be inaccurate  ***'	
                     goto 490
                  endif
                  if(iram.eq.1.or.iram.eq.nram) then
                     if(niter.eq.1.or.niter/20*20.eq.niter) then 
                        write(6,'(a11,i4,2x,e15.7)') 
     +                        'iteration #',niter,temp
                     endif
                  endif
                  B2=0.d0
                  do 626 i=1,nL 	   
                     if(ind(i).gt.0) goto 626
                     B2i(i)=0.d0
                     do imn=1,uvmax(i)
                        B2i(i)=B2i(i)
     +                         +dconjg(as0(i,imn))*asp(i,imn)
                        B2i(i)=B2i(i)
     +                         +dconjg(bs0(i,imn))*bsp(i,imn)
                     enddo
                     B2=B2+B2i(i)
 626              continue
                  A0=B0*Bj
                  if(cdabs(A0).lt.1.d-100) goto 490
                  A0=-Aj*B2/A0
                  do 627 i=1,nL
                     if(ind(i).gt.0) goto 627	   
                     do imn=1,uvmax(i)
                        as1(i,imn)=as1(i,imn)-Bj*ast(i,imn)
                        bs1(i,imn)=bs1(i,imn)-Bj*bst(i,imn)
                        as1(i,imn)=asp(i,imn)-A0*as1(i,imn)
                        bs1(i,imn)=bsp(i,imn)-A0*bs1(i,imn)
                     enddo
 627              continue
                  B0=0.d0
                  do 628 i=1,nL 	   
                     if(ind(i).gt.0) goto 628
                     cext0=c1i(i)/c0i(i)
                     if(cext0.lt.small) then
                        ind(i)=1
                        goto 628
                     endif
                     B0=B0+B2i(i)
 628              continue
                  call transs(k,nL,r0,nmax0,nmax,uvmax,fint,ek,
     +               besj,besy,drot,as1,bs1,ast,bst,ind,icase)
                  A0=0.d0
                  do 629 i=1,nL
                     if(ind(i).gt.0) goto 629
                     do imn=1,uvmax(i)
                        n=dsqrt(dble(imn))
                        ast(i,imn)=aMie(i,n)*ast(i,imn)
     +                             +as1(i,imn)
                        bst(i,imn)=bMie(i,n)*bst(i,imn)
     +                             +bs1(i,imn)
                        A0=A0+dconjg(as0(i,imn))*ast(i,imn)
                        A0=A0+dconjg(bs0(i,imn))*bst(i,imn)
                     enddo
 629              continue
                  if(cdabs(A0).lt.1.d-100) goto 490
                  Aj=B0/A0
                  niter=niter+1
                  goto 62
C----------------------------------------------------------------
c  computing total and differential scattering cross sections and 
c  efficiencies for radiation pressure, using the formulas given  
c  by Xu, Physics Letters A 249, 30 (1998)
C----------------------------------------------------------------
 490              do i=1,nL
                     ind(i)=0
                  enddo

c                  open(17,file='scacof.dat',status='unknown')
c                  do 282 i=1,nL
c                     write(17,'(a2,i3)') 'i=',i
c                     do 281 j=1,uvmax(i)
c                        n=dsqrt(dble(j))
c                        m=j-n*n-n
c                        write(17,'(2i5,4e17.8)') 
c     +                     m,n,dble(as(i,j)),dimag(as(i,j)),
c     +                     dble(bs(i,j)),dimag(bs(i,j))
c 281                 continue
c 282              continue
c                  close(17)

                  icase=2
                  call transs(k,nL,r0,nmax0,nmax,uvmax,fint,ek,
     +               besj,besy,drot,as,bs,as1,bs1,ind,icase)
                  do i=1,nL
                     do imn=1,uvmax(i)
                        at(imn)=as(i,imn)+as1(i,imn)
                        bt(imn)=bs(i,imn)+bs1(i,imn)
                     enddo
                     do n=1,nmax(i)
                        n1=n+1
                        n2=2*n
                        rn=1.0d0/dble(n*n1)
                        p=fnr(n)*fnr(n+2)/fnr(n2+1)
     +                    /fnr(n2+3)/dble(n1)
                        t=fnr(n-1)*fnr(n+1)/fnr(n2-1)
     +                    /fnr(n2+1)/dble(n)
                        sc=0.d0
                        temp=0.d0
                        do m=-n,n
                           iL=n*(n+1)+m
                           sc=sc+dble(dconjg(as(i,iL))*at(iL))
                           sc=sc+dble(dconjg(bs(i,iL))*bt(iL))
                           rm=dble(m)*rn
                           A0=rm*bt(iL)
                           B0=rm*at(iL)              
                           if(n.eq.nmax(i)) goto 51
                           u=(n+1)*(n+2)+m
                           fnp=fnr(n+m+1)*fnr(n-m+1)*p
                           A0=A0+fnp*at(u)
                           B0=B0+fnp*bt(u)
 51                        if(n.eq.1.or.iabs(m).gt.n-1) goto 52
                           u=(n-1)*n+m
                           fn=fnr(n+m)*fnr(n-m)*t
                           A0=A0+fn*at(u)
                           B0=B0+fn*bt(u)
 52                        temp=temp+dble(dconjg(as(i,iL))*A0)
                           temp=temp+dble(dconjg(bs(i,iL))*B0)
                        enddo
                        if(indpol.lt.1) then
                           cscaxi(i)=cscaxi(i)+sc
                           cscax=cscax+sc
                           cprxi(i)=cprxi(i)+temp
                           cprx=cprx+temp
                        else
                           cscayi(i)=cscayi(i)+sc
                           cscay=cscay+sc
                           cpryi(i)=cpryi(i)+temp
                           cpry=cpry+temp
                        endif
                     enddo
                  enddo
C-------------------------------------------------------------------
c  computing total and differential extinction and absorption 
c  cross-sections [see, for example, Xu, Appl. Opt. 36, 9496 (1997)] 
C-------------------------------------------------------------------
                  do j=1,nL
                     cz=dcos(k*r0(3,j))
                     sz=dsin(k*r0(3,j))
                     cmz=dcmplx(cz,-sz)
                     A=0.d0
                     B=0.d0
                     do n=1,nmax(j)
                        rn=fnr(2*n+1)
                        m0=n*n+n+1
                        u0=n*n+n-1
                        A=A+rn*(as(j,m0)+bs(j,m0))
                        B=B+rn*(as(j,u0)-bs(j,u0))
                     enddo
                     if(indpol.lt.1) then
                        cextxi(j)=cextxi(j)+dble((A-B)*cmz)
                        cextx=cextx+dble((A-B)*cmz)
                     else
                        cextyi(j)=cextyi(j)-dimag((A+B)*cmz)
                        cexty=cexty-dimag((A+B)*cmz)
                     endif
                  enddo
                  do j=1,nL
                     do n=1,nmax(j)
                        A=ref(j)*dcmplx(rsr(n,j),-rsi(n,j))
                        temp1=-dimag(A)
                        A=px(n,j)*(ref(j)*rsx(n,j)
     +                    -dcmplx(rsr(n,j),rsi(n,j)))
                        temp=cdabs(A)*cdabs(A)
                        if(temp.eq.0.d0) then 
                           dn=0.d0
                        else
                           dn=temp1/temp
                        endif
                        A=dcmplx(r0(5,j),-r0(6,j))
     +                    *dcmplx(rsr(n,j),-rsi(n,j))
                        temp1=-dimag(A)
                        A=px(n,j)*(rsx(n,j)
     +                    -ref(j)*dcmplx(rsr(n,j),rsi(n,j)))
                        temp=cdabs(A)*cdabs(A)
                        if(temp.eq.0.d0) then 
                           cn=0.d0
                        else
                           cn=temp1/temp
                        endif
                        do m=-n,n
                           i=n*n+n+m      
                           temp1=dn*cdabs(as(j,i))*cdabs(as(j,i))
     +	                        +cn*cdabs(bs(j,i))*cdabs(bs(j,i))
                           if(indpol.lt.1) then
                              cabsxi(j)=cabsxi(j)+temp1
                              cabsx=cabsx+temp1
                           else
                              cabsyi(j)=cabsyi(j)+temp1
                              cabsy=cabsy+temp1
                           endif
                        enddo
                     enddo
                  enddo
C-----------------------------------------------------------------
c  computing two-dimensional angular distribution of the total 
c  scattered field in an array of (npng+1)*nang2 specified by
c  the input of (sang,pang)
c  Note that the definition used here for the scattering amplitude  
c  matrix is slightly different from that given by van de Hulst 
c  and by Bohren & Huffman. The incident field components refer to 
c  the x-z plane instead of the scattering plane defined by the
c  z axis and the scattering direction.      10/20/02
C-----------------------------------------------------------------
                  do i=1,nang
                     iang=2*nang-i
                     dang(i)=sang*dble(i-1)
                     dang(iang)=180.0d0-dang(i)
                     theta=dang(i)*pione/180.0d0
                     xt=dcos(theta)
                     st=dsin(theta)
                     call tipitaud(nmax0,xt)
                     do jc=1,npng
                        azphi=pang*pih*dble(jc-1)/90.0d0
                        sphi=dsin(azphi)
                        cphi=dcos(azphi)
                        do imn=1,nmp
                           at(imn)=0.d0
                           bt(imn)=0.d0
                           atj(imn)=0.d0
                           btj(imn)=0.d0
                        enddo
                        do 31 j=1,nL
                           sb=r0(1,j)*cphi+r0(2,j)*sphi
                           sb=sb*st
                           cb=r0(3,j)*xt
                           cz=k*(sb+cb)
                           sz=k*(sb-cb)
                           A=dcmplx(dcos(cz),-dsin(cz))
                           B=dcmplx(dcos(sz),-dsin(sz))
                           do 32 imn=1,uvmax(j)
                              n=dsqrt(dble(imn))
                              if(n.gt.nmax(j)) goto 31
                              if(idMie.gt.0) then 
                                 m=imn-n*n-n
                                 if(iabs(m).ne.1) goto 32
                              endif
                              at(imn)=at(imn)+A*as(j,imn)
                              bt(imn)=bt(imn)+A*bs(j,imn)
                              atj(imn)=atj(imn)+B*as(j,imn)
                              btj(imn)=btj(imn)+B*bs(j,imn)
 32                        continue
 31                     continue
                        if(indpol.lt.1) then
                           s2x(jc,i)=0.d0
                           s4x(jc,i)=0.d0
                           s2x(jc,iang)=0.d0
                           s4x(jc,iang)=0.d0
                        else
                           s3y(jc,i)=0.d0
                           s1y(jc,i)=0.d0
                           s3y(jc,iang)=0.d0
                           s1y(jc,iang)=0.d0
                        endif
                        A=0.d0
                        B=0.d0
                        Aj=0.d0
                        Bj=0.d0
                        do j=1,nmax0
                           imn=(j-1)*(j+2)/2+1
                           u=j*j+j
                           A=A+at(u)*tau(imn)
                           B=B+bt(u)*tau(imn)
                           if(i.ne.iang) then
                              t=(-1)**(j+1)
                              Aj=Aj+atj(u)*tau(imn)*t
                              Bj=Bj+btj(u)*tau(imn)*t
                           endif
                        enddo
                        if(indpol.lt.1) then
                           s2x(jc,i)=s2x(jc,i)+A
                           s4x(jc,i)=s4x(jc,i)
     +                               +B*dcmplx(0.d0,-1.d0)
                           if(i.ne.iang) then
                              s2x(jc,iang)=s2x(jc,iang)+Aj
                              s4x(jc,iang)=s4x(jc,iang)
     +                              +Bj*dcmplx(0.d0,-1.d0)
                           endif
                        else
                           s3y(jc,i)=s3y(jc,i)-A
                           s1y(jc,i)=s1y(jc,i)
     +                               +B*dcmplx(0.d0,1.d0)
                           if(i.ne.iang) then
                              s3y(jc,iang)=s3y(jc,iang)-Aj
                              s1y(jc,iang)=s1y(jc,iang)+Bj
     +                                  *dcmplx(0.d0,1.d0)
                           endif
                        endif
                        rm=1.d0
                        do 302 m=1,nmax0
                           A=0.d0
                           B=0.d0
                           A2=0.d0
                           B2=0.d0
                           Aj=0.d0
                           Bj=0.d0
                           Aj2=0.d0
                           Bj2=0.d0
                           rm=-rm
                           do 303 j=m,nmax0
                              imn=(j-1)*(j+2)/2+m+1 
                              u=j*j+j+m             
                              v=u-2*m               
                              A0=at(u)*tau(imn)
     +                           +bt(u)*pi(imn)
                              B0=rm*(at(v)*tau(imn)
     +                           -bt(v)*pi(imn))
                              A=A+A0+B0
                              A2=A2+A0-B0
                              A0=at(u)*pi(imn)
     +                           +bt(u)*tau(imn)
                              B0=rm*(at(v)*pi(imn)
     +                           -bt(v)*tau(imn))
                              B=B+A0-B0
                              B2=B2+A0+B0
                              if(i.ne.iang) then
                                 t=(-1)**(j+m+1)
                                 p=-t
                                 A0=atj(u)*tau(imn)*t
     +                              +btj(u)*pi(imn)*p
                                 B0=rm*(atj(v)*tau(imn)*t
     +                              -btj(v)*pi(imn)*p)
                                 Aj=Aj+A0+B0
                                 Aj2=Aj2+A0-B0
                                 A0=atj(u)*pi(imn)*p
     +                              +btj(u)*tau(imn)*t
                                 B0=rm*(atj(v)*pi(imn)*p
     +                              -btj(v)*tau(imn)*t)
                                 Bj=Bj+A0-B0
                                 Bj2=Bj2+A0+B0
                              endif
 303                       continue	         
                           temp=dble(m)*azphi
                           sb=dsin(temp)
                           cb=dcos(temp)
                           if(indpol.lt.1) then
                              s2x(jc,i)=s2x(jc,i)+A*cb
     +                                  +A2*dcmplx(0.d0,1.d0)*sb
                              s4x(jc,i)=s4x(jc,i)
     +                                  +B*dcmplx(0.d0,-1.d0)*cb
     +                                  +B2*sb
                              if(i.ne.iang) then
                                 s2x(jc,iang)=s2x(jc,iang)
     +                              +Aj*cb
     +                              +Aj2*dcmplx(0.d0,1.d0)*sb
                                 s4x(jc,iang)=s4x(jc,iang)
     +                              +Bj*dcmplx(0.d0,-1.d0)*cb
     +                              +Bj2*sb
                              endif
                           else
                              s3y(jc,i)=s3y(jc,i)-A*cb
     +                                  +A2*dcmplx(0.d0,-1.d0)*sb
                              s1y(jc,i)=s1y(jc,i)
     +                                  +B*dcmplx(0.d0,1.d0)*cb
     +                                  -B2*sb
                              if(i.ne.iang) then
                                 s3y(jc,iang)=s3y(jc,iang)-Aj*cb
     +                              +Aj2*dcmplx(0.d0,-1.d0)*sb
                                 s1y(jc,iang)=s1y(jc,iang)
     +                              +Bj*dcmplx(0.d0,1.d0)*cb
     +                              -Bj2*sb
                              endif
                           endif  
 302                    continue
                     enddo
                  enddo

C-----------------------------------------------------------------------
c  computing heat-source functions (i.e., the three-dimensional internal 
c  electric field distributions) for each of the component spheres [see 
c  Xu et al. Physical Review E 60 (1999)]
c  in spherical coordinates when idphoto=1 or in Cartesian coordinates 
c  when idphoto=2
c  idphoto=0 for not calculating  
C-----------------------------------------------------------------------
                  if(iram.eq.1) then
                     call internd(nL,idphoto,nphoto,dphi,istart,
     +                   iend,istep,indpol,idMie,x,nmax,ref0,px,
     +                   rsx,rsr,rsi,NXMAX,py0,py,dpy,as,bs)
                  endif

                  indpol=indpol+2
                  factor=factor2
                  if(indpol.lt.3) then
                     write(6,'(a8,i4,a32)') ' orien.#',iram,
     +                  '  Solving for y-pol. inci. state'
                     goto 18
                  endif
                  do i=1,nang2
                     i22(i)=i22(i)+cdabs(s2x(1,i))*cdabs(s2x(1,i))
                     i21(i)=i21(i)+cdabs(s4x(1,i))*cdabs(s4x(1,i))
                     i11(i)=i11(i)+cdabs(s1y(1,i))*cdabs(s1y(1,i))
                     i12(i)=i12(i)+cdabs(s3y(1,i))*cdabs(s3y(1,i))
                     do jc=1,npng
                        call mueller(s1y(jc,i),s2x(jc,i),s3y(jc,i),
     +                               s4x(jc,i),smue)
                        do j=1,4
                           do m=1,4
                              mue(j,m,jc,i)=mue(j,m,jc,i)+smue(j,m)
                           enddo
                        enddo
                     enddo
                  enddo
                  cbakx=cbakx
     +                  +cdabs(s2x(1,nang2))*cdabs(s2x(1,nang2))
                  cbaky=cbaky
     +                  +cdabs(s1y(1,nang2))*cdabs(s1y(1,nang2))
                  cz=4.0d0/(gcs*k*k)
                  if(idpq.eq.1) then
                     temp1=s2x(1,1)*cz
                     temp2=s1y(1,1)*cz
                     write(12,'(i6,1x,2f6.1,1x,4e14.6)') iram,
     +                  thet/pih*90.d0,phai/pih*90.0d0,temp1,
     +                  dimag(s2x(1,1))*cz,temp2,
     +                  dimag(s1y(1,1))*cz
                     if(idc.gt.0) then 
                        write(6,'(i7,3f7.1)') iram,dang(1),
     +                     dcos(thet),phai*90.d0/pih
                     else 
                        write(6,'(i7,3f7.1)') iram,dang(1),
     +	                     thet*90.d0/pih,phai*90.d0/pih
                     endif
                  endif
                  if(idMie.eq.1) goto 19
                  if(iram.lt.10) then
                     write(cnr1,'(i1)') iram
                     tailn='00'//cnr1
                  else
                     if(iram.lt.100) then
                        write(cnr2,'(i2)') iram
                        tailn='0'//cnr2
                     else
                        write(cnr3,'(i3)') iram
                        tailn=cnr3
                     endif
                  endif
 19            continue
            enddo
         enddo
      enddo
C
C  ending doloop for nbeta,nphi,ntheta
C
      if(nram.eq.1) then
	 write(6,*) 'amplitude scattering ',
     +      'matrix is in the output file amp.out' 
         open(23,file='amp.out',status='unknown')
         write(23,'(a)') 
     +      'amp.out   (amplitude scattering matrix)'
         write(23,'(a,f7.4,3x,a,a20)')  
     +      'wavelength: ',w,'input filename: ',FLNAME
         write(23,'(a)') 
     +      'sphere#,x,y,z,r,complex refractive index:'
	 do i=1,nL,nL
            write(23,'(i5,6f11.4)') i,r0(1,i),r0(2,i),
     +                r0(3,i),r0(4,i),r0(5,i),r0(6,i)
         enddo
         write(23,'(a)') 
     +   'scattering angle, s2x(complex), s3y(complex)'
         write(23,'(a)') 
     +   '                  s4x(complex), s1y(complex)'
         do jc=1,npng
            temp=pang*dble(jc-1)
            write(23,'(a18,f8.3)') 
     +         'phi (in degrees): ',temp
            do i=1,nang2  
               write(23,'(f8.2,4e15.6)') dang(i),
     +            dble(s2x(jc,i)),dimag(s2x(jc,i)),
     +            dble(s3y(jc,i)),dimag(s3y(jc,i))
               write(23,'(8x,4e15.6)') 
     +            dble(s4x(jc,i)),dimag(s4x(jc,i)),
     +            dble(s1y(jc,i)),dimag(s1y(jc,i))
            enddo
         enddo
      endif  
      if(iram.ne.nram) then
         write(6,*) 'Note: iram is not equal to nram!'
         write(6,'(a,i6,a,i6)') 'iram: ',iram,'nram: ',nram
      endif
      if(idpq.eq.1) close(12)
      cz=nram
      do i=1,nang2
         i11(i)=i11(i)/cz
         i21(i)=i21(i)/cz
         i22(i)=i22(i)/cz
         i12(i)=i12(i)/cz
         inat(i)=i11(i)+i22(i)+i12(i)+i21(i)
         pol(i)=(i11(i)+i12(i)-i22(i)-i21(i))/inat(i)
         do jc=1,npng
            do j=1,4
               do m=1,4
                  mue(j,m,jc,i)=mue(j,m,jc,i)/cz
               enddo
            enddo
         enddo
      enddo
      cz=cz*k*k
      cscax=2.0d0*twopi*cscax/cz
      cscay=2.0d0*twopi*cscay/cz
      csca=0.5d0*(cscax+cscay)
      cextx=twopi*cextx/cz
      cexty=twopi*cexty/cz
      cext=0.5d0*(cextx+cexty)
      cabsx=2.0d0*twopi*cabsx/cz
      cabsy=2.0d0*twopi*cabsy/cz
      cabs=0.5d0*(cabsx+cabsy)
      cprx=2.0d0*twopi*cprx/cz
      cpry=2.0d0*twopi*cpry/cz
      cpr=0.5d0*(cprx+cpry)
      assym=(cprx+cpry)/(cscax+cscay)
      assym0=0.5d0*(cprx/cscax+cpry/cscay)
      cbakx=2.0d0*twopi*cbakx/cz
      cbaky=2.0d0*twopi*cbaky/cz
      cbak=0.5d0*(cbakx+cbaky)
      write(6,'(5x,2e15.6)') assym,assym0
      write(6,'(6x,a4,9x,a4,9x,a4,9x,a4,9x,a3,8x,a12)')
     +   'Cext','Cabs','Csca','Cbak','Cpr','<cos(theta)>'
      write(6,'(2x,6e13.5)') cext,cabs,csca,cbak,cext-cpr,assym
      cscax=0.d0
      cscay=0.d0
      cextx=0.d0
      cexty=0.d0
      cabsx=0.d0
      cabsy=0.d0
      cprx=0.d0
      cpry=0.d0
      do i=1,nL
         cscax=cscax+cscaxi(i)
         cscay=cscay+cscayi(i)
         cextx=cextx+cextxi(i)
         cabsx=cabsx+cabsxi(i)
         cexty=cexty+cextyi(i)
         cabsy=cabsy+cabsyi(i)
         cprx=cprx+cprxi(i)
         cpry=cpry+cpryi(i)
      enddo
      assymx=cprx/cscax
      assymy=cpry/cscay
      assym0=0.5d0*(assymx+assymy)
      cscax=2.0d0*twopi*cscax/cz
      cscay=2.0d0*twopi*cscay/cz
      csca=0.5d0*(cscax+cscay)
      cextx=twopi*cextx/cz
      cexty=twopi*cexty/cz
      cext=0.5d0*(cextx+cexty)
      cabsx=2.0d0*twopi*cabsx/cz
      cabsy=2.0d0*twopi*cabsy/cz
      cabs=0.5d0*(cabsx+cabsy)
      cprx=2.0d0*twopi*cprx/cz
      cpry=2.0d0*twopi*cpry/cz
      cpr=0.5d0*(cprx+cpry)
      assym=cpr/csca	
      write(6,'(2x,6e13.5)') cext,cabs,csca,cbak,cext-cpr,assym
      do i=1,nL
         cabsxi(i)=4.0d0*pione*cabsxi(i)/cz
         cabsyi(i)=4.0d0*pione*cabsyi(i)/cz
         cextxi(i)=2.0d0*pione*cextxi(i)/cz
         cextyi(i)=2.0d0*pione*cextyi(i)/cz
         cscaxi(i)=4.0d0*pione*cscaxi(i)/cz
         cscayi(i)=4.0d0*pione*cscayi(i)/cz
         cprxi(i)=4.0d0*pione*cprxi(i)/cz
         cpryi(i)=4.0d0*pione*cpryi(i)/cz
         cscai(i)=0.5d0*(cscaxi(i)+cscayi(i))
         cexti(i)=0.5d0*(cextxi(i)+cextyi(i))
         cabsi(i)=0.5d0*(cabsxi(i)+cabsyi(i))
         cpri(i)=0.5d0*(cprxi(i)+cpryi(i))
         cpri(i)=cscai(i)+cabsi(i)-cpri(i)
         assymxi(i)=cprxi(i)/cscaxi(i)
         assymyi(i)=cpryi(i)/cscayi(i)
         assymi(i)=0.5d0*(cprxi(i)+cpryi(i))/csca
         cprxi(i)=cscaxi(i)+cabsxi(i)-cprxi(i)
         cpryi(i)=cscayi(i)+cabsyi(i)-cpryi(i)
         write(6,'(i5,5e15.6)') i,cexti(i),cabsi(i),
     +                   cscai(i),cpri(i),assymi(i)
      enddo
      write(6,'(/,a)') 'efficiencies for radiation pressure'
      write(6,'(5x,6e13.5)') assym,assym0
      write(6,'(5x,6e13.5)') assym,cext-cpr,assymx,cextx-cprx,
     +                                      assymy,cexty-cpry
      do i=1,nL
         write(6,'(i5,6e13.5)') i,assymi(i),cpri(i),
     +       assymxi(i),cprxi(i),assymyi(i),cpryi(i)  	   
      enddo
      betami=betami*90.d0/pih
      betamx=betamx*90.d0/pih
      thetmi=thetmi*90.d0/pih
      thetmx=thetmx*90.d0/pih
      phaimi=phaimi*90.d0/pih
      phaimx=phaimx*90.d0/pih
      flout='cr'//fileout
      open(33,file=flout,status='unknown')
      write(33,'(a20,a47)') 
     +   flout,'(Total and individual-particle cross sections)'
      write(33,'(a32,2x,a22)') 
     +   'input sphere-aggregate filename:',FLNAME
      write(33,'(a19,3i5)') 'nbeta,nthet,nphai: ',
     +                       nbeta,nthet,nphai
      write(33,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +          betami,betamx,thetmi,thetmx,phaimi,phaimx
      write(33,'(a28,i5)') '# of orientations averaged: ',nram
      write(33,'(12x,a4,11x,a4,11x,a4,11x,a3,8x,a12)')
     +    'Cext','Cabs','Csca','Cpr','<cos(theta)>'
      write(33,'(a5,5e15.6)') 'total',cext,cabs,csca,cext-cpr,
     +     assym
      do i=1,nL
         write(33,'(i5,5e15.6)') i,cexti(i),cabsi(i),cscai(i),
     +           cpri(i),assymi(i)
      enddo
      close(33)       
      cz=pione*gcvr*gcvr
      assym=cpr/csca
      assymx=cprx/cscax
      assymy=cpry/cscay
      cabsxv=cabsx/cz
      cabsyv=cabsy/cz
      cextxv=cextx/cz
      cextyv=cexty/cz
      cscaxv=cscax/cz
      cscayv=cscay/cz
      cprxv=cprx/cz
      cprxv=cextxv-cprxv
      cpryv=cpry/cz
      cpryv=cextyv-cpryv
      cscav=0.5d0*(cscaxv+cscayv)
      cextv=0.5d0*(cextxv+cextyv)
      cabsv=0.5d0*(cabsxv+cabsyv)
      cprv=0.5d0*(cprxv+cpryv)
      cbakxv=cbakx/cz
      cbakyv=cbaky/cz
      cbakv=0.5d0*(cbakxv+cbakyv)
      temp=gcvr*gcvr/gcs
      cabsxs=cabsxv*temp
      cabsys=cabsyv*temp
      cextxs=cextxv*temp
      cextys=cextyv*temp
      cscaxs=cscaxv*temp
      cscays=cscayv*temp
      cprxs=cprxv*temp
      cprys=cpryv*temp
      cscas=cscav*temp
      cexts=cextv*temp
      cabss=cabsv*temp
      cprs=cprv*temp
      cbakxs=cbakxv*temp
      cbakys=cbakyv*temp
      cbaks=cbakv*temp
 222  format(1x,a1,6e13.5)
 221  format(6x,a5,8x,a5,8x,a5,8x,a5,8x,a4,5x,a12)        
      write(6,221)
     +   'Qextv','Qabsv','Qscav','Qbakv','Qprv','<cos(theta)>'	
      write(6,222) 't',cextv,cabsv,cscav,cbakv,cprv,assym
      write(6,222) 'x',cextxv,cabsxv,cscaxv,cbakxv,cprxv,assymx
      write(6,222) 'y',cextyv,cabsyv,cscayv,cbakyv,cpryv,assymy
      write(6,221)
     +   'Qexts','Qabss','Qscas','Qbaks','Qprs','<cos(theta)>'
      write(6,222) 't',cexts,cabss,cscas,cbaks,cprs,assym
      write(6,222) 'x',cextxs,cabsxs,cscaxs,cbakxs,cprxs,assymx
      write(6,222) 'y',cextys,cabsys,cscays,cbakys,cprys,assymy
      temp=-(cabs+csca-cext)/cext
      write(6,'(/,a37,e14.5)') 
     +   'Accuracy of this numerical solution: ',temp
      open(12,file=fileout,status='unknown')
      write(12,'(a20,a16,a18,a4,f8.3,a5,f8.3)') fileout, 
     +   '--- input file: ',FLNAME,' xv:',xv,'  xs:',xs
      write(12,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +          betami,betamx,thetmi,thetmx,phaimi,phaimx        
      write(12,'(a19,3i5,6x,a28,i5)') 
     +   'nbeta,nthet,nphai: ',nbeta,nthet,nphai,
     +   '# of orientations averaged: ',nram
      write(12,221)
     +   'Cext','Cabs','Csca','Cbak','Cpr','<cos(theta)>'
      write(12,222) 
     +   't',cext,cabs,csca,cbak,cext-cpr,assym
      write(12,222)
     +   'x',cextx,cabsx,cscax,cbakx,cextx-cprx,assymx
      write(12,222)
     +   'y',cexty,cabsy,cscay,cbaky,cexty-cpry,assymy
      write(12,221)
     +   'Qextv','Qabsv','Qscav','Qbakv','Qprv','<cos(theta)>'
      write(12,222) 't',cextv,cabsv,cscav,cbakv,cprv,assym
      write(12,222) 'x',cextxv,cabsxv,cscaxv,cbakxv,cprxv,assymx
      write(12,222) 'y',cextyv,cabsyv,cscayv,cbakyv,cpryv,assymy
      write(12,221)
     +   'Qexts','Qabss','Qscas','Qbaks','Qprs','<cos(theta)>'
      write(12,222) 't',cexts,cabss,cscas,cbaks,cprs,assym
      write(12,222) 'x',cextxs,cabsxs,cscaxs,cbakxs,cprxs,assymx
      write(12,222) 'y',cextys,cabsys,cscays,cbakys,cprys,assymy
      write(12,'(/)')
      write(12,*) 
     +'For the scattering plane of phi=0 (i.e., the x-z plane) only:'
      write(12,'(2x,a4,5x,a6,6x,a4,6x,a5,8x,a5,8x,a5,8x,a5)')
     +    's.a.','itotal','pol.','S1*S1','S4*S4','S3*S3','S2*S2'
      do i=1,nang2
         write(12,'(f7.1,e13.5,f8.4,4e13.5)') 
     +      dang(i),inat(i),pol(i),i11(i),i21(i),i12(i),i22(i) 
      enddo
      close(12)
      open(11,file='mueller.out',status='unknown')
      write(11,'(a11,6x,a16)') 'mueller.out','(Mueller matrix)'
      write(11,'(a16,a20)') 'Input filename: ',FLNAME
      write(11,'(a19,3i5)') 'nbeta,nthet,nphai: ',nbeta,nthet,nphai
      write(11,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +                betami,betamx,thetmi,thetmx,phaimi,phaimx
      write(11,'(a28,i5)') '# of orientations averaged: ',nram
      write(11,'(a30,f9.3)') 'interval of scattering angle: ',sang
      write(11,'(a30,f9.3)') 'interval of azimuth angle:    ',pang
      do jc=1,npng
         t=pang*dble(jc-1)
         write(11,'(a18,f8.3)') 'phi (in degrees): ',t
         do i=1,nang2
            write(11,'(f7.1,4e16.7)') 
     +         dang(i),mue(1,1,jc,i),mue(1,2,jc,i),mue(1,3,jc,i),
     +                 mue(1,4,jc,i)
            write(11,'(7x,4e16.7)') 
     +                 mue(2,1,jc,i),mue(2,2,jc,i),mue(2,3,jc,i),
     +                 mue(2,4,jc,i)
            write(11,'(7x,4e16.7)') 
     +                 mue(3,1,jc,i),mue(3,2,jc,i),mue(3,3,jc,i),
     +                 mue(3,4,jc,i)
            write(11,'(7x,4e16.7)') 
     +                 mue(4,1,jc,i),mue(4,2,jc,i),mue(4,3,jc,i),
     +                 mue(4,4,jc,i)
         enddo
      enddo
      write(11,*) 'phi (in degrees): 360'
      do i=1,nang2
         write(11,'(f7.1,4e16.7)') 
     +      dang(i),mue(1,1,1,i),mue(1,2,1,i),mue(1,3,1,i),
     +              mue(1,4,1,i)
         write(11,'(7x,4e16.7)') 
     +              mue(2,1,1,i),mue(2,2,1,i),mue(2,3,1,i),
     +              mue(2,4,1,i)
         write(11,'(7x,4e16.7)') 
     +              mue(3,1,1,i),mue(3,2,1,i),mue(3,3,1,i),
     +              mue(3,4,1,i)
         write(11,'(7x,4e16.7)') 
     +              mue(4,1,1,i),mue(4,2,1,i),mue(4,3,1,i),
     +              mue(4,4,1,i)
      enddo
      close(11)
      STOP
      END

      SUBROUTINE orientcd(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,
     &                  MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI,
     &                  BETA,THETA,PHI)
C Arguments:
      INTEGER MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI
      double precision BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX
      double precision BETA(MXBETA),THETA(MXTHET),PHI(MXPHI)
C Local variables:
      INTEGER J
      double precision DELTA
C***********************************************************************
C Given: BETAMI=minimum value of beta (radians)
C        BETAMX=maximum value of beta (radians)
C        THETMI=minimum value of theta (radians)
C        THETMX=maximum value of theta (radians)
C        PHIMIN=minimum value of phi (radians)
C        PHIMAX=maximum value of phi (radians)
C        MXBETA,MXTHET,MXPHI=dimensions of the arrays BETA,THETA,PHI
C        NBETA=number of values of beta
C        NTHETA=number of values of theta
C        NPHI=number of values of PHI
C Returns: BETA(1-NBETA)=beta values (radians)
C          THETA(1-NTHETA)=theta values (radians)
C          PHI(1-NPHI)=phi values (radians)
C Purpose: to generate a sequence of desired target orientations
C Present version assumes:
C        beta to be uniformly distributed between BETAMI and BETAMX
C        cos(theta) to be uniformly distributed between cos(THETMI) and
C                   cos(THETMX)
C        phi to be uniformly distributed between PHIMIN and PHIMAX
C        If NPHI=1, first angle is THETMI, last angle is THETMX
C        If NPHI>1, then angles are midpoints of intervals of equal
C            range in theta subdividing range from THETMI to THETMX
C***********************************************************************
      BETA(1)=BETAMI            
      IF(NBETA.GT.1)THEN
         DELTA=(BETAMX-BETAMI)/DBLE(NBETA-1)
         DO 1000 J=2,NBETA
           BETA(J)=BETA(1)+DELTA*DBLE(J-1)
 1000    CONTINUE
      ENDIF
      IF(NPHI.EQ.1.AND.NTHETA.GT.1)THEN
         DELTA=(DCOS(THETMX)-DCOS(THETMI))/DBLE(NTHETA-1)
         THETA(1)=THETMI
      ELSE
         DELTA=(DCOS(THETMX)-DCOS(THETMI))/DBLE(NTHETA)
         THETA(1)=DACOS(DCOS(THETMI)+0.5d0*DELTA)
      ENDIF
      IF(NTHETA.GT.1)THEN
         DO 2000 J=2,NTHETA
            THETA(J)=DACOS(DCOS(THETA(1))+DELTA*DBLE(J-1))
 2000    CONTINUE
      ENDIF
c      DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI)
c      PHI(1)=PHIMIN+0.5D0*DELTA
      PHI(1)=PHIMIN
      IF(NPHI.GT.1)THEN
         DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI-1)
         DO 3000 J=2,NPHI
            PHI(J)=PHI(1)+DELTA*DBLE(J-1)
 3000    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE orientud(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,
     &                  MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI,
     &                  BETA,THETA,PHI)
C Arguments:
      INTEGER MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI
      double precision BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX
      double precision BETA(MXBETA),THETA(MXTHET),PHI(MXPHI)
C Local variables:
      INTEGER J
      double precision DELTA
C***********************************************************************
C Given: BETAMI=minimum value of beta (radians)
C        BETAMX=maximum value of beta (radians)
C        THETMI=minimum value of theta (radians)
C        THETMX=maximum value of theta (radians)
C        PHIMIN=minimum value of phi (radians)
C        PHIMAX=maximum value of phi (radians)
C        MXBETA,MXTHET,MXPHI=dimensions of the arrays BETA,THETA,PHI
C        NBETA=number of values of beta
C        NTHETA=number of values of theta
C        NPHI=number of values of PHI
C Returns:  BETA(1-NBETA)=beta values (radians)
C           THETA(1-NTHETA)=theta values (radians)
C           PHI(1-NPHI)=phi values (radians)
C Note: it is assumed that target orientation weight function
C       can be factored into WGTA*WGTB -- i.e., that rotations
C       around a1 are decoupled from orientation of a1.
C Purpose: to generate a sequence of desired target orientations
C Present version assumes:
C        beta to be uniformly distributed between BETAMI and BETAMX
C        theta to be uniformly distributed between THETMI and THETMX
C        phi to be uniformly distributed between PHIMIN and PHIMAX
C            first angle is THETMI, last angle is THETMX
C***********************************************************************
      BETA(1)=BETAMI
      IF(NBETA.GT.1)THEN
         DELTA=(BETAMX-BETAMI)/DBLE(NBETA-1)
         DO 1000 J=2,NBETA
           BETA(J)=BETA(1)+DELTA*DBLE(J-1)
 1000    CONTINUE
      ENDIF
      THETA(1)=THETMI
      IF(NTHETA.GT.1)THEN
         DELTA=(THETMX-THETMI)/DBLE(NTHETA-1)
         DO 2000 J=2,NTHETA
            THETA(J)=THETA(1)+DELTA*DBLE(J-1)
 2000    CONTINUE
      ENDIF
      PHI(1)=PHIMIN
      IF(NPHI.GT.1)THEN
         DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI-1)
         DO 3000 J=2,NPHI
            PHI(J)=PHI(1)+DELTA*DBLE(J-1)
 3000    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE abMiexud(X,REFREL,NP,NMAX,NM,AN,BN,NADD,
     +                    RSR,RSI,RSX,PX,AR,AI,BR,BI,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NMAX,NM,NADD,NSTOP,NX,K,N
      COMPLEX*16 REFREL,AN(NP),BN(NP)
      DOUBLE PRECISION AR(NP),AI(NP),BR(NP),BI(NP),
     +   RSR(NMAX),RSI(NMAX),RSX(NMAX),PX(NMAX)
      IF(EPS.GT.1.D0.OR.EPS.LT.0.D0) EPS=1.0D-20
      IF(NADD.NE.0) EPS=0.D0
      CTC=EPS
      XM=DBLE(REFREL)
      YM=DIMAG(REFREL)
      XMX=X*XM
      YMX=X*YM
      RP2=XMX*XMX+YMX*YMX
      NSTOP=X+4.D0*X**.3333D0
      NSTOP=NSTOP+2
      NM=NSTOP+NADD
      XN=DSQRT(XM**2+YM**2)*X
      NX=1.1D0*XN+10.D0
      if(NX-NM.lt.10) NX=NM+10
      write(6,*) 'Wiscombe criterion:',NSTOP 
      write(6,*) 'NADD:',NADD
      write(6,*) 'NX:',NX
      IF(NX.GT.NMAX) THEN
         WRITE(6,*) 'Parameter NXMAX too small'
         WRITE(6,*) '  NXMAX must be greater than', NX
         WRITE(6,*) 'Please correct NXMAX in main code,' 
         WRITE(6,*) '  recompile, then try again'
         STOP
      ENDIF
      IF(NM.GT.NP) THEN
         WRITE(6,*) 'Parameter np too small'
         WRITE(6,*) '  np must be greater than', NM
         WRITE(6,*) 'Please correct np in gmm01f.par,' 
         WRITE(6,*) '  recompile the code, then try again'
         STOP
      ENDIF
C  DOWN RECURSION FOR RATIOS RSR,RSI,RSX,PNR,PNI,PX
      PNX=X/DBLE(2*NX+3)
      PNR=XMX/DBLE(2*NX+3)
      PNI=YMX/DBLE(2*NX+3)
      DO 5 K=1,NX
         N=NX-K+1
         CN=DBLE(N)
         ALN=(2.D0*CN+1.D0)*XMX/RP2-PNR
         BEN=(2.D0*CN+1.D0)*YMX/RP2+PNI
         RSR(N)=-CN*XMX/RP2+ALN
         RSI(N)=CN*YMX/RP2-BEN
         PZD=ALN*ALN+BEN*BEN
         PNR=ALN/PZD
         PNI=BEN/PZD
         RSX(N)=(CN+1.D0)/X-PNX
         IF(N.EQ.1) GO TO 20
         PNX=X/(2.D0*CN+1.D0-PNX*X)
         PX(N)=PNX
    5    CONTINUE
   20 SNM1X=DSIN(X)
      CNM1X=DCOS(X)
      IF(X-0.1D0)21,22,22
   21 SNX=X**2./3.D0-X**4./30.D0+X**6./840.D0-X**8./45360.D0
      GO TO 23
   22 SNX=SNM1X/X-CNM1X
   23 CNX=CNM1X/X+SNM1X
      DO 10 N=1,NX
         PX(N)=SNX      !preparing for the calculation of Cabs
         C=DBLE(N)
         DCNX=CNM1X-C*CNX/X
         DSNX=RSX(N)*SNX
C  CALCULATION OF EXTERIOR COEFFICIENTS AN AND BN
         ANNR=RSR(N)*SNX-XM*DSNX
         ANNI=RSI(N)*SNX-YM*DSNX
         TA1=RSR(N)*SNX-RSI(N)*CNX
         TA2=RSI(N)*SNX+RSR(N)*CNX
         ANDR=TA1-XM*DSNX+YM*DCNX
         ANDI=TA2-XM*DCNX-YM*DSNX
         AND=ANDR*ANDR+ANDI*ANDI
         BNNR=(XM*RSR(N)-YM*RSI(N))*SNX-DSNX
         BNNI=(XM*RSI(N)+YM*RSR(N))*SNX
         TB1=RSR(N)*SNX-RSI(N)*CNX
         TB2=RSR(N)*CNX+RSI(N)*SNX
         BNDR=XM*TB1-YM*TB2-DSNX
         BNDI=XM*TB2+YM*TB1-DCNX
         BND=BNDR*BNDR+BNDI*BNDI
         AR(N)=(ANNR*ANDR+ANNI*ANDI)/AND
         AI(N)=(ANNI*ANDR-ANNR*ANDI)/AND
         BR(N)=(BNNR*BNDR+BNNI*BNDI)/BND
         BI(N)=(BNNI*BNDR-BNNR*BNDI)/BND
C  MIE SERIES CONVERGENCE TEST IS MADE BY TESTING AN'S AND BN'S
         TI=AR(N)*AR(N)+AI(N)*AI(N)+BR(N)*BR(N)+BI(N)*BI(N)
         TI=TI/(AR(1)*AR(1)+AI(1)*AI(1)+BR(1)*BR(1)+BI(1)*BI(1))
         IF(TI-CTC) 16,18,18
   18    IF(NM-N) 15,15,6
    6    IF(N-NX)7,15,15
    7    M=N+1
         SNX=PX(M)*SNX
         CNM2X=CNM1X
         CNM1X=CNX
         CNX=(2.D0*C+1.D0)*CNM1X/X-CNM2X
   10 CONTINUE
      GO TO 15
   16 WRITE(6,*) '*** NOTE THAT THE FIELD-EXPANSION TRANCATION'
      WRITE(6,*) '*** IS DETERMINED BY eps GIVEN IN THE INPUT'
      WRITE(6,*) '*** FILE gmm01f.in'
      WRITE(6,*) '*** IN CASE YOU NEED A HIGHER ORDER, eps MUST'
      WRITE(6,'(a,e9.1)')
     +          ' *** BE SMALLER THAN THE CURRENT VALUE',EPS
   15 NM=N
      DO I=1,NM
         AN(I)=DCMPLX(AR(I),-AI(I))
         BN(I)=DCMPLX(BR(I),-BI(I))
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION ran1d(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,
     *IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,
     *EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1d=dmin1(AM*iy,RNMX)
      return
      END

	subroutine cofsrd(nmax)
	include 'gmm01f.par'
	parameter (nmp=np*(np+2))	
	double precision cofsr(nmp),lnfacd,c
	common/crot/cofsr
        i=0
	do n=1,nmax
	   do m=-n,n
	      i=i+1       
	      c=lnfacd(dble(n-m))-lnfacd(dble(n+m))
	      cofsr(i)=0.5d0*c
c	      c=0.5d0*c
c	      cofsr(i)=dexp(c)
           enddo
        enddo
	return
	end

	subroutine cofd0(nmax)
	implicit double precision (a-h,o-z)	
	include 'gmm01f.par'
	parameter (nmp=np*(np+2))	
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)	
	integer v
	double precision lnfacd
	common/cofmnv0/cof0(ni0)
	common/crot/cofsr(nmp)	
        common/fnr/fnr(0:2*(np+2))
	i=0
	sm=-0.5d0*dble((-1)**nmax)
	do m=-nmax,nmax
	   ns=max(1,iabs(m))
	   sm=-sm
	   do n=ns,nmax
	      inm=n*(n+1)-m
	      do v=ns,nmax
	         i=i+1
	         ivm=v*(v+1)+m
	         c=cofsr(inm)+cofsr(ivm)	            
	         c=sm*dexp(c)
	         c0=fnr(2*n+1)*fnr(2*v+1)
	         c1=fnr(n)*fnr(v)*fnr(n+1)*fnr(v+1)
	         c0=c0/c1
                 cof0(i)=c*c0
              enddo
           enddo
        enddo
	return
	end

	subroutine cofnv0(nmax)
	include 'gmm01f.par'
	integer n,v
	double precision c1,lnfacd,cnv(np,np)
	common/cfnv/cnv
	do n=1,nmax
	   do v=n,nmax
	      c1=lnfacd(dble(2*n))+lnfacd(dble(2*v))
	      c1=c1-lnfacd(dble(2*n+2*v))
	      c1=c1+2.d0*lnfacd(dble(n+v))
	      c1=c1-lnfacd(dble(n))-lnfacd(dble(v))
	      cnv(n,v)=c1
	   enddo
	enddo
	return
	end

C  subroutine gau0.f generates tabulated values for  
C  Gaunt coefficients up to n=v=n_max
	subroutine gau0(nmax)
	include 'gmm01f.par'
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
	integer v,qmax,uvmax,iga0(ni0)
	double precision ga0(ng0)
	common/g0/ga0
	common/ig0/iga0
        na=0
        uvmax=nmax*(nmax+2)
        i=0
	do m=-nmax,nmax
	   ns=max(1,iabs(m))
	   do n=ns,nmax
	      do v=ns,nmax
	         call gxurcd0(-m,n,v,qmax,na)
	         i=i+1
                 iga0(i)=na
	         na=na+qmax+1
              enddo
	   enddo
        enddo
	return
	end

c  transforms the rectangular coordinates (x,y,z)
c  to spherical coordinates (r,theta,phi)
	subroutine carsphd(x,y,z,r,xt,sphi,cphi)
	double precision x,y,z,r,xt,sphi,cphi
	r=dsqrt(x*x+y*y+z*z)
	if(r.eq.0.d0) then
	   xt=1.d0
	   sphi=0.d0
	   cphi=1.d0
	   return
	endif
	xt=z/r
	if(y.eq.0.d0.and.x.eq.0.d0) then
	   sphi=0.d0
	   cphi=1.d0
	   return
	endif
	sphi=dsqrt(x*x+y*y)
	cphi=x/sphi
	sphi=y/sphi
	return
	end

c  subroutine besseljd.f  (in double precision arithmetic)
c  returns an array of the spherical Bessel function of the  
c  first kind with a real argument: j_0,j_1,j_2,...,j_{NC}
c  uses Ru Wang's ratio method for the downward recursive 
c  calculation of the Riccati-Bessel function Psi_n(z)=z j_n(z) 
c  [see Xu et al., Physical Review E, v.60, 2347-2365 (1999)]
      SUBROUTINE besseljd(NC,X,BESJ)
      INTEGER NC,NX,K,N
      DOUBLE PRECISION X,BESJ(0:NC),PN,CN,X2
      DO K=1,NC
         BESJ(K)=0.D0
      ENDDO
      IF(DABS(X).LT.1.D-12) THEN
         BESJ(0)=1.D0
         BESJ(1)=1.D0/3.D0*X
         RETURN
      ENDIF
c  down-recursively calculating an array of the ratio functions  
c  P_{NC},P_{NC-1},...,P(2) stored in the same array for j_n, 
c  starting with an asymptotic value P_{NX+1}=X/(2 NX+3) at the 
c  highest order NX+1, where NX=NC+1.1X+10 
      NX=1.1D0*X+10.D0
      NX=NC+NX   
      PN=X/DBLE(2*NX+3)
      DO 5 K=1,NX-1
         N=NX-K+1
         CN=DBLE(N)
         PN=X/(DBLE(2*N+1)-PN*X)
         IF(N.GT.NC) GOTO 5
         BESJ(N)=PN
    5 CONTINUE
C  calculating j_0(x) and j_1(x)
      IF(DABS(X)-0.1D0) 10,11,11
   10 X2=X*X
      BESJ(0)=1.D0-X2/72.D0
      BESJ(0)=1.D0-X2/42.D0*BESJ(0)
      BESJ(0)=1.D0-X2/20.D0*BESJ(0)
      BESJ(0)=1.D0-X2/6.D0*BESJ(0)
      BESJ(1)=1.D0/45360.D0-X2/3991680.D0
      BESJ(1)=1.D0/840.D0-X2*BESJ(1)
      BESJ(1)=1.D0/30.0d0-X2*BESJ(1)
      BESJ(1)=X*(1.D0/3.0d0-X2*BESJ(1))
      GOTO 12
   11 BESJ(0)=DSIN(X)/X
      BESJ(1)=(DSIN(X)/X-DCOS(X))/X      
c  calculating j_2,...,j_{NC} 
   12 DO 20 N=2,NC
         BESJ(N)=BESJ(N)*BESJ(N-1)
   20 CONTINUE
      RETURN
      END

c  sub. besselyd.f  (in double precision arithmetic)
c  returns an array of the spherical Bessel function of
c  the second kind with a real argument: y_0,y_1,...,y_n
       subroutine besselyd(n,x,besy)
       integer i,n
       double precision x,besy(0:n),besyn,x2
       if(x.eq.0.d0) then
         write(6,*) 'bad argument in sub. besselyd'
         stop
       endif
       if(dabs(x)-0.1d0)10,11,11
  10   x2=x*x
       besyn=1.d0-x2/72.d0
       besyn=1.d0-x2/42.d0*besyn
       besyn=1.d0-x2/20.d0*besyn
       besyn=1.d0-x2/6.d0*besyn
       besy(0)=1.d0-x2/56.d0
       besy(0)=1.d0-x2/30.d0*besy(0)
       besy(0)=1.d0-x2/12.d0*besy(0)
       besy(0)=1.d0-0.5d0*x2*besy(0)
       besy(0)=-besy(0)/x
       goto 12
  11   besyn=dsin(x)/x
       besy(0)=-dcos(x)/x
  12   besy(1)=besy(0)/x-besyn
       do i=2,n
         besy(i)=dble(2*i-1)/x*besy(i-1)-besy(i-2)
       enddo
       return
       end

c  subroutine rotcoef is originally written by Mackowski, Fuller, and 
c  Mishchenko (taken from the code scsmfo1b.for released to public by 
c  the authors at 8/1999)
c  the rotational coefficients are required for the implementation of 
c  Mackowski's three-step numerical technique in subroutine rtr.f for 
c  decomposition of translation matrix into rotational and axial 
c  translational parts [see Mackowski, Proc. R. Soc. Lond. A 433, 599 
c  (1991)]
c  Yu-lin Xu   12/2000         

      subroutine rotcoef(cbe,nmax)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2))
      implicit double precision (a-h,o-z)
      double precision dk0(-2*np:2*np),dk01(-2*np:2*np)
      common/rot/bcof(0:np+2),dc(-np:np,0:nmp)
      common/fnr/fnr(0:2*(np+2))
      sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
      cbe2=.5d0*(1.d0+cbe)
      sbe2=.5d0*(1.d0-cbe)
      in=1
      dk0(0)=1.d0
      sben=1.d0
      dc(0,0)=1.d0
      dk01(0)=0.d0
      do n=1,nmax
         nn1=n*(n+1)
         in=-in
         sben=sben*sbe/2.d0
         dk0(n)=dble(in)*sben*bcof(n)
         dk0(-n)=dble(in)*dk0(n)
         dk01(n)=0.d0
         dk01(-n)=0.d0
         dc(0,nn1+n)=dk0(n)
         dc(0,nn1-n)=dk0(-n)
         do k=-n+1,n-1
            kn=nn1+k
            dkt=dk01(k)
            dk01(k)=dk0(k)
            dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)
     1             *fnr(n+k-1)*dkt)/(fnr(n+k)*fnr(n-k))
            dc(0,kn)=dk0(k)
         enddo
         im=1
         do m=1,n
            im=-im
            fmn=1.d0/fnr(n-m+1)/fnr(n+m)
            m1=m-1
            dkm0=0.d0
            do k=-n,n
               kn=nn1+k
               dkm1=dkm0
               dkm0=dc(m1,kn)
               if(k.eq.n) then
                  dkn1=0.d0
               else
                  dkn1=dc(m1,kn+1)
               endif
               dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1
     1           -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1
     1              -dble(k)*sbe*dc(m1,kn))*fmn
               dc(-m,nn1-k)=dble((-1)**(k)*im)*dc(m,kn)
            enddo
         enddo
      enddo
      return
      end

c  subroutine cofxuds0.f returns the two classes of vector 
c  (axial) translation coefficients for a given combination of 
c  (m,n,m,v) and a given dimensionless translation distance kd 
cu uses subroutine gid0.f 
	subroutine cofxuds0(nmax,m,n,v,sja,sya,A,B,Aj,Bj)
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
	integer v,p,qmax
	double precision sja(0:n+v+1),sya(0:n+v+1)
	complex*16 A,B,Aj,Bj,signz
	common/ig0/iga0(ni0)
	common/g0/ga0(ng0)
	common/cofmnv0/cof0(ni0)
	fa(m,p)=dble(-2*m*p*(p-1))
	fb(n,v,p)=dble(p*p-(n+v+1)*(n+v+1))
     +            *dble(p*p-(n-v)*(n-v))/dble(4*p*p-1)
	if(iabs(m).gt.n.or.iabs(m).gt.v) then 
	   write(6,*) '|m|>n or v in subroutine cofxuds0.f'
	   stop
	endif
	A=0.d0
	B=0.d0
	Aj=0.d0
	Bj=0.d0
	call gid0(nmax,m,n,v,id)
	c=cof0(id)
	ig=iga0(id)
	nv2=v*(v+1)+n*(n+1)
	signz=dcmplx(0.d0,1.d0)**(n+v)
	p=n+v+2
	qmax=min(n,v)
	do i=1,qmax+1
	   p=p-2
	   cp=dble(nv2-p*(p+1))*ga0(ig+i)
	   sj=sja(p)
	   sy=sya(p)
	   A=A+dcmplx(sj,sy)*signz*cp
	   Aj=Aj+sj*signz*cp
	   signz=-signz
	enddo
	A=A*c
	Aj=Aj*c
	if(m.eq.0) return
 	signz=dcmplx(0.d0,1.d0)**(n+v+1)
 	p=n+v
	do 20 i=1,qmax
	   p=p-2
	   signz=-signz
	   if(i.eq.1) then
              cp=dble(2*p+3)*fa(m,p+3)      
              cp=cp*ga0(ig+1)/dble((p+3)*(p+2))
              goto 21
           endif
           if(i.eq.qmax) then
              if(p.eq.0) goto 22
              nv2=p*(p+1)
              cp=dble(2*p+3)*fa(m,p+1)
              cp=-cp*ga0(ig+i+1)/dble(nv2)
              goto 21
           endif
 22	   c4=fa(m,p+2)
           cp=-dble((p+1)*(p+2))*fb(n,v,p+2)*ga0(ig+i)
	   cp=cp+dble((p+2)*(p+1))*fb(n,v,p+1)*ga0(ig+i+1)
 	   cp=cp*dble(2*p+3)/c4
 21	   sj=sja(p+1)
	   sy=sya(p+1)
 	   B=B+dcmplx(sj,sy)*signz*cp
 	   Bj=Bj+sj*signz*cp
 20	continue
	B=B*c
	Bj=Bj*c
	return
	end

c  subroutine tipitaud.f 
c  calculates pi(m,n) & tau(m,n) up to a specified nmax for all 
c  m=0,1,...n at a given x
c  pi(m,n) & tau(m,n) calculated are normalized by  
c         C_mn=[(2n+1)(n-m)!/n/(n+1)/(n+m)!]^(1/2)
c  Yu-lin Xu    12/2000

	subroutine tipitaud(nmax,x)
	include 'gmm01f.par'
	parameter (nmp0=(np+1)*(np+4)/2)
	implicit double precision (a-h,o-z)
	common/fnr/fnr(0:2*(np+2))
	common/pitau/pi(nmp0),tau(nmp0)	

	nt=(nmax+1)*(nmax+4)/2          ! calculates pi up to nmax+1
	if(nt.gt.nmp0.or.dabs(x).gt.1.d0) then
	   write(6,*) 'dimension or argument wrong in sub. tipitaud'
	   write(6,*) 'argument: ',x
	   stop
	endif
	sx=dsqrt(1.d0-x*x)
	pi(1)=0.d0                     ! pi(0,1)  pi(0,n)=0 when m=0
	pi(2)=dsqrt(.75d0)             ! pi(1,1)
	pi(3)=0.d0                     ! pi(0,2)
	t125=dsqrt(1.25d0) 
	pi(4)=t125*x                      ! pi(1,2)  
	pi(5)=t125*sx                     ! pi(2,2)
	imn=5
	do i=3,nmax+1
	   n=i
	   imn=imn+1
	   pi(imn)=0.d0                ! pi(0,n)=0
	   do 11 j=2,n
	      m=j-1
	      imn=imn+1
	      i1=(n-2)*(n+1)/2+m+1
	      if(m.eq.n-1) then
	         pi(imn)=fnr(n-1)*fnr(2*n+1)/fnr(n+1)*x*pi(i1)
	         goto 11
	      endif
	      i2=(n-3)*n/2+m+1
	      t=fnr(n)*fnr(2*n-3)
	      t=fnr(n-2)*fnr(n-m-1)*fnr(n+m-1)/t
	      pi(imn)=fnr(2*n-1)*x*pi(i1)-t*pi(i2)
	      t=fnr(n+1)*fnr(n-m)*fnr(n+m)
	      t=fnr(n-1)*fnr(2*n+1)/t
	      pi(imn)=t*pi(imn)
 11	   continue
	   imn=imn+1
	   i1=(n-2)*(n+1)/2+n
	   t=fnr(n-1)*fnr(n+1)
	   t=dsqrt(.5d0)*fnr(n)*fnr(2*n+1)/t
	   pi(imn)=t*sx*pi(i1)
	enddo
	tx=x*sx
	tau(1)=-dsqrt(1.5d0)*sx          ! tau(0,1)
	tau(2)=pi(2)*x                   ! tau(1,1)
	tau(3)=-dsqrt(7.5d0)*tx          ! tau(0,2)   
	tau(4)=t125*(2.d0*x*x-1.d0)      ! tau(1,2)
	tau(5)=t125*tx                   ! tau(2,2)
	imn=5
	do i=3,nmax
	   n=i
	   do 21 j=1,n+1
	      m=j-1
	      imn=imn+1
	      if(m.eq.0) then
	         i1=(n-2)*(n+1)/2+1
	         i2=(n-3)*n/2+1
	         t=fnr(2*n-3)
	         t=fnr(n-2)*fnr(n)/t
	         tau(imn)=fnr(2*n-1)*x*tau(i1)-t*tau(i2)
	         t=fnr(n-1)*fnr(n+1)
	         t=fnr(2*n+1)/t
	         tau(imn)=t*tau(imn)
	         goto 21
	      endif
	      i1=n*(n+3)/2+m+1
	      t=fnr(n)*fnr(2*n+3)
	      t=fnr(n+2)*fnr(2*n+1)*fnr(n-m+1)*fnr(n+m+1)/t
	      tau(imn)=t*pi(i1)-dble(n+1)*x*pi(imn)
	      tau(imn)=tau(imn)/dble(m)
 21	   continue
	enddo
	return
	end

c  subroutine mueller.f
c  returns the values of 16 elements of the 4x4 
c  Mueller matrix calculated from the known 2x2 
c  amplitude scattering matrix 
c  using Bohren & Huffman's formulas (p.65)
	subroutine mueller(s1,s2,s3,s4,s)
	double precision s(4,4),s1s,s2s,s3s,s4s
	complex*16 s1,s2,s3,s4
	complex*16 s2s3c,s1s4c,s2s4c,s1s3c,s1s2c
	complex*16 s3s4c,s2s1c,s4s3c,s2cs4,s3cs1
	s1s=cdabs(s1)**2
	s2s=cdabs(s2)**2
	s3s=cdabs(s3)**2
	s4s=cdabs(s4)**2
	s2s3c=s2*dconjg(s3)
	s1s4c=s1*dconjg(s4)
	s2s4c=s2*dconjg(s4)
	s1s3c=s1*dconjg(s3)
	s1s2c=s1*dconjg(s2)
	s3s4c=s3*dconjg(s4)
	s2s1c=s2*dconjg(s1)
	s4s3c=s4*dconjg(s3)
	s2cs4=dconjg(s2)*s4
	s3cs1=dconjg(s3)*s1
	s(1,1)=0.5d0*(s1s+s2s+s3s+s4s)
	s(1,2)=0.5d0*(s2s-s1s+s4s-s3s)
	s(1,3)=s2s3c+s1s4c
	s(1,4)=dimag(s2s3c-s1s4c)
	s(2,1)=0.5d0*(s2s-s1s-s4s+s3s)
	s(2,2)=0.5d0*(s2s+s1s-s4s-s3s)
	s(2,3)=s2s3c-s1s4c
	s(2,4)=dimag(s2s3c+s1s4c)
	s(3,1)=s2s4c+s1s3c
	s(3,2)=s2s4c-s1s3c
	s(3,3)=s1s2c+s3s4c
	s(3,4)=dimag(s2s1c+s4s3c)
	s(4,1)=dimag(s2cs4+s3cs1)
	s(4,2)=dimag(s2cs4-s3cs1)
	s(4,3)=dimag(s1s2c-s3s4c)
	s(4,4)=s1s2c-s3s4c
	return
	end

c  lnfacd.f  (double precision function)
c  returns ln(z!)  z>-1.0
c  based on Lanczos' method [see Xu, Journal of Computational 
c  Physics, v.139, 137-165 (1998)]
      double precision function lnfacd(z)
      integer i
      double precision z,a,b,cp,c0(11)
      data c0/0.16427423239836267d5, -0.48589401600331902d5,
     +        0.55557391003815523d5, -0.30964901015912058d5,
     +        0.87287202992571788d4, -0.11714474574532352d4,
     +        0.63103078123601037d2, -0.93060589791758878d0,
     +        0.13919002438227877d-2,-0.45006835613027859d-8,
     +        0.13069587914063262d-9/ 
      a=1.d0
      cp=2.5066282746310005d0
      b=z+10.5d0
      b=(z+0.5d0)*dlog(b)-b
      do i=1,11
        z=z+1.d0
        a=a+c0(i)/z
      enddo
      lnfacd=b+dlog(cp*a)
      return
      end
 
c  gxurcd0.f to compute Gaunt coefficients a(-m,n,m,v,p)
cu uses lnfacd.f to compute ln(z!)
	subroutine gxurcd0(m,n,v,qmax,na)
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
	integer v,qmax,p
	double precision lnfacd,cnv(np,np),ga0(ng0)
	common/cfnv/cnv
	common/g0/ga0
	fb(n,v,p)=dble(p-(n+v+1))*dble(p+(n+v+1))
     +           *dble(p-(n-v))*dble(p+(n-v))
     +           /(dble(2*p+1)*dble(2*p-1))	
	if(iabs(m).gt.n.or.iabs(m).gt.v) then
	   write(6,*) 'warning: |m|>n or v in gxurcd0'
	   qmax=-1	   
	   return
	endif
	qmax=min(n,v)
	nq=qmax+1
	if(n.le.v) then 
	   c1=cnv(n,v)
	else
	   c1=cnv(v,n)
	endif
	c1=c1-lnfacd(dble(n-m))-lnfacd(dble(v+m))
	ga0(na+1)=dexp(c1)
	if(qmax.lt.1) return	
 	p=n+v
	do 8 i=2,nq
	   p=p-2
	   if(m.eq.0) then
	      c1=fb(n,v,p+1)
	      c2=fb(n,v,p+2)
	      goto 2
	   endif
	   c1=fb(n,v,p+1)
	   c2=dble(4*m*m)+fb(n,v,p+2)+fb(n,v,p+3)
	   if(i.eq.2) goto 2
	   c3=-fb(n,v,p+4)
	   goto 4
  2	   ga0(na+i)=c2*ga0(na+i-1)/c1
	   goto 8
  4	   ga0(na+i)=(c2*ga0(na+i-1)+c3*ga0(na+i-2))/c1
  8     continue
	return
	end

	subroutine gid0(nmax,m,n,iv,id)
	nt=nmax*(nmax+1)*(2*nmax+1)/3+nmax*nmax
	ns=max(1,iabs(m))
 	nc0=nmax-iabs(m)
 	id=nc0*(nc0+1)*(2*nc0+1)/6
	if(m) 10,11,12
 10	id=id+(n-ns)*(nc0+1)+iv-ns+1  
	return
 11	id=id+(n-ns)*nmax+iv
	return
 12	id=id+(nmax-n)*(nc0+1)+nmax-iv
	id=nt-id
	return
	end

      subroutine rtr(anpt,nodrj,nodri,ekt,drot)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter(nmp=np*(np+2))
      parameter (nrc=4*np*(np+1)*(np+2)/3+np)      
      double precision drot(nrc)
      complex*16 anpt(2,nmp),ant(2,2*np),amt(2,-np:np),a,b,
     +           ekt(np),eks(-np:np),atr(2,np,nmp)
      common/tran/atr
c
      eks(0)=1.d0
      nmax=max(nodrj,nodri)
      do m=1,nmax
         eks(m)=ekt(m)
         eks(-m)=dconjg(eks(m))
      enddo
      irc=0
      do n=1,nodrj
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.d0
            amt(2,m)=0.d0
         enddo
         do k=-n,n
            kn=n1+k
            a=eks(k)*anpt(1,kn)
            b=eks(k)*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
         enddo
         do m=-n,n
            imn=n1+m
            anpt(1,imn)=amt(1,m)
            anpt(2,imn)=amt(2,m)
         enddo
      enddo
      mmax=min(nodrj,nodri)
      do m=-mmax,mmax
         n1=max(1,iabs(m))
         do n=n1,nodrj
            imn=n*(n+1)+m
            do ip=1,2
               ant(ip,n)=anpt(ip,imn)
            enddo
         enddo
         do n=n1,nodri
            imn=n*(n+1)+m
            a=0.d0
            b=0.d0
            do l=n1,nodrj
               ml=l*(l+1)+m
               a=a+atr(1,n,ml)*ant(1,l)
     1            +atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)
     1            +atr(2,n,ml)*ant(1,l)
            enddo
            anpt(1,imn) = a
            anpt(2,imn) = b
         enddo
      enddo
      in=1
      irc=0
      do n=1,nodri
         in=-in
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.d0
            amt(2,m)=0.d0
         enddo
         sik=-in
         do k=-n,n
	    sik=-sik
	    kn=n1+k
	    a=sik*anpt(1,kn)
            b=sik*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
         enddo
         sik=-in
         do m=-n,n
            sik=-sik
            imn=n1+m
            anpt(1,imn)=amt(1,m)*eks(-m)*sik
            anpt(2,imn)=amt(2,m)*eks(-m)*sik
         enddo
      enddo
      return
      end

      subroutine trv(anpt,nodrj,nodri)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter(nmp=np*(np+2))
      complex*16 anpt(2,nmp),ant(2,2*np),a,b,atr(2,np,nmp)
      common/tran/atr
      mmax=min(nodrj,nodri)
      do m=-mmax,mmax
         n1=max(1,iabs(m))
         do n=n1,nodrj
            imn=n*(n+1)+m
            do ip=1,2
               ant(ip,n)=anpt(ip,imn)
            enddo
         enddo
         do n=n1,nodri
            imn=n*(n+1)+m
            a=0.d0
            b=0.d0
            do l=n1,nodrj
               ml=l*(l+1)+m
               a=a+atr(1,n,ml)*ant(1,l)
     1            +atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)
     1            +atr(2,n,ml)*ant(1,l)
            enddo
            anpt(1,imn) = a
            anpt(2,imn) = b
         enddo
      enddo
      return
      end

c subroutine psiy.f to evaluate Psi_n(y)/y^2 needed for a sphere in  
c the calculation of its source function of photophoretic force, 
c where y is complex and Psi_n(y)=y j_n(y) with j_n(y) being the 
c Bessel function of the first kind

      SUBROUTINE psiy(FRAC,X,REFREL,NC,NPY,PY,DPY)

c FRAC ---- the fraction of the sphere radius ranging [-1,1]
c X ------- the size parameter of the sphere
c REFREL -- the complex refractive index of the sphere
c NC ------ the dimension of PY that needs to be calculated
c NPY ----- the maximum dimension of PY and DPY
c PY ------ the final output of this subroutine, Psi_n(y)/y^2
c DPY ----- the ratio of Psi'_n(y) to Psi_n(y)
c NC must be smaller than NPY

      INTEGER NC,NPY,NX,K,N
      COMPLEX*16 REFREL,Y,PY(NPY),DPY(NPY),PNY,SNY,Y2
      DOUBLE PRECISION FRAC,X,XN,CN
      DOUBLE PRECISION U,V,SU,CU,SHV,CHV,YM2,PY1R,PY1I

      Y=REFREL*X*FRAC
      Y2=Y*Y
      XN=CDABS(REFREL)*X
      NX=1.1D0*XN+10.D0
      IF(NC.GT.NPY) THEN
         WRITE(6,*) 'The dimension of PY too small'
         STOP
      ENDIF
      IF(NC.GT.NX) THEN
         WRITE(6,*) 'Sub. psiy: the dimension NC of PY too large'
         STOP
      ENDIF

C Down recursion for calculating ratios PY(NX) to PY(2)

      IF(CDABS(Y).LT.1.D-10) THEN
         PY(1)=1.D0/3.D0
         RETURN
      ENDIF   
      PNY=Y/DBLE(2*NX+3)
      DO 5 K=1,NX
         N=NX-K+1
         CN=DBLE(N)
         DPY(N)=(CN+1.0D0)/Y-PNY
         IF(N.EQ.1) GO TO 20
         PNY=Y/(2.D0*CN+1.D0-PNY*Y)
         PY(N)=PNY
    5 CONTINUE

c Calculating PY_1(y)=Psi_1(y)/y^2

   20 IF(CDABS(Y).LT.0.1D0) THEN
         SNY=1.D0/45360.D0-Y2/3991680.D0
         SNY=1.D0/840.D0-Y2*SNY
         SNY=1.D0/30.0d0-Y2*SNY
         SNY=1.D0/3.0d0-Y2*SNY
      ELSE
         U=DBLE(Y)
         V=DIMAG(Y)
         CU=DCOS(U)
         SU=DSIN(U)
         CHV=DCOSH(V)
         SHV=DSINH(V)
         YM2=CDABS(Y)*CDABS(Y)
         PY1R=U*SU*CHV+V*CU*SHV
         PY1R=PY1R/YM2-CU*CHV
         PY1I=V*SU*CHV-U*CU*SHV
         PY1I=-(PY1I/YM2-SU*SHV)
         SNY=DCMPLX(PY1R,PY1I)/Y2
      ENDIF    
      
c calculating PY_n from the ratios and PY_1

      PY(1)=SNY
      DO 10 N=2,NC
         PY(N)=PY(N)*PY(N-1)  
   10 CONTINUE
      RETURN
      END

      function plgndrd(l,mr,x)
      integer l,mr,m,index,i,ll
      double precision x,fact,pll,pmm,pmmp1,somx2
      double precision plgndrd,lnfacd
      m=mr
      index=0
      if(m.lt.0) then
         index=1
         m=-m
      endif
      if(dabs(x).gt.1.d0) then
	 write(6,*) 'n,m,x: ', l,(-2*index+1)*m,x
         pause 'bad arguments in plgndrd'
      end if
      if(m.gt.l) then
	 plgndrd=0.d0
         return
      end if
      pmm=1.d0
      if(m.gt.0) then
        somx2=dsqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndrd=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndrd=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndrd=pll
        endif
      endif
      plgndrd=-plgndrd
      if(m/2*2.eq.m) plgndrd=-plgndrd
      if(index.gt.0) then
         fact=lnfacd(dble(l-m))-lnfacd(dble(l+m))
         fact=dexp(fact)
         plgndrd=-plgndrd*fact
         if(m/2*2.eq.m) plgndrd=-plgndrd
      endif 
      return
      end
C
C  subroutine internd.f 
C  computing three-dimensional internal electric field distributions for 
C  each component spheres [see Xu et al. Physical Review E 60:2347 (1999)]
C  idphoto=1: in spherical coordinates
C  idphoto=2: in Cartesian coordinates 
C
	subroutine internd(nL,idphoto,nphoto,dphi,istart,iend,istep,
     +     indpol,idMie,x,nmax,ref,px,rsx,rsr,rsi,NXMAX,py0,py,dpy,
     +     as,bs)
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
	parameter (nmp=np*(np+2),nphimax=1800)
	parameter (nmp0=(np+1)*(np+4)/2)
	integer nmax(nLp)
	double precision x(nLp),phi(nphimax),tsf(nphimax),cs(nphimax),
     +     px(np,nLp),rsr(np,nLp),rsi(np,nLp),rsx(np,nLp),lnfacd
	complex*16 A,B,cmz,A2,B2,A0,B0,ci,cphm,Aj(2*nphimax),
     +     Bj(2*nphimax),Aj2(2*nphimax),cph(2*nphimax),ref(nLp),
     +     py0(NXMAX),py(NXMAX),dpy(NXMAX),as(nLp,nmp),bs(nLp,nmp)
	CHARACTER cnr1*1,cnr2*2,cnr3*3,cnr4*4,fout*7
	common/pitau/pi(nmp0),tau(nmp0)
       
        if(idphoto.lt.1) return
        if(nphoto.lt.1) return
	pih=dacos(0.d0)
	pd=pih/90.0d0
	ci=dcmplx(0.d0,1.d0)
	if(dphi.lt.0.1d0.and.dphi.ge.1.d-5) then
	   dphi=0.1d0
	   write(6,*) 'Warning: dphi on the last line of the input'
	   write(6,*) '     file gmm01f.in has been changed to 0.1'
	   write(6,*) 
     +       'If you really need a smaller dphi, please increase'
	   write(6,*) 
     +       'the parameter nphimax in the subroutine internd.f'
	   write(6,*) ' accordingly, recompile, and run again'
	endif
	if(nphoto/2*2.ne.nphoto) then
	   nphoto=nphoto-1
	   write(6,*) 
     +     'Warning: nphoto (on the last line of gmm01f.in) must be'
	   write(6,*) 
     +     '  an even number and has been changed to', nphoto
	endif
	if(dphi.lt.1.d-5) then
	  nphi=1
	else
	  nphi=180.0d0/dphi
	endif
	do i=1,nphi
	   phi(i)=dble(i-1)*dphi
	   temp=pd*phi(i)
	   t1=dcos(temp)
	   cs(i)=t1
	   t2=dsin(temp)
	   cph(i)=dcmplx(t1,t2)
	enddo 	
	if(idphoto.gt.1) goto 772

	sz=360.d0/dble(nphoto)  	
        do j=istart,iend,istep
           if(j.lt.10) then
	      write(cnr1,'(i1)') j
	      cnr4='000'//cnr1
	   else
	      if(j.lt.100) then
	         write(cnr2,'(i2)') j
	         cnr4='00'//cnr2
	      else
	         if(j.lt.1000) then
	            write(cnr3,'(i3)') j
	            cnr4='0'//cnr3
	         else
	            write(cnr4,'(i4)') j
	         endif
	      endif
	   endif
           if(indpol.lt.1) then
	      fout='xsf'//cnr4
	   else
	      fout='ysf'//cnr4
	   endif
	   open(16,file=fout,status='unknown')
	   write(16,'(2x,a7)') fout
	   write(16,'(a18,f6.3,a20,f8.5,a1,f7.5,a1)') 
     +        '  Size parameter: ',x(j),'  Refractive index: ',
     +            dble(ref(j)),'+',dabs(dimag(ref(j))),'i'
           if(indpol.lt.1) then
	      write(16,*) 'incident wave x-polarized'
	   else
	      write(16,*) 'incident wave y-polarized'
	   endif
	   write(16,*) ''
	   write(16,'(a9,i5)') '  nphoto:',nphoto
           write(16,*)
     +     'lists (nphoto+1)*(nphoto+1) points for each phi in'
           write(16,*)
     +     'spherical coordinates r/a [0,1] and theta [0,360]' 
           write(16,*)
     +     '(a is sphere radius and theta is counted counterclockwise' 
           write(16,*)
     +     'from the positive z axis along the incident direction)' 
	   write(16,*) ''
	   write(16,'(2a9,4x,a29,f5.1,a1)')
     +     'r/a','theta(deg)','|E|^2/E0^2 at phi(deg) (dphi:',dphi,')'
	   write(16,'(22x,1800(f5.1,9x))') (phi(ip),ip=1,nphi)
           d=1.d0
           call psiy(d,x(j),ref(j),nmax(j),NXMAX,py0,dpy)
           cmz=ref(j)*x(j)
           do ii=1,nmax(j)
              py0(ii)=py0(ii)*cmz*cmz
           enddo
           do ii=1,nphoto+1
	      d=dble(ii-1)/dble(nphoto)
	      d=dsqrt(d)
	      z0=-sz
              do 78 jj=1,nphoto+1
                 z0=z0+sz
                 st=z0*pd
                 xt=dcos(st)
                 st=dsin(st)
                 A2=d*x(j)*ref(j)
                 call psiy(d,x(j),ref(j),nmax(j),NXMAX,py,dpy)
                 call tipitaud(nmax(j),xt)
	         do ip=1,nphi	            
                    Aj(ip)=0.d0
                    Bj(ip)=0.d0
                    Aj2(ip)=0.d0
                    tsf(ip)=0.d0
                 enddo
                 do 87 n=1,nmax(j)
                    if(cdabs(A2).lt.1.d-10.and.n.gt.1) goto 87
                    A0=px(n,j)*py0(n)
                    A0=dcmplx(0.d0,1.d0)*ref(j)/A0
                    A=ref(j)*rsx(n,j)-dcmplx(rsr(n,j),-rsi(n,j))
                    A=A0/A
                    B=rsx(n,j)-ref(j)*dcmplx(rsr(n,j),-rsi(n,j))
                    B=A0/B
                    A0=ci**n*py(n)
                    pn=dble(2*n+1)/dble(n*(n+1))
                    pn=dsqrt(pn)
                    do 88 m=-n,n
                       m0=iabs(m)
                       if(nL.eq.1.or.idMie.eq.1) then
                          if(m0.ne.1) goto 88
                       endif
	               i=n*n+n+m    !sequence # for scattering coeff.
                       i0=(n-1)*(n+2)/2+m0+1 !sequence # for pi & tau
                       if(m.lt.0) then
                          p=(-1.d0)**m0
                          t=p*tau(i0)
                          p=-p*pi(i0)
                       else
                          t=tau(i0)
                          p=pi(i0)
                       endif
                       if(cdabs(A2).lt.1.d-10) then
                         B2=-ci*A*as(j,i)*t*2.d0/3.d0
                         B0=A*as(j,i)*p*2.d0/3.d0
                         do ip=1,nphi
                            cphm=cph(ip)**m
                            temp=(-1)**m
                            if(st.gt.0.d0.and.cs(ip).lt.0.d0) 
     +                         cphm=cphm*temp
                            if(st.lt.0.d0.and.cs(ip).gt.0.d0)
     +                         cphm=cphm*temp
                            Aj(ip)=Aj(ip)+A0*B2*cphm
                            Bj(ip)=Bj(ip)+A0*B0*cphm
                         enddo
                         goto 98
                       endif
                       B2=-ci*A*dpy(n)*as(j,i)*t
                       B2=B2+B*bs(j,i)*p
                       B0=ci*B*bs(j,i)*t 
                       B0=B0+A*as(j,i)*dpy(n)*p
                       do ip=1,nphi
                          cphm=cph(ip)**m
                          temp=(-1)**m
                          if(st.gt.0.d0.and.cs(ip).lt.0.d0) 
     +                       cphm=cphm*temp
                          if(st.lt.0.d0.and.cs(ip).gt.0.d0) 
     +                       cphm=cphm*temp
                          Aj(ip)=Aj(ip)+A0*B2*A2*cphm
                          Bj(ip)=Bj(ip)+A0*B0*A2*cphm
                       enddo
  98	               B2=-ci*dble(n*(n+1))*A*as(j,i)
	               t=plgndrd(n,m,xt)
	               rn=lnfacd(dble(n-m))-lnfacd(dble(n+m))
                       rn=0.5d0*rn
                       rn=dexp(rn)
                       t=t*pn*rn
	               do ip=1,nphi
	                  cphm=cph(ip)**m
                          temp=(-1)**m
                          if(st.gt.0.d0.and.cs(ip).lt.0.d0) 
     +                       cphm=cphm*temp
                          if(st.lt.0.d0.and.cs(ip).gt.0.d0) 
     +                       cphm=cphm*temp
	                  Aj2(ip)=Aj2(ip)+A0*B2*t*cphm
	               enddo             
  88	            continue
  87	         continue
                 do ip=1,nphi
	            t=cdabs(Aj(ip))**2+cdabs(Bj(ip))**2
	            t=t+cdabs(Aj2(ip))**2
	            tsf(ip)=t
	         enddo
	         write(16,'(2f10.4,1800e14.6)') 
     +                 d,z0,(tsf(ip),ip=1,nphi)
  78	      continue
	   enddo
	enddo
	close(16)
	return

 772	cz=2.d0/dble(nphoto)		
        do j=istart,iend,istep
           if(j.lt.10) then
	      write(cnr1,'(i1)') j
	      cnr4='000'//cnr1
	   else
	      if(j.lt.100) then
	         write(cnr2,'(i2)') j
	         cnr4='00'//cnr2
	      else
	         if(j.lt.1000) then
	            write(cnr3,'(i3)') j
	            cnr4='0'//cnr3
	         else
	            write(cnr4,'(i4)') j
	         endif
	      endif
	   endif
           if(indpol.lt.1) then
	      fout='xsf'//cnr4
	   else
	      fout='ysf'//cnr4
	   endif
	   open(16,file=fout,status='unknown')
	   write(16,'(2x,a7)') fout
	   write(16,'(a18,f6.3,a20,f8.5,a1,f7.5,a1)') 
     +        '  Size parameter: ',x(j),'  Refractive index: ',
     +            dble(ref(j)),'+',dabs(dimag(ref(j))),'i'
           if(indpol.lt.1) then
	      write(16,*) 'incident wave x-polarized'
	   else
	      write(16,*) 'incident wave y-polarized'
	   endif
           write(16,*) ''
           write(16,'(a9,i5)') '  nphoto:',nphoto
           write(16,*)
     +     'lists (nphoto+1)*(nphoto+1) points for each phi'
           write(16,*)
     +     'in Cartesian coordinates x0/a [-1,1] and z/a [-1,1]' 
           write(16,*)
     +     '(a is sphere radius, the positive z axis is along the'
           write(16,*)
     +     'incident direction, and x0 always has the same sign as x)'
           write(16,*) ''
           write(16,'(2a9,4x,a29,f5.1,a1)')
     +        'x0/a','z/a','|E|^2/E0^2 at phi(deg) (dphi:',dphi,')'
	   write(16,'(22x,1800(f5.1,9x))') (phi(ip),ip=1,nphi)
           d=1.d0
           call psiy(d,x(j),ref(j),nmax(j),NXMAX,py0,dpy)
           cmz=ref(j)*x(j)
           do ii=1,nmax(j)
              py0(ii)=py0(ii)*cmz*cmz
           enddo
           do ii=1,nphoto+1
	      x0=-1.d0+dble(ii-1)*cz
              do 782 jj=1,nphoto+1
                 z0=-1.d0+dble(jj-1)*cz
                 d=dsqrt(z0*z0+x0*x0)
                 if(d.gt.0.99999d0) then
                    do ip=1,nphi
                       tsf(ip)=0.d0
                    enddo
                    goto 862
                 endif
                 if(d.lt.1.d-10) then
                    xt=1.d0
                 else
                    xt=z0/d
                 endif
                 A2=d*x(j)*ref(j)
                 call psiy(d,x(j),ref(j),nmax(j),NXMAX,py,dpy)
                 call tipitaud(nmax(j),xt)
                 do ip=1,nphi	            
                    Aj(ip)=0.d0
                    Bj(ip)=0.d0
                    Aj2(ip)=0.d0
                    tsf(ip)=0.d0
                 enddo
                 do 872 n=1,nmax(j)
                    if(cdabs(A2).lt.1.d-10.and.n.gt.1) goto 872
                    A0=px(n,j)*py0(n)
                    A0=dcmplx(0.d0,1.d0)*ref(j)/A0
                    A=ref(j)*rsx(n,j)-dcmplx(rsr(n,j),-rsi(n,j))
                    A=A0/A
                    B=rsx(n,j)-ref(j)*dcmplx(rsr(n,j),-rsi(n,j))
                    B=A0/B
                    A0=ci**n*py(n)
                    pn=dble(2*n+1)/dble(n*(n+1))
                    pn=dsqrt(pn)
                    do 882 m=-n,n
                       m0=iabs(m)
                       if(nL.eq.1.or.idMie.eq.1) then
                          if(m0.ne.1) goto 882
                       endif
	               i=n*n+n+m    !sequence # for scattering coeff.
                       i0=(n-1)*(n+2)/2+m0+1 !sequence # for pi & tau
                       if(m.lt.0) then
                          p=(-1.d0)**m0
                          t=p*tau(i0)
                          p=-p*pi(i0)
                       else
                          t=tau(i0)
                          p=pi(i0)
                       endif
                       if(cdabs(A2).lt.1.d-10) then
                         B2=-ci*A*as(j,i)*t*2.d0/3.d0
                         B0=A*as(j,i)*p*2.d0/3.d0
                         do ip=1,nphi
                            cphm=cph(ip)**m
                            temp=(-1)**m
                            if(x0.gt.0.d0.and.cs(ip).lt.0.d0) 
     +                         cphm=cphm*temp
                            if(x0.lt.0.d0.and.cs(ip).gt.0.d0) 
     +                         cphm=cphm*temp
                            Aj(ip)=Aj(ip)+A0*B2*cphm
                            Bj(ip)=Bj(ip)+A0*B0*cphm
                         enddo
                         goto 982
                       endif
                       B2=-ci*A*dpy(n)*as(j,i)*t
                       B2=B2+B*bs(j,i)*p
                       B0=dcmplx(0.d0,1.d0)*B*bs(j,i)*t 
                       B0=B0+A*as(j,i)*dpy(n)*p 
                       do ip=1,nphi
                          cphm=cph(ip)**m
                          temp=(-1)**m
                          if(x0.gt.0.d0.and.cs(ip).lt.0.d0) 
     +                       cphm=cphm*temp
                          if(x0.lt.0.d0.and.cs(ip).gt.0.d0) 
     +                       cphm=cphm*temp
                          Aj(ip)=Aj(ip)+A0*B2*A2*cphm
                          Bj(ip)=Bj(ip)+A0*B0*A2*cphm
                       enddo
 982	               B2=-ci*dble(n*(n+1))*A*as(j,i)
	               t=plgndrd(n,m,xt)
	               rn=lnfacd(dble(n-m))-lnfacd(dble(n+m))
                       rn=0.5d0*rn
                       rn=dexp(rn)
                       t=t*pn*rn
	               do ip=1,nphi
	                  cphm=cph(ip)**m
                          temp=(-1)**m
                          if(x0.gt.0.d0.and.cs(ip).lt.0.d0) 
     +                       cphm=cphm*temp
                          if(x0.lt.0.d0.and.cs(ip).gt.0.d0) 
     +                       cphm=cphm*temp
	                  Aj2(ip)=Aj2(ip)+A0*B2*t*cphm
	               enddo            
 882	            continue
 872	         continue
                 do ip=1,nphi
	            t=cdabs(Aj(ip))**2+cdabs(Bj(ip))**2
	            t=t+cdabs(Aj2(ip))**2
	            tsf(ip)=t
	         enddo
 862	         write(16,'(2f10.4,1800e14.6)') 
     +                x0,z0,(tsf(ip),ip=1,nphi)
 782	      continue
	   enddo
	enddo
	close(16)
	return
	end

      subroutine transs(k,nL,r0,nmax0,nmax,uvmax,fint,ek,
     +             besj,besy,drot,as,bs,as1,bs1,ind,icase)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2),nmp0=(np+1)*(np+4)/2)
      parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)	
      parameter (nrc=4*np*(np+1)*(np+2)/3+np)
      parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)	
      integer u,v,nmax(nLp),uvmax(nLp),ind(nLp)
      double precision r0(6,nLp),drot(nrc),k,
     +          besj(0:2*np+1),besy(0:2*np+1)
      complex*16 atr(2,np,nmp),
     +   ek(np),as(nLp,nmp),bs(nLp,nmp),as1(nLp,nmp),
     +   bs1(nLp,nmp),at1(2,nmp),t1,t2,ephi
      common/tran/atr
      common/rot/bcof(0:np+2),dc(-np:np,0:nmp)
      common/fnr/fnr(0:2*(np+2))
      common/pitau/pi(nmp0),tau(nmp0)
      common/ig0/iga0(ni0)
      common/g0/ga0(ng0)
      common/cofmnv0/cof0(ni0)
      common/crot/cofsr(nmp)
      do i=1,nL
         do imn=1,uvmax(i)
            as1(i,imn)=dcmplx(0.d0,0.d0)
            bs1(i,imn)=dcmplx(0.d0,0.d0)
         enddo
      enddo
      do 11 i=1,nL-1
         if(ind(i).gt.0) goto 11
         do 10 j=i+1,nL
            if(ind(j).gt.0) goto 10
            x0=r0(1,i)-r0(1,j)
            y0=r0(2,i)-r0(2,j)
            z0=r0(3,i)-r0(3,j)
            call carsphd(x0,y0,z0,d,xt,sphi,cphi)
            temp=(r0(4,i)+r0(4,j))/d
            if(temp.le.1.d0) goto 1
            temp1=r0(4,i)
            if(r0(4,j).lt.r0(4,i)) temp1=r0(4,j)
            temp1=(r0(4,i)+r0(4,j)-d)/temp1
            if(temp1.gt.0.001d0) then
               write(6,'(/)')
               write(6,*) 'Fatal error in aggregate configuration:'
               write(6,*) '        SPHERES OVERLAPPED' 
               write(6,'(e12.4)') temp1
               write(6,'(i6,4e12.4)') i,(r0(ilist,i),ilist=1,4)
               write(6,'(i6,4e12.4)') j,(r0(ilist,j),ilist=1,4)
               write(6,*) 'This code calculates aggregate-scattering'
               write(6,*) 'of non-intersecting spheres. The ratio of'
               write(6,*) '(r1+r2)/d must not exceed 1, where r1+r2'
               write(6,*) 'is the sum of the radii of any pair of'
               write(6,*) 'spheres and d is the separation distance' 
               write(6,*) 'between the two sphere centers. This ratio'
               write(6,*) 'for the pair of spheres listed is greater'
               write(6,*) 'than 1 and the overlapped portion exceeds'
               write(6,*) '0.001 of the (smaller) sphere-radius as '
               write(6,*) 'shown above. The sequence numbers and the'
               write(6,*) 'Cartesian coordinates of the two sphere' 
               write(6,*) 'centers, the radii of the two spheres are'
               write(6,*) 'also shown. Please correct the input file' 
               write(6,*) 'for the aggregate configuration and then'
               write(6,*) 'run again.'
               stop
            endif
  1         if(temp.le.fint) goto 10
            ephi=dcmplx(cphi,sphi)
            nlarge=max(nmax(i),nmax(j))
            do m=1,nlarge
               ek(m)=ephi**m
            enddo
            xd=k*d
            nbes=2*nlarge+1
            call besseljd(nbes,xd,besj)
            call besselyd(nbes,xd,besy)
            call rotcoef(xt,nlarge)
            irc=0
            do n=1,nlarge
               n1=n*(n+1)
               do u=-n,n
                  do m=-n,n
                     imn=n1+m
                     irc=irc+1
                     drot(irc)=dc(u,imn)
                  enddo
               enddo
            enddo
            nsmall=min(nmax(i),nmax(j))
C  the formulation used here for the calculation of vector translation 
C  coefficients are from Cruzan, Q. Appl. Math. 20, 33 (1962) and  
C  Xu, J. Comput. Phys. 139, 137 (1998)
            do m=-nsmall,nsmall
               n1=max(1,iabs(m))
               do n=n1,nlarge
                  do v=n1,nlarge
                     iuv=v*(v+1)+m
                     if(icase.lt.2) then
                        call cofxuds0(nmax0,m,n,v,besj,besy,
     +                      atr(1,n,iuv),atr(2,n,iuv),t1,t2)
                     else 
                        call cofxuds0(nmax0,m,n,v,besj,besy,
     +                      t1,t2,atr(1,n,iuv),atr(2,n,iuv))
                     endif
                     if(x0.eq.0.d0.and.y0.eq.0.d0) then
                        if(z0.lt.0.d0) then
                           sic=dble((-1)**(n+v))
                           atr(1,n,iuv)=sic*atr(1,n,iuv)
                           atr(2,n,iuv)=-sic*atr(2,n,iuv)
                        endif
                     endif
                  enddo
               enddo
            enddo
            do iuv=1,uvmax(j)
               at1(1,iuv)=as(j,iuv)
               at1(2,iuv)=bs(j,iuv)
            enddo
            If(nmax(j).lt.nlarge) then
               do n=nmax(j)+1,nlarge
                  do m=-n,n
                     iuv=n*(n+1)+m
                     at1(1,iuv)=dcmplx(0.d0,0.d0)
                     at1(2,iuv)=dcmplx(0.d0,0.d0)
                  enddo
               enddo
            endif              
            if(x0.eq.0.d0.and.y0.eq.0.d0) then
               call trv(at1,nmax(j),nmax(i))
            else 
               call rtr(at1,nmax(j),nmax(i),ek,drot)
            endif
            do imn=1,uvmax(i)
               as1(i,imn)=as1(i,imn)+at1(1,imn)
               bs1(i,imn)=bs1(i,imn)+at1(2,imn)
            enddo
            do m=-nsmall,nsmall
               n1=max(1,iabs(m))
               do n=n1,nlarge
                  do v=n1,nlarge
                     iuv=v*(v+1)+m
                     sic=dble((-1)**(n+v))
                     atr(1,n,iuv)=sic*atr(1,n,iuv)
                     atr(2,n,iuv)=-sic*atr(2,n,iuv)
                  enddo
               enddo
            enddo
            do iuv=1,uvmax(i)
               at1(1,iuv)=as(i,iuv)
               at1(2,iuv)=bs(i,iuv)
            enddo
            If(nmax(i).lt.nlarge) then
               do n=nmax(i)+1,nlarge
                  do m=-n,n
                     iuv=n*(n+1)+m
                     at1(1,iuv)=dcmplx(0.d0,0.d0)
                     at1(2,iuv)=dcmplx(0.d0,0.d0)
                  enddo
               enddo
            endif 
            if(x0.eq.0.d0.and.y0.eq.0.d0) then
               call trv(at1,nmax(i),nmax(j))
            else
               call rtr(at1,nmax(i),nmax(j),ek,drot)
            endif
            do imn=1,uvmax(j)
               as1(j,imn)=as1(j,imn)+at1(1,imn)
               bs1(j,imn)=bs1(j,imn)+at1(2,imn)
            enddo
 10      continue
 11   continue
      return
      end
