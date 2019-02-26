      program electronic_band
	  implicit none
	  integer, parameter :: na=2,qqqqq=100,nk=qqqqq**3,nb=5*na
	  real, DIMENSION(3) :: b1 = (/ -1.,  1.  ,  1. /)
	  real, DIMENSION(3) :: b2 = (/  1., -1.  ,  1. /)
	  real, DIMENSION(3) :: b3 = (/  1.,  1.  , -1. /)
	  real*8 PI,gamma_
	  real*8 qx,qy,qz,up,ur,us,omega,summation1 !,ltcst
	  integer is,iiq,p,r,sss,ib,iq,j !M,mtl,qqqqq,
	  integer nkband,nkn,i,ik
	  REAL*8 PAR(13),EB(nk,nb),eb_band(801,10),VK(3) !,C(3),S(3)
      REAL*8 VR(10,10),VI(10,10),EN(10) ,NKL(801,3)
	  real*8 HR(10,10),HI(10,10)
	  OPEN(unit=20,file='tbparameters.dat',status='old')
	  open(unit=21,file='nkl_801.dat',status='old')
	  OPEN(unit=22,file='eb.dat',status='replace')
	  open(unit=23,file='DOS.dat',status='replace')
c	  open(unit=24,file='HamiltonianR.dat',status='replace')
c	  open(unit=25,file='HamiltonianI.dat',status='replace')

	  PI = dacos(-1.d0)
	  gamma_= 0.001
c	  mtl=2 !Read the parameters of different sets for every materials
	        !Each set consists of 13 parameters.

c	  DO M=1,MTL
	  READ(20,*)PAR
c	  END DO

      WRITE(*,*) 'PAR =  ', PAR
c     ltcst=PAR(14)    ! lattice constant

	  IS=0
	  QX=0.
	  QY=0.
	  QZ=0.

c	  qqqqq=qqqq

	  IIQ=0
	  do p=1,qqqqq
	     do r=1,qqqqq
	        do sss=1,qqqqq
	           IIQ=IIQ+1
		       up=(2.*float(p)-float(qqqqq)-1.)/(2.*float(qqqqq))
		       ur=(2.*float(r)-float(qqqqq)-1.)/(2.*float(qqqqq))
	           us=(2.*float(sss)-float(qqqqq)-1.)/(2.*float(qqqqq))

	           QX=up * b1(1) + ur * b2(1) + us * b3(1)
	           QY=up * b1(2) + ur * b2(2) + us * b3(2)
	           QZ=up * b1(3) + ur * b2(3) + us * b3(3)

	           VK(1)=QX !* (2*PI/ltcst)
	           VK(2)=QY !* (2*PI/ltcst)
	           VK(3)=QZ !* (2*PI/ltcst)

               CALL TBHF(PAR,VK,EN,VR,VI,IS,Hr,Hi)
	           DO IB=1,10
	              EB(iiq,IB)=EN(IB)
	           END DO
	        end do ! for up
	     end do ! for ur
	  end do ! for us

	  WRITE(*,*) 'iiq= ,nk= ', IIQ,nk
!!!!!!!!!!! calculate density of states !!!!!!!!!!!!!

	  do omega= -20.0, 20.0, 0.1
	     summation1=0.0
	     do IQ=1,nk
	        do J=1,10
	           summation1= summation1 +
     &         gamma_/(2.*PI)/( (omega-EB(iq,j))**2 + gamma_**2/4. )
	        end do
	     end do
	     write (23,"(2F20.6)") omega, summation1
	  end do
!!!!!!!!!!! calculate band structure!!!!!!!!!!!!!!!!!	
	  nkband=801

	  do nkn=1,nkband,1
	     read(21,*) ( NKL(nkn,i), i=1,3 )
	  end do

	  DO IK=1,nkband
	     DO I=1,3
	        VK(I)=NKL(IK,I)/20.
	     end do
            CALL TBHF(PAR,VK,EN,VR,VI,IS,Hr,Hi)
	     DO IB=1,10
	        EB_band(IK,IB)=EN(IB)
	     end do
	   
c	     write(24,"(4F16.5)") (vk(i),i=1,3)
c	     write(25,"(4F16.5)") (vk(i),i=1,3)
	   	
c	     do i=1,10
c 	        write(24,7001) (Hr(i,j),j=1,10)
c	     end do

c	     do i=1,10
c	        write(25,7001) (Hi(i,j),j=1,10)
c	     end do

c	     write(24,*)
c	     write(25,*)
	  end do
 	
	  do ik=1,nkband
	     write(22,7001) float( ik-1 ), (EB_band(IK,IB),Ib=1,10)
	  end do
 7001 FORMAT(12F12.5)
	  END program electronic_band
