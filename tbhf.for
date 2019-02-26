C*************************************************
C     This routine computes the eigenvalues (EN) and
C     eigenvectors (VR and VI) at K-point vk given
C     LCAO matrix matrix elements the tbparameters.dat file.
C     IS=0 (EIGENVALUE ONLY) , 1(EIGENVECTORS AND EIGENVALUES)
C     Order of parameters:
C       E(s,a), E(p,a), E(s,c), E(p,c), V(s,s), V(x,x)
C       V(x,y), V(sa,pc), V(sc,pa), E(s*,a), E(s*,c)
C       V(s*a,pc), V(pa,s*c)
      SUBROUTINE TBHF(PAR,VK,EN,VR,VI,IS,Hr,Hi)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HR(10,10),HI(10,10),EN(10),FV1(10),FV2(10),FM1(2,10)
     &,PAR(13),C(3),S(3),VK(3)
      DIMENSION VR(10,10),VI(10,10)
      COMPLEX*16 H(10,10),G0,G1,G2,G3
      COMPLEX*16 CG0,CG1,CG2,CG3
      COMPLEX*16 CQ1,CQ0
      CQ1=CMPLX(1.0,0.0)
      CQ0=CMPLX(0.0,0.0)
      PI=dacos(-1.d0)
      IG=1
      NSTATE=10
      MMX=10
      DO J=1,3
         C(J)=COS(VK(J)*PI/2.)
         S(J)=SIN(VK(J)*PI/2.)
      END DO

      G0=CMPLX(C(1)*C(2)*C(3),-S(1)*S(2)*S(3))
      G1=CMPLX(-C(1)*S(2)*S(3),S(1)*C(2)*C(3))
      G2=CMPLX(-S(1)*C(2)*S(3),C(1)*S(2)*C(3))
      G3=CMPLX(-S(1)*S(2)*C(3),C(1)*C(2)*S(3))
      CG0=CONJG(G0)
      CG1=CONJG(G1)
      CG2=CONJG(G2)
      CG3=CONJG(G3)
      DO I=1,NSTATE
         DO J=1,I
            H(I,J)=CQ0
         END DO
      END DO

C     DIAGONAL TERMS
      H(1,1)=CQ1*PAR(1)
      H(2,2)=CQ1*PAR(3)
      H(3,3)=CQ1*PAR(2)
      H(4,4)=CQ1*PAR(2)
      H(5,5)=CQ1*PAR(2)
      H(6,6)=CQ1*PAR(4)
      H(7,7)=CQ1*PAR(4)
      H(8,8)=CQ1*PAR(4)
C     OFF-DIAGONAL TERMS
      H(2,1)=PAR(5)*CG0
      H(3,2)=-PAR(9)*G1
      H(4,2)=-PAR(9)*G2
      H(5,2)=-PAR(9)*G3
      H(6,1)=PAR(8)*CG1
      H(7,1)=PAR(8)*CG2
      H(8,1)=PAR(8)*CG3
C
      H(6,3)=PAR(6)*CG0
      H(7,4)=PAR(6)*CG0
      H(8,5)=PAR(6)*CG0
      H(6,4)=PAR(7)*CG3
      H(6,5)=PAR(7)*CG2
      H(7,5)=PAR(7)*CG1
      H(7,3)=PAR(7)*CG3
      H(8,3)=PAR(7)*CG2
      H(8,4)=PAR(7)*CG1
C     EXCITED S STATES
      H(9,9)=PAR(10)
      H(10,10)=PAR(11)
      H(9,6)=PAR(12)*G1
      H(9,7)=PAR(12)*G2
      H(9,8)=PAR(12)*G3
      H(10,3)=-PAR(13)*CG1
      H(10,4)=-PAR(13)*CG2
      H(10,5)=-PAR(13)*CG3
      DO I=1,NSTATE
         DO J=1,I
            HR(I,J)=REAL(H(I,J))
            HI(I,J)=DIMAG(H(I,J))
         END DO
      END DO

      CALL CH(MMX,NSTATE,HR,HI,EN,IS,VR,VI,FV1,FV2,FM1,IERR)
c      IF(IERR.NE.0) then
c      WRITE(*,*) '  EISPACK ERROR, IER=',IER
c      stop
c      end if
      RETURN
      END SUBROUTINE TBHF

