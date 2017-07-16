c-----------------------------------------------------------------------------c
      PROGRAM RbMain

c     TO RUN YOU NEED:
c     HyperFine.o
c     angmom.o
c     lapack.a
c     blas.a    ! from lapack or another source

      IMPLICIT REAL*8 (a-b,d,p,q,k,w)
      IMPLICIT INTEGER (n)
      IMPLICIT COMPLEX*16 (c)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/params/ipar

      qmin=1.d-6
      qmax=50d0
      Nq=512+1
      Nphi=30+1

      angle=20d0

      Nk=10
      kmax=5d-1

      Bmin=163d0
      Bmax=165d0
      NB=1

c     Here is where you choose the parameters with which to run the calculation
c      ipar=1                    ! Jose singlet and triplet
c      ipar=2                    ! |b> channel params
      ipar=3                    ! From McNeil

      PRINT*
      IF(ipar.EQ.1) THEN
         PRINT*,'Using Jose Singlet and Triplet Parameters'
      END IF
      IF(ipar.EQ.2) THEN
         PRINT*,'Using Parameters for |b> Channel'
      END IF
      IF(ipar.EQ.3) THEN
         PRINT*,'Using Parameters from McNeil'
      END IF
      PRINT*

      CALL Setup(qmin,qmax,angle)
      CALL Print_Threebody(kmax,Nk,Bmin,Bmax,NB)

      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------------c
      SUBROUTINE Setup(qmin,qmax,angle)

      IMPLICIT REAL*8 (a-h,j-m,o-z)
      IMPLICIT INTEGER (i,n)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/dat/pi,hbarc,mRb,Ehf

c     Here are the constants for the problem
      pi=2d0*asin(1d0)
      amu=931.494028d6
      hbarc=197.326968d0
      clight=2.99792458d17
      mu_Bohr=5.788381804d-5
      a_Bohr=0.052917720859d0
      
      mRb=84.911789738d0*amu/hbarc
      mu_Elec=1.0011596521859d0*mu_Bohr
      hPlanck=2*pi*hbarc/clight
      Ehf=3.036d9*hPlanck	! eV
      Bhf=1d4*Ehf/mu_elec	! Gauss
      Hartree=alpha*hbarc/a_Bohr

c     q-grid
      dt=LOG(qmax/qmin)/(Nq-1)  ! This is the different iteration spacing

      DO i=1,Nq
         qlist(i)=qmin*EXP((i-1)*dt) ! These are the q's
      END DO

c     Angle grid
      phi0=Pi/angle              ! Angle to lower path by
      dphi=phi0/(Nphi-1)

      DO i=1,Nphi
         phi(i)=phi0-(i-1)*dphi
      END DO

c     Now I calculate the weights for the integration
c     I use the Bode 8=point type defined in McNeil's code
      CALL Weights

      RETURN
      END
c-----------------------------------------------------------------------------c

c-----------------------------------------------------------------------------c
      SUBROUTINE Print_Threebody(kmax,Nk,Bmin,Bmax,NB)

      IMPLICIT REAL*8 (a-b,d-h,k,m,o-z)
      IMPLICIT COMPLEX*16 (c-c)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/dat/Pi,hbarc,mRb,Ehf
      COMMON/const/Wfnorm,g1,g2,b1,b2
      DATA cI/(0d0,1d0)/

      dk=kmax/Nk
      dB=(Bmax-Bmin)/NB

      DO j=1,NB

         B=(j-1)*dB+Bmin
         ik_last=0

         PRINT*,'-----------------------------------------------------'
         PRINT*,'  Magnetic Field = ',B,' gauss'
         PRINT*,'-----------------------------------------------------'

         CALL TwoBodyHyperfine(B)

         PRINT 13
         PRINT 12

         DO i=1,Nk

            k=i*dk
            ik=INT(LOG(k/qlist(1))/dt+1)

c--------------------------------------------------------
            fk_Boltz=8.617343d-5 ! eV/K
            vlight=2.99792458d17

            EBreakUp=(3*qlist(ik)**2/4-g1**2)/mRb
            EBrkUpmK=EBreakUp*hbarc/fk_Boltz*1d6
c--------------------------------------------------------

            IF(ik.EQ.ik_last) GOTO 1
            ik_last=ik

            CALL ThreeBody(ik,cfscatt)

            cCot=1d0/cfscatt/qlist(ik)

c-------------------------------------------------------------
            cdel=-cI*LOG(1d0+2*cI*qlist(ik)*cfscatt)/2d0 !in radians

            mu=mRb/SQRT(3d0)
            IF(EBreakUp.LE.0d0) THEN
               RecomRate=0d0
            ELSE
               kBrkUp=SQRT(2*mu*EBreakUp)
               RecomFactor=0.75d0*192d0*pi**2/(mu*kBrkUp**4)*vlight*
     +              1d-42
               RecomRate=RecomFactor*(1d0-EXP(-4d0*DIMAG(cdel)))
            END IF
c-------------------------------------------------------------

            PRINT 10,EBrkUpmK,RecomRate,cCot,DIMAG(cdel)

 13         FORMAT(6X,'E (micro-K)',8X,'Recom Rate (cm^6/s)',10X,
     +           'Real kcot/k',10X,'Imaginary kcot/k')
 12         FORMAT(1X,4(20('='),4X))

 10         FORMAT(3X,f15.8,9X,e15.8,7X,f15.8,11X,e15.8,6X,e15.8)
 11         FORMAT(f15.8,2X,e15.8)

 1          CONTINUE

         END DO                 ! k loop
         
      END DO                    ! B loop

      RETURN
      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------------c
      SUBROUTINE ThreeBody(ik,cfscatt)

      IMPLICIT REAL*8 (a-b,d-h,k-m,o-z)
      IMPLICIT INTEGER (i,j,n)
      IMPLICIT COMPLEX*16 (c)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/dat/pi,hbarc,mRb,Ehf
      COMMON/const/Wfnorm,g1,g2,b1,b2
      DATA cI/(0d0,1d0)/,c1/(1d0,0d0)/
      DATA nrank/2/
      DIMENSION cKern(nmax,nmax),IPIVOT(nmax),cXsoln(nmax)

      epsln=1d-25
      ck=qlist(ik)
      cEk=3d0/4d0*ck**2/mRb-g1**2/mRb+cI*epsln/mRb !Bound state energy

      Ndim=Nq+Nphi              !Max # of pts for one "nrank"

      cphz0=EXP(-cI*phi(1))

c     Now for the complex X

      DO ip=1,Ndim
         IF(ip.LE.Nq) THEN
            cp=qlist(ip)*cphz0
         ELSE
            cp=qlist(Nq)*EXP(-cI*phi(ip-Nq))
         END IF

c     Now for the same with q (integration variable)
               DO iq=1,Ndim
                  IF(iq.LE.Nq) THEN
                     cq=qlist(iq)*cphz0
                     cdq_dphi=cq*dt*wt(iq)
                  ELSE
                     cq=qlist(Nq)*EXP(-cI*phi(iq-Nq))
                     cdq_dphi=cq*cI*dphi*wt(iq)
                  END IF

c     Now we calculate the M matrix.  This is going to have nrank in it
                  DO ir1=1,nrank
                     DO ir2=1,nrank
                        cKern(ip+(ir1-1)*Ndim,iq+(ir2-1)*Ndim)=
     +                       -cKernel(cp,cq,cEk,ir1,ir2)
     +                       *cq**2*cdq_dphi
                     END DO     !end of ir2
                  END DO        !end of ir1
               END DO           !end of q loop

c     Now add the delta function
               DO ir1=1,nrank
                  cKern(ip+(ir1-1)*Ndim,ip+(ir1-1)*Ndim)=c1+
     +                 cKern(ip+(ir1-1)*Ndim,ip+(ir1-1)*Ndim)
               END DO           !end of ir1

c     Now we calc the 2*Z
               DO ir=1,nrank
                  cXsoln(ip+(ir-1)*Ndim)=2*cZfun(cp,ck,cEk,ir,1)
               END DO

      END DO                    !end of p loop

c     Now we calc X=M^-1*2*Z
      CALL ZGESV(Ndim*nrank,1,cKern,nmax,IPIVOT,cXsoln,nmax,INFO)

c     Now we get the real X matrix

      DO iq=1,Ndim
         IF(iq.LE.Nq) THEN
            cq=qlist(iq)*cphz0
            cdq_dphi=cq*dt*wt(iq)
         ELSE
            cq=qlist(Nq)*EXP(-cI*phi(iq-Nq))
            cdq_dphi=cq*cI*dphi*wt(iq)
         END IF

         DO ir=1,nrank
            cKern(ik,iq+(ir-1)*Ndim)=cKernel(ck,cq,cEk,1,ir)*
     +           cq**2*cdq_dphi
         END DO                 !end of ir loop

      END DO                    !end of q loop

      cTmatrix=2d0*cZfun(ck,ck,cEk,1,1)

      DO iq=1,Ndim
         DO ir=1,nrank
            cTmatrix=cTmatrix+cKern(ik,iq+(ir-1)*Ndim)*
     +           cXsoln(iq+(ir-1)*Ndim)
         END DO                 !end of ir loop
      END DO                    !end of q loop

      cfscatt=-mRb*cTmatrix/(3*pi*Wfnorm)

      RETURN
      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------------c
      COMPLEX *16 FUNCTION cKernel(cp,cq,cEk,ir1,ir2)

      IMPLICIT REAL*8 (a-b,d-h,k-m,o-z)
      IMPLICIT INTEGER (i,j,n)
      IMPLICIT COMPLEX*16 (c)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/dat/pi,hbarc,mRb,Ehf
      COMMON/const/Wfnorm,g1,g2,b1,b2
      DATA cI/(0d0,1d0)/,c1/(1d0,0d0)/
      DATA nrank/2/
      DIMENSION cKern(nmax,nmax),cXsoln(nmax)

c     Here is the energy for the tau matrix
      cEtau=cEk-3d0/4d0*cq**2/mRb

      cKernel=(0d0,0d0)

      DO i=1,2
         cKernel=cKernel+2d0*cZfun(cp,cq,cEk,ir1,i)*
     +        cTaufun(cEtau,i,ir2)/(2*pi**2)
      END DO

      RETURN
      END
c-----------------------------------------------------------------------------c

c-----------------------------------------------------------------------------c
      COMPLEX*16 FUNCTION cZfun(cp,cq,cEk,ir1,ir2)

c     This function returns Z (uses gaussian integration)
c     integral over the angle theta

      IMPLICIT REAL*8 (a-b,d-h,l,m,o-z)
      IMPLICIT COMPLEX*16 (c-c)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/dat/Pi,hbarc,mRb,Ehf
      COMMON/const/Wfnorm,g1,g2,beta1,beta2
      COMMON/hf/lmat,Echan,Ps,Pt,cStr
      DIMENSION cZmath(2,2),cZmat(2,2),lmat(2,2)
      DIMENSION clmat(2,2)
      DIMENSION Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      REAL*8 Xgauss(ngauss),Wgauss(ngauss)
      DATA c1/(1d0,0d0)/

      CALL GAUSSLEG(-1d0,1d0,Xgauss,Wgauss,ngauss)

c     here are the cp^2 and cq^2 factors (less room in actual calc)
      cp2=cp**2
      cq2=cq**2
c     here is the "dot" term
      cpq=cp*cq

c     here the initial term in the integration
      cangle=(0d0,0d0)

c     Here is the actual integral using Gaussian Integration

c     First we do element (1,1)
      DO i=1, ngauss

         cpqX=cpq*Xgauss(i)     ! Dot product with Xgauss weight

c     Here is the p+1/2q squared term
         cksq1=cp2+1d0/4d0*cq2-cpqX
c     Here is the q+1/2p squared term
         cksq2=cq2+1d0/4d0*cp2-cpqX

         cG0=1d0/(cp2+cq2-mRb*cEk+cpqX)*
     +        beta1**4/((cksq1+beta1**2)*(cksq2+beta1**2)) ! Term to integrate
         cangle=cangle+Wgauss(i)*cG0 ! Integral Answer
      END DO

      cZmath(1,1)=-mRb*cangle/2d0 ! 2 comes from the angular integration

c     Now for the off diagonal bits
      cZmath(1,2)=(0d0,0d0)
      cZmath(2,1)=(0d0,0d0)

c     Now for the (2,2) element
      cangle=(0d0,0d0)

      DO i=1, ngauss

         cpqX=cpq*Xgauss(i)     ! Dot product with Xgauss weight

c     Here is the p+1/2q squared term
         cksq1=cp2+1d0/4d0*cq2-cpqX
c     Here is the q+1/2p squared term
         cksq2=cq2+1d0/4d0*cp2-cpqX

         cG0=1d0/(cp2+cq2-mRb*cEk+cpqX+Echan(2)**2)*
     +        beta2**4/((cksq1+beta2**2)*(cksq2+beta2**2)) ! Term to integrate
         cangle=cangle+Wgauss(i)*cG0 ! Integral Answer
      END DO

      cZmath(2,2)=-mRb*cangle/2d0

      DO i=1,2
         DO j=1,2
            clmat(i,j)=c1*lmat(i,j)
         END DO
      END DO

c     Now with the Lambda matrix
      DO i=1,2
         DO j=1,2
            chold=0d0
            DO k=1,2
               chold=chold+clmat(i,k)*cZmath(k,j)
            END DO
            cZmat(i,j)=chold
         END DO
      END DO

      cZfun=cZmat(ir1,ir2)

      RETURN
      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------------c
      COMPLEX*16 FUNCTION cTaufun(cEtau,ir1,ir2)

      IMPLICIT REAL*8 (p,a,h,m,e,b,l,g,w)
      IMPLICIT COMPLEX*16 (c)
      IMPLICIT INTEGER (i,j)
      COMMON/dat/pi,hbarc,mRb,Ehf
      COMMON/hf/LamMat,Echan,Ps,Pt,cStr
      COMMON/const/Wfnorm,g1,g2,b1,b2
      DIMENSION cTaumat(2,2)
      DIMENSION LamMat(2,2),Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      DATA cI/(0d0,1d0)/

      i1=1
      IF(i1.EQ.1) THEN
c-------------------FROM MATHEMATICA----------------------------------------------------------------------
c-------------------FROM YamawithVaryB_CHANNELB.nb----ON 2-24-2010----------------------------------------

         cTaumat(1,1)=(8*(b1**2 + cEtau*mRb)**2*Pi**2*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -      8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))))/
     -  (Pi*(64*cEtau**3*mRb**3*Pi**2 + 8*cEtau**2*mRb**2*Pi*
     -        (-(mRb*(b1**3*cStr(1,1) + b2**3*cStr(2,2))) + 
     -          8*Pi*(2*b1**2 - b2**2 - Echan(2)**2 - 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))) - 
     -       b1**4*(b1*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -          64*Pi**2*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)) + 
     -          8*mRb*Pi*(b1*b2**2*cStr(1,1) + b2**3*cStr(2,2) + b1*cStr(1,1)*Echan(2)**2 + 
     -             2*b1*b2*cStr(1,1)*Sqrt(-(cEtau*mRb) + Echan(2)**2))) + 
     -       b1**2*cEtau*mRb*(b1*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -          8*mRb*Pi*(b1*b2**2*cStr(1,1) - 2*b2**3*cStr(2,2) + b1*cStr(1,1)*(b1**2 + Echan(2)**2) + 
     -             2*b1*b2*cStr(1,1)*Sqrt(-(cEtau*mRb) + Echan(2)**2)) + 
     -          64*Pi**2*(b1**2 - 2*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))))) + 
     -    b1**4*mRb*Sqrt(cEtau*mRb)*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -       8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))*Log(cEtau*mRb) + 
     -    2*b1**4*mRb*Sqrt(cEtau*mRb)*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -       8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))*Log(-(1/Sqrt(cEtau*mRb))))

         cTaumat(1,2)=(64*(b1**2 + cEtau*mRb)**2*Pi**3*cStr(1,2)*(-b2**2 + cEtau*mRb - Echan(2)**2 - 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))/
     -  (Pi*(64*cEtau**3*mRb**3*Pi**2 + 8*cEtau**2*mRb**2*Pi*
     -        (-(mRb*(b1**3*cStr(1,1) + b2**3*cStr(2,2))) + 
     -          8*Pi*(2*b1**2 - b2**2 - Echan(2)**2 - 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))) - 
     -       b1**4*(b1*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -          64*Pi**2*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)) + 
     -          8*mRb*Pi*(b1*b2**2*cStr(1,1) + b2**3*cStr(2,2) + b1*cStr(1,1)*Echan(2)**2 + 
     -             2*b1*b2*cStr(1,1)*Sqrt(-(cEtau*mRb) + Echan(2)**2))) + 
     -       b1**2*cEtau*mRb*(b1*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -          8*mRb*Pi*(b1*b2**2*cStr(1,1) - 2*b2**3*cStr(2,2) + b1*cStr(1,1)*(b1**2 + Echan(2)**2) + 
     -             2*b1*b2*cStr(1,1)*Sqrt(-(cEtau*mRb) + Echan(2)**2)) + 
     -          64*Pi**2*(b1**2 - 2*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))))) + 
     -    b1**4*mRb*Sqrt(cEtau*mRb)*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -       8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))*Log(cEtau*mRb) + 
     -    2*b1**4*mRb*Sqrt(cEtau*mRb)*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -       8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))*Log(-(1/Sqrt(cEtau*mRb))))

         cTaumat(2,1)=(64*(b1**2 + cEtau*mRb)**2*Pi**3*cStr(1,2)*(-b2**2 + cEtau*mRb - Echan(2)**2 - 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))/
     -  (Pi*(64*cEtau**3*mRb**3*Pi**2 + 8*cEtau**2*mRb**2*Pi*
     -        (-(mRb*(b1**3*cStr(1,1) + b2**3*cStr(2,2))) + 
     -          8*Pi*(2*b1**2 - b2**2 - Echan(2)**2 - 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))) - 
     -       b1**4*(b1*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -          64*Pi**2*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)) + 
     -          8*mRb*Pi*(b1*b2**2*cStr(1,1) + b2**3*cStr(2,2) + b1*cStr(1,1)*Echan(2)**2 + 
     -             2*b1*b2*cStr(1,1)*Sqrt(-(cEtau*mRb) + Echan(2)**2))) + 
     -       b1**2*cEtau*mRb*(b1*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -          8*mRb*Pi*(b1*b2**2*cStr(1,1) - 2*b2**3*cStr(2,2) + b1*cStr(1,1)*(b1**2 + Echan(2)**2) + 
     -             2*b1*b2*cStr(1,1)*Sqrt(-(cEtau*mRb) + Echan(2)**2)) + 
     -          64*Pi**2*(b1**2 - 2*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))))) + 
     -    b1**4*mRb*Sqrt(cEtau*mRb)*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -       8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))*Log(cEtau*mRb) + 
     -    2*b1**4*mRb*Sqrt(cEtau*mRb)*(8*cEtau*mRb*Pi*cStr(1,1) + b2**3*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2)) - 
     -       8*Pi*cStr(1,1)*(b2**2 + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2)))*Log(-(1/Sqrt(cEtau*mRb))))

         cTaumat(2,2)=(Pi*(8*cEtau**2*mRb**2*Pi*cStr(2,2) + b1**2*cEtau*mRb*(16*Pi*cStr(2,2) + b1*mRb*(cStr(1,2)**2 - cStr(1,1)*cStr(2,2))) + 
     -       b1**4*(8*Pi*cStr(2,2) + b1*mRb*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)))) + 
     -    b1**4*mRb*Sqrt(cEtau*mRb)*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2))*Log(cEtau*mRb) + 
     -    2*b1**4*mRb*Sqrt(cEtau*mRb)*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2))*Log(-(1/Sqrt(cEtau*mRb))))/
     -  (8.*(b1**2 + cEtau*mRb)**2*Pi**2*((b1**3*b2**3*mRb**2*cStr(1,2)**2*
     -         (-(b1**2*Pi) + cEtau*mRb*Pi - b1*Sqrt(cEtau*mRb)*Log(cEtau*mRb) - 2*b1*Sqrt(cEtau*mRb)*Log(-(1/Sqrt(cEtau*mRb))))
     -         )/(64.*(b1**2 + cEtau*mRb)**2*Pi**3*(b2**2 - cEtau*mRb + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))) + 
     -      (1 + (b2**3*mRb*cStr(2,2))/(8.*Pi*(b2**2 - cEtau*mRb + Echan(2)**2 + 2*b2*Sqrt(-(cEtau*mRb) + Echan(2)**2))))*
     -       (1 - (b1**3*mRb*cStr(1,1)*(-(b1**2*Pi) + cEtau*mRb*Pi - b1*Sqrt(cEtau*mRb)*Log(cEtau*mRb) - 
     -              2*b1*Sqrt(cEtau*mRb)*Log(-(1/Sqrt(cEtau*mRb)))))/(8.*(b1**2 + cEtau*mRb)**2*Pi**2))))

         ctauDen=-1d0
c---------------------------------------------------------------------------------------------------------
      END IF

      cThold=-cTaumat(ir1,ir2)*1d0/ctauDen

      cTaufun=cThold

      RETURN
      END
c-----------------------------------------------------------------------------c
c-----------------------------------------------------------------------------c
      COMPLEX*16 FUNCTION cDetTau(cEtau)

      IMPLICIT REAL*8 (p,a,h,m,e,b,l,g,w)
      IMPLICIT COMPLEX*16 (c)
      IMPLICIT INTEGER (i,j)
      COMMON/dat/pi,hbarc,mRb,Ehf
      COMMON/hf/LamMat,Echan,Ps,Pt,cStr
      COMMON/const/Wfnorm,g1,g2,b1,b2
      DIMENSION cTaumat(2,2)
      DIMENSION LamMat(2,2),Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      DATA cI/(0d0,1d0)/

c---------------------------------------------------------------------
      cFac=(1d0,0d0)*mRb/(8d0*pi)
      cDetLam=cStr(1,1)*cStr(2,2)-cStr(1,2)**2
      cEnDel=SQRT(-cEtau*mRb+Echan(2)**2)

      c1=-(SQRT(-cEtau*mRb)+b1)**2
      c2=-(cEnDel+b2)**2
      c1c2=c1*c2

      ch11=cFac*b1**3/c1
      ch22=cFac*b2**3/c2

      cDetJ=(1d0/cDetLam)**2*((cStr(2,2)-cDetLam*ch11)*
     +     (cStr(1,1)-cDetLam*ch22)-cStr(1,2)**2)

      ctauDen=cDetJ*cDetLam

      cTaumat(1,1)=cStr(1,1)+ch22*cDetLam
      cTaumat(1,2)=cStr(1,2)
      cTaumat(2,1)=cStr(2,1)
      cTaumat(2,2)=cStr(2,2)+ch11*cDetLam
c---------------------------------------------------------------------

      cT1122=cTaumat(1,1)*cTaumat(2,2)
      cT1212=cTaumat(1,2)*cTaumat(2,1)

      cDT=(cT1122-cT1212)*1d0/ctauDen**2

      cDetTau=cDT

      RETURN
      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------------c
      SUBROUTINE TwoBodyHyperfine(B)

      IMPLICIT REAL*8 (a,b,g,l,e,p,r,v,w,m,h)
      IMPLICIT COMPLEX*16 (c)
      COMMON/hf/LamMat,Echan,Ps,Pt,cStr
      COMMON/const/Wfnorm,gamma1,gamma2,b1,b2
      COMMON/dat/Pi,hbarc,mRb,Ehf
      COMMON/params/ipar
      DIMENSION LamMat(2,2),Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      DIMENSION lmat(2,2),lscale(2,2),lhold(5,5)
      cI=(0d0,1d0)
      c1=(1d0,0d0)

c------CLOSED CHANNEL---------
      ich=3

      CALL HyperSub(B,lhold,Echan,Ps,Pt)

      DO i=1,5
         Echan(i)=SQRT(Echan(i)*mRb*Ehf/hbarc)
      END DO

c----------------------
      Echan(2)=Echan(ich)
c----------------------

c---------INITIALIZING SCALES-------
         lscale(1,1)=1d0
         lscale(1,2)=1.2d0
         lscale(2,1)=lscale(1,2)
         lscale(2,2)=0.75d0

         beta1scale=1d0
         beta2scale=18.9761d0

      IF(ipar.EQ.1) THEN
c----------------------------------------------------------------------
c----------------------------JOSE PARAMETERS---------------------------
c     FROM recomb.f (McNeil's one channel code)

c----------SINGLET------------------
         as=126.27305582869954d0 ! JAM fit March 26, 2009 fit to pCot(delta)
         r0=17.933315071786954d0 ! RbRb_scattering_3.26.nb
         v4=52.17034593703047d0
         
         ls=2*pi*as/(mRb*(as-2*r0))*(3*as-4*r0-sqrt(9*as**2-16*r0*as))
         beta1=(3+sqrt(9d0-16*r0/as))/(2*r0)
c-----------------------------------

c----------TRIPLET------------------
         as=-19.71763203438355d0 !  JAM fit March 26, 2009 fit to pCot(delta)
         r0=11.303238103698735d0 !  RbRb_scattering_3.26.nb
         v4=133.69226101569353d0

         lt=2*pi*as/(mRb*(as-2*r0))*(3*as-4*r0-sqrt(9*as**2-16*r0*as))
         beta2=(3+sqrt(9d0-16*r0/as))/(2*r0)
c-----------------------------------
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      END IF
      IF(ipar.EQ.2) THEN
c----------SCALES------------------------c
c     Calculated in YamawithVaryB_CHANNELB.nb Mathematica Notebook
c     Last Edited on 2-12-10
         lscale(1,1)=1d0
         lscale(1,2)=1.2d0
         lscale(2,1)=lscale(1,2)
         lscale(2,2)=0.75d0

         beta2scale=18.9761d0
c----------------------------------------c
c--------------------FROM MATHEMATICA FROM LAST YEAR----------------
c     In YamawithVaryB_CHANNELB.nb Mathematica Notebook
c     Last Edited on 2-9-10
         gs=.00823028d0
         gt=-.03814310d0
         bs=.252917d0
         bt=.217046
         
         ls=-4d0*pi/(mRb*bs)*(1d0+gs/bs)**2
         lt=-4d0*pi/(mRb*bt)*(1d0+gt/bt)**2
         beta1=bs
         beta2=bt
c-------------------------------------------------------------------
      END IF

      lmat(1,1)=Ps(1,1)*ls+Pt(1,1)*lt
      lmat(1,2)=Ps(1,ich)*ls+Pt(1,ich)*lt
      lmat(2,1)=Ps(ich,1)*ls+Pt(ich,1)*lt
      lmat(2,2)=Ps(ich,ich)*ls+Pt(ich,ich)*lt

      DO i=1,2
         DO j=1,2
            cStr(i,j)=lmat(i,j)*lscale(i,j)
         END DO
      END DO


      bs=SQRT(cStr(1,1)/(ls*Ps(1,1)/beta1**2+lt*Pt(1,1)/beta2**2))
      bt=SQRT(cStr(2,2)/(ls*Ps(ich,ich)/beta1**2+lt*Pt(ich,ich)/
     +     beta2**2))

      IF(ipar.EQ.2) THEN
         b1=(bs+bt)/2d0
         b2=b1*beta2scale
      ELSE
         b1=bs*beta1scale
         b2=bt*beta2scale
      END IF

      LamMat(1,1)=lhold(1,1)
      LamMat(1,2)=lhold(1,3)
      LamMat(2,1)=lhold(3,1)
      LamMat(2,2)=lhold(3,3)

      IF(ipar.EQ.3) THEN
c----------------FROM McNeil----------------------------------
         DATA betaS/ 0.2572039403130236d0/,betaT/ 0.22506901222498846d0/
         DATA LamS/ -2.594093887416936d-7/,LamT/ -1.910039188030657d-7/
         DATA LamScale/4.183d0/,betaScale/5.9d0/ ! Fits from McNeil-checked in YamawithVaryB_CHANNELB.nb
                                ! on 2-24-2010
         
         cStr(1,1)=Ps(1,1)*LamS+Pt(1,1)*LamT
         cStr(1,2)=-(Ps(1,ich)*LamS+Pt(1,ich)*LamT)
         cStr(2,1)=-(Ps(ich,1)*LamS+Pt(ich,1)*LamT)
         cStr(2,2)=(Ps(ich,ich)*LamS+Pt(ich,ich)*LamT)*LamScale
         
         b1=(betaS+betaT)/2d0
         b2=b1*betaScale
c-------------------------------------------------------------
      END IF

      CALL CalcGamma()

      RETURN
      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------------c
      SUBROUTINE CalcGamma()

      IMPLICIT REAL*8 (a-b,d-h,l-m,o-z)
      IMPLICIT COMPLEX*16 (c)
      PARAMETER (nphimax=257,nmax=8193+nphimax,ngauss=32)
      COMMON/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      COMMON/hf/LamMat,Echan,Ps,Pt,cStr
      COMMON/const/Wfnorm,g1,g2,b1,b2
      COMMON/dat/Pi,hbarc,mRb,Ehf
      DIMENSION LamMat(2,2),Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      DATA c1/(1d0,0d0)/,cI/(0d0,1d0)/

      gnew=1d0
      gold=gnew

      DO i=1,1000
         f=Denom(gnew)
         fp=DDenom(gnew)

         gnew=gnew-f/fp

         IF(gold.LE.gnew) THEN
            ghold=gnew-gold
         ELSE
            ghold=gold-gnew
         END IF

         IF(ghold.LE.1d-15) GOTO 1

         gold=gnew
      END DO

 1    CONTINUE

      g1=gnew
      g2=g1

      PRINT*
      PRINT*,'Gamma = ',g1

      cEn=-c1*g1**2/mRb
      cDethold1=cDetTau(cEn)
      cDethold2=cTaufun(cEn,1,1)*cTaufun(cEn,2,2)-
     +     cTaufun(cEn,1,2)*cTaufun(cEn,2,1)

      PRINT*
      PRINT*,'Elements of Tau matrix for Energy = ',REAL(cEn)
      PRINT*,'Det(Tau)1 =  ',cDethold1,' from subroutine'
      PRINT*,'Det(Tau)2 =  ',cDethold2

      PRINT*
      PRINT*,'Strength Matrix = '
      DO i=1,2
         PRINT*,(REAL(cStr(i,j)),j=1,2)
      END DO

      PRINT*
      PRINT*,'Beta Open = ',b1
      PRINT*,'Beta Closed = ',b2

      Wfnorm=1d0*mRb**2*b1**3/(8d0*pi*g1*(b1+g1)**3)

c     Here are the corrections to the normalization
c     155.2G = 1.9375
c     156G = 2.01928

      PRINT*
      PRINT*,'Waveform Normalization = ',1d0/Wfnorm
      PRINT*

c--------------------------------------------------------
c     ALL INPUT PARAMETERS
c      PRINT*
c      PRINT*,'m = ',mRb
c      PRINT*,'gamma = ',g1
c      PRINT*,'b1 = ',b1
c      PRINT*,'b2 = ',b2
c      PRINT*,'delta = ',Echan(2)
c      PRINT*,'lambda o = ',REAL(cStr(1,1))
c      PRINT*,'lambda x = ',REAL(cStr(1,2))
c      PRINT*,'lambda c = ',REAL(cStr(2,2))
c      PRINT*
c--------------------------------------------------------

      PRINT*,'Det Tau Testing'
      DO i1=-2,2
         cEk=-g1**2/mRb*(1d0+i1*0.01d0)

         cDethold2=cTaufun(cEk,1,1)*cTaufun(cEk,2,2)-
     +        cTaufun(cEk,1,2)*cTaufun(cEk,2,1)

         PRINT*,'Subroutine ',cDetTau(cEk),' Function ',cDethold2
      END DO
      PRINT*


      RETURN
      END
c-----------------------------------------------------------------------------c

c-----------------------------------------------------------------------------c
      REAL*8 FUNCTION Denom(En)

      IMPLICIT REAL*8 (a-b,d-h,l-m,o-z)
      IMPLICIT COMPLEX*16 (c)
      COMMON/hf/LamMat,Echan,Ps,Pt,cStr
      COMMON/const/Wfnorm,ga1,g2,b1,b2
      COMMON/dat/Pi,hbarc,mRb,Ehf
      DIMENSION LamMat(2,2),Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      DIMENSION Str(2,2)
      DATA c1/(1d0,0d0)/

      DO i=1,2
         DO j=1,2
            Str(i,j)=REAL(cStr(i,j))
         END DO
      END DO

      g1=En

c     From hmatrix.nb on 2-16-2010

      Denom=b1**3*b2**3*mRb**2*(-cStr(1,2)**2 + cStr(1,1)*cStr(2,2)) + 
     -  64*(b1 + g1)**2*Pi**2*
     -   (b2**2 + g1**2 + Echan(2)**2 + 2*b2*Sqrt(g1**2 + Echan(2)**2))
     -    + 8*mRb*Pi*(b1**2*b2**3*cStr(2,2) + 
     -     2*b1*b2**3*g1*cStr(2,2) + b2**3*g1**2*cStr(2,2) + 
     -     b1**3*cStr(1,1)*(b2**2 + g1**2 + Echan(2)**2 + 
     -        2*b2*Sqrt(g1**2 + Echan(2)**2)))

      RETURN
      END
c-----------------------------------------------------------------------------c

c-----------------------------------------------------------------------------c
      REAL*8 FUNCTION DDenom(En)

      IMPLICIT REAL*8 (a-b,d-h,l-m,o-z)
      IMPLICIT COMPLEX*16 (c)
      COMMON/hf/LamMat,Echan,Ps,Pt,cStr
      COMMON/const/Wfnorm,ga1,g2,b1,b2
      COMMON/dat/Pi,hbarc,mRb,Ehf
      DIMENSION LamMat(2,2),Echan(5),Ps(5,5),Pt(5,5),cStr(2,2)
      DIMENSION Str(2,2)

      DO i=1,2
         DO j=1,2
            Str(i,j)=REAL(cStr(i,j))
         END DO
      END DO

      g1=En

c     From hmatrix.nb on 2-16-2010

      DDenom=64*(b1 + g1)**2*Pi**2*
     -   (2*g1 + (2*b2*g1)/Sqrt(g1**2 + Echan(2)**2)) + 
     -  128*(b1 + g1)*Pi**2*(b2**2 + g1**2 + Echan(2)**2 + 
     -     2*b2*Sqrt(g1**2 + Echan(2)**2)) + 
     -  8*mRb*Pi*(2*b1*b2**3*cStr(2,2) + 2*b2**3*g1*cStr(2,2) + 
     -     b1**3*cStr(1,1)*(2*g1 + (2*b2*g1)/Sqrt(g1**2 + Echan(2)**2))
     -     )

      RETURN
      END
c-----------------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GAUSSLEG(X1,X2,X,W,N)

c     Taken from McNeil's code

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      pi=2.d0*asin(1.d0)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
         Z=COS(pi*(I-.25D0)/(N+.5D0))
 1       CONTINUE
         P1=1.D0
         P2=0.D0
         DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
 11      CONTINUE
         PP=N*(Z*P1-P2)/(Z*Z-1.D0)
         Z1=Z
         Z=Z1-P1/PP
         IF(ABS(Z-Z1).GT.EPS)GO TO 1
         X(I)=XM-XL*Z
         X(N+1-I)=XM+XL*Z
         W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
         W(N+1-I)=W(I)
 12   CONTINUE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Weights

c     From McNeil's code

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax)
      common/q_grid/ dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      dimension bode8(8)
      data bode8/1978d0,5888d0,-928d0,10496d0,-4540d0,10496d0,-928d0
     +     ,5888d0/

      ndim=Nq+Nphi
      
      
c     Bode's 8-point	  
      do i=1,Nq
         wt(i)=4*bode8(mod(i-1,8)+1)/14175
      enddo
      wt(1)=wt(1)/2
      wt(Nq)=wt(Nq)/2
      do i=1,Nphi
         wt(i+Nq)=4*bode8(mod(i-1,8)+1)/14175
      enddo
      wt(Nq+1)=wt(Nq+1)/2
      wt(Nq+Nphi)=wt(Nq+Nphi)/2
      
      return
      end	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
