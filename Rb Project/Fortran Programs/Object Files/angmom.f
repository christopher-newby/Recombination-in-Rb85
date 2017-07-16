c	Angular momentum package 
c
c		VCC   = Vector coupling coefficients
c               Clebsch = Clebsch-Gordon coefficients
c		Racah = Racah 6-j coefficients
c		Winej = Wigner 9-j coefficients
c
c   CAUTION:  ALL ANGULAR MOMENTA ARE TWICE THEIR ACTUAL VALUES !!!!
c
      FUNCTION WIGNER3J(JX1,JX2,JX3,MX1,MX2,MX3)
	  implicit double precision (a-h,o-z)
      WIGNER3J=0.0
      IF(MX1+MX2+MX3.NE.0)RETURN
      PHZ=(JX2+MX3-JX1)/2
      ROOT=SQRT(FLOAT(JX3+1))
      WIGNER3J= (-1.)**PHZ*VCC(JX1,JX2,JX3,MX1,MX2)/ROOT
      RETURN
      END
C
      Function CG(jx1,jx2,jx3,mx1,mx2,mx3)
      implicit double precision (a-h,o-z)
	  cg=0.d0
	  if(mx1+mx2.ne.mx3)return
	  cg=vcc(jx1,jx2,jx3,mx1,mx2)
	  return
	  end
c      	
C
      Function Clebsch(jx1,mx1,jx2,mx2,jx3,mx3) ! re-order input
      implicit double precision (a-h,o-z)
	  clebsch=0.d0
	  if(mx1+mx2.ne.mx3)return
	  clebsch=vcc(jx1,jx2,jx3,mx1,mx2)
	  return
	  end
c
      FUNCTION VCC(JX1,JX2,JX3,MX1,MX2)
	  implicit double precision (a-h,o-z)
      DIMENSION F(34),FACT(33)
      EQUIVALENCE (FACT(1),F(2))
      DATA F   /1.0D+00,1.0D+00,2.0D+00,6.0D+00,2.4D+01,1.20D+02,
     17.20D+02,5.040D+03,4.0320D+04,3.62880D+05,3.628800D+06,
     23.9916800D+07,4.79001600D+08,6.22702080D+09,8.71782912D+10,
     31.30767437D+12,2.09227899D+13,3.55687428D+14,6.40237371D+15,
     41.21645100D+17,2.43290200D+18,5.10909422D+19,1.12400073D+21,
     52.58520167D+22,6.20448401D+23,1.55112100D+25,4.03291461D+26,
     61.08888694D+28,3.04888344D+29,8.84176198D+30,2.65252859D+32,
     78.22283864D+33,2.63130836D+35,8.68331760D+36/
      J1=JX1
      J2=JX2
      J3=JX3
      M1=MX1
      M2=MX2
      IF(J1-J2)20,10,10
   10 IF(J3-J2)30,15,15
   15 CNTRL=0.0
      GO TO 40
   20 IF(J3-J1)30,25,25
   25 CNTRL=-1.0
      IT=J1
      J1=J2
      J2=IT
      IT=M1
      M1=M2
      M2=IT
      GO TO 40
   30 CNTRL=1.0
      IT=J2
      J2=J3
      J3=IT
      M2=-M1-M2
   40 JZ1=(J1+J2-J3)/2
      IF(JZ1)70,45,45
   45 JZ2=(J1+J3-J2)/2
      IF(JZ2)70,50,50
   50 JZ3=(J2+J3-J1)/2
      IF(JZ3)70,55,55
   55 IF(J1-IABS(M1))70,60,60
   60 IF(J2-IABS(M2))70,65,65
   65 IF(J3-IABS(M1+M2))70,75,75
   70 VCC=0.0
      GO TO 150
   75 JT1=(J1-J3+M2)/2
      JT2=(J2-J3-M1)/2
      NUMIN=MAX0 (JT1,JT2,0)
      JT3=(J1-M1)/2
      JT4=(J2+M2)/2
      NUMAX=MIN0 (JT3,JT4,JZ1)
      JT5=(J2-M2)/2
      IF(NUMAX-NUMIN)70,80,80
   80 J4=J1/2
      J5=J3/2
      FNUSM=0.0
      DO 100 NU=NUMIN,NUMAX
      JY1=JT4-NU
      JY2=NU-JT1
      JY3=JZ1-NU
      FCTOR=YXFCT(JT3-NU,J4)*YXFCT(NU-JT2,J5)
      C2=FACT(JY1)*FACT(JY2)*FACT(JY3)*FACT(NU)
  100 FNUSM=FNUSM+PHASEF(NU)*FCTOR/C2
      FCTOR=YXFCT((J1+J2+J3)/2+1,JZ2)*YXFCT(J4,(J1+M1)/2)*YXFCT(J4,JT3)*
     1YXFCT(J5,(J3+M1+M2)/2)*YXFCT(J5,(J3-M1-M2)/2)*FACT(JZ1)*FACT(JZ3)*
     2FACT(JT4)*FACT(JT5)*FLOAT(J3+1)
      Z5=SQRT(FCTOR)*FNUSM
      IF(CNTRL)120,130,110
  110 VCC=SQRT(FLOAT(J2+1)/FLOAT(J3+1))*PHASEF(JT3)*Z5
      GO TO 150
  120 VCC=Z5*PHASEF(JZ1)
      GO TO 150
  130 VCC=Z5
  150 RETURN
      END
c
      FUNCTION PHASEF(N)
	implicit double precision (a-h,o-z)
      PHASEF=FLOAT(1-2*IABS(N-2*(N/2)))
      RETURN
      END
c
      FUNCTION YXFCT(M,N)
	implicit double precision (a-h,o-z)
C     COMPUTES NFACT/MFACT
      YXFCT=1.0
      NUMAX=M-N
      IF(NUMAX)30,100,20
   20 CNTRL=-1.0
      FCTOR=FLOAT(N+1)
      GO TO 40
   30 CNTRL=1.0
      NUMAX=-NUMAX
      FCTOR=FLOAT(M+1)
   40 CONTINUE
      DO 50 NU=1,NUMAX
      YXFCT=YXFCT*FCTOR
   50 FCTOR=FCTOR+1.0
      IF(CNTRL)70,100,100
   70 YXFCT=1.0/YXFCT
  100 CONTINUE
      RETURN
      END
c
      FUNCTION WIGNER6J(J1,J2,J3,J4,J5,J6)
	implicit double precision (a-h,o-z)
        jph=(j1+j2+j4+j5)/2
    	phz=(-1.)**jph
	WIGNER6J=PHZ*RACAH(J1,J2,J5,J4,J3,J6)
        return
        end
c
      FUNCTION RACAH(J1,J2,J3,J4,J5,J6)
	implicit double precision (a-h,o-z)
      DIMENSION F(34),FACT(33)
      EQUIVALENCE (FACT(1),F(2))
      DATA F   /1.0D+00,1.0D+00,2.0D+00,6.0D+00,2.4D+01,1.20D+02,
     17.20D+02,5.040D+03,4.0320D+04,3.62880D+05,3.628800D+06,
     23.9916800D+07,4.79001600D+08,6.22702080D+09,8.71782912D+10,
     31.30767437D+12,2.09227899D+13,3.55687428D+14,6.40237371D+15,
     41.21645100D+17,2.43290200D+18,5.10909422D+19,1.12400073D+21,
     52.58520167D+22,6.20448401D+23,1.55112100D+25,4.03291461D+26,
     61.08888694D+28,3.04888344D+29,8.84176198D+30,2.65252859D+32,
     78.22283864D+33,2.63130836D+35,8.68331760D+36/
	RACAH=0.0
	Z1=DELR(J1,J2,J5)
	IF (Z1.EQ.0.0) GO TO 90
   	Z1=DELR(J3,J4,J5)*Z1
	IF (Z1.EQ.0.0) GO TO 90
	Z2=DELR(J1,J3,J6)
	IF (Z2.EQ.0.0) GO TO 90
	Z2=DELR(J2,J4,J6)*Z2
	IF (Z2.EQ.0.0) GO TO 90
	Z1=SQRT(Z1/Z2)*Z2
      JT1=(J1+J2+J5)/2
      JT2=(J3+J4+J5)/2
      JT3=(J1+J3+J6)/2
      JT4=(J2+J4+J6)/2
      JZ1=(J1+J2+J3+J4)/2
      JZ2=(J1+J4+J5+J6)/2
      JZ3=(J2+J3+J5+J6)/2
      NUMIN=MAX0 (JT1,JT2,JT3,JT4)
      NUMAX=MIN0 (JZ1,JZ2,JZ3)
	IF (NUMAX.LT.NUMIN) GO TO 90
	PHASE=PHASEF(NUMIN+JZ1)*Z1
      DO 80 NU=NUMIN,NUMAX
      JY1=NU-JT1
      JY2=NU-JT2
      JY3=NU-JT3
      JY4=JZ1-NU
      JY5=JZ2-NU
      JY6=JZ3-NU
      FCTOR=FACT(JY1)*FACT(JY2)*FACT(JY3)*YXFCT(NU+1,NU-JT4)
     1*FACT(JY4)*FACT(JY5)*FACT(JY6)
	RACAH=RACAH+PHASE/FCTOR
	PHASE=-PHASE
80	CONTINUE
90     RETURN
      END
c
      FUNCTION DELR(J1,J2,J3)
	implicit double precision (a-h,o-z)
      DIMENSION F(34),FACT(33)
      EQUIVALENCE (FACT(1),F(2))
      DATA F   /1.0D+00,1.0D+00,2.0D+00,6.0D+00,2.4D+01,1.20D+02,
     17.20D+02,5.040D+03,4.0320D+04,3.62880D+05,3.628800D+06,
     23.9916800D+07,4.79001600D+08,6.22702080D+09,8.71782912D+10,
     31.30767437D+12,2.09227899D+13,3.55687428D+14,6.40237371D+15,
     41.21645100D+17,2.43290200D+18,5.10909422D+19,1.12400073D+21,
     52.58520167D+22,6.20448401D+23,1.55112100D+25,4.03291461D+26,
     61.08888694D+28,3.04888344D+29,8.84176198D+30,2.65252859D+32,
     78.22283864D+33,2.63130836D+35,8.68331760D+36/
      JZ1=(J1+J2-J3)/2
      IF(JZ1)30,10,10
   10 JZ2=(J1-J2+J3)/2
      IF(JZ2)30,20,20
   20 JZ3=(J2+J3-J1)/2
      IF(JZ3)30,40,40
   30 DELR=0.0
      GO TO 100
   40 JZ4=(J1+J2+J3)/2+1
      IF(JZ3-JZ2)80,50,50
   50 IF(JZ3-JZ1)70,60,60
   60 DELR=YXFCT(JZ4,JZ3)*FACT(JZ1)*FACT(JZ2)
      GO TO 100
   70 DELR=YXFCT(JZ4,JZ1)*FACT(JZ2)*FACT(JZ3)
      GO TO 100
   80 IF(JZ2-JZ1)70,90,90
   90 DELR=YXFCT(JZ4,JZ2)*FACT(JZ1)*FACT(JZ3)
  100 CONTINUE
      RETURN
      END
c
      FUNCTION WIGNER9J(J1,J2,J3,J4,J5,J6,J7,J8,J9)
	implicit double precision (a-h,o-z)
      JT1=IABS(J1-J9)
      JT2=IABS(J2-J6)
      JT3=IABS(J4-J8)
      MUMIN=MAX0 (JT1,JT2,JT3)
      JZ1=J1+J9
      JZ2=J2+J6
      JZ3=J4+J8
      MUMAX=MIN0 (JZ1,JZ2,JZ3)
      IF(MUMAX-MUMIN)40,10,10
   10 FMUSM=0.0
       DO 20 MU=MUMIN,MUMAX,2
	T1=RACAH(J1,J4,J9,J8,J7,MU)
	T2=RACAH(J2,J5,MU,J4,J8,J6)
	T3=RACAH(J3,J6,J1,MU,J9,J2)
	PROD=T1*T2*T3*FLOAT(MU+1)
   20 FMUSM=FMUSM+PROD
      WIGNER9J=FMUSM*PHASEF((J2+J3+J5+J6+J8+J9)/2+J1+J4)
   30 CONTINUE
      RETURN
   40 WIGNER9J=0.0
      GO TO 30
      END
c
	subroutine rationalize(decimal,inum,idenom)
        implicit double precision (a-h,o-z)
c
c	by David Mercer, 1989 (modified 1991 by J.A.McNeil)
c	When given a decimal number, this program finds close
c	rational numbers.
c
	integer whole,sign
        tol=.00000000001d0
	sign=1
	if(decimal.lt.0)sign=-1
	whole=int(decimal*sign)
	x=abs(decimal)-float(whole)
	if(abs(x).lt.tol)then
		inum=0
		idenom=1
		guess=0.d0
		goto 20
	endif
        errorsv=1.d0
	inum=1
	idenom=2
10	rnum=float(inum)
	rden=float(idenom)
	guess=rnum/rden
        error=abs(x-guess)
        if(error.lt.errorsv)then
           inumsv=inum
           idnomsv=idenom
        endif
	if(abs(x-guess).le.tol)goto 20
        if(idenom.gt.999999) goto 30
	idenom=idenom+1
	if(x.gt.guess)inum=inum+1
	goto 10
 30     inum=inumsv
        idenom=idnomsv
 20	inum=sign*(inum+whole*idenom)
        return
	end
c================================================================
