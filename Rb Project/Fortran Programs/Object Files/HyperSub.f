c----------------------------------------------------------------------------------------------------------------------c

c       This program can be used to calculate the Projection matrices as well as the Deltas based upon a specific B-field
c       The last two subroutines are all that are needed, if the constants are used properly.
c       Will probably make whole file into an object file, and then send testHyperfine the Psing,Ptrip matrices, as well
c       as a matrix to hold the Deltas in, since this program stores them in a matrix.  Need to edit testHyperfine
c       so that it takes in a magnetic field and then outputs the specific values for that specific strength.
c       The code is basically already set up for this, it just needs to be tweaked

c       Haven't done this yet...

c----------------------------------------------------------------------------------------------------------------------c


cccccccccccccccccccccccccccccc		
	SUBROUTINE HyperSub(Bval,LamMat,Ech,Ps,Pt)
	implicit none
c       Constants
	double precision pi,amu,hbarc,clight,mu_Bohr,a_Bohr,Hartree,
     +    mu_Elec
	double precision hPlanck,alpha
	common/constants/pi,hbarc,clight,mu_Bohr,a_Bohr,mu_Elec,Hartree

	double precision M_Rb85
	double precision M_red,Ehf,Bhf
	common/Rbvalues/M_red,Ehf,Bhf
c       Working variables
	integer i,j,ipot
	
	double precision dB,Bmax,Bfield,Bval
	integer nB
	common/Bgrid2/ dB,Bfield(200),nB
	
	double precision V

	REAL*8 LamMat(5,5)

	double precision Echan
	common/E_channel/Echan(5)
	double precision Psmatrix(5,5),Ptmatrix(5,5)
	COMMON/ProjMat/Psmatrix,Ptmatrix

	REAL*8 LM(2,2),Ech(5),Ps(5,5),Pt(5,5)
	
C       
c       Constants
C       
	pi=2*asin(1d0)
	amu=931.49456d0
	hbarc=197.326968d0
	clight=2.99792458d17
	mu_Bohr=5.788381804d-5
	a_Bohr=0.052917720859d0
	
	M_Rb85=84.911789738d0*amu
	M_red=M_Rb85/2.d0
	mu_Elec=1.0011596521859d0*mu_Bohr
	hPlanck=2*pi*hbarc/clight
	Ehf=3.036d9*hPlanck	! eV
	Bhf=1d4*Ehf/mu_elec	! Gauss
	Hartree=alpha*hbarc/a_Bohr
c       
c       B-grid
c       
	Bmax=200d0		!Gauss
	nB=11
	dB=Bmax/(nB-1)
	do i=1,nB
	   Bfield(i)=(i-1)*dB/Bhf ! B in units of hyperfine equivalent B-field
	enddo
c       
	Bfield(1)=Bval/Bhf
	call Hyperfine(1)

	CALL CalcLambda(LamMat)

	DO i=1,5
	   Ech(i)=Echan(i)
	   DO j=1,5
	      Ps(i,j)=Psmatrix(i,j)
	      Pt(i,j)=Ptmatrix(i,j)
	   END DO
	END DO
	
	end
c       
C-----------------------------------------------------------------------
c	
	Subroutine Hyperfine(iB)
	implicit none
	integer nchan
	parameter(nchan=5)
	
c       Working variables
	integer i,j,iB,nB,im
	double precision ps,pt,dB,Bfield
	common/Bgrid2/ dB,Bfield(200),nB
	
	double precision Echan
	common/E_channel/Echan(nchan)
	
	integer f(2),mf(3),i_pow
	double precision H_hf_11,H_hf_22,H_1(2,2),H_2(2,2)
	double precision Eplus(2),Eminus(2),Eminus3,Egs
	double precision a,b,d,cos2theta,sign,sine(2),cosine(2)
	double precision sz_1,sz_2,sz_3
	double precision racah,clebsch,Srme,wigner6j
	double precision Psing,Psvector(nchan,2)
	double precision Psmatrix(nchan,nchan),Ptmatrix(nchan,nchan)
	data f/2,3/,mf/-1,-2,-3/


	COMMON/angle/sine,cosine
	COMMON/ProjMat/Psmatrix,Ptmatrix
	
	H_hf_11=(f(1)*(f(1)+1d0)-5d0/2*(5d0/2+1)-3d0/4)/6d0 ! f = 2
	H_hf_22=(f(2)*(f(2)+1d0)-5d0/2*(5d0/2+1)-3d0/4)/6d0 ! f = 3
	Srme=0.5d0*sqrt(3d0)	!  1/(2*clebsch(1,1,2,0,1,1))
	i_pow=5d0/2+3d0/2	!  i + 1 + 1/2
	do i=1,2		! f = (2,3)
	   do j=1,2
	      racah=(-1d0)**(i_pow+f(j))*sqrt(2*(2*f(i)+1d0))*
     +	   wigner6j(5,1,2*f(i),2,2*f(j),1) ! independent of mf
	      sz_1=clebsch(2*f(i),2*mf(1),2,0,2*f(j),2*mf(1))*racah*Srme
	      sz_2=clebsch(2*f(i),2*mf(2),2,0,2*f(j),2*mf(2))*racah*Srme
	      H_1(i,j)=2*Bfield(iB)*sz_1 ! mf = -1
	      H_2(i,j)=2*Bfield(iB)*sz_2 ! mf = -2
	      if(i.eq.j)then
		 if(i.eq.1)then
		    H_1(i,i)=H_1(i,i)+H_hf_11
		    H_2(i,i)=H_2(i,i)+H_hf_11
		 else
		    H_1(i,i)=H_1(i,i)+H_hf_22
		    H_2(i,i)=H_2(i,i)+H_hf_22
		 endif
	      endif
	      if((i.eq.2).and.(j.eq.2))then ! mf = -3
		 sz_3=clebsch(2*f(i),2*mf(3),2,0,2*f(j),2*mf(3))*racah*
     +	      Srme
		 Eminus3=H_hf_22+2*Bfield(iB)*sz_3
	      endif
	   enddo
	enddo
	
	a=(H_1(1,1)+H_1(2,2))/2
	d=(H_1(1,1)-H_1(2,2))/2
	b=H_1(1,2)
	cos2theta=d/sqrt(b**2+d**2)
	cosine(1)=sqrt(0.5d0*(1+cos2theta))
	sign=1d0
	if(b.le.0d0)sign=-1d0
	sine(1)=sign*sqrt(0.5d0*(1-cos2theta))
	Eplus(1)=a+sqrt(b**2+d**2)
	Eminus(1)=a-sqrt(b**2+d**2)
	
	a=(H_2(1,1)+H_2(2,2))/2
	d=(H_2(1,1)-H_2(2,2))/2
	b=H_2(1,2)
	cos2theta=d/sqrt(b**2+d**2)
	cosine(2)=sqrt(0.5d0*(1+cos2theta))
	sign=1d0
	if(b.le.0d0)sign=-1d0
	sine(2)=sign*sqrt(0.5d0*(1-cos2theta))
	Eplus(2)=a+sqrt(b**2+d**2)
	Eminus(2)=a-sqrt(b**2+d**2)

	Egs=2*Eminus(2)
	Echan(1)=2*Eminus(2)-Egs          !  |2,-2> x |2,-2>
	Echan(2)=Eminus(1)+Eminus3-Egs    !  |2,-1> x |3,-3> 
	Echan(3)=Eminus(2)+Eplus(2)-Egs   !  |2,-2> x |3,-2> 
	Echan(4)=Eplus(1)+Eminus3-Egs     !  |3,-1> x |3,-3>
	Echan(5)=Eplus(2)+Eplus(2)-Egs    !  |3,-2> x |3,-2> 
c       
c       Psvector(ichan=1,5; iF=1(F=4), iF=2(F=5) )
	Psvector(1,1)=sine(2)**2*Psing(2,2,-2,4)+cosine(2)**2*
     +    Psing(3,3,-2,4)
     +    -sine(2)*cosine(2)*(Psing(2,3,-2,4)+Psing(3,2,-2,4))
	Psvector(1,2)=cosine(2)**2*Psing(3,3,-2,5)
     +    -sine(2)*cosine(2)*(Psing(2,3,-2,5)+Psing(3,2,-2,5))
	Psvector(2,1)=(-sine(1)*(Psing(2,3,-1,4)+Psing(3,2,-3,4))
     +    +cosine(1)*(Psing(3,3,-1,4)+Psing(3,3,-3,4)))/sqrt(2d0)
	Psvector(2,2)=(-sine(1)*(Psing(2,3,-1,5)+Psing(3,2,-3,5))
     +    +cosine(1)*(Psing(3,3,-1,5)+Psing(3,3,-3,5)))/sqrt(2d0)
	Psvector(3,1)=(-sine(2)**2*(Psing(3,2,-2,4)+Psing(2,3,-2,4))
     +    +cosine(2)**2*(Psing(3,2,-2,4)+Psing(2,3,-2,4))
     +    -2*sine(2)*cosine(2)*(Psing(2,2,-2,4)-Psing(3,3,-2,4))
     +    )/sqrt(2d0)
	Psvector(3,2)=(-sine(2)**2*(Psing(3,2,-2,5)+Psing(2,3,-2,5))
     +    +cosine(2)**2*(Psing(3,2,-2,5)+Psing(2,3,-2,5))
     +    -2*sine(2)*cosine(2)*(Psing(2,2,-2,5)-Psing(3,3,-2,5))
     +    )/sqrt(2d0)	 
	Psvector(4,1)=(cosine(1)*(Psing(2,3,-1,4)+Psing(3,2,-3,4))
     +    +sine(1)*(Psing(3,3,-1,4)+Psing(3,3,-3,4)))/sqrt(2d0)
	Psvector(4,2)=(cosine(1)*(Psing(2,3,-1,5)+Psing(3,2,-3,5))
     +    +sine(1)*(Psing(3,3,-1,5)+Psing(3,3,-3,5)))/sqrt(2d0)
	Psvector(5,1)=sine(2)**2*Psing(3,3,-2,4)+cosine(2)**2*
     +    Psing(2,2,-2,4)
     +    +sine(2)*cosine(2)*(Psing(2,3,-2,4)+Psing(3,2,-2,4))
	Psvector(5,2)=sine(2)**2*Psing(3,3,-2,5)+cosine(2)**2*
     +    Psing(2,2,-2,5)
     +    +sine(2)*cosine(2)*(Psing(2,3,-2,5)+Psing(3,2,-2,5))
	
	do i=1,nchan
	   do j=1,nchan
	      Psmatrix(i,j)=Psvector(i,1)*Psvector(j,1)+Psvector(i,2)*
     +	   Psvector(i,2)
	      Ptmatrix(i,j)=-Psmatrix(i,j)
	      if(i.eq.j)Ptmatrix(i,j)=1d0-Psmatrix(i,j)
	   enddo
	enddo	
	
	return
	end
c       
C-----------------------------------------------------------------------
c
	function Psing(f1,f2,m1,F)
	implicit none
	integer i,f1,f2,m1,F,ipow
	double precision Psing,nineJ0,wigner6j,clebsch
	ipow=.5d0 + F-5d0-2d0*f2-f1-5d0/2
	nineJ0=(-1d0)**ipow*wigner6j(5,2*f1,1,2*f2,5,2*F)/sqrt(2d0*
     +    (2d0*F+1d0))
	ipow=5d0+f1+f2+2*F+1
	Psing=(-1d0)**ipow*clebsch(2*f1,2*m1,2*f2,2*(-4-m1),2*F,-8)*
     +    nineJ0*
     +    sqrt((2d0*f1+1)*(2d0*f2+1)*(2d0*F+1))
	  
	return
	end
c
c-------------------------------------------------------------------------

c---------------------------------------------------------------------------c
	SUBROUTINE CalcLambda(LamMat)

	IMPLICIT REAL*8 (f,l,o,a,w,b,c)
	IMPLICIT INTEGER (i)
	REAL*8 sine(2),cosine(2)
	COMMON/angle/sine,cosine
	DIMENSION LamMat(5,5),alpha(2),beta(2)

	alpha(1)=cosine(2)
	alpha(2)=sine(2)
	beta(1)=-sine(2)
	beta(2)=cosine(2)

c       Here are the constant elements
	l11=0d0
	l13=0d0
	l33=0d0

c       Here is the calc for the h's

	DO if1=2,3
	   DO if2=2,3
	      DO if3=2,3

		 IF((if1-if2).LE.0d0) THEN
		    if12low=if2-if1
		 ELSE
		    if12low=if1-if2
		 END IF   

		 if12high=if1+if2

		 IF((if2-if3).LE.0d0) THEN
		    if23low=if3-if2
		 ELSE
		    if23low=if2-if3
		 END IF   

		 if23high=if2+if3

		 DO if12=if12low,if12high
		    DO if23=if23low,if23high

		       IF((if23-if1).LE.(if12-if3)) THEN
			  iFlow=if12-if3
		       ELSE
			  iFlow=if23-if1
		       END IF

		       IF((if23+if1).LE.(if12+if3)) THEN
			  iFhigh=if23+if1
		       ELSE
			  iFhigh=if12+if3
		       END IF

		       DO iFtot=iFlow,iFHigh

			  f1=if1*1d0
			  f2=if2*1d0
			  f3=if3*1d0
			  f12=if12*1d0
			  f23=if23*1d0
			  F=iFtot*1d0

			  l11=l11+R6J(f1,f2,f3,f12,f23,F)*
     +	                      (alpha(if1-1)*alpha(if2-1)*alpha(if3-1))**2
			  l13=l13+R6J(f1,f2,f3,f12,f23,F)*
     +	                        (alpha(if1-1)*alpha(if2-1)*alpha(if3-1))*
     +                          (alpha(if1-1)*beta(if2-1)*alpha(if3-1))
			  l33=l33+R6J(f1,f2,f3,f12,f23,F)*
     +                        (alpha(if1-1)*beta(if2-1)*alpha(if3-1))**2
		       END DO
		    END DO
		 END DO
	      END DO
	   END DO
	END DO

c       Here are the 1x elements from the lambda matrix
	l12=0d0
	l14=0d0
	l15=0d0

c       Here are the 2x elements from the lambda matrix
	l21=l12
	l22=0d0
	l23=0d0
	l24=0d0
	l25=0d0

c       Here are the 3x elements from the lambda matrix
	l31=l13
	l32=l23
	l34=0d0
	l35=0d0

c       Here are the 4x elements from the lambda matrix
	l41=l14
	l42=l24
	l43=l34
	l44=0d0
	l45=0d0
c       Here are the 5x elements from the lambda matrix
	l51=l15
	l52=l25
	l53=l35
	l54=l45
	l55=0d0

c       Now we place these values in the lambda matrix
	LamMat(1,1)=l11
	LamMat(1,2)=l12
	LamMat(1,3)=l13
	LamMat(1,4)=l14
	LamMat(1,5)=l15

	LamMat(2,1)=l21
	LamMat(2,2)=l22
	LamMat(2,3)=l23
	LamMat(2,4)=l24
	LamMat(2,5)=l25

	LamMat(3,1)=l31
	LamMat(3,2)=l32
	LamMat(3,3)=l33
	LamMat(3,4)=l34
	LamMat(3,5)=l35

	LamMat(4,1)=l41
	LamMat(4,2)=l42
	LamMat(4,3)=l43
	LamMat(4,4)=l44
	LamMat(4,5)=l45

	LamMat(5,1)=l51
	LamMat(5,2)=l52
	LamMat(5,3)=l53
	LamMat(5,4)=l54
	LamMat(5,5)=l55

	RETURN
	END
c---------------------------------------------------------------------------c

c---------------------------------------------------------------------------c
	REAL*8 FUNCTION R6j(f1,f2,f3,f12,f23,F)
	
	IMPLICIT REAL*8 (f,r,p,c)
	IMPLICIT INTEGER (j,m)

	fhat(ft)=SQRT(2*ft+1)

	j1=INT(2d0*f1)
	j2=INT(2d0*f2)
	j3=INT(2d0*f12)
	j4=INT(2d0*f3)
	j5=INT(2d0*F)
	j6=INT(2d0*f23)
	m2=INT(2d0*-2d0)
	m4=INT(2d0*-4d0)
	m6=INT(2d0*-6d0)

	phz=(-1)**(f1+f2+f3+F)
	CG1=Clebsch(j1,m2,j2,m2,j3,m4)
	CG2=Clebsch(j3,m4,j4,m2,j5,m6)
	CG3=Clebsch(j4,m2,j2,m2,j6,m4)
	CG4=Clebsch(j6,m4,j1,m2,j5,m6)
	SIXJ=phz*fhat(f12)*fhat(f23)*
     +        WIGNER6J(j1,j2,j3,j4,j5,j6)
	R6j=CG1*CG2*CG3*CG4*SIXJ

	RETURN
	END
c---------------------------------------------------------------------------c
