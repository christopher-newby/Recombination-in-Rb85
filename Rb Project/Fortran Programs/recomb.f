       program recomb
c     TO RUN
c     NEEDS:
c     lapack and some blas.
c     also needs HyperSub.o and angmom.o

c
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax)
	  common/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
          COMMON/Bpoint/iB
          COMMON/Blambda/flam2
          COMMON/HSang/fHSang

c  Parameters of the problem...
	  Nq=1024+1     ! number of points on q-grid
	  qmin=1.d-8
	  qmax=20d0
	  Nphi=64+1    ! number of points on phi-grid

          fHSang=25d0           ! HS angle is pi/fHSang

          DO i=8,8
             iB=i+1
c     18 is the highest we can go for the numbers as of 3-10-2010 for the dipole form factors
c     1 is beginning of the positive scattering lengths for these numbers
c     ALL THESE ARE FOR icalc=5!!!!!
c     20 and 21 for 157.25 and 157.75 respectively
c     22-30 for 157.1-157.9 including 157.5

             Lwave=0            ! partial wave, currently implemented for icalc = 2,3,40 only
             Nk=300              ! number of k-values to print out 
             pkmax=5d-2        ! maximum k-value
	  
c  Model choices...
c	  m Etau = m Ek-3*q**2/4 = 3 (k^2-q^2)/4 - gamma^2
c
c             icalc=0            ! ! LO:   g[q]=1, tau=tau0=1/(-gamma + sqrt(-m Etau))
c             icalc=1            !(10)	! NLO  - about zero:  g[q] = 1, tau->tau0*(1-r0/2 3(q^2-k^2)/4*tau0)
c             icalc=2            !(11)	! NNLO - about zero:  g[q] = 1, tau->tau0*(1-r0/2 3(q^2-k^2)/4*tau0 + (r0/2 3(q^2-k^2)/4*tau0)^2 )
c try not to use	  icalc=3 !(21)  ! NLO  - about pole: tau-->tau0*(1+r0*(gamma+Sqrt[])/2)
c try not to use	  icalc=4 !(22)  ! NNLO - about pole: tau-->tau0*(1+r0*(gamma+Sqrt[])/2+(r0*(gamma+Sqrt[])/2)**2)
             icalc=5            !(2)	! Rank 1 - dipole form factor: V[p,q]=lambda g[p] g[q];  g[q] = 1/(1+q^2/beta^2)
c             icalc=7            !		! delta-shell (DO NOT USE!!!!)

c  2-body scattering choices...
c             iPot=1             !  He-He  (LM2M2)
c             iPot=2             !  Rb-Rb  (triplet -Jose D'Incao)
c             iPot=3             !  Rb-Rb  (singlet- Jose)
c             iPot=4             !  Rb-Rb  scattering length for fixed B-field
             iPot=5             ! Rb-Rb scattering length for fixed B from Newby
                                ! Data files start from 155G and go up in increments of 0.5G to 165G

c
c  Atom selection
c             iAtom=1            ! He  (Needs qmax=20)
             iAtom=2            ! Rb85
c  Set up...	  
             call setup(qmin,qmax)
c             call diagnosis
             call PrintRecomb(pkmax,Nk,qmin,qmax,iB)

          END DO

      call CPUTime

	  stop
	  end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine PrintRecomb(pkmax,Nk,qmin,qmax,iB)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax)
	  character*10 Atom_label(2)
          character*5 Bfieldhold
          character*1 icalchold
          character*4 HSang
          character*40 fileloc1,fileloc2,fileloc3
	  common/q_grid/dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/hyper/ vlight,fmu_Bohr,a_Bohr,fmu_Elec,Hartree,Ehf
          COMMON/HSang/fHSang
	  data cI/(0d0,1d0)/,fk_Boltz/8.617343d-5/  ! eV/K
	  data Atom_label/' He-He    ',' Rb85-Rb85'/
c type of output
c	  iK3=0  !  prints phase shifts for linear span of momenta
	  iK3=1  !  prints K3 factors for logrithmic span of energies

          Bf=.5d0*(iB-1)+155d0  !1-19
c          Bf=157.25d0           !20
c          Bf=157.75d0           !21
c          Bf=(iB-22)*.1d0+157.0d0 !22-30
          write(Bfieldhold,1000)Bf
          write(icalchold,1001)icalc
          write(HSang,1002)fHSang
          IF(iAtom.EQ.1) THEN
             fileloc1="Data/He/K3vsE_icalc"//icalchold//"_He.dat"
             fileloc2="Data/He/Delta_icalc"//icalchold//"_He.dat"
             fileloc3="Data/He/Cot_icalc"//icalchold//"_He.dat"
          ELSE
             fileloc1="Data/Real/B"//Bfieldhold//".dat"
c             fileloc1="Data/icalc"//icalchold//"/K3vsE_B"//Bfieldhold//"_1.dat"
             fileloc2="Data/icalc"//icalchold//"/Delta_B"//Bfieldhold//"_1.dat"
             fileloc3="Data/icalc"//icalchold//"/Cot_B"//Bfieldhold//"_1.dat"
          END IF
 1000     FORMAT(f5.1)
 1001     FORMAT(i1)
 1002     FORMAT(f4.1)
	   	  
	  open(unit=22,status='unknown',file=fileloc1)
c	  open(unit=23,status='unknown',file=fileloc2)
c	  open(unit=24,status='unknown',file=fileloc3)

	  write(6,*)' ' 
	  write(6,*)'             ',Atom_label(iAtom),' Scattering: L = ',Lwave
	  write(6,*)' '

c  Header for printing out phase shifts...
	  if(iK3.eq.0)then
	    write(6,*)'      k(eV)                 k Cot(delta)                      delta                           1 / (k f)                        K3 (cm^6/s)              Ebreakup (mK)'
	    write(6,*)'  =========        ==========================       ===========================      ============================      =========================     =================='

c     pkmin here!!!!!!!!!!!!!!!!!!!!!!!!!!
	    pkmin=1.3d-2
	    dpk=(pkmax-pkmin)/(Nk-1)
	  endif

c  Header for printing out recombination rates...
	  if(iK3.eq.1)then
	     EBrkUpmin=1d-10*fk_Boltz/hbarc  ! Start at 1 nK
	     EBrkUpmax=1d-3*fk_Boltz/hbarc  ! End at 1 milli-K  (Also good to end at 1 microK)
	     dlogEBrkUp=log(EBrkUpmax/EBrkUpmin)/(Nk-1)
	     write(6,*)'  E(micro-K)       K3 (cm^6/s)                 1 / (k f)'
	     write(6,*)'  =========      ==============       ==========================='	
	  endif
	  
	  ik_last=0
	  do i=1,Nk
		if(iK3.eq.0)then
		    pk=(i-1)*dpk+pkmin
		else
			EBrkUp=EBrkUpmin*exp((i-1)*dlogEBrkUp)
			pk=sqrt(4*(fm*EBrkUp+gamma**2)/3)
		endif
			ik=int(log(pk/qlist(1))/dt+1)
			if(ik.eq.ik_last)goto 900  ! skip over any repeat values
			ik_last=ik
			EBreakUp=(3*qlist(ik)**2/4-gamma**2)/fm
			
	        call ThreeBody(ik,cfscatt)
			cpCot=1d0/cfscatt + cI*qlist(ik)
			cdel=-cI*log(1d0+2*cI*qlist(ik)*cfscatt)/2d0*180/pi
			
c  save first two values to extrapolate scattering length to 0
		    if(i.eq.1)then
			   q1=qlist(ik)
			   Cot1=1d0/(q1*cfscatt) + cI
			   a1=-1d0/(q1**(2*Lwave+1)*Real(Cot1))
			endif
			if(i.eq.2)then
			   q2=qlist(ik)
			   Cot2=1d0/(q2*cfscatt) + cI
			   a2=-1d0/(q2**(2*Lwave+1)*Real(Cot2))
			endif
			
			fmu=fm/sqrt(3d0)
			fkBrkUp=sqrt(2*fmu*EBreakUp)
			recomfactor=0.75d0*(2*Lwave+1)*192*pi**2/(fmu*fkBrkUp**4)  *vlight*1d-42  ! .75? cm-->lab ?			
			EBrkUpmK=EBreakUp*hbarc/fk_Boltz*1d6
			RecomRate=recomfactor*(1d0-exp(-4*dimag(cdel)*pi/180)) 
			
			if((real(cdel).gt.0d0).and.(Lwave.eq.0))cdel=cdel-180d0
			if((real(cdel).lt.0d0).and.(Lwave.eq.1))cdel=cdel+180d0
			if(iK3.eq.0)then
			   if(EBrkUpmK.le.0d0)RecomRate=0d0
                           write(6,101)qlist(ik)*hbarc,cpCot,cdel,1d0/cfscatt/qlist(ik),RecomRate,EBrkUpmK
			else
			   if(EBrkUpmK.le.0d0)RecomRate=0d0
		       write(6,103)EBrkUpmK,RecomRate,1d0/cfscatt/qlist(ik),cdel
			endif
                        IF(RecomRate.GT.0) THEN
                           write(22,102)EBrkUpmK,RecomRate
                        END IF
c                        write(23,102)EBrkUpmK,DIMAG(cdel)
c                        write(24,102)EBrkUpmK,DIMAG(cpCot)

 900	    continue
	  enddo
	  if(iK3.eq.0)then
	     as_0=(q1**2-q2**2)/(q2**2*q1*Cot1-q1**2*q2*Cot2)
	     print*,' '
	     print*,'  Scattering length(volume) = ',as_0,' nm^(',2*Lwave+1,')'
	  endif

 	  close(22)
c          close(23)
c          close(24)
	  	  
 100  format(f12.5,2x,2f18.8)
 101  format(f12.5,2x,2f16.8,2x,2f16.8,2x,2f16.8,6X,e18.10,8X,f16.8)
 102  format(f12.5,2x,e18.10)
 103  format(f12.5,2x,e18.10,2x,2f16.9,2x,2f16.8)
 104  format(f12.5,2x,e18.10,2x,e18.10)
	  return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setup(qmin,qmax)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax,Ngauss=32,nBmax=500)
      common/data/ pi,hbarc,fm,gammaExp
	  common/hyper/ vlight,fmu_Bohr,a_Bohr,fmu_Elec,Hartree,Ehf
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/ iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/q_grid/ dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
	  common/gauss/Xgauss(Ngauss),Wgauss(Ngauss)
	  common/ThreeBodyForce/H0,BigLam
	  common/Bgrid/ dB,Bfield(nBmax),nB
          COMMON/Bpoint/iB
          COMMON/Blambda/flam2
          COMMON/HSang/fHSang
	  character*22 ID(9),IDPot(6),Atoms(2)

	  data    ID/'   LO  (g=1)          ','   NLO (p=0)          ','   NNLO (p=0)         ',
     &           '   NLO (p=gamma)      ','   NNLO (p=gamma)     ','   Rank 1 (dipole)    ',
     &           '   Rank 2 (dipole)    ','   Delta Shell        ','   undefined          '/
	  data IDPot/'   LM2M2 (He-He)      ','   DIncao (RbRb-sing) ','   DIncao (RbRb-trip) ',
     &           '   DIncao (as@B fixed)','   NEWBY Rb85-Rb85    ','  DIncao (Rb-Rb)      '/
	  data Atoms/'   He - He            ','   Rb85 - Rb85        '/
          DIMENSION a0list(500),r0list(500),v4list(500),Bhold(500)
          DIMENSION flhold(5,5),Echan(5),Ps(5,5),Pt(5,5)

c Universal constants...
      pi=2.d0*asin(1.d0)
      hbarc=197.3269631d0  !  2006 CODATA from NIST website
	  amu=931.494028d6
          fkB=1.3806504d-23  ! Boltzmann's constant in J/K
          eCharge=1.602176487d-19  ! Coulombs
	  fk_Boltz=fkB/eCharge  ! eV/K	    
	  vlight=2.99792458d17
	  hPlanck=2*pi*hbarc/vlight
	  fmu_Bohr=5.7883817555d-5
	  fmu_Elec= 1.00115965218111d0*fmu_Bohr  !  2008 CODATA from NIST
	  a_Bohr=0.052917720859d0
	  alpha=7.2973525376d-03  ! Jose D'Incao
	  Hartree=alpha*hbarc/a_Bohr
c Atom-specific constants...
	  if(iAtom.eq.1)then
		A_atom=4.002602d0
		BE_2K=1.3034834d-3 ! 2-body BE in K (this value is for He-He)
		BE_2=-BE_2K*fk_Boltz/hbarc  !2-body BE in nm^-1
		BE_3K=2.28d-3  ! 3-body BE in K
		Ek=-BE_3K*fk_Boltz/hbarc  ! 3-body BE in nm^-1
		Ehf=0d0
	  endif
	  
	  if(iAtom.eq.2)then
		A_atom=84.911789738d0
		BE_2K=1.3034834d-3 ! 2-body BE in K (this value is for He-He)
		BE_2=-BE_2K*fk_Boltz/hbarc  !2-body BE in nm^-1
		Ehf=3.035732439d9*hPlanck  ! eV
	  endif
	  
	  fm=A_atom*amu/hbarc  ! atomic mass in nm^-1
	  gammaExp=sqrt(-fm*BE_2)  ! gamma in nm^-1
	  Bhf=1d4*Ehf/fmu_elec  ! Gauss

c
c  B-grid and B-value (without looping)
c
          Bmin=155d0
          Bmax=165d0
          nB=20
          dB=(Bmax-Bmin)/nB
          DO i=1,nB
             Bfield(i)=(i-1)*dB+Bmin
          END DO

          Bfield(nB+1)=157.25d0
          Bfield(nB+2)=157.75d0

          DO i=1,10
             Bfield(i+nB+2)=157d0+0.1d0*i
          END DO

c  gauss grid
	  call gaussleg(-1d0,1d0,Xgauss,Wgauss,Ngauss)
c

c  2-body scattering length and range parameters
	  if(iAtom.eq.1)then  ! LM2M2: 
             as=10.036600731540563d0 ! JAM fit April 16, 2008 fit to pCot(delta)
             r0=0.7283732465957433d0
             v4=0.008678647255100813d0
	  endif
	  
	  if(iAtom.eq.2) then   ! defaults to Triplet potential
		   as=-19.71763203438355d0  !  JAM fit March 26, 2009 fit to pCot(delta)
		   r0=11.303238103698735d0  !  RbRb_scattering_3.26.nb
		   v4=133.69226101569353d0
	     if(ipot.eq.2)then  ! Jose D'Incao: Rb-Rb Triplet
		   as=-19.71763203438355d0  !  JAM fit March 26, 2009 fit to pCot(delta)
		   r0=11.303238103698735d0  !  RbRb_scattering_3.26.nb
		   v4=133.69226101569353d0
	     endif
	     if(ipot.eq.3)then  ! Jose D'Incao: Rb-Rb Singlet
		   as=126.27305582869954d0  ! JAM fit March 26, 2009 fit to pCot(delta)
		   r0=17.933315071786954d0  ! RbRb_scattering_3.26.nb
		   v4=52.17034593703047d0
	     endif
		 if(ipot.eq.4)then   ! Scattering length fixed by Rb-Rb code for fixed B-field
c			as=-0.13215410d4*a_Bohr  ! B=150 G
c			as=0.24954132d6*a_Bohr   ! B=155.4 G
c			as=0.74205232d4*a_Bohr	 ! B=156 G
c 			as=0.58761518d3*a_Bohr   ! B=160 G
			as=4.43772d0  ! 0.84268977d2*a_Bohr	 ! B=165 G
			r0=-4.68892d0  !
			v4=0d0
		 endif


c-----------------------NEWBY's params------------------------------
                 IF(ipot.EQ.5) THEN
c     Data from CalcParams.nb on 2-26-2010

                    OPEN(unit=13,status='unknown',file='Params.dat')

                    DO i=1,nB
                       READ(13,*) Bhold(i),fhold,a0list(i),r0list(i)
                    END DO

                    CLOSE(13)

                    a0list(nB+1)=1958.3276537d0
                    r0list(nB+1)=10.8340347d0
                    a0list(nB+2)=1491.1076223d0
                    r0list(nB+2)=10.5032384d0

c                  --------------------------------------
                    a0list(23)=2146.8208146d0
                    r0list(23)=10.9169065d0
                    a0list(24)=2012.7168600d0
                    r0list(24)=10.9821916d0
                    a0list(25)=1901.5289253d0
                    r0list(25)=10.8017312d0
                    a0list(26)=1795.7340899d0
                    r0list(26)=10.7429498d0
                    a0list(27)=1699.2638359d0
                    r0list(27)=10.6734233d0
                    a0list(28)=1610.7848356d0
                    r0list(28)=10.6055348d0
                    a0list(29)=1530.0444138d0
                    r0list(29)=10.5071223d0
                    a0list(30)=1454.3430645d0
                    r0list(30)=10.4638492d0
                    a0list(31)=1384.7215258d0
                    r0list(31)=10.3958969d0

                    as=a0list(iB)*a_Bohr
                    r0=r0list(iB)
                    v4=0d0

                 END IF
c-------------------------------------------------------------------

	  endif

          IF(iAtom.EQ.1) THEN
             flam2=1d0
          ELSE
             CALL HyperSub(Bfield(iB),flhold,Echan,Ps,Pt)
             flam2=flhold(1,1)
          END IF

c  3-body parameters: H0 determined by fit to 3-body binding energy in boson_bnd.f
	  if(icalc.eq.0)then  !  ERE-LO
		BigLam=qmax
		H0=0d0
		if(iAtom.eq.1)H0=-3.422841902725d0 !  LO
	  endif
	  if(icalc.eq.1)then  !  ERE-NLO
		BigLam=qmax
		H0=0d0
		if(iAtom.eq.1)H0 = -.695989355339d0  !  NLO !old H0=-.69598932d0
	  endif
	  if(icalc.eq.2)then  !  ERE-NNLO
		BigLam=qmax
		H0=0d0
		if(iAtom.eq.1)H0=-.0887321325d0
	  endif	  
	  if(icalc.eq.3)then  !  ERE-NLO-pole
		BigLam=qmax
		H0=0d0
		if(iAtom.eq.1)H0=-.695536691799d0
	  endif
	  if(icalc.eq.4)then  !  ERE-NNLO-pole
		   BigLam=qmax
		H0=0d0
		   if(iAtom.eq.1)H0=-.087197503787d0
	  endif

c	Form Factor parameters
	  if(icalc.eq.5)then   !  Rank 1
		BigLam=qmax
		H0=0d0
			if(iAtom.eq.1)then  ! FF parameters fit which in turn give these ERE parameters
                           as=10.036600731540563d0 ! April 16, 2008 fit of dipole FF to cCot(delta)
                           r0=0.7575939136251458d0
                           v4=0.008486943422048065d0
                           H0 = 0d0!  .345238941775d0 !  Rank 1
			endif
	  endif

	  if(icalc.eq.7)then  ! delta shell
		BigLam=qmax
		H0=0d0
		if(ipot.eq.1)then  ! delta shell parameters fit which in turn give these ERE parameters
		   H0 =-.144977424522d0
		   a_shell=0.5619497440709951d0
		   gamma=0.10351768834489634d0
		   flam=-8*pi*a_shell**2*gamma/fm/(1-exp(-2*a_shell*gamma))
		   as=a_shell*fm*flam/(4*pi*a_shell+fm*flam)
		   r0=2*a_shell*(1-4*a_shell*pi/(fm*flam))/3
		   v4=a_shell**3/45*(1-12*a_shell*pi/(fm*flam))
		   beta=a_shell ! for printing out
		endif
	  endif

c  exponential q-grid
	  dt=log(qmax/qmin)/(nq-1)
	  do i=1,nq
		 qlist(i)=qmin*exp((i-1)*dt)
	  enddo
c  angle grid for Hetherington & Schick integration
	  phi0=pi/fHSang
	  dphi=phi0/(nphi-1)
	  do iphi=1,nphi
	     phi(iphi)=phi0-(iphi-1)*dphi
	  enddo
c
c  integration weights
c	  int_type=1  ! Euler
c	  int_type=2  ! Simpson's 1/3
c	  int_type=3  ! Bode 4-point
	  int_type=4  ! Bode 8-point
	  call weights(int_type)

	  write(6,*)' '
	  write(6,*)'          Atom-Atom Recombination Scattering'
	  write(6,*)' '

      write(6,100) fm,BE_2*hbarc,gammaExp
	  write(6,101) Nq,Nphi,qlist(nq),dt
c  Two body parameters...

	  call TwoBody(iID)
	  alpha=4*pi/(fm*flam)
	  
	  write(6,*)' '
	  write(6,*)'                           Two-body Parameters for ',IDPot(ipot)
	  write(6,*)'     flam           beta          alpha          gamma            as            r0             v4            r0_eff        WFnorm*'
	  write(6,*)'   ===========   ============  =============  =============  =============  =============  =============  =============  ============='
	  write(6,102)flam,beta,alpha,gamma,as,r0,v4,r0-4*v4*gamma**2,WFnorm
	  write(6,*)'                                                                                                     (* WFnorm = 1 uses T-matrix form)'
	  write(6,*)' '
	  write(6,*)'  icalc = ',icalc,'  --> ',ID(iID)
	  write(6,*)'  iPot  = ',iPot, '  --> ',IDPot(iPot)
	  write(6,*)'  iAtom = ',iAtom,'  --> ',Atoms(iAtom)
	  write(6,*)'  Ehf = ',Ehf,' eV'
          write(6,*)' '
          write(6,*)'  Bfield = ',Bfield(iB), ' gauss'
      return
 100	format(' m = ',e15.9,'  nm^-1,  BE = ',f14.12,' eV,  gamma(Exp) = ',f14.12,' nm^-1')
 101	format(' Nq = ',i4,', Nphi = ',i4,',  qmax = ',f12.6,', dt = ',f14.12)
 102	format(8(1x,f14.8),1x,e14.8)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine TwoBody(iID)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/h_matrix/cH(2,2),ctau_det
	  common/q_grid/ dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi

	  data gammaRank/0.1d0/ !  <----- LM2M2
	  
c  Rank 1: {flam,beta} fit to {as,r0}, e.g. {as=5.403,r0=1.7494} gives {flam=-5.05005,beta=1.41601}
	  flam=2*pi*as/(fm*(as-2*r0))*(3*as-4*r0-sqrt(9*as**2-16*r0*as))
	  beta=(3+sqrt(9d0-16*r0/as))/(2*r0)
	  nrank=1

	  if(icalc.eq.0)then
	     gamma=(1.d0-sqrt(1.d0-2.d0*r0/as))/r0    !  LO
		 WFnorm=fm**2/(8*pi*gamma)
		 iID=1
		 return
	  else
	  	     gammaERE=(1.d0-sqrt(1.d0-2.d0*r0/as))/r0   !  start at LO
		 do i=1,1000                             ! iterate to get gamma exactly including the v4 term
		    gammanew=1/as+r0*gammaERE**2/2-v4*gammaERE**4
			if(abs(gammanew-gammaERE).lt.1d-14)goto 10
			gammaERE=gammanew
		 enddo
		 write(6,*)' STOP: No convergence in gamma (ERE)'
		 stop
 10		 r0_eff=r0-4*v4*gammaERE**2
	     r0gam=r0_eff*gammaERE
	  endif
	  
	  if(icalc.eq.1)then  !  NLO in range expansion
	     gamma=(1.d0-sqrt(1.d0-2.d0*r0/as))/r0   !  ERE 
		 WFnorm=fm**2/(8*pi*gamma*(1d0+r0*gamma))
		 iID=2
		 return
	  endif

	  if(icalc.eq.2)then  !  NNLO in range expansion
	     gamma=(1.d0-sqrt(1.d0-2.d0*r0/as))/r0    !  ERE
		 WFnorm=fm**2/(8*pi*gamma*(1d0+r0*gamma+(r0*gamma)**2)) 
		 iID=3
		 return
	  endif
	  
	  if(icalc.eq.3)then  !  ERE-NLO - expansion about pole
	     gamma=gammaERE
		 WFnorm=fm**2/(8*pi*gamma)/(1d0+r0gam)
		 iID=4
		 return
	  endif

	  if(icalc.eq.4)then  !  ERE-NNLO - expansion about pole
	     gamma=gammaERE
		 WFnorm=fm**2/(8*pi*gamma)/(1d0 + r0gam + r0gam**2)
		 iID=5
		 return
	  endif
	  
	  if(icalc.eq.5)then
		 gamma=beta*(sqrt(-fm*beta*flam/(8*pi))-1d0)
		 WFnorm=1d0
		 iID=6
		 return
	  endif

	  if(icalc.eq.7)then  ! delta shell
		   WFnorm=1d0
		   beta=a_shell ! for printing out only
		   iID=8
		   return
	  endif

	  write(6,*)'  Invalid calculation option, icalc = ',icalc
	  stop
	  end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ThreeBody(ik,cfscatt)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/q_grid/ dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
	  common/h_matrix/cH(2,2),ctau_det
	  dimension cKern(nmax,nmax)
	  dimension cXsoln(nmax),Ipivot(nmax)
	  data cI/(0d0,1d0)/,c1/(1d0,0d0)/
c  setup factors
	  eps=1d-10
	  fac=1d0/(2*pi**2)
	  ck=qlist(ik)
	  cEk=3*ck**2/(4*fm)-gamma**2/fm + cI*eps/fm

c Load the kernel for both p and q complex
	  ndim=nq+nphi
	  cphz0=exp(-cI*phi(1))
	  do ip=1,ndim
	     if(ip.le.Nq)then
	        cp=qlist(ip)*cphz0
		 else
		    cp=qlist(Nq)*exp(-cI*phi(ip-Nq))
		 endif
	  do iq=1,ndim
	     if(iq.le.Nq)then
	        cq=qlist(iq)*cphz0
			cdq_dphi=cq*dt*wt(iq)
		 else
		    cq=qlist(Nq)*exp(-cI*phi(iq-Nq))
			cdq_dphi=cq*cI*dphi*wt(iq)
		 endif
		 cEtau=cEk-3*cq**2/(4*fm)
		do ir=1,nrank
			ipstart=(ir-1)*ndim
		do is=1,nrank
			iqstart=(is-1)*ndim
			cKern(ipstart+ip,iqstart+iq)=-cKernel(ir,is,cp,cq,cEk)*cq**2*cdq_dphi
		enddo  ! is
		enddo  ! ir

	  enddo  ! iq
		 do ir=1,nrank
			ipstart=(ir-1)*ndim
			cKern(ipstart+ip,ipstart+ip)=c1+cKern(ipstart+ip,ipstart+ip)
			cXsoln(ipstart+ip)=cBorn(ir,1,cp,ck,cEk) ! cXsoln fed into the linear equation solver; on exit, is solution
c			if((ip.lt.2))write(6,666)ir,cp,ck,cBorn(ir,1,cp,ck,cEk)
		 enddo

	  enddo  ! ip
c 666	  format('ir, cp,ck,cBorn =',i5,2x,2e18.12,2x,2e18.12,2x,2e18.12)

c  general linear equation solver gives the Xmatrix for q complex

	  call ZGESV(nrank*ndim,1,cKern,nmax,IPIVOT,cXsoln,nmax,INFO)
	  
		do is=1,nrank
			iqstart=(is-1)*ndim		 
		do iq=1,ndim
	     if(iq.le.Nq)then
	        cq=qlist(iq)*cphz0
			cdq_dphi=cq*dt*wt(iq)
		 else
		    cq=qlist(Nq)*exp(-cI*phi(iq-Nq))
			cdq_dphi=cq*cI*dphi*wt(iq)
		 endif
		 cEtau=cEk-3*cq**2/(4*fm)
			cKern(ik,iqstart+iq)=cKernel(1,is,ck,cq,cEk)*cq**2*cdq_dphi
		enddo  ! iq

		enddo  ! is

c  operate with the LS equation to get Tmatrix for q=k real

	  cTmatrix=cBorn(1,1,ck,ck,cEk) ! Born term
	  do ir=1,nrank
		 iqstart=(ir-1)*ndim
	  do iq=1,ndim
	     cTmatrix=cTmatrix+cKern(ik,iqstart+iq)*cXsoln(iqstart+iq)
	  enddo
	  enddo
	  
	  cfscatt=-fm*cTmatrix/(3*pi*WFnorm)

	  return
	  end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function cKernel(ir,is,cp,cq,cEk)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (ngauss=32)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/gauss/Xgauss(Ngauss),Wgauss(Ngauss)
	  common/ThreeBodyForce/H0,BigLam
          COMMON/Blambda/flam2

	  
	   cEtau=cEk-3*cq**2/(4*fm)  ! cEk = 3 k^2/4m - gamma^2/m

	   if(icalc.eq.5)then
		  beta2=beta**2
		  cp2=cp**2
		  cq2=cq**2
		  cpq=cp*cq
		  cangle=(0d0,0d0)
		  do i=1,Ngauss
		    cpqx=cpq*Xgauss(i)
			cden1=beta2+cp2+cq2/4+cpqx
			cden2=beta2+cq2+cp2/4+cpqx
			cden3=cp2+cq2-fm*cEk+cpqx
			cden=cden1*cden2*cden3
			cangle=cangle+Wgauss(i)/cden*FLeg(Xgauss(i),Lwave)
		  enddo
		  cKernel=4/(beta*pi)*(1+gamma/beta)**2*(beta**4*cangle+2*H0/BigLam**2)
     &            /(1-((beta+gamma)/(beta+sqrt(-fm*cEtau)))**2)*flam2
		  return
	   endif

c  .... all other cases:

		cKernel=dcmplx(0d0,0d0)
		do it=1,nrank
			cKernel=cKernel+2*cZfun(ir,it,cp,cq,cEk)*ctaufun(it,is,cEtau)/(2*pi**2)
		enddo
	
	  return
	  end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function cBorn(ir,is,cp,cq,cEk)  !  Born term used only for Harms rank 2 calculation, icalc=6
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (ngauss=32)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/gauss/Xgauss(Ngauss),Wgauss(Ngauss)
	  common/ThreeBodyForce/H0,BigLam
          COMMON/Blambda/flam2

		  
	    if(icalc.eq.5)then
	      cBorn=dcmplx(0d0,0d0)
		  cp2=cp**2
		  cq2=cq**2
		  cpq=cp*cq
	      cangle=(0d0,0d0)
		  beta2=beta**2
		  do i=1,Ngauss
		    cpqx=cpq*Xgauss(i)
			cden1=beta2+cp2+cq2/4+cpqx
			cden2=beta2+cq2+cp2/4+cpqx
			cden3=cp2+cq2-fm*cEk+cpqx
			cden=cden1*cden2*cden3
			cangle=cangle+Wgauss(i)/cden*FLeg(Xgauss(i),Lwave)
		  enddo
c		  cBorn=-8*pi*fIS*gamma/fm*(1+gamma/beta)**3*beta**4*cangle  ! from quartet.f
		  ZB=gamma*(beta+gamma)/fm
		  fnorm=8*pi*(1+gamma/beta)**2/(beta*fm)
		  cBorn=ZB*fnorm*2*(-fm)*(beta**4*cangle/2+H0/BigLam**2)*flam2
		  return
	    endif

		 if(icalc.eq.7)then
			agam=a_shell*gamma
			ZB=(exp(2*agam)-1)*(2*gamma**2)/fm/(exp(2*agam)-2*agam-1)
			cBorn=2*ZB*cZfun(ir,1,cp,cq,cEk)
			return
		 endif

c  All other cases...

		cBorn=2*cZfun(ir,1,cp,cq,cEk)
		
	  return
	  end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function cZfun(ir,is,cp,cq,cEk)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (Ngauss=32,nBmax=500)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/gauss/Xgauss(Ngauss),Wgauss(Ngauss)
	  common/h_matrix/cH(2,2),ctau_det
	  common/ThreeBodyForce/H0,BigLam
c	  DATA flamH0/1D0,0.159192027205D0/  !  Rank 2
          COMMON/Blambda/flam2

		  cp2=cp**2
		  cq2=cq**2
		  cpq=cp*cq
	      cangle=(0d0,0d0)
	  
	  if(icalc.eq.5)then
		  beta2=beta**2
		  do i=1,Ngauss
		    cpqx=cpq*Xgauss(i)
			cden1=beta2+cp2+cq2/4+cpqx
			cden2=beta2+cq2+cp2/4+cpqx
			cden3=cp2+cq2-fm*cEk+cpqx
			cden=cden1*cden2*cden3
			cangle=cangle+Wgauss(i)/cden*FLeg(Xgauss(i),Lwave)
		  enddo
		  cZfun=-fm*beta**4*(cangle/2+H0/BigLam**2)*flam2
		  return
	  endif



	  if(icalc.eq.7)then  ! delta-shell
		  do i=1,Ngauss
		    cpqx=cp*cq*Xgauss(i)
			capq4=a_shell*sqrt(cp2+cq2/4+cpqx)
			caqp4=a_shell*sqrt(cp2/4+cq2+cpqx)
			cg1=sin(capq4)/capq4
			cg2=sin(caqp4)/caqp4
			cden=cp2+cq2-fm*cEk+cpqx
			cangle=cangle+Wgauss(i)*cg1*cg2/cden
		  enddo
		  fnorm=(8*pi*a_shell**2*gamma/fm)/(1-exp(-2*a_shell*gamma))
		  cZfun=-fm*(cangle/2+H0/BigLam**2)*fnorm*flam2
		  return
	  endif

c  All other cases:

			do i=1,Ngauss
				cpqx=cpq*Xgauss(i)
				cG0=1d0/(cp2+cq2-fm*cEk+cpqx)
				cangle=cangle+Wgauss(i)*cG0*FLeg(Xgauss(i),Lwave)
			enddo
			cZfun=-fm*(cangle/2+H0/BigLam**2)*flam2

	  return
	  end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function ctaufun(ir,is,cEtau)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/h_matrix/cH(2,2),ctau_det

		   cIp=-sqrt(-fm*cEtau)
		   ctau0=1d0/(-gamma-cIp)

	  if(icalc.eq.0)then   ! LO
	     ctaufun=(-4*pi/fm)*ctau0
		 return
	  endif
	  	  
	  if(icalc.eq.1)then   ! NLO - ERE
		   carg=-fm*cEtau-gamma**2
		   ctaufun=(-4*pi/fm)*ctau0*(1d0+r0*carg*ctau0/2)
		   return
	  endif

	  if(icalc.eq.2)then   ! NNLO - ERE
		   carg=-fm*cEtau-gamma**2
		   ctaufun=(-4*pi/fm)*ctau0*(1d0+r0*carg*ctau0/2+(r0*carg*ctau0/2)**2)  !+v4*(carg*ctau0)**2)
		   return
	  endif
	  
	  if(icalc.eq.3)then   ! NLO - ERE - expanded about the pole
		   r0_eff=r0-4*v4*gamma**2
		   cr0rt= r0_eff*(sqrt(-fm*cEtau)+gamma)/2 
		   ctaufun=(-4*pi/fm)*ctau0*(1d0 + cr0rt)
		   return
	  endif

	  if(icalc.eq.4)then   ! NNLO - ERE - expanded about the pole
		   r0_eff=r0-4*v4*gamma**2
		   cr0rt= r0_eff*(sqrt(-fm*cEtau)+gamma)/2 
		   ctaufun=(-4*pi/fm)*ctau0*(1d0 + cr0rt + cr0rt**2)
		   return
	  endif
	  
	  if(icalc.eq.5)then  ! Rank 1 with dipole form factor
		 ctaufun=-8*pi*(1+gamma/beta)**2/(fm*beta)/(1-((beta+gamma)/(beta+sqrt(-fm*cEtau)))**2  )
		 print*,' Error: ctaufun should not be called for FF option'
	     return
	  endif
	  	  	  


	  if(icalc.eq.7) then
		  crt=sqrt(-fm*cEtau)
		  cH_shell=-gamma*(1-exp(-2*a_shell*crt))/(1-exp(-2*a_shell*gamma))/crt
	      ctaufun=-1d0/(1d0+cH_shell)
	      return
	  endif
	  
	  write(6,*)' Error in ctaufun: STOP'
	  stop
	  end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine weights(int)
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax)
	  common/q_grid/ dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
	  dimension bode4(4),bode8(8)
	  data bode4/14d0,32d0,12d0,32d0/
	  data bode8/1978d0,5888d0,-928d0,10496d0,-4540d0,10496d0,-928d0,5888d0/
	  ndim=nq+nphi
c  Euler is the default weight
	      do i=1,ndim
		    wt(i)=1d0
		  enddo
		    wt(1)=.5d0
			wt(Nq)=.5d0
			wt(Nq+1)=.5d0
			wt(ndim)=.5d0
c  Simpson's 1/8
	  if(int.eq.2)then
	      do i=1,Nq+1,2
		    wt(i)=4d0/3
			wt(i+1)=2d0/3
		  enddo
		  wt(1)=1d0/3
		  wt(Nq)=1d0/3
	      do i=Nq+2,ndim-1,2
		    wt(i)=4d0/3
			wt(i+1)=2d0/3
		  enddo
		  wt(Nq+1)=1d0/3
		  wt(ndim)=1d0/3
	  endif
c  Bode's 4-point	  
	  if(int.eq.3)then
	      do i=1,Nq+Nphi
		    wt(i)=2*bode4(mod(i-1,4)+1)/45
		  enddo
		  wt(1)=wt(1)/2
		  wt(Nq+Nphi)=wt(Nq+Nphi)/2
	  endif

c  Bode's 8-point	  
	  if(int.eq.4)then
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
	  endif

	  return
	  end
c	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GAUSSLEG(X1,X2,X,W,N)
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
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diagnosis
      implicit real*8 (a-b,d-h,o-z)
	  implicit complex*16(c-c)
      parameter (nphimax=257,nmax=8193+nphimax,ngauss=32)
      common/data/ pi,hbarc,fm,gammaExp
	  common/range/ as,r0,v4,v6,flam,beta,gamma,a_shell,WFnorm
	  common/control/iAtom,ipot,icalc,nrank,Lwave,Lmax
	  common/q_grid/ dt,dphi,qlist(nmax),phi(nphimax),wt(nmax),Nq,Nphi
	  common/gauss/Xgauss(Ngauss),Wgauss(Ngauss)
	  common/ThreeBodyForce/H0,BigLam
	  data cI/(0d0,1d0)/,c1/(1d0,0d0)/

	    icalcsave=icalc
	    fac=4*pi/(2*pi)**3
	    cphz0=exp(-cI*phi(1))

	    print*,'  '
	    print*,'  >>>>>Begin Diagnosis<<<<< '
	    print*,'  '
	    eps=1d-10
	    ik=2
	    ck=c1*qlist(ik)
	    cEk=3*ck**2/(4*fm)-gamma**2/fm+cI*eps

	     ip=3
		 iq=4
	     cp=qlist(ip)*cphz0
	     cq=qlist(iq)*cphz0
		 cdq_phi=cq*dt*wt(iq)
		 cEtau=cEk-3*cq**2/(4*fm)
		 carg=-fm*cEtau-gamma**2
		 ctau0=1d0/(gamma-sqrt(-fm*cEtau))
c		   ctau=ctau0*(1d0-r0*carg*ctau0/2)
		 ctau=ctaufun(1,1,cEtau)
		 cM1=fac*(2*cZfun(1,1,cp,cq,cEk))*ctau
		 cKern=-cM1*cq**2*cdq_phi
		 cM1alt=cKernel(1,1,cp,cq,cEk)

		print*,' m     = ',fm,'  H0 = ',H0,'   BigLam = ',BigLam
	    print*,' ik = ',ik,' fk = ',ck
		print*,' cp = ',cp,'  cq = ',cq
	  print*,' cEk              = ',cEk
	  print*,' 2*Zfun           = ',2*cZfun(1,1,cp,cq,cEk)
	  print*,' 3/4(q^2-k^2)     = ',carg
	  print*,' 1/(gamma-sqrt()) = ',ctau0
	  print*,' r0*carg/2        = ',r0*carg/2
	  print*,' tau              = ',ctau*fm/(4*pi)
	  print*,' M1               = ',cM1
	  print*,' M1 (alt)         = ',cM1alt
	  print*,' kernel           = ',cKern
	  print*,' Born             = ',cBorn(1,1,cp,cq,cEk)

	   icalc=icalcsave
	   print*,'  '
	   print*,'  >>>>>End Diagnosis<<<<<'
	   stop
c 100   format(i4,i4,5(2x,2e18.12))
	   end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function FLEG(X,NL)
      IMPLICIT REAL*8 (A-H,O-Z)

      p0=1.d0
      p1=x
      if(nl.eq.0) then
         fleg=p0
      elseif(nl.eq.1) then
         fleg=p1
      else
         do l=2,nl
            fleg=(dfloat(2*l-1)*x*p1-dfloat(l-1)*p0)/dfloat(l)
            p0=p1
            p1=fleg
         end do
      end if

      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CPUTime
      real*4 tarray(2),time

      call ETime(tarray,time)
      write(*,1) tarray(1),tarray(2),time,time/60.
 1    format(' Elapsed CPU time:  user',f10.2,',  sys',
     x     f7.2,',  tot',f10.2,' sec',' = ',f7.2,' min')

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


