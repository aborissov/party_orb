PROGRAM reltest		! Relativistic Particle Code: Main Program (JT Dec 2015)
!+ sets up initial field conditions, particle grid and hands over to rkdriver

USE GLOBAL
!USE lare_functions
USE mpi_routines
USE M_DRIVERR, ONLY: RKDRIVE
USE bourdin_fields, ONLY: bour_ini, bour_fini
USE NLFF_fields, ONLY: NLFFF_ini, NLFFF_fini
USE MHDp_fields, ONLY: MHDp_ini, MHDp_fini
USE M_products, ONLY: DOT, CROSS
USE M_fields, ONLY: FIELDS
USE gammadist_mod, ONLY: random_gamma
USE omp_lib

IMPLICIT NONE

 INTEGER :: NOK, NBAD
 INTEGER :: NKEEP,time_no
 INTEGER :: pos_no_x,pos_no_y,pos_no_z,pos_no_alpha,pos_no_ekin
 INTEGER :: pnmax, fcount
 INTEGER, DIMENSION(3) :: pos_no_r
 
 LOGICAL :: file_exists, reset_flag
 
 REAL(num), DIMENSION(3) :: gds, lbox
 !REAL(num), DIMENSION(NKEEPMAX) :: TT 
 !REAL(num), DIMENSION(NKEEPMAX,3) :: S, TOTAL
 
 CHARACTER(LEN=35)	 :: finfile, sumname
 
 !welcome screen
 WRITE(*,*) "====RELATIVISTIC particle code====="
 WRITE(*,*) "using ", FMOD, " fields module."
 bourdinflag=.FALSE.
 l3dflag=.FALSE.
 l2dflag=.FALSE.
 analyticalflag=.FALSE.
 NLFFflag=.FALSE.
 MHDpflag=.FALSE.
 !reset_flag=.FALSE.

 ! read in input parameters
 CALL read_param
    
 ! initial setup options depend on chosen environment
  IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d"))) THEN
   c_ndims=3
   xee=xe
   yee=ye
   zee=ze
   l3dflag=.TRUE.
   CALL MPI_INIT(errcode)
   CALL mpi_initialise      ! mpi_routines.f90
   
   ltimes=ltimes*l3dtimenorm
   print*, 'loaded times:', ltimes
   IF (((nframes.GT.1)).AND.((T1/tscl.lt.ltimes(1)).OR.(T2/tscl.gt.ltimes(nframes)))) THEN
    PRINT*, 'FATAL ERROR!' 
    PRINT*, '(normalised) start/end times of particle range go beyond Lare grid of times'
    PRINT*, '-> RETHINK normalisation, ADD in more snapshots, or LIMIT orbit lifetime.'
    STOP
   ENDIF
   PRINT*, '..evaluating particle array against lare grid..' 
  ELSE IF ((str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN
   l2dflag=.TRUE.
   c_ndims=2   
   xee=xe
   yee=ye
   zee=ze
   CALL MPI_INIT(errcode)
   CALL mpi_initialise_2d
   IF ((R1(3).NE.R2(3)).OR.(RSTEPS(3).GT.1)) THEN
    PRINT*, '-FATAL ERROR-' 
    PRINT*, 'Lare2D data requires single z position value in initial grid!!'
    STOP
   ENDIF
   IF (((nframes.GT.1)).AND.((T1/tscl.lt.ltimes(1)).OR.(T2/tscl.gt.ltimes(nframes)))) THEN
    PRINT*, '-FATAL ERROR-' 
    PRINT*, '(normalised start/end times of particle range go beyond Lare grid of times)'
    PRINT*, '-> RETHINK normalisation, ADD in more snapshots, or LIMIT orbit lifetime.'
    STOP
   ENDIF
  ELSE IF ((str_cmp(FMOD, "SEP")).OR.(str_cmp(FMOD, "sep"))) THEN
    analyticalflag=.TRUE.
    PRINT*, '..evaluating particle array against analytical field bounds..'
  ELSE IF ((str_cmp(FMOD, "CMT")).OR.(str_cmp(FMOD, "cmt"))) THEN
      !CMT  setup?
  ELSE IF ((str_cmp(FMOD, "TEST")).OR.(str_cmp(FMOD, "test"))) THEN
      !test setup?
   ! ALEXEI: temporary hardcoded limits added to make things run. Make this
   ! whole limit selector work better somehow.
   xee=(/-1., 1./)
   yee=(/-1., 1./)
   zee=(/-1., 1./)
  ELSE IF ((str_cmp(FMOD, "NLFF")).OR.(str_cmp(FMOD, "nlff"))) THEN
   NLFFflag=.TRUE.
   CALL NLFFF_ini
   xee=(/myx(6),myx(nx-5)/)
   yee=(/myy(6),myy(ny-5)/)
   zee=(/myz(6),myz(nz-5)/)     
   print*, '----'
   PRINT*, '..evaluating particle array against NLFF grid..' 
  ELSE IF ((str_cmp(FMOD, "MHDp")).OR.(str_cmp(FMOD, "mhdp"))) THEN
   MHDpflag=.TRUE.
   CALL MHDp_ini
   xee=(/myx(6),myx(nx-5)/)
   yee=(/myy(6),myy(ny-5)/)
   zee=(/myz(6),myz(nz-5)/)     
   print*, '----'
   PRINT*, '..evaluating particle array against Paolos MHD grid..' 
  ELSE IF ((str_cmp(FMOD, "BOR")).OR.(str_cmp(FMOD, "bor"))) THEN
   bourdinflag=.TRUE.
   CALL bour_ini      ! read in data
   xee=(/myx(6),myx(nx-5)/)
   yee=(/myy(6),myy(ny-5)/)
   zee=(/myz(6),myz(nz-5)/)
   print*, '----'
   PRINT*, '..evaluating particle array against BOURDIN grid..'
  ELSE
   PRINT*, "incorrect module selection, choose from:"
   PRINT*, "['l3d','l2d','sep','CMT','test','bour', 'nlff']"
   STOP
  END IF
  
  ! chose particle range xe/ye/ze -> particles must not start outside this range!
  IF (((R1(1)/lscl).le.xee(1)).OR.((R2(1)/lscl).ge.xee(2)))  THEN
    WRITE(*,*) '..particles not within x extent '
    print*, R1(1)/lscl
    print*, xee(1) 
    print*, R2(1)/lscl
    print*, xee(2) 
    gflag=.true.
   ENDIF
   IF (((R1(2)/lscl).le.yee(1)).OR.((R2(2)/lscl).ge.yee(2)))  THEN
    WRITE(*,*) '..particles not within y extent '
    gflag=.true.
   ENDIF
   IF (((R1(3)/lscl).le.zee(1)).OR.((R2(3)/lscl).ge.zee(2)))  THEN
    WRITE(*,*) '..particles not within z extent '
    gflag=.true.
   ENDIF
   IF (gflag) THEN
    PRINT*, '-FATAL ERROR-'
    WRITE(*,*) '(particle grid out of bounds)'
    STOP
   ELSE 
    WRITE(*,*) 'fine!'
   ENDIF
   
  DO i=1,3
   IF (RSTEPS(i).EQ.1) THEN 	;! if rsteps=1, 1/(rsteps-1)=1/0!!
    gds(i)=1.0_num 
   ELSE 
    gds(i)=1.0_num/REAL(RSTEPS(i)-1)
   ENDIF
  ENDDO
  lbox=(/R2(1)-R1(1),R2(2)-R1(2),R2(3)-R1(3)/)

  nparticles=RSTEPS(1)*RSTEPS(2)*RSTEPS(3)*(AlphaSteps-1)*EkinSteps
  dalpha = (AlphaMax-Alphamin)/(Alphasteps - 1.0d0) !added by S.Oskoui

  T1Keep=T1
  T2Keep=T2
 
  IF (JTo4) open(49,file=dlocR//'finishr.tmp' ,recl=1024,status='unknown')
  
  PRINT*, ''
  ! restart our calculation for certain particles within a given grid?
  ! can use to divvy up the same grid to different CPUs.
  IF (p_restart) THEN	! p_restart_no is the place a calculation starts again at
   PRINT*, '--RESTARTING p_grid--'
   pn=p_restart_no-1
   pos_no_x = pn / (RSTEPS(2)*RSTEPS(3))  
   pos_no_y = MOD(pn / RSTEPS(3), RSTEPS(2))
   pos_no_z = MOD(pn,RSTEPS(3))
   fcount=1
  ELSE
   PRINT*, '--starting particle grid from beginning--'   
   pn=0
   pos_no_x=0
   pos_no_y=0
   pos_no_z=0
   fcount=0
  ENDIF
  
  IF (p_stop) THEN	!p_stop_no is the particle number where we stop calculating
   pnmax=p_stop_no-1
  ELSE			! otherwise iterate to total number of expected particles
   pnmax=nparticles
  ENDIF
  
  IF (writesum) THEN
   WRITE(sumname,"(A,'Rsum',i1'.dat')"),dlocR, fcount
   inquire(file=sumname, exist=file_exists)
   DO WHILE (file_exists)
    print*, 'summary file called ', sumname, 'encountered! moving up one..'
    fcount=fcount+1
    WRITE(sumname,"(A,'Rsum',i1'.dat')"),dlocR, fcount
    inquire(file=sumname, exist=file_exists)
    IF (fcount.gt.9) THEN
      PRINT*, 'TEN summary files encountered: delete or move some!'
      EXIT 
    ENDIF
   ENDDO
   print*, 'dumping start and end times, positions and energies to ', sumname
   open(39,file=sumname,recl=1024,status='unknown')
  ENDIF

  maxwellEfirst=.TRUE.
  !$OMP PARALLEL
  DO pos_no_x = 0, RSTEPS(1)-1, 1
   DO pos_no_y = 0, RSTEPS(2)-1, 1
    1066 DO pos_no_z = 0, RSTEPS(3)-1, 1
    ! ALEXEI: do we really need to check that we're not exceeding the max number of particles?
    IF (pn .GT. pnmax) THEN
      EXIT
    END IF

     DO pos_no_alpha = 2, AlphaSteps,1
      DO pos_no_ekin = 1, EkinSteps,1
       pos_no_r = (/pos_no_x,pos_no_y, pos_no_z/)
       
       IF (RANDOMISE_R) THEN
        CALL init_random_seed()
        CALL RANDOM_NUMBER(tempr)
        RSTART   = R1+tempr*lbox	!randomise position in bounds set in input
       ELSE
        RSTART   = R1+lbox*(pos_no_r*1.0_num)*gds
       ENDIF
       RSTARTKEEP=RSTART     !remember where we started
       
       pn= pn + 1
       !call progress(pn,nparticles) ! generate the progress bar.
       
       !IF (JTo4) write(49,"(I4)",advance='no'), pn	   
	   
       1000 format ("particle no. ",i7,"/",i7, ", R=(",ES9.2,",",ES9.2,",",ES9.2,")") 
       print 1000, pn,nparticles, RSTART

       T1=T1Keep
       T2=T2Keep

       ! calculate pitch angle
       IF (RANDOMISE_A) THEN
        CALL init_random_seed()
        CALL RANDOM_NUMBER(tempa)
        alpha = Alphamin+tempa*(Alphamax-Alphamin)
       ELSE
        alpha = Alphamin+dalpha*(pos_no_alpha -1)	!added by S.Oskoui
       ENDIF
       alpha = alpha*Pi/180.0d0				! RADEG: added by S.Oskoui

       ! calculate kinetic energy - needed to define parallel and perp velocity based on alpha
       !pos_no_ekin starts from 0, if started from 1 then (stepekin-1)
       IF (RANDOMISE_E) THEN
        !Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)*tempe   
	Ekin= random_gamma(1.5_num, kb*maxwellpeaktemp, maxwellEfirst)
	EKin=Ekin*6.242e18  !(convert to eV)
	maxwellEfirst = .FALSE.
       ELSE
        Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)
       ENDIF
       
       !PRINT*,'Normalising:'
       RSTART=RSTART/Lscl
       RSTARTKEEP=RSTARTKEEP/Lscl
       T1=T1/Tscl
       T2=T2/Tscl
       
       ! WARNING passing in dimensional Ekin into mu calc
       CALL JTMUcalc(MU,USTART,GAMMASTART,Ekin,Alpha,RSTART,T1,T2, reset_flag)
       
       !print*, pos_no_r
       
       IF (reset_flag) THEN 
        ! very small initial B detected
        IF (RANDOMISE_R) THEN 
	 ! if positions are random, then go back and generate a new random position	 
	 print*, 'initial |B| is too small: trying a new position'
	 pn=pn-1
	 !pos_no_z=pos_no_z-1
	 !CYCLE
	 GO TO 1066
	 ! remember, random variable seed based on system clock, so will take a second or two to generate a new seed!
	ELSE   
         ! if the position is SPECIFIED then skip this one	 
	 print*, 'initial |B| is too small: skipping this particle'
         CYCLE
	ENDIF
       ENDIF
       
       IF (JTo4) write(49,"(I4)",advance='no'), pn
       
       USTARTKEEP=USTART  
       GAMMASTARTKEEP=GAMMASTART       
       Ekin = Ekin*AQ/Ekscl        
       
       Rlost=.FALSE.
       !Call the rk sophisticated driver, which then works out the arrays for the
       !time steps and positions.
       !CALL RKDRIVE(RSTART,USTART,GAMMASTART,MU,T1,T2,EPS,H1,NOK,NBAD,TT,S,TOTAL) ! Unnecessary TT, S, Total?
       CALL RKDRIVE(RSTART,USTART,GAMMASTART,MU,T1,T2,EPS,H1,NOK,NBAD)
       
       IF (writesum) write(39,*) Tscl*(T1), Lscl*RSTARTKEEP, oneuponAQ*(GAMMASTARTKEEP-1)*m*c*c, &
       Tscl*(T2), Lscl*RSTART, oneuponAQ*(GAMMASTART-1)*m*c*c
        
       NKEEP = (NOK +NBAD)/NSTORE

      ! CALL WRITE_ENDTIME(RSTART,T2,MU,VPARSTART)
      
      END DO
     END DO
    END DO
   END DO
  END DO 
  !$OMP END PARALLEL
  IF (writesum) CLOSE(39)
  IF (JTo4) CLOSE(49)
 !CALL MAKEFILE(time_no)
  
 !END DO
 
 IF ((str_cmp(FMOD, "LARE")).OR.(str_cmp(FMOD, "lare"))) THEN	!forget arrays at end
  CALL mpi_close                     ! mpi_routines.f90
  CALL MPI_FINALIZE(errcode)
 ENDIF
 IF ((str_cmp(FMOD, "BOUR")).OR.(str_cmp(FMOD, "bour"))) THEN
  CALL bour_fini      		! deallocate stuff, leave everything nice and tidy
 ENDIF 
 IF ((str_cmp(FMOD, "NLFF")).OR.(str_cmp(FMOD, "nlff"))) THEN
  CALL NLFFF_fini      		
 ENDIF 
 IF ((str_cmp(FMOD, "MHDp")).OR.(str_cmp(FMOD, "mhdp"))) THEN
  CALL MHDp_fini      		
 ENDIF 
!------------------------------------------------------------------------------!
 Contains
!------------------------------------------------------------------------------!
SUBROUTINE JTMUcalc(mu,USTART,GAMMASTART, Ekin,alpha,RSTART,T1,T2, resetflag)

  REAL(num), DIMENSION(3),INTENT(IN) :: RSTART
  REAL(num), INTENT(IN) :: T1,T2, Ekin, Alpha				! 
  REAL(num), INTENT(OUT) :: mu, gammastart, Ustart
  REAL(num), DIMENSION(3) :: B,El,a2,a3,a4,a5,a6,a7,a8,a9,a10,ue, RT
  !REAL(num) :: magB,vtot,vperp, vparstart!,Erest
  REAL(num) :: modB,vtot, gamma, Etemp
  LOGICAL, INTENT(OUT)   :: resetflag
 
 resetflag=.FALSE.
 
 !calculate B, E, V at this point/time:
 CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2)

 !calculate magnitude of B
 modB=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
 modB=sqrt(modB)
 RT=RSTART

 IF (modb.le.lowbthresh) THEN 
  resetflag=.TRUE.
 ENDIF

 !print*, 'RT=', RT
 !PRINT*, 'B=', B
 !STOP
 !print*, modB

 !Erest = (M*c*c)*oneuponAQ
 
 ! E X B drift
 ue=cross(El,B)/dot(B,B)  !*0.5
 Etemp=Ekin
 
! print*, 'Etemp=', Ekin
 
 !need to check if 1/2mUE^2 is covered by the initial KE
 IF (Ekin.lt.0.5d0*m*vscl*vscl*dot(ue,ue)*6.242e18) THEN
  PRINT*, 'WARNING: raising Initial KE to account for local UE drift'
  Etemp=0.5*m*vscl*vscl*dot(ue,ue)*6.242e18
 ENDIF
 
 gamma=Etemp/m/c/c*AQ+1.0d0
 
 ! vtot^2=vpar^2+vperp^2+UE^2!!, UE should NOT come out of Gamma!
 !vtot=c/gamma*sqrt(gamma*gamma-1)			! I *think* vtot is dimensional
 !vtot=sqrt(c*c/gamma/gamma*(gamma*gamma-1)-vscl*vscl*dot(ue,ue))
 !vtot=sqrt(c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue))			! THRELFALL ET AL 2015
 vtot=sqrt(c*c*(1.0d0-1.0d0/gamma/gamma)-vscl*vscl*dot(ue,ue))
 ! vtot=sqrt(c*c-c/gamma*c/gamma-vce*vce*dot(ue,ue))			! THRELFALL ET AL 2015
 
 !print*, 'gamma ', gamma
 !print*, 'uE^2:', vscl*vscl*dot(ue,ue)
 !print*, 'c*c-c/gamma*c/gamma', c*c-c/gamma*c/gamma
 !print*, 'c*c*(1.0d0-1.0d0/gamma/gamma)', c*c*(1.0d0-1.0d0/gamma/gamma)
 !PRINT*, c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue)
 !print*, 'vtot ', vtot 
 
 IF (vtot/=vtot) THEN
  PRINT*, VTOT
  PRINT*, 'Vtot NaN ENCOUNTERED'
  STOP
 ENDIF
 
 !Ustart=(c*sqrt(gamma*gamma-1))*cos(alpha)
 Ustart=gamma*vtot*cos(alpha)
 
 
 !print*, vtot*cos(alpha)
 
 mu=0.5_num*m*vtot*vtot*sin(alpha)*sin(alpha)*gamma*gamma/modB		
 
 USTART=Ustart/vscl					! hence, non-dimensionalising..
 GAMMASTART=gamma
 mu=mu/m/vscl/vscl					! no bscl normalising factor - using normalised B's already!
 !PRINT*, 'vtot=', vtot
 !PRINT*, 'gamma=', gamma
 !PRINT*, 'mu=', mu
 !PRINT*, 'modB=', modB
  !calculate mu
  !mu=0.5_num*vperp*vperp/magB*gammastart*gammastart
 !STOP
 ! WRITE (19,*) RStart,T1,Ekin,Alpha, mu*magB, 0.5_num*vparstart*vparstart
  !WRITE (19,*) vtot,vperp,vparstart,El,B,magB,mu

END SUBROUTINE

END PROGRAM reltest
