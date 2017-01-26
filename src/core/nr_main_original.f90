PROGRAM NREL	! Non-Relativistic Particle Code: Main Program (JT Dec 2015)
!+ sets up initial field conditions, particle grid and hands over to rkdriver

USE GLOBAL
USE M_DRIVERN, ONLY: RKDRIVE
USE lare_functions
USE mpi_routines
USE bourdin_fields, ONLY: bour_ini, bour_fini
USE M_products, ONLY: DOT, CROSS
USE M_fields, ONLY: FIELDS
USE gammadist_mod, ONLY: random_gamma

IMPLICIT NONE

 INTEGER :: NOK, NBAD
 INTEGER :: NKEEP,time_no
 INTEGER :: pos_no_x,pos_no_y,pos_no_z,pos_no_alpha,pos_no_ekin
 INTEGER :: pnmax
 INTEGER, DIMENSION(3) :: pos_no_r
 
 REAL(num), DIMENSION(3) :: gds, lbox
 REAL(num), DIMENSION(3) :: RSTART,RSTARTKEEP!, R1,R2
 REAL(num), DIMENSION(NKEEPMAX) :: TT 
 REAL(num), DIMENSION(NKEEPMAX,3) :: S, TOTAL
 REAL(num), DIMENSION(3) :: Et,Bt,DBDXt,DBDYt,DBDZt,DBDTt,DEDXt,DEDYt,DEDZt,DEDTt,Vft
 CHARACTER(LEN=35)	 :: finfile


! CALL MPI_INIT(errcode)		! mpi initialise
 !welcome screen:
 WRITE(*,*) "=====NON-RELATIVISTIC particle code===="
 WRITE(*,*) "using ", FMOD, " fields module."
 bourdinflag=.FALSE.
 l3dflag=.FALSE.
 l2dflag=.FALSE.
 analyticalflag=.FALSE.

 !read in the parameters from the input file
 CALL read_param

 ! initial setup options depend on chosen environment
  IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d"))) THEN
   c_ndims=3
   !allocate(dims(ndims))
   l3dflag=.TRUE.
   CALL MPI_INIT(errcode)
   CALL mpi_initialise      ! mpi_routines.f90
   IF (((nframes.GT.1)).AND.((T1/tscl.lt.ltimes(0)).OR.(T2/tscl.gt.ltimes(nframes-1)))) THEN
    PRINT*, 'FATAL ERROR!' 
    PRINT*, '(normalised) start/end times of particle range go beyond Lare grid of times'
    PRINT*, '-> RETHINK normalisation, ADD in more snapshots, or LIMIT orbit lifetime.'
    STOP
   ENDIF
   PRINT*, '..evaluating particle array against lare grid..' 
  ELSE IF ((str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN
   l2dflag=.TRUE.
   c_ndims=2
   !ndims=ndims
   !allocate(dims(ndims))
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
   PRINT*, "['l3d','l2d','sep','CMT','test','bour']"
   STOP
  END IF
  
  ! chose particle range xe/ye/ze -> particles must not start outside this range!
   IF (((R1(1)/lscl).le.xe(1)).OR.((R2(1)/lscl).ge.xe(2)))  THEN
    WRITE(*,*) '..particles not within x extent '
    gflag=.true.
   ENDIF
   IF (((R1(2)/lscl).le.ye(1)).OR.((R2(2)/lscl).ge.ye(2)))  THEN
    WRITE(*,*) '..particles not within y extent '
    gflag=.true.
   ENDIF
   IF (((R1(3)/lscl).le.ze(1)).OR.((R2(3)/lscl).ge.ze(2)))  THEN
    WRITE(*,*) '..particles not within z extent '
    gflag=.true.
   ENDIF
   IF (gflag) THEN
    WRITE(*,*) 'terminating: particle grid out of bounds.'
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

  !pn=0
  nparticles=RSTEPS(1)*RSTEPS(2)*RSTEPS(3)*(AlphaSteps-1)*EkinSteps
  dalpha = (AlphaMax-Alphamin)/(Alphasteps - 1.0d0) !added by S.Oskoui
  
  !Adjust T2 to use loop value.
  !T2=time_no*1.0_num

  T1Keep=T1
  T2Keep=T2

  IF (JTo4)  open(48,file=dlocR//'finishnr.tmp' ,recl=1024,status='unknown')
  PRINT*, ''
  IF (p_restart) THEN
   PRINT*, '--RESTARTING p_grid--'
   pn=p_restart_no-1
   pos_no_x = pn / (RSTEPS(2)*RSTEPS(3))  
   pos_no_y = MOD(pn / RSTEPS(3), RSTEPS(2))
   pos_no_z = MOD(pn,RSTEPS(3))
  ELSE
   PRINT*, '--starting particle grid from beginning--'   
   pn=0
   pos_no_x=0
   pos_no_y=0
   pos_no_z=0
  ENDIF 

  IF (p_stop) THEN	!p_stop_no is the particle number where we stop calculating
   pnmax=p_stop_no-1
  ELSE			! otherwise iterate to total number of expected particles
   pnmax=nparticles
  ENDIF

  maxwellEfirst=.TRUE.
  DO WHILE (pos_no_x .LE. RSTEPS(1)-1)
   DO WHILE (pos_no_y .LE. RSTEPS(2)-1)
    DO WHILE ((pos_no_z .LE. RSTEPS(3)-1).AND.(pn .LE. pnmax))
     DO pos_no_alpha = 2, AlphaSteps,1
      DO pos_no_ekin = 1, EkinSteps,1

      !redefine rstart so that it uses the value from before, not the position of last particle!
       
       pos_no_r = (/pos_no_x,pos_no_y, pos_no_z/)
       IF (RANDOMISE_R) THEN
        CALL init_random_seed()
        CALL RANDOM_NUMBER(tempr)
        RSTART   = R1+tempr*lbox	!randomise position in bounds set in input
       ELSE
        RSTART   = R1+lbox*(pos_no_r*1.0_num)*gds
       ENDIF
       pn= pn + 1
       
       !call progress(pn,nparticles) ! generate the progress bar.
       
       IF (JTo4) write(48,"(I4)",advance='no'), pn	   
           
       IF (nparticles.gt.1000) THEN 
        print 1111, pn,nparticles, RSTART     
        1111 format ("particle no. ",i4,"/",i4, ", R=(",ES9.2,",",ES9.2,",",ES9.2,")")
       ELSE 
        print 1112, pn,nparticles, RSTART     
        1112 format ("particle no. ",i3,"/",i3, ", R=(",ES9.2,",",ES9.2,",",ES9.2,")")
       ENDIF    

       T1=T1Keep
       T2=T2Keep

       IF (RANDOMISE_A) THEN
        alpha = Alphamin+dalpha*(pos_no_alpha -1)*tempa	!added by S.Oskoui
       ELSE
        alpha = Alphamin+dalpha*(pos_no_alpha -1)	!added by S.Oskoui
       ENDIF
       alpha = alpha*Pi/180.0d0				! RADEG: added by S.Oskoui

       IF (RANDOMISE_E) THEN
        !Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)*tempe   
	Ekin= random_gamma(1.5_num, kb*maxwellpeaktemp, maxwellEfirst)
	maxwellEfirst = .FALSE.
       ELSE
        Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)
       ENDIF

   !pos_no_ekin starts from 0, if started from 1 then (stepekin-1)
 
   !alpha = pi/(no of steps+1) if fullangle is 1 (ie, steps from >=0 to >Pi (but not including Pi))
   !alpha = pi/2/(no of steps) if fullangle is 0 (steps from 0 to Pi/2 inclusive)
  
       !VPARSTARTKEEP=VPARSTART
       RSTARTKEEP=RSTART
       !PRINT*,'Normalising:'
       RSTART=RSTART/Lscl
       RSTARTKEEP=RSTARTKEEP/Lscl
       T1=T1/Tscl
       T2=T2/Tscl

       !PRINT*,'T1=',T1
       !PRINT*,'T2=',T2
       !PRINT*,'VPARSTART=',VPARSTART
       !PRINT*,'RSTART=',RSTART
       !PRINT*, '***********************************************************'
       
       !call sub to calculate mu from Ekinetic (initial), vpar (to get Epar),
       !initial position and initial time.
       CALL CALC2_MU(MU,vparstart,Ekin,Alpha,RSTART,T1,T2)		! calculates vparstart
       
       
       Ekin = Ekin*AQ/Ekscl !convert energy from ev to joules (inc. in non-rel)
       
       !Call the rk sophisticated driver, which then works out the arrays for the
       !time steps and positions.
       CALL RKDRIVE(RSTART,VPARSTART,MU,T1,T2,EPS,H1,NOK,NBAD,TT,S,TOTAL)

       NKEEP = (NOK +NBAD)/NSTORE

       !CALL WRITE_ENDTIME(RSTART,T2,MU,VPARSTART)

      END DO
     END DO
      pos_no_z=pos_no_z+1
    END DO
    pos_no_y=pos_no_y+1
    pos_no_z=0
   END DO
   pos_no_x=pos_no_x+1
   pos_no_y=0
  END DO
  IF (JTo4) CLOSE(48)
 
 
 !END DO
 
IF ((str_cmp(FMOD, "LARE")).OR.(str_cmp(FMOD, "lare"))) THEN	!forget arrays at end
  CALL mpi_close                     ! mpi_routines.f90
  CALL MPI_FINALIZE(errcode)
 ENDIF

 IF ((str_cmp(FMOD, "BOUR")).OR.(str_cmp(FMOD, "bour"))) THEN
  CALL bour_fini      		! deallocate stuff, leave everything nice and tidy
 ENDIF 

!------------------------------------------------------------------------------!
 Contains
!------------------------------------------------------------------------------!
SUBROUTINE CALC2_MU(mu,vparstart,Ekin,alpha,RSTART,T1,T2)

  REAL(num), DIMENSION(3),INTENT(IN) :: RSTART
  REAL(num), INTENT(IN) :: T1,T2, Ekin, Alpha
  REAL(num), INTENT(OUT) :: mu, vparstart
  REAL(num), DIMENSION(3) :: B,El,a2,a3,a4,a5,a6,a7,a8,a9,a10,ue
  REAL(num) :: magB,vtot,vperp,Erest
 
 !calculate B, E, V at this point/time:
 CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2)

 !calculate magnitude of B
 magB=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
 magB=sqrt(magB)

 Erest = (M*c*c)*oneuponAQ
 ! E X B drift
 ue=cross(El,B)/dot(B,B)  !*0.5
 
 vtot=sqrt(2.d0*Ekin-dot(ue,ue))

 !vtot= sqrt(((2.d0*Ekin)/Erest)*c**2 - dot(ue,ue))   ! Normalised vtot 

 !VPERP = sqrt((((2.d0*Ekin))/511769.6629)*c**2 - dot(ue,ue))* sin(alpha)
 !VPARSTART = sqrt((((2.d0*Ekin))/511769.6629)*c**2 - dot(ue,ue))* cos(alpha)

  vparstart=vtot*cos(alpha) !I think vparstart IS normalised - 
  !PRINT*, 'Vtot',vtot
  !PRINT *, 'Alpha', Alpha
  !PRINT*, 'cosalpha',cos(Alpha)

  vperp=vtot*sin(alpha)

 !calculate mu
  mu=0.5_num*vperp*vperp/magB
  !mu=vperp**2/magB/2.0_num
  !mu  = (M*(Vscl*Vperp)**2)/(abs(Q)*magB*2d0)                             !This dont work ??
  !mu = (511769.6629*(Vscl*Vperp)**2)/(c**2*B0*magB*2.0)
  !mu = (((((2.d0*Ekin))/511769.6629)*c**2 - dot(ue,ue))* sin(alpha)**2 )/ magB/ 2d0
  !PRINT*, 'El', El
  !PRINT*, 'B', B
  !PRINT*, 'magB', magB
  !PRINT*,"vtot",vtot
  !PRINT*, 'Vparstart', vparstart
  !PRINT*,'ue', ue
  !PRINT*, 'dot(ue,ue)',dot(ue,ue) 
  !!PRINT*, 'c',c
  !!PRINT*, 'M',M
  !!PRINT*, 'Q',q  
  !PRINT*,'mu',mu
  !!PRINT*, 'E0', E0
  !!PRINT*, 'B0', B0 
  !PRINT*, 'Erest', Erest
  !!PRINT*, '(L/Tscl/Vscl)', (L/Tscl/Vscl)
  !!PRINT*, '(L/Vscl/Tscl)', (L/Vscl/Tscl)
  !!PRINT*, '(1.0_num/L)', (1.0_num/L)

 !some data dumps for checking
 !  PRINT*, "In Calc2_mu"
  !PRINT*,"alpha",alpha

  !stop
  !PRINT*, "Ekin:",Ekin,"; vparstart:",vparstart
  !PRINT*, "B",B
  !PRINT*, "magB", magB
  !PRINT*, "mu:",mu
  !PRINT*,"muB+0.5m*vpar^2",mu*magB+0.5*vparstart**2
  !PRINT*,"pitch angle",acos(vparstart/(sqrt(2*Ekin))) 
  !PRINT*, "leaving calc2_mu"

  !output Ekin,alpha,Eperp,Epar
  !vtot,vperp,vpar,Efield,Bfield,mu
  !WRITE (19,*) RStart,T1,Ekin,Alpha, mu*magB, 0.5_num*vparstart*vparstart
  ! WRITE (19,*) "RStart,T1,Ekin,alpha, mu*magB, 0.5*vparstart**2"
  !WRITE (19,*) vtot,vperp,vparstart,El,B,magB,mu
 !WRITE (19,*) "vtot,vperp,vparstart,El,B,magB,mu"
END SUBROUTINE
!------------------------------------------------------------------------------!
END PROGRAM NREL
