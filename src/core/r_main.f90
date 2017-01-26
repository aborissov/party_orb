PROGRAM reltest		! Relativistic Particle Code: Main Program (JT Dec 2015)
!+ sets up initial field conditions, particle grid and hands over to rkdriver

USE GLOBAL
USE lare_functions
USE mpi_routines
USE M_DRIVERR, ONLY: RKDRIVE
USE bourdin_fields, ONLY: bour_ini, bour_fini
USE M_products, ONLY: DOT, CROSS
USE M_fields, ONLY: FIELDS
USE gammadist_mod, ONLY: random_gamma

IMPLICIT NONE

 INTEGER :: NOK, NBAD, particle_tstart, particle_tend,j
 INTEGER :: NKEEP,time_no,size_r
 INTEGER :: pos_no_x,pos_no_y,pos_no_z,pos_no_alpha,pos_no_ekin
 INTEGER :: pnmax
 INTEGER, DIMENSION(3) :: pos_no_r
 integer, parameter    :: print_data = 10
 
 REAL(num), DIMENSION(3) :: gds, lbox
 REAL(num), DIMENSION(3) :: RSTART,RSTARTKEEP!, R1,R2
 !REAL(num), DIMENSION(NKEEPMAX) :: TT 
 !REAL(num), DIMENSION(NKEEPMAX,3) :: S, TOTAL
 REAL(num)                       :: TT 
 REAL(num), DIMENSION(3)         :: S, TOTAL
 REAL(num), DIMENSION(3) :: Et,Bt,DBDXt,DBDYt,DBDZt,DBDTt,DEDXt,DEDYt,DEDZt,DEDTt,Vft,ue
 REAL(num)               :: rrr1,ttt1,etaa1
 CHARACTER(LEN=35)	 :: finfile
 logical                :: message_flag

 
 bourdinflag=.FALSE.
 l3dflag=.FALSE.
 l2dflag=.FALSE.
 analyticalflag=.FALSE.

 ! read in input parameters
 CALL read_param
    
 ! initial setup options depend on chosen environment
  IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d"))) THEN
   c_ndims=3
   !allocate(dims(ndims))
   l3dflag=.TRUE.
   CALL MPI_INIT(errcode)
   CALL mpi_initialise      ! mpi_routines.f90
   IF (((nframes.GT.1)).AND.((T1/tscl.lt.ltimes(0)).OR.(T2/tscl.gt.ltimes(nframes-1)))) THEN
    print *, 'main particle initial time = ', T1/tscl, ' first snapshot time = ', ltimes(1), ' particle end time = ', T2/tscl, ' last snapshot time = ', ltimes(nframes-1)
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
    print *, 'main particle initial time = ', T1/tscl, ' first snapshot time = ', ltimes(1), ' particle end time = ', T2/tscl, ' last snapshot time = ', ltimes(nframes-1)
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
    print *, 'r_main R1(1)/lscl = ', R1(1)/lscl, ' R2(1)/lscl = ', R2(1)/lscl, 'grid bounds ', xe(1), xe(2)
    WRITE(*,*) '..particles not within x extent '
    gflag=.true.
   ENDIF
   IF (((R1(2)/lscl).le.ye(1)).OR.((R2(2)/lscl).ge.ye(2)))  THEN
    print *, 'r_main R1(2)/lscl = ', R1(2)/lscl, ' R2(2)/lscl = ', R2(2)/lscl, 'grid bounds ', ye(1), ye(2)
    WRITE(*,*) '..particles not within y extent '
    gflag=.true.
   ENDIF
   IF (((R1(3)/lscl).le.ze(1)).OR.((R2(3)/lscl).ge.ze(2)))  THEN
    WRITE(*,*) '..particles not within z extent '
    gflag=.true.
   ENDIF
   IF (gflag) THEN
    PRINT*, '-FATAL ERROR-'
    WRITE(*,*) '(particle grid out of bounds)'
    STOP
   ELSE 
    if (rank .eq. 0) WRITE(*,*) 'fine!'
   ENDIF
   
  DO i=1,3
   IF (RSTEPS(i).EQ.1) THEN 	;! if rsteps=1, 1/(rsteps-1)=1/0!!
    gds(i)=1.0_num 
   ELSE 
    gds(i)=1.0_num/REAL(RSTEPS(i)-1)
   ENDIF
  ENDDO
  lbox=(/R2(1)-R1(1),R2(2)-R1(2),R2(3)-R1(3)/)

  nparticles=RSTEPS(1)*RSTEPS(2)*RSTEPS(3)*AlphaSteps*EkinSteps
  if (Alphasteps .eq. 1) then
    dalpha = (AlphaMax-Alphamin)/(Alphasteps) 
  else
    dalpha = (AlphaMax-Alphamin)/(Alphasteps - 1.0d0) 
  endif
  
  !Adjust T2 to use loop value.
  !T2=time_no*1.0_num

  T1Keep=T1
  T2Keep=T2
 
  if (p_restart) then
    IF (JTo4) open(49,file=dlocR//'finishr.tmp',recl=1024,status='unknown',position = 'append')
    IF (JTo4) open(50,file=dlocR//'startr.tmp' ,recl=1024,status='unknown',position = 'append')
    IF (JTo4) open(51,file=dlocR//'failedr.tmp' ,recl=1024,status='unknown',position = 'append')
    IF (JTo4) open(52,file=dlocR//'exitr.tmp' ,recl=1024,status='unknown',position = 'append')
  else
    IF (JTo4) open(49,file=dlocR//'finishr.tmp',recl=1024,status='unknown')
    IF (JTo4) open(50,file=dlocR//'startr.tmp' ,recl=1024,status='unknown')
    IF (JTo4) open(51,file=dlocR//'failedr.tmp' ,recl=1024,status='unknown')
    IF (JTo4) open(52,file=dlocR//'exitr.tmp' ,recl=1024,status='unknown')
  endif
  
  PRINT*, ''
  ! restart our calculation for certain particles within a given grid?
  ! can use to divvy up the same grid to different CPUs.
  IF (p_restart) THEN	! p_restart_no is the place a calculation starts again at
   PRINT*, '--RESTARTING p_grid at ', p_restart_no, '--'
   pn=p_restart_no-1
   pos_no_x = pn / (RSTEPS(2)*RSTEPS(3)*EkinSteps*AlphaSteps)  
   pos_no_y = MOD(pn / (RSTEPS(3)*EkinSteps*AlphaSteps), RSTEPS(2))
   pos_no_z = MOD(pn/(EkinSteps*AlphaSteps),RSTEPS(3))
   pos_no_alpha = mod(pn/EkinSteps, AlphaSteps)
   pos_no_ekin = mod(pn,EkinSteps)
  ELSE
   if (rank .eq. 0) PRINT*, '--starting particle grid from beginning--'   
   pn=0
   pos_no_x=0
   pos_no_y=0
   pos_no_z=0
   pos_no_alpha=0
   pos_no_ekin=0
  ENDIF
  
  IF (p_stop) THEN	!p_stop_no is the particle number where we stop calculating
   pnmax=p_stop_no-1
  ELSE			! otherwise iterate to total number of expected particles
   pnmax=nparticles
  ENDIF

  !print *, 'main pos_no_x = ', pos_no_x, 'RSTEPS(1)', RSTEPS(1)
  !print *, 'main pos_no_y = ', pos_no_y, 'RSTEPS(2)', RSTEPS(2)
  !print *, 'main pos_no_z = ', pos_no_z, 'RSTEPS(3)', RSTEPS(3)
  !print *, 'main pn', pn, ' pnmax ', pnmax
  maxwellEfirst=.TRUE.
  H1=H1/Tscl
  dt_min = H1


  !if (rank .eq. 1) then
  !      finish_data%pos = (/1.0_num,2.0_num,3.0_num/)
  !      finish_data%par_vel = 2.5_num
  !      finish_data%perp_vel = 3.2_num
  !      finish_data%kin_e = 45.2E22_num
  !      finish_data%p_angle = 21.0_num
  !      finish_data%p_failed = 1
  !      finish_data%p_exit = -1
  !      finish_data%par_num = 50
  !      print *, 'finish data = ', finish_data
  !      call mpi_send(finish_data,1,mpi_finish_data2,0,0,mpi_comm_world,errcode)
  !else if (rank .eq. 0) then
  !      print *, 'recv data1 = ', recv_data
  !      call mpi_recv(recv_data,1,mpi_finish_data2,1,0,mpi_comm_world,errcode)
  !      print *, 'recv data2 = ', recv_data
  !      stop
  !endif
 
 !welcome screen
 if (rank .eq. 0) then
   WRITE(*,*) "====RELATIVISTIC particle code====="
   WRITE(*,*) "using ", FMOD, " fields module."
 endif

  
  if (rank .eq. 0) then
 DO WHILE (pos_no_x .LE. RSTEPS(1)-1)
   DO WHILE (pos_no_y .LE. RSTEPS(2)-1)
    DO WHILE (pos_no_z .LE. RSTEPS(3)-1)
     DO WHILE (pos_no_alpha .LE. AlphaSteps-1)
      DO WHILE ((pos_no_ekin .LE. EkinSteps-1) .AND. (pn .LE. pnmax))
       particle_tstart = time8()
       
       pos_no_r = (/pos_no_x,pos_no_y, pos_no_z/)
       !print*, tempr
       IF (RANDOMISE_R) THEN
        CALL init_random_seed()
        CALL RANDOM_NUMBER(tempr)
        RSTART   = R1+tempr*lbox	!randomise position in bounds set in input
       ELSE
        RSTART   = R1+lbox*(pos_no_r*1.0_num)*gds
       ENDIF
       
       pn= pn + 1
       
       !call progress(pn,nparticles) ! generate the progress bar.
       
       !IF (JTo4) write(49,"(I4)",advance='no'), pn	   
	   

       T1=T1Keep
       T2=T2Keep

       !IF (RANDOMISE_A) THEN
       ! alpha = Alphamin+dalpha*(pos_no_alpha -1)*tempa	!added by S.Oskoui
       !ELSE
       ! alpha = Alphamin+dalpha*(pos_no_alpha -1)	!added by S.Oskoui
       !ENDIF
       IF (RANDOMISE_A) THEN
        alpha = Alphamin+dalpha*(pos_no_alpha)*tempa	!added by S.Oskoui
       ELSE
        alpha = Alphamin+dalpha*(pos_no_alpha)	!added by S.Oskoui
       ENDIF
       alpha = alpha*Pi/180.0d0				! RADEG: added by S.Oskoui
       
       if ((print_number .eq. 1) .and. (nproc .eq. 1)) then
       IF (nparticles.gt.1000) THEN 
        print 1000, pn,nparticles, RSTART,alpha*57.3
        1000 format ("particle no. ",i4,"/",i4, ", R=(",ES9.2,",",ES9.2,",",ES9.2,"), alpha = ",ES9.2,",")
       ELSE 
        print 1001, pn,nparticles, RSTART ,alpha*57.3    
        1001 format ("particle no. ",i3,"/",i3, ", R=(",ES9.2,",",ES9.2,",",ES9.2,"), alpha = ",ES9.2,",")
       ENDIF    
       endif

       
       !pos_no_ekin starts from 0, if started from 1 then (stepekin-1)
       IF (RANDOMISE_E) THEN
        !Ekin=EKinLow+(EKinHigh-EKinLow)*pos_no_ekin/(EkinSteps*1.0d0)*tempe   
        do while (.TRUE.)
	  Ekin= random_gamma(1.5_num, kb*maxwellpeaktemp, maxwellEfirst)
          CALL FIELDS(RSTART/Lscl,T1/Tscl,Et,Bt,DBDXt,DBDYt,DBDZt,DBDTt,DEDXt,DEDYt,DEDZt,DEDTt,Vft,T1,T2,rrr1,ttt1,etaa1)
          ue=cross(Et,Bt)/dot(Bt,Bt)  !*0.5
          if (c*c-c/(Ekin/m/c/c+1.0d0)*c/(Ekin/m/c/c+1.0d0)-vscl*vscl*dot(ue,ue) .gt. 0.0_num) exit
        enddo
	EKin=Ekin*6.242e18  !(convert to eV)
	!maxwellEfirst = .FALSE.
       ELSE
        if (EkinSteps .eq. 1) then
          Ekin=EKinLow+(EKinHigh-EKinLow)*(pos_no_ekin+1)/(EkinSteps*1.0d0)
        else
          Ekin=EKinLow+(EKinHigh-EKinLow)*(pos_no_ekin)/((EkinSteps-1)*1.0d0)        ! Switch to this if want to have energy inclusive of endpoints
        endif
        !print *, 'main Ekin = ', Ekin, EKinLow,(EKinHigh-EKinLow),(pos_no_ekin+1),(EkinSteps*1.0d0)
       ENDIF
       if ((print_number .eq. 1) .and. (nproc .eq. 1)) print*, 'kinetic energy (in eV) = ', EKIN
       !alpha = pi/(no of steps+1) if fullangle is 1 (ie, steps from >=0 to >Pi (but not including Pi))
       !alpha = pi/2/(no of steps) if fullangle is 0 (steps from 0 to Pi/2 inclusive)
  
       RSTARTKEEP=RSTART
       USTARTKEEP=USTART
       GAMMASTARTKEEP=GAMMASTART
       
       !PRINT*,'Normalising:'
       RSTART=RSTART/Lscl
       RSTARTKEEP=RSTARTKEEP/Lscl
       T1=T1/Tscl
       T2=T2/Tscl
       
       ! WARNING passing in dimensional Ekin into mu calc
       CALL JTMUcalc(MU,USTART,GAMMASTART,Ekin,Alpha,RSTART,T1,T2)
              
       Ekin = Ekin*AQ/Ekscl        
       
       !Call the rk sophisticated driver, which then works out the arrays for the
       !time steps and positions.
       if (nproc .eq. 1) then
          CALL RKDRIVE(pn,RSTART,USTART,GAMMASTART,MU,T1,T2,EPS,H1,NOK,NBAD,TT,S,TOTAL)
          NKEEP = (NOK +NBAD)/NSTORE
       else
             !print *, 'r_main recv = ', recv_data
          call mpi_recv(recv_data,1,mpi_finish_data2,mpi_any_source,mpi_any_tag,mpi_comm_world,status,errcode)
             !print *, 'r_main recv = ', recv_data
             !print *, 'r_main received tag', status(mpi_tag)
          start_data%R_s = RSTART
          start_data%u_s = USTART
          start_data%gamma_s = GAMMASTART
          start_data%mu_s = MU
          start_data%t1_s = T1
          start_data%t2_s = T2
          start_data%eps_s = EPS
          start_data%h1_s = H1
          start_data%nok_s = NOK
          start_data%nbad_s = NBAD
          start_data%tt_s = TT
          start_data%s_s = S
          start_data%total_s = TOTAL
          start_data%par_num = pn
          sender = status(mpi_source)
          !print *, 'rank ', rank ,' particle ', pn, ' send_data = ', start_data
          call mpi_send(start_data,1,mpi_start_data,sender,0,mpi_comm_world,errcode) ! tag 0 means there are still particles left to run
          !print *, 'particle ', pn, ' running on rank', sender

          ! if not on first step then write the data received from the worker
          ! (if first step tag should be 0)
          if (status(mpi_tag) .eq. 1) then
             !print *, 'writing particle', recv_data%par_num, ' data'
             IF (JTo4) write(49,*),recv_data%par_num,recv_data%pos,recv_data%par_vel,recv_data%perp_vel,recv_data%kin_e,recv_data%p_angle,recv_data%rho,recv_data%tfinal
             IF (JTo4) write(50,*),recv_data%par_num,recv_data%s_pos,recv_data%s_par_vel,recv_data%perp_vel,recv_data%s_kin_e,recv_data%s_p_angle,recv_data%s_rho
             IF ((JTo4) .and. (recv_data%p_exit .eq. 5)) write(51,*),recv_data%par_num,recv_data%s_pos,0,0   ! note that when running on multiple processes the failed file doesn't get any parallel velocity or gamma info written
             IF (JTo4) write(52,*),recv_data%par_num,recv_data%p_exit
             NKEEP = (NOK +NBAD)/NSTORE
             if (modulo(recv_data%par_num,print_data) .eq. 0) print *, 'particle ', recv_data%par_num, ' running on rank ', status(mpi_source), ' written to output files'
          endif
       endif
        

      ! CALL WRITE_ENDTIME(RSTART,T2,MU,VPARSTART)

        if ((print_number .eq. 1) .and. (nproc .eq. 1)) then
          particle_tend = time8()
          if (particle_tend .ne. particle_tstart) print *,'particle ',pn,' time elapsed = ',particle_tend - particle_tstart
        endif
      
       pos_no_ekin=pos_no_ekin+1
      END DO
      pos_no_alpha=pos_no_alpha+1
      pos_no_ekin=0
     END DO
     pos_no_z=pos_no_z+1
     pos_no_alpha=0
    END DO
    pos_no_y=pos_no_y+1
    pos_no_z=0
   END DO
   pos_no_x=pos_no_x+1
   pos_no_y=0
  END DO

  if (nproc .ne. 1) then
     do j = 1,nproc-1
        !call mpi_probe(j,mpi_any_tag,mpi_comm_world,status,errcode)
        call mpi_recv(recv_data,1,mpi_finish_data2,mpi_any_source,mpi_any_tag,mpi_comm_world,status,errcode)
        if ((status(mpi_tag) .eq. 1)) then
             !print *, 'writing particle', recv_data%par_num, ' data'
             IF (JTo4) write(49,*),recv_data%par_num,recv_data%pos,recv_data%par_vel,recv_data%perp_vel,recv_data%kin_e,recv_data%p_angle,recv_data%rho,recv_data%tfinal
             IF (JTo4) write(50,*),recv_data%par_num,recv_data%s_pos,recv_data%s_par_vel,recv_data%perp_vel,recv_data%s_kin_e,recv_data%s_p_angle,recv_data%s_rho
             IF ((JTo4) .and. (recv_data%p_exit .eq. 5)) write(51,*),recv_data%par_num,recv_data%s_pos,0,0   ! note that when running on multiple processes the failed file doesn't get any parallel velocity or gamma info written
             IF (JTo4) write(52,*),recv_data%par_num,recv_data%p_exit
             NKEEP = (NOK +NBAD)/NSTORE
             if (modulo(recv_data%par_num,print_data) .eq. 0) print *, 'particle ', recv_data%par_num, ' running on rank ', status(mpi_source), ' written to output files'
          endif
     enddo
  endif
  if (nproc .ne. 1) then
     do j = 1,nproc-1
        call mpi_send(start_data,1,mpi_start_data,j,1,mpi_comm_world,errcode) ! tag 1 means time to stop, don't do anything with data
     enddo
  endif
  
  IF (JTo4) CLOSE(49)
  IF (JTo4) CLOSE(50)
  IF (JTo4) CLOSE(51)
  IF (JTo4) CLOSE(52)

  else
     call mpi_send(finish_data,1,mpi_finish_data2,0,0,mpi_comm_world,errcode) ! don't write this data (tag 0)
     do
        call mpi_recv(start_data,1,mpi_start_data,0,mpi_any_tag,mpi_comm_world,status,errcode)
        if (status(mpi_tag) .ne. 1) then
          CALL RKDRIVE(start_data%par_num,start_data%R_s,start_data%u_s,start_data%gamma_s,start_data%mu_s,&
            start_data%t1_s,start_data%t2_s,start_data%eps_s,start_data%h1_s,start_data%nok_s,&
            start_data%nbad_s,start_data%tt_s,start_data%s_s,start_data%total_s)
            !if (rank .eq. 6) call sleep(2)
        else
          exit
        endif
     enddo
  endif

 !CALL MAKEFILE(time_no)
  
 !END DO
 
 IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d")).OR.(str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN	!forget arrays at end
  CALL mpi_close                     ! mpi_routines.f90
  CALL MPI_FINALIZE(errcode)  
 ENDIF

 IF ((str_cmp(FMOD, "BOUR")).OR.(str_cmp(FMOD, "bour"))) THEN
  CALL bour_fini      		! deallocate stuff, leave everything nice and tidy
 ENDIF 

!------------------------------------------------------------------------------!
 Contains
!------------------------------------------------------------------------------!
SUBROUTINE JTMUcalc(mu,USTART,GAMMASTART, Ekin,alpha,RSTART,T1,T2)

  REAL(num), DIMENSION(3),INTENT(IN) :: RSTART
  REAL(num), INTENT(IN) :: T1,T2, Ekin, Alpha				! 
  REAL(num), INTENT(OUT) :: mu, gammastart, Ustart
  REAL(num), DIMENSION(3) :: B,El,a2,a3,a4,a5,a6,a7,a8,a9,a10,ue, RT
  !REAL(num) :: magB,vtot,vperp, vparstart!,Erest
  REAL(num) :: modB,vtot, gamma,rho,temperature,eta
 
 !calculate B, E, V at this point/time:
 !CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2)
 CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2,rho,temperature,eta)

 !calculate magnitude of B
 modB=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
 modB=sqrt(modB)
 RT=RSTART

 !print*, 'RT=', RT
 !PRINT*, 'B=', B
 !STOP
 !print*, modB

 !Erest = (M*c*c)*oneuponAQ
 
 ! E X B drift
 ue=cross(El,B)/dot(B,B)  !*0.5
 
 gamma=Ekin/m/c/c*AQ+1.0d0
 
 if (c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue) .le. 0.0_num) then
   print *, 'Warning: vtot calculation in JTMUcalc gives unphysical answer. This is most likely due to E cross B drift being faster than total energy of the particle. Need to rethink initial energy and/or location. Stopping.'
   stop
 endif
 ! vtot^2=vpar^2+vperp^2+UE^2!!, UE should NOT come out of Gamma!
 !vtot=c/gamma*sqrt(gamma*gamma-1)			! I *think* vtot is dimensional
 !vtot=sqrt(c*c/gamma/gamma*(gamma*gamma-1)-vscl*vscl*dot(ue,ue))
 vtot=sqrt(c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue))			! THRELFALL ET AL 2015
 ! vtot=sqrt(c*c-c/gamma*c/gamma-vce*vce*dot(ue,ue))			! THRELFALL ET AL 2015
 
 !Ustart=(c*sqrt(gamma*gamma-1))*cos(alpha)
 Ustart=gamma*vtot*cos(alpha)

 
 mu=0.5_num*m*vtot*vtot*sin(alpha)*sin(alpha)*gamma*gamma/modB		
 
 USTART=Ustart/vscl					! hence, non-dimensionalising..
 GAMMASTART=gamma
 mu=mu/m/vscl/vscl					! no bscl normalising factor - using normalised B's already!


END SUBROUTINE

END PROGRAM reltest
