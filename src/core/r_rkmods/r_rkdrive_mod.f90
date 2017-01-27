MODULE M_driverR

    USE global
    USE M_derivsR, ONLY: DERIVS 
    USE M_rkqsR, ONLY: RKQS
    USE M_fields, ONLY: FIELDS
    USE M_rkcolR
    USE M_products, ONLY: DOT, CROSS

IMPLICIT NONE

  PRIVATE
  PUBLIC :: RKDRIVE

  CONTAINS

SUBROUTINE RKDRIVE(PAR_NUM,RSTART,USTART,GAMMASTART,MU,T1,T2,EPS,H1,NOK,NBAD,TT,S,TOTAL)
 !##################################################################
 !Driver routine with adaptive stepsize control. It goes from T1 to
 !T2 with accuracy eps. Hmin is the minimum allowed stepsize. nok and 
 !nbad are the number of good and bad (i.e. retried) steps. RSTART is 
 !replaced by the end values.
!##################################################################
  
 IMPLICIT NONE
 INTEGER				:: NOK, NBAD
 INTEGER				:: I 
 INTEGER				:: UNDERFLOW
 INTEGER, INTENT(IN)                    :: par_num
 REAL(num), PARAMETER			:: TINY=1E-20
 REAL(num), INTENT(IN)			:: EPS, H1
 REAL(num), INTENT(INOUT)		:: T1,T2,MU
 REAL(num), INTENT(INOUT), DIMENSION(3) :: RSTART
 REAL(num), DIMENSION(3)		:: DRDT, R
 REAL(num), DIMENSION(3)		:: E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf
 REAL(num), DIMENSION(3)		:: bb, GRADB
 REAL(num), DIMENSION(3)		:: ENERGY
 REAL(num), DIMENSION(3)		:: UE
 REAL(num), DIMENSION(5)		:: RSCAL, RERR
 !REAL(num), DIMENSION(NKEEPMAX)		:: TT
 !REAL(num), DIMENSION(NKEEPMAX,3)	:: S, TOTAL
 REAL(num)              		:: TT           ! what does this variable do?
 REAL(num), DIMENSION(3)        	:: S, TOTAL     ! what do these variables do?
 REAL(num) 				:: H, HDID, HNEXT, T, H0, H_old
 REAL(num) 				:: U, USTART, DUDT, vpar
 REAL(num) 				:: GAMMA, GAMMASTART, DGAMMADT
 REAL(num) 				:: efct,e1,e2,e3, vtot_non, ek,rho,temperature,eta
 REAL(num) 				:: gyrofreq, gyroperiod, gyrorad, Epar, theta
 REAL(num)				:: MODB, oMODB,  DMODBDS, MODGRADB
 CHARACTER(LEN=30)			:: rvfilename
 CHARACTER(LEN=35)			:: tempfile, tempfile2, tempfile3
 integer                                :: collisions_initial
 real(num)                              :: U_old,mu_old,gamma_old,tmpnum

 call mpi_comm_size(mpi_comm_world,nproc,errcode)

  T=T1
  !TT(1) = T1
  TT = T1
  H=SIGN(H1,T2-T1)
  H0 = H
  NOK = 0
  NBAD = 0
  DO I = 1,3
   R(I) = RSTART(I)
   !S(1,I) = RSTART(I)
  ENDDO
  U = USTART
  GAMMA = GAMMASTART
  !DO I=1,3
  ! TOTAL(1,I) = 0.
  !END DO
UNDERFLOW=0

 efct=oneuponAQ
 IF (writervs) WRITE(rvfilename,"(A,'RV',I8.8,'.dat')"),dlocR,par_num    ! particleno in old code
 IF (writervs)  open(29,file=rvfilename,recl=1024,status='unknown')     	 
 IF ((JTo2).AND.(q.gt.0)) WRITE(tempfile2,"(A,'d',I8.8,'p.tmp')"),dlocR,par_num    !
 IF ((JTo2).AND.(q.lt.0)) WRITE(tempfile2,"(A,'d',I8.8,'e.tmp')"),dlocR,par_num    !
 IF (JTo2)  open(56,file=tempfile2,recl=1024,status='unknown')
 IF ((JTo3).AND.(q.gt.0)) WRITE(tempfile3,"(A,'f',I8.8,'p.tmp')"),dlocR,par_num    !
 IF ((JTo3).AND.(q.lt.0)) WRITE(tempfile3,"(A,'f',I8.8,'e.tmp')"),dlocR,par_num    !
 IF (JTo3)  open(57,file=tempfile3,recl=1024,status='unknown')
 
!print*, "R=", R
 CALL DERIVS (T, R, DRDT, U, DUDT,GAMMA,DGAMMADT,MU,T1,T2)
 !CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
 CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
 bb=B/sqrt(dot(B,B))
 
 UE=cross(E,B)/dot(B,B)
 MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
 oMODB=1.0_num/MODB
 GRADB(1) = DOT(B,DBDX)*oMODB
 GRADB(2) = DOT(B,DBDY)*oMODB
 GRADB(3) = DOT(B,DBDZ)*oMODB
 MODGRADB=SQRT(GRADB(1)*GRADB(1)+GRADB(2)*GRADB(2)+GRADB(3)*GRADB(3))
 
 !PRINT*, 'ini E-field:', E
 
 vpar=U/GAMMA
 Epar=dot(bb,E)

 e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
 e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))!*gamma*gamma	! don't think this gamma^2 should be here.
 e3=efct*0.5_num*M*dot(UE,UE)*vscl*vscl

 ek=efct*(gamma-1)*m*c*c

! vperpsq=sum((DRDT-VPAR*bb)**2)
! rL=MoQ*vscl*vperpsq/bscl/sqrt(dot(B,B))
 
 !JTmu=0.5_num*m*
! rL=MoQ*2.0_num*MU			! <- I don't think mu defined here has a mass in it (JT)
 !vperpsq=2.0_num*MU*sqrt(dot(B,B))*M
 !vtot_non = sqrt(Vpar*Vpar + 2.0_num*MU*sqrt(dot(B,B)))
 !rL=MoQ/sqrt(dot(B,B))*sqrt(2.0_num*((e1+e2+e3)-0.5_num*M*Vscl*Vpar*Vscl*Vpar)*oM)
 

 gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
 gyroperiod=1.0_num/gyrofreq
 !gyrorad=sqrt((e2+e3)*AQoM*2.0_num)/gyrofreq
 gyrorad=MoAQ*Vscl/Bscl/gamma*sqrt(2.0_num*MU/sqrt(dot(B,B)))
   theta = atan(sqrt(2.0_num*mu*modB)/U)
   if (theta .lt. 0) theta = pi + theta
 
  CALL DERIVS (T, R, DRDT, U, DUDT,GAMMA,DGAMMADT, MU, T1, T2)
  
  IF (writervs)  write(29,*)Tscl*(T-T1),	&   !1
  Lscl*R,						&   !2,3,4
  Vscl*VPAR,						&   !5
  MU*sqrt(dot(B,B)),					&   !6
  vscl*sqrt(sum((DRDT-U/gamma*bb)**2)),			&   !7 !d|r_perp|/dt ?	
  Escl*E,						&   !8,9,10
  Bscl*B,						&   !11,12,13
  vscl*vtot_non,					&   !14 vtot_non
  H*tscl,							&   !15
  e3,							&   !16
  ek,							&   !17  
  sqrt(dot(Vf,Vf))*Vscl,			        &   !18
  Escl*Epar,						&   !19
  vscl*UE,						&   !20,21,22
  vscl*U,						&   !23
  Vscl*Vf,						&   !24,25,26
  gyrofreq,gyroperiod,gyrorad,				&   !27,28,29
  theta, temperature*tempscl, eta*Vscl*Lscl*mu0_si          !30,31,32
    
  if (nproc .eq. 1) then
    IF (JTo4) write(50,*), par_num,R*Lscl, EPAR*Escl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl
  else
    finish_data%s_pos = R*Lscl
    finish_data%s_par_vel = VPAR*Lscl/Tscl
    finish_data%s_perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
    finish_data%s_kin_e = (gamma - 1.0)*m*c*c/qp
    finish_data%s_p_angle = theta
    finish_data%s_rho = rho*rhoscl
    finish_data%par_num = par_num
  endif
  !print *, 'rkdrive epsilon = ', m*Escl/(q*Bscl*Bscl*Lscl)

!****************************** Main Time-Loop Starts **************
  !PRINT *, "R=",R
  !PRINT *, "VPARSTART=",VPARSTART
  !PRINT *, "Tscl*(T-T1)=",Tscl*(T-T1)
  !PRINT *, "H=",H
  !PRINT *, "rkdrive B=",B
  
  H_old = H
  DO NSTP = 1, NSTPMAX
   if ((scattering .eq. 1) .and. (eta .ne. 0.0_num)) then
        H = H_old
        call rkcol(mu,gamma,U,B,DUDT,DBDT,UE,H,T,eta,temperature,rho)
   endif
     if (isnan(U)) then
        if (JTo4 .and. (nproc .eq. 1)) write(51,*) par_num,RSTART,USTART,GAMMASTART
        exit
     endif
   CALL DERIVS (T, R, DRDT, U, DUDT,GAMMA,DGAMMADT,MU, T1, T2)
   CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
   vpar=U/GAMMA
   ek=efct*(gamma-1)*m*c*c
   bb=B/sqrt(dot(B,B))
   UE=cross(E,B)/dot(B,B)
   Epar=dot(bb,E)
  !PRINT*, 'r_rkdrive Epar = ', Epar*Escl
  !print *, 'r_rkdrive epar 2= ', (E(1)*B(1) + E(2)*B(2) + E(3)*B(3))/sqrt(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
  !print *, 'rkdrive fields E = ', E, ' B = ', B


   MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
   oMODB=1.0_num/MODB

   
   e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
   e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))*gamma*gamma
   e3=efct*0.5_num*M*dot(UE,UE)*vscl*vscl
   
   gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
   gyroperiod=1.0_num/gyrofreq
   gyrorad=MoAQ*Vscl/Bscl/gamma*sqrt(2.0_num*MU/sqrt(dot(B,B)))
   theta = atan(sqrt(2.0_num*mu*modB)/U)
   if (theta .lt. 0) theta = pi + theta
   
   !PRINT*, '--------------------------'
   !print 667, NSTP,NSTPMAX, R, B, E
   !667 format (I9,'/',I9,' R=[',ES9.2,',',ES9.2,',',ES9.2,'], B=[',ES9.2,',',ES9.2,',',ES9.2,'], E=[',ES9.2,',',ES9.2,',',ES9.2,']')
   
    IF ((writervs).AND.(everystepswitch)) write(29,*)Tscl*(T-T1),	&   !1
    Lscl*R,						&   !2,3,4
    Vscl*VPAR,						&   !5
    MU*sqrt(dot(B,B)),					&   !6
    vscl*sqrt(sum((DRDT-U/gamma*bb)**2)),			&   !7 !d|r_perp|/dt ?	
    Escl*E,						&   !8,9,10
    Bscl*B,						&   !11,12,13
    vscl*vtot_non,					&   !14 vtot_non
    H*tscl,							&   !15
    e3,							&   !16
    ek,							&   !17  
    sqrt(dot(Vf,Vf))*Vscl,				&   !18
    Escl*Epar,						&   !19
    vscl*UE,						&   !20,21,22
    vscl*U,						&   !23
    Vscl*Vf,						&   !24,25,26
    gyrofreq,gyroperiod,gyrorad,			&   !27,28,29
    theta, temperature*tempscl, eta*Vscl*Lscl*mu0_si          !30,31,32
   
   DO I = 1,3       !Scaling used to monitor accuracy
    RSCAL(I) = ABS(R(I))+ABS(H*DRDT(I)) + TINY
   ENDDO
   RSCAL(4)=ABS(U)+ABS(H*DUDT) + TINY
   RSCAL(5)=ABS(GAMMA)+ABS(H*DGAMMADT) + TINY
   RSCAL =1

   IF((T+H-T2)*(T+H-T1) > 0.) THEN
    H=T2-T  			 !if stepsize can overshoot, decrease
   END IF 

   !print 668, R
   !668 format ('R2=[',ES9.2,',',ES9.2,',',ES9.2,']')

   !CALL RKQS(R,DRDT,U,DUDT,GAMMA,DGAMMADT,T,H,MU,EPS,RSCAL,HDID,HNEXT,T1,T2, UNDERFLOW)	! T modified here.
   if ((scattering .eq. 0) .or. (eta .eq. 0.0_num)) then
      collisions = 0
      CALL RKQS(R,DRDT,U,DUDT,GAMMA,DGAMMADT,T,H,MU,EPS,RSCAL,HDID,HNEXT,T1,T2,B,UNDERFLOW,RERR) ! T modified here.
   else if ((scattering .eq. 1) .and. (eta .ne. 0.0_num)) then
      collisions = 1
      CALL RKQS(R,DRDT,U,DUDT,GAMMA,DGAMMADT,T,H,MU,EPS,RSCAL,HDID,HNEXT,T1,T2,B,UNDERFLOW,RERR) ! T modified here.
   else
      print *, 'rkdrive: make up your mind on what timestep you want. stopping'
      stop
   endif
   !print *, 'rkdrive 246 H = ', H*tscl, HDID*tscl, HNEXT*tscl


   
   !print 680, R
   !680 format ('R4=[',ES9.2,',',ES9.2,',',ES9.2,']')
    
    MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
    oMODB=1.0_num/MODB
    GRADB(1) = DOT(B,DBDX)*oMODB
    GRADB(2) = DOT(B,DBDY)*oMODB
    GRADB(3) = DOT(B,DBDZ)*oMODB
    MODGRADB=SQRT(GRADB(1)*GRADB(1)+GRADB(2)*GRADB(2)+GRADB(3)*GRADB(3))
    
    DMODBDS=dot(B,B(1)*DBDX+B(2)*DBDY+B(3)*DBDZ)*oMODB*oMODB
    gyrorad=MoAQ*Vscl/Bscl/gamma*sqrt(2.0_num*MU/sqrt(dot(B,B)))/lscl
    
   ! print*, 'NSTP:', nstp, 'T=', t
   ! print*, '|b|=', modb
    
   IF (JTo2) write(56,*) NSTP, T, H, DUDT, DRDT, DGAMMADT
   IF (JTo3) write(57,*)	NSTP, B, DBDX, DBDY, DBDZ, E, &
   			   	DEDX, DEDY,DEDZ, DRDT, DUDT, DGAMMADT, &
   			   	Epar, T, R, U, MODGRADB, gyrorad   
   
   IF (HDID == H) THEN
    NOK = NOK+1
   ELSE
    NBAD = NBAD+1
   ENDIF

   IF ((MOD(NSTP,NSTORE)==0).AND.((NSTP/NSTORE)+1.GE.NKEEPMAX)) THEN	! JT fix to array overallocation in l.143
    !print 1001, T
    !1001 format ("ERROR: small timestep, not enough points in T array, exiting at T=",ES9.2)
    print*, 'ERROR: small timestep, not enough points in T array, EXITING..'
    !IF (JTo4) write(49,*), 'S'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, -2
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 0
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = -2
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2=T
    USTART = U
    GAMMASTART = GAMMA
 !print *, 'rkdrive finish_data = ', finish_data
    RETURN 
   ENDIF
    
  
!This is for storing every NSTORE step
   IF (MOD(NSTP,NSTORE)==0) THEN

    !TT((NSTP/NSTORE)+1) = T	
    TT = T	
    
    CALL DERIVS (T, R, DRDT, U, DUDT,GAMMA,DGAMMADT,MU,T1,T2)
    !CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
    CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
    
 
    IF (((bourdinflag).OR.(l3dflag).OR.(l2dflag)).AND.(SUM(E).EQ.0.0_num).AND.(SUM(B).EQ.0.0_num) &
 			      .AND.(SUM(DBDX).EQ.0.0_num).AND.(SUM(DBDY).EQ.0.0_num) &
			      .AND.(SUM(DBDZ).EQ.0.0_num).AND.(SUM(DEDX).EQ.0.0_num) &
			      .AND.(SUM(DEDY).EQ.0.0_num).AND.(SUM(DEDZ).EQ.0.0_num) &
!			      .AND.(SUM(DBDT).EQ.0.0_num).AND.(SUM(DEDT).EQ.0.0_num) &
    			      .AND.(SUM(Vf).EQ.0.0_num)) THEN	! not technically beyond bourdin range
    if ((print_number .eq. 1).and. (nproc .eq. 1)) print *, 'special box extent exit'
    !IF (JTo4) write(49,*), 'B'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, -1
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 0
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = -1
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    USTART = U
    GAMMASTART = GAMMA
 !print *, 'rkdrive finish_data = ', finish_data
    RETURN
   ENDIF
    
    
    
    bb=B/sqrt(dot(B,B))
    UE=cross(E,B)/dot(B,B)
    VPAR=U/GAMMA
    Epar=dot(bb,E)
    
    ENERGY(1)=VPAR*VPAR
    ENERGY(2)=MU*sqrt(dot(B,B))
    ENERGY(3)=sum((DRDT-VPAR*bb)**2)

    !DO I = 1,3
    !  S((NSTP/NSTORE)+1,I) = R(I)
    !  TOTAL((NSTP/NSTORE)+1,I) = ENERGY(I)
    !ENDDO

    e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
    e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))*gamma*gamma
    e3=efct*0.5_num*M*dot(UE,UE)*vscl*lscl
    ek=efct*(gamma-1)*m*c*c 

    gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
    gyroperiod=1.0_num/gyrofreq
    gyrorad=MoAQ*Vscl/Bscl/gamma*sqrt(2.0_num*MU/sqrt(dot(B,B)))
   theta = atan(sqrt(2.0_num*mu*modB)/U)
   if (theta .lt. 0) theta = pi + theta
    
      IF (writervs)  write(29,*)Tscl*(T-T1),	&   !1
      Lscl*R,						&   !2,3,4
      Vscl*VPAR,					&   !5
      MU*sqrt(dot(B,B)),				&   !6
      vscl*sqrt(sum((DRDT-U/gamma*bb)**2)),		&   !7 !d|r_perp|/dt ?	
      Escl*E,						&   !8,9,10
      Bscl*B,						&   !11,12,13
      vscl*vtot_non,					&   !14 vtot_non
      H*tscl,						&   !15
      e3,						&   !16
      ek,						&   !17  
      sqrt(dot(Vf,Vf))*Vscl,						&   !18
      Escl*Epar,					&   !19
      vscl*UE,						&   !20,21,22
      vscl*U,						&   !23
      Vscl*Vf,					        &   !24,25,26
      gyrofreq,gyroperiod,gyrorad,			&   !27,28,29
      theta, temperature*tempscl, eta*Vscl*Lscl*mu0_si          !30,31,32
    ENDIF	!(ends mod(nstp,nstore) loop)

    ! EXIT CRITERIA:
    ! normal exit (if time is up):
    IF((T-T2)*(T2-T1) >= 0.) THEN         !Are we done?
      if ((print_number .eq. 1).and. (nproc .eq. 1))  PRINT *, 'time exit'
      !IF (JTo4) write(49,*), 'T'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, 1 
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 0
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = 1
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
      !print*, Tscl*T
      DO I = 1,3
        RSTART(I)=R(I)
      ENDDO
    !  T2=T
      USTART = U
      GAMMASTART=GAMMA
 !print *, 'rkdrive finish_data = ', finish_data
      RETURN                            !normal exit
    ENDIF

! JT exit conditions:
   IF (((analyticalflag).OR.(l3dflag).OR.(l2dflag).OR.(bourdinflag)) &
   			    .AND.((R(1).GE.xe(2)).OR.(R(1).LE.xe(1)) &	! beyond simulation range
    			      .OR.(R(2).GE.ye(2)).OR.(R(2).LE.ye(1)) &
			      .OR.(R(3).GE.ze(2)).OR.(R(3).LE.ze(1)))) THEN
    if ((print_number .eq. 1).and. (nproc .eq. 1)) print *, 'box extent exit'
    !IF (JTo4) write(49,*), 'B'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, 2
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 0
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = 2
        !print *, 'rkdrive finish_data = ', finish_data
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    USTART = U
    GAMMASTART = GAMMA
    RETURN
   ENDIF
   IF ((abs(H).lt.EPS).AND.(abs(HNEXT).lt.EPS).AND.(scattering .eq. 0)) THEN ! both this and the next step are unbelievably small so quit before we get stuck!
    if ((print_number .eq. 1).and. (nproc .eq. 1)) print *, 'timestep shrink'
    !IF (JTo4) write(49,*), 'H'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, 3
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 0
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = 3
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    USTART = U
    GAMMASTART = GAMMA
 !print *, 'rkdrive finish_data = ', finish_data
    RETURN
   ENDIF
   IF (UNDERFLOW.EQ.1) THEN	! JT fix to array overallocation in l.143
    if ((print_number .eq. 1).and. (nproc .eq. 1)) print*, 'ERROR: timestep UNDERFLOW in RKQS, EXITING..'
    !IF (JTo4) write(49,*), 'U'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, 4
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 0
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = 4
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
    IF (JTo3) THEN
      WRITE(tempfile,"(A,'O',I8.8,'.ufl')"),dlocR,par_num    !
      open(55,file=tempfile,recl=1024,status='unknown')
    
      !CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
      CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
    
      vpar=U/GAMMA
      ek=efct*(gamma-1)*m*c*c
      bb=B/sqrt(dot(B,B))
      UE=cross(E,B)/dot(B,B)
      Epar=dot(bb,E)

      e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
      e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))*gamma*gamma
      e3=efct*0.5_num*M*dot(UE,UE)*vscl*vscl
   
      gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
      gyroperiod=1.0_num/gyrofreq
      gyrorad=MoAQ*Vscl/Bscl/gamma*sqrt(2.0_num*MU/sqrt(dot(B,B)))
   
      write(55,*)Tscl*(T-T1),				&   !1
      Lscl*R,						&   !2,3,4
      Vscl*VPAR,						&   !5
      MU*sqrt(dot(B,B)),					&   !6
      vscl*sqrt(sum((DRDT-U/gamma*bb)**2)),			&   !7 !d|r_perp|/dt ?	
      Escl*E,						&   !8,9,10
      Bscl*B,						&   !11,12,13
      vscl*vtot_non,					&   !14 vtot_non
      H,							&   !15
      e3,							&   !16
      ek,							&   !17  
      gamma,						&   !18
      Escl*Epar,						&   !19
      vscl*UE,						&   !20,21,22
      vscl*U,						&   !23
      Vscl*DRDT!,						&   !24,25,26
      !gyrofreq,gyroperiod,gyrorad				   !27,28,29
      CLOSE(55)
    ENDIF
    
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2=T
    USTART = U
    GAMMASTART = GAMMA
 !print *, 'rkdrive finish_data = ', finish_data
    RETURN 
   ENDIF
   
   ! timestep limiter
   IF (HNEXT.lt.0.001*(T2-T1)) THEN	
    H=HNEXT
   ENDIF
   
  ENDDO                      !if we get to nstpmax...
 if ((print_number .eq. 1).and. (nproc .eq. 1)) PRINT *, 'NSTP Loop ended'
 !IF (JTo4) 	write(49,*), 'N'
    IF (nproc .eq. 1) then
        IF (JTo4) write(49,*), par_num, R*Lscl, VPAR*Lscl/Tscl, 2.0_num*MU*modB*Bscl*vscl*vscl,(gamma - 1.0)*m*c*c/qp,theta,rho*rhoscl,T*Tscl
        IF (JTo4) write(52,*), par_num, 5
    ELSE
        finish_data%pos = R*Lscl
        finish_data%par_vel = VPAR*Lscl/Tscl
        finish_data%perp_vel = 2.0_num*MU*modB*Bscl*vscl*vscl
        finish_data%kin_e = (gamma - 1.0)*m*c*c/qp
        finish_data%p_angle = theta
        finish_data%p_failed = 1
        finish_data%rho = rho*rhoscl
        finish_data%tfinal = T*Tscl
        finish_data%p_exit = 5
        call mpi_send(finish_data,1,mpi_finish_data2,0,1,mpi_comm_world,errcode)
    ENDIF
 IF (JTo3)  	CLOSE(57)
 IF (JTo2)  	CLOSE(56)
 IF (writervs)	CLOSE(29)
 !print *, 'rkdrive finish_data = ', finish_data
 RETURN

END SUBROUTINE RKDRIVE

END MODULE M_driverR
