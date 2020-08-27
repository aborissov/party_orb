MODULE M_driverN

    USE global
    USE M_derivsN, ONLY: DERIVS 
    USE M_rkqsN, ONLY: RKQS
    USE M_fields, ONLY: FIELDS
    USE M_products, ONLY: DOT, CROSS

IMPLICIT NONE

  PRIVATE
  PUBLIC :: RKDRIVE

  CONTAINS

SUBROUTINE RKDRIVE(RSTART,VPARSTART,MU,T1,T2,EPS,H1,NOK,NBAD,TT,S,TOTAL)
 !##################################################################
 !Driver routine with adaptive stepsize control. It goes from T1 to
 !T2 with accuracy eps. Hmin is the minimum allowed stepsize. nok and 
 !nbad are the number of good and bad (i.e. retried) steps. RSTART is 
 !replaced by the end values.
!##################################################################
  
 IMPLICIT NONE
 !INTEGER, PARAMETER :: num = KIND(1.0)

 !max. no of functions and max no. of steps to be stored
 INTEGER				:: NOK, NBAD, I, UNDERFLOW 
 INTEGER, DIMENSION(6), SAVE 		:: bcroute = (/-1,-1,-1,-1,-1,-1/)
 REAL(num), PARAMETER			:: TINY=1E-20
 REAL(num), INTENT(IN)			:: EPS, H1,MU
 REAL(num), INTENT(INOUT)		:: T1,T2
 REAL(num), INTENT(INOUT), DIMENSION(3) :: RSTART
 REAL(num), DIMENSION(3)		:: DRDT, R
 REAL(num),  DIMENSION(3)		:: RLAST
 REAL(num), DIMENSION(3)		:: E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf
 REAL(num), DIMENSION(3)		:: bb, GRADB, ENERGY, UE, UGB, UC
 REAL(num), DIMENSION(4)		:: RSCAL 
 REAL(num), DIMENSION(NKEEPMAX)		:: TT
 REAL(num), DIMENSION(NKEEPMAX,3)	:: S, TOTAL
 REAL(num) 				:: H, HDID, HNEXT, T
 REAL(num) 				:: VPAR, VPARSTART, DVPARDT
 REAL(num) 				:: efct,e1,e2,e3, vtot_non
 REAL(num) 				:: gyrofreq, gyroperiod, gyrorad, Epar
 REAL(num)				:: MODB, oMODB,  DMODBDS, MODGRADB
 CHARACTER(LEN=30)			:: rvfilename
 CHARACTER(LEN=35)			:: tempf, tempf2, tempf3

  T=T1
  TT(1) = T1
  H=SIGN(H1,T2-T1)
  NOK = 0
  NBAD = 0
  DO I = 1,3
   R(I) = RSTART(I)
   S(1,I) = RSTART(I)
  ENDDO
  VPAR = VPARSTART
  DO I=1,3
   TOTAL(1,I) = 0.
  END DO
UNDERFLOW=0

 efct=oneuponAQ
 IF (writervs) WRITE(rvfilename,"(A,'RV',I8.8,'.dat')"),dlocN,pn    !
 IF (writervs)  open(29,file=rvfilename,recl=1024,status='unknown')     	 
 IF ((JTo2).AND.(q.gt.0)) WRITE(tempf2,"(A,'d',I8.8,'p.tmp')"),dlocN,pn    !
 IF ((JTo2).AND.(q.lt.0)) WRITE(tempf2,"(A,'d',I8.8,'e.tmp')"),dlocN,pn    !
 IF (JTo2)  open(56,file=tempf2,recl=1024,status='unknown')
 IF ((JTo3).AND.(q.gt.0)) WRITE(tempf3,"(A,'f',I8.8,'p.tmp')"),dlocN,pn    !
 IF ((JTo3).AND.(q.lt.0)) WRITE(tempf3,"(A,'f',I8.8,'e.tmp')"),dlocN,pn    !
 IF (JTo3)  open(57,file=tempf3,recl=1024,status='unknown')

!print*, "R=", R
 CALL DERIVS (T, R, DRDT, VPAR, DVPARDT,MU,T1,T2, UGB, UC)
 CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
 bb=B/sqrt(dot(B,B))
 UE=Vscl*cross(E,B)/dot(B,B)
 Epar=dot(E,bb)

 e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
 e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))
 e3=efct*0.5_num*M*dot(UE,UE)*Vscl*Vscl

! vperpsq=sum((DRDT-VPAR*bb)**2)
! rL=MoQ*vscl*vperpsq/bscl/sqrt(dot(B,B))
 
 !JTmu=0.5_num*m*
! rL=MoQ*2.0_num*MU			! <- I don't think mu defined here has a mass in it (JT)
 !vperpsq=2.0_num*MU*sqrt(dot(B,B))*M
 vtot_non = sqrt(Vpar*Vpar + 2.0_num*MU*sqrt(dot(B,B)))
 !rL=MoQ/sqrt(dot(B,B))*sqrt(2.0_num*((e1+e2+e3)-0.5_num*M*Vscl*Vpar*Vscl*Vpar)*oM)

    gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
    gyroperiod=1.0_num/gyrofreq
    gyrorad=MoAQ*Vscl/Bscl*sqrt(2.0_num*MU/sqrt(dot(B,B)))
 
  CALL DERIVS (T, R, DRDT, VPAR, DVPARDT,MU, T1, T2, UGB, UC)
! Initial dump of variables to file
  IF (writervs) write(29,*)Tscl*(T-T1),				&  !1
      Lscl*R,							&  !2,3,4
      Vscl*VPAR,						&  !5
      MU*sqrt(dot(B,B)),					&  !6
      vscl*sqrt(sum((DRDT-VPAR*bb)**2)),			&  !7
      Escl*E,							&  !8,9,10
      Bscl*B,							&  !11,12,13
      vscl*vtot_non,						&  !14   vtot_non
      e2,							&  !15
      e3,							&  !16
      0.5_num*vtot_non*vtot_non+dot(ue,ue),			&  !17  
      efct*Q*Escl*Vscl*dot(DRDT,E),				&  !18
      efct*M*Vscl*Lscl*MU*dot(B,DBDT)/sqrt(dot(B,B)),		&  !19 
      Vscl*Vf,							&  !20,21,22
      Escl*Epar,						&  !23
      vscl*UE,							&  !24,25,26
      gyrofreq,gyroperiod,gyrorad				   !27,28,29

!****************************** Main Time-Loop Starts **************
  DO NSTP = 1, NSTPMAX
   !PRINT*, NSTP
   CALL DERIVS (T, R, DRDT, VPAR, DVPARDT,MU, T1, T2, UGB, UC)
     
    vtot_non = sqrt(Vpar*Vpar + 2.0_num*MU*sqrt(dot(B,B)))
    gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
    gyroperiod=1.0_num/gyrofreq
    gyrorad=MoAQ*Vscl/Bscl*sqrt(2.0_num*MU/sqrt(dot(B,B)))
    e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
    e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))
    e3=efct*0.5_num*M*dot(UE,UE)*Vscl*Vscl
   
    IF ((writervs).AND.(everystepswitch)) write(29,*)Tscl*(T-T1),&  !1
      Lscl*R,							&  !2,3,4
      Vscl*VPAR,						&  !5
      MU*sqrt(dot(B,B)),					&  !6
      vscl*sqrt(sum((DRDT-VPAR*bb)**2)),			&  !7
      Escl*E,							&  !8,9,10
      Bscl*B,							&  !11,12,13
      vscl*vtot_non,							&  !14   vtot_non
      e2,							&  !15
      e3,							&  !16
      0.5_num*vtot_non*vtot_non+dot(ue,ue),			&  !17  
      efct*Q*Escl*Vscl*dot(DRDT,E),				&  !18
      efct*M*Vscl*Lscl*MU*dot(B,DBDT)/sqrt(dot(B,B)),		&  !19 
      Vscl*Vf,							&  !20,21,22
      !H,							&  !23
      Escl*Epar,						&  !23
      vscl*UE,							&  !24,25,26
      gyrofreq,gyroperiod,gyrorad				   !30,31,32
   

   DO I = 1,3       !Scaling used to monitor accuracy
    RSCAL(I) = ABS(R(I))+ABS(H*DRDT(I)) + TINY
   ENDDO
   RSCAL(4)=ABS(VPAR)+ABS(H*DVPARDT) + TINY
   RSCAL =1

   IF((T+H-T2)*(T+H-T1) > 0.) THEN
    H=T2-T  			 !if stepsize can overshoot, decrease
   END IF 

   
   RLAST=R
   CALL RKQS(R,DRDT,VPAR,DVPARDT,T,H,MU,EPS,RSCAL,HDID,HNEXT,T1,T2, UNDERFLOW)	! T modified here.
   
   IF ((analyticalflag).OR.(l3dflag).OR.(l2dflag).OR.(NLFFflag).OR. &
   (bourdinflag).OR.(testflag).OR.(FREflag).OR.(CMTflag).OR.(NLFFflag).OR.(MHDpflag))  THEN 
   !organising boundary conditions:
   ! CHECK X bottom and top, Y bottom and top, Z bottom and top.
   
    IF ((R(1).LE.spatial_extent_lower_bound(1)).AND.(DRDT(1).LT.0)) THEN	! if bottom bound AND heading downwards 
       101 SELECT CASE(bcroute(1))
       CASE(-1) ! setup on first use of subroutine
        IF (str_cmp(xlowbc, bcchoices(1))) THEN 	!transparent
         bcroute(1)=1
        ELSE IF (str_cmp(xlowbc, bcchoices(2))) THEN 	!partial reflective
         bcroute(1)=2  	
        ELSE IF (str_cmp(xlowbc, bcchoices(3))) THEN 	!fully reflective
         bcroute(1)=3   
        ELSE
         PRINT*, "confused about boundary condition choice"
         STOP
        END IF
        GO TO 101 ! now actually head back and select case we want!
       CASE(1)
        print *, 'box extent exit'
        IF (JTo4) write(49,*), 'X'
        DO I = 1,3
         RSTART(I)=R(I)
        ENDDO
        T2 = T
        VPARSTART=VPAR
        RETURN
       CASE(2)
        IF (sqrt(sum((DRDT-VPAR*bb)**2))/sqrt(VPAR*VPAR).ge.tanthetathresh) THEN
         VPAR=-VPAR	
	 Rlost=.FALSE.	 
	 R=RLAST
         print *, 'bounce' 
         !H=SIGN(H1,T2-T1)	! if there is a reflection, do we need to reset the step size?
         CYCLE
        ELSE
         print *, 'box extent exit'
         IF (JTo4) write(49,*), 'X'
         DO I = 1,3
          RSTART(I)=R(I)
         ENDDO
         T2 = T
	 VPARSTART=VPAR
         RETURN
        ENDIF
       CASE(3)
        VPAR=-VPAR
	Rlost=.FALSE.
	R=RLAST
        print *, 'bounce'
        !H=SIGN(H1,T2-T1)
        CYCLE
       END SELECT
    ENDIF
    IF ((R(1).GE.spatial_extent_upper_bound(1)).AND.(DRDT(1).GT.0)) THEN	! if top bound AND heading upwards
       102 SELECT CASE(bcroute(2))
       CASE(-1)
        IF (str_cmp(xupbc, bcchoices(1))) THEN 
         bcroute(2)=1
        ELSE IF (str_cmp(xupbc, bcchoices(2))) THEN
         bcroute(2)=2  
        ELSE IF (str_cmp(xupbc, bcchoices(3))) THEN
         bcroute(2)=3   
        ELSE
         PRINT*, "confused about boundary condition choice"
         STOP
        END IF
        GO TO 102 
       CASE(1)
        print *, 'box extent exit'
        IF (JTo4) write(49,*), 'X'
        DO I = 1,3
         RSTART(I)=R(I)
        ENDDO
        T2 = T	
	VPARSTART=VPAR
        RETURN
       CASE(2)
        IF (sqrt(sum((DRDT-VPAR*bb)**2))/sqrt(VPAR*VPAR).ge.tanthetathresh) THEN
         VPAR=-VPAR
	 Rlost=.FALSE.
	 R=RLAST
         print *, 'bounce' 
         !H=SIGN(H1,T2-T1)	! if there is a reflection, do we need to reset the step size?
         CYCLE
        ELSE
         print *, 'box extent exit'
         IF (JTo4) write(49,*), 'X'
         DO I = 1,3
          RSTART(I)=R(I)
         ENDDO
         T2 = T
	 VPARSTART=VPAR
         RETURN
        ENDIF
       CASE(3)
        VPAR=-VPAR	
	Rlost=.FALSE.
	R=RLAST
        print *, 'bounce'
        !H=SIGN(H1,T2-T1)
        CYCLE
       END SELECT
    ENDIF
    IF ((R(2).LE.spatial_extent_lower_bound(2)).AND.(DRDT(2).LT.0)) THEN	! if bottom bound AND heading downwards 
       103 SELECT CASE(bcroute(3))
       CASE(-1) ! setup on first use of subroutine
        IF (str_cmp(ylowbc, bcchoices(1))) THEN !transparent
         bcroute(3)=1
        ELSE IF (str_cmp(ylowbc, bcchoices(2))) THEN !partial reflective
         bcroute(3)=2  
        ELSE IF (str_cmp(ylowbc, bcchoices(3))) THEN !fully reflective
         bcroute(3)=3   
        ELSE
         PRINT*, "confused about boundary condition choice"
         STOP
        END IF
        GO TO 103 ! now actually head back and select case we want!
       CASE(1)
        print *, 'box extent exit'
        IF (JTo4) write(49,*), 'X'
        DO I = 1,3
         RSTART(I)=R(I)
        ENDDO
        T2 = T	 
	VPARSTART=VPAR
        RETURN
       CASE(2)
        IF (sqrt(sum((DRDT-VPAR*bb)**2))/sqrt(VPAR*VPAR).ge.tanthetathresh) THEN
         VPAR=-VPAR
	 Rlost=.FALSE.	
	 R=RLAST
         print *, 'bounce'
         !H=SIGN(H1,T2-T1)	
         CYCLE
        ELSE
         print *, 'box extent exit'
         IF (JTo4) write(49,*), 'X'
         DO I = 1,3
          RSTART(I)=R(I)
         ENDDO
         T2 = T	
	 VPARSTART=VPAR
         RETURN
        ENDIF
       CASE(3)
        VPAR=-VPAR	 
	R=RLAST
	Rlost=.FALSE.
        print *, 'bounce'
        !H=SIGN(H1,T2-T1)
        CYCLE
       END SELECT
    ENDIF
    IF ((R(2).GE.spatial_extent_upper_bound(2)).AND.(DRDT(2).GT.0)) THEN	! if top bound AND heading upwards
       104 SELECT CASE(bcroute(4))
       CASE(-1) 
        IF (str_cmp(yupbc, bcchoices(1))) THEN 
         bcroute(4)=1
        ELSE IF (str_cmp(yupbc, bcchoices(2))) THEN 
         bcroute(4)=2  
        ELSE IF (str_cmp(yupbc, bcchoices(3))) THEN
         bcroute(4)=3   
        ELSE
         PRINT*, "confused about boundary condition choice"
         STOP
        END IF
        GO TO 104 
       CASE(1)
        print *, 'box extent exit'
        IF (JTo4) write(49,*), 'X'
        DO I = 1,3
         RSTART(I)=R(I)
        ENDDO
        T2 = T
	VPARSTART=VPAR
        RETURN
       CASE(2)        
        IF (sqrt(sum((DRDT-VPAR*bb)**2))/sqrt(VPAR*VPAR).ge.tanthetathresh) THEN
         VPAR=-VPAR
	 Rlost=.FALSE.
	 R=RLAST
         print *, 'bounce' 
         !H=SIGN(H1,T2-T1)	! if there is a reflection, do we need to reset the step size?
         CYCLE
        ELSE
         print *, 'box extent exit'
         IF (JTo4) write(49,*), 'X'
         DO I = 1,3
          RSTART(I)=R(I)
         ENDDO
         T2 = T
 	 VPARSTART=VPAR
         RETURN
        ENDIF
       CASE(3)	
        R=RLAST
        VPAR=-VPAR
	Rlost=.FALSE.
        print *, 'bounce'
        !H=SIGN(H1,T2-T1)
        CYCLE
       END SELECT
    ENDIF   
    IF ((R(3).LE.spatial_extent_lower_bound(3)).AND.(DRDT(3).LT.0)) THEN	! if bottom bound AND heading downwards 
       105 SELECT CASE(bcroute(5))
       CASE(-1) ! setup on first use of subroutine
        IF (str_cmp(zlowbc, bcchoices(1))) THEN
         bcroute(5)=1
        ELSE IF (str_cmp(zlowbc, bcchoices(2))) THEN 
         bcroute(5)=2  
        ELSE IF (str_cmp(zlowbc, bcchoices(3))) THEN 
         bcroute(5)=3   
        ELSE
         PRINT*, "confused about boundary condition choice"
         STOP
        END IF
        GO TO 105
       CASE(1)
        print *, 'box extent exit'
        IF (JTo4) write(49,*), 'X'
        DO I = 1,3
         RSTART(I)=R(I)
        ENDDO
        T2 = T        
	VPARSTART=VPAR
        RETURN
       CASE(2)        
        IF (sqrt(sum((DRDT-VPAR*bb)**2))/sqrt(VPAR*VPAR).ge.tanthetathresh) THEN
         VPAR=-VPAR
	 R=RLAST
	 Rlost=.FALSE.
         print *, 'bounce' 
         !H=SIGN(H1,T2-T1)	! if there is a reflection, do we need to reset the step size?
         CYCLE
        ELSE
         print *, 'box extent exit'
         IF (JTo4) write(49,*), 'X'
         DO I = 1,3
          RSTART(I)=R(I)
         ENDDO
         T2 = T
	 VPARSTART=VPAR
         RETURN
        ENDIF
       CASE(3)
        VPAR=-VPAR
	R=RLAST
	Rlost=.FALSE.
        print*, 'bounce'
        !print*, vpar, R(3), DRDT(3)
	!H=SIGN(H1,T2-T1)
        CYCLE
       END SELECT
    ENDIF
    IF ((R(3).GE.spatial_extent_upper_bound(3)).AND.(DRDT(3).GT.0)) THEN	! if TOP zbound AND heading upwards 
       106 SELECT CASE(bcroute(6))
       CASE(-1)
        IF (str_cmp(zupbc, bcchoices(1))) THEN 
         bcroute(6)=1
        ELSE IF (str_cmp(zupbc, bcchoices(2))) THEN 
         bcroute(6)=2  
        ELSE IF (str_cmp(zupbc, bcchoices(3))) THEN
         bcroute(6)=3   
        ELSE
         PRINT*, "confused about boundary condition choice"
         STOP
        END IF
        GO TO 106
       CASE(1)
        print *, 'box extent exit'
        IF (JTo4) write(49,*), 'X'
        DO I = 1,3
         RSTART(I)=R(I)
        ENDDO
        T2 = T
        VPARSTART=VPAR
        RETURN
       CASE(2)
        IF (sqrt(sum((DRDT-VPAR*bb)**2))/sqrt(VPAR*VPAR).ge.tanthetathresh) THEN
        VPAR=-VPAR	
	R=RLAST
	Rlost=.FALSE. 
        print *, 'bounce'
         !H=SIGN(H1,T2-T1)	! if there is a reflection, do we need to reset the step size?
         CYCLE
        ELSE
         print *, 'box extent exit'
         IF (JTo4) write(49,*), 'X'
         DO I = 1,3
          RSTART(I)=R(I)
         ENDDO
         T2 = T
         VPARSTART=VPAR
         RETURN
        ENDIF
       CASE(3)
        VPAR=-VPAR	
	R=RLAST
	Rlost=.FALSE. 
        print *, 'bounce'
        !H=SIGN(H1,T2-T1)
        CYCLE
       END SELECT
    ENDIF      
   ENDIF
   
   
   
    MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
    oMODB=1.0_num/MODB
    GRADB(1) = DOT(B,DBDX)*oMODB
    GRADB(2) = DOT(B,DBDY)*oMODB
    GRADB(3) = DOT(B,DBDZ)*oMODB
    MODGRADB=SQRT(GRADB(1)*GRADB(1)+GRADB(2)*GRADB(2)+GRADB(3)*GRADB(3))
    
    DMODBDS=dot(B,B(1)*DBDX+B(2)*DBDY+B(3)*DBDZ)*oMODB*oMODB
    
   IF (JTo2) write(56,*) NSTP, T, H, DVPARDT, DRDT
   IF (JTo3) write(57,*) NSTP, B, DBDX, DBDY, DBDZ, E, &
   			   DEDX, DEDY,DEDZ, DRDT, DVPARDT, &
   			   	    Epar, T, R, VPAR, MODGRADB, gyrorad   
   
   IF (HDID == H) THEN
    NOK = NOK+1
   ELSE
    NBAD = NBAD+1
   ENDIF

   IF ((MOD(NSTP,NSTORE)==0).AND.((NSTP/NSTORE)+1.GE.NKEEPMAX)) THEN	! JT fix to array overallocation in l.143
    print*, 'ERROR: small timestep, not enough points in T array, exiting..'
    IF (JTo4) write(49,*), 'S'
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2=T
    VPARSTART = VPAR
    RETURN 
   ENDIF
!This is for storing every NSTORE step

   IF (MOD(NSTP,NSTORE)==0) THEN

    TT((NSTP/NSTORE)+1) = T	
    
    CALL DERIVS (T, R, DRDT, VPAR, DVPARDT,MU,T1,T2, UGB, UC)
    CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
    
    IF (((bourdinflag).OR.(l3dflag).OR.(l2dflag)).AND.(SUM(E).EQ.0.0_num).AND.(SUM(B).EQ.0.0_num) &
 			      .AND.(SUM(DBDX).EQ.0.0_num).AND.(SUM(DBDY).EQ.0.0_num) &
			      .AND.(SUM(DBDZ).EQ.0.0_num).AND.(SUM(DEDX).EQ.0.0_num) &
			      .AND.(SUM(DEDY).EQ.0.0_num).AND.(SUM(DEDZ).EQ.0.0_num) &
!			      .AND.(SUM(DBDT).EQ.0.0_num).AND.(SUM(DEDT).EQ.0.0_num) &
    			      .AND.(SUM(Vf).EQ.0.0_num)) THEN	! not technically beyond bourdin range
    print *, 'special box extent exit'
    IF (JTo4) write(49,*), 'X'
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    VPARSTART=VPAR
    RETURN
   ENDIF
    
    bb=B/sqrt(dot(B,B))
    UE=Vscl*cross(E,B)/dot(B,B)

    ENERGY(1)=VPAR*VPAR
    ENERGY(2)=MU*sqrt(dot(B,B))
    ENERGY(3)=sum((DRDT-VPAR*bb)**2)

    DO I = 1,3
      S((NSTP/NSTORE)+1,I) = R(I)
      TOTAL((NSTP/NSTORE)+1,I) = ENERGY(I)
    ENDDO

    e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
    e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))
    e3=efct*0.5_num*M*dot(UE,UE)*Vscl*Vscl 
    vtot_non = sqrt(Vpar*Vpar + 2.0_num*MU*sqrt(dot(B,B)))    
    gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
    gyroperiod=1.0_num/gyrofreq
    gyrorad=MoAQ*Vscl/Bscl*sqrt(2.0_num*MU/sqrt(dot(B,B)))
    
    IF (writervs) write(29,*)Tscl*(T-T1),		&  !1
      Lscl*R,							&  !2,3,4
      Vscl*VPAR,						&  !5
      MU*sqrt(dot(B,B)),					&  !6
      vscl*sqrt(sum((DRDT-VPAR*bb)**2)),			&  !7
      Escl*E,							&  !8,9,10
      Bscl*B,							&  !11,12,13
      vscl*vtot_non,							&  !14   vtot_non
      e2,							&  !15
      e3,							&  !16
      0.5_num*vtot_non*vtot_non+dot(ue,ue),			&  !17  
      efct*Q*Escl*Vscl*dot(DRDT,E),				&  !18
      efct*M*Vscl*Lscl*MU*dot(B,DBDT)/sqrt(dot(B,B)),		&  !19 
      Vscl*Vf,							&  !20,21,22
      !H,							&  !23
      Escl*Epar,							&  !23
      vscl*UE,							&  !24,25,26
      gyrofreq,gyroperiod,gyrorad				   !30,31,32
      !WARNING THE PREVIOUS IF IS A SINGLE LINE - THE FOLLOWING ENDIF ENDS THE MOD(NSTP,NSTORE) IF
    ENDIF

    ! EXIT CRITERIA:
    ! normal exit (if time is up):
    IF((T-T2)*(T2-T1) >= 0.) THEN         !Are we done?
      PRINT *, 'time exit'
      IF (JTo4) write(49,*), 'T'
      !print*, Tscl*T
      DO I = 1,3
        RSTART(I)=R(I)
      ENDDO
    !  T2=T
	VPARSTART=VPAR
      RETURN                            !normal exit
    ENDIF

   !IF (((analyticalflag).OR.(l3dflag).OR.(l2dflag).OR.(bourdinflag)) &
   !			    .AND.((R(1).GE.xe(2)).OR.(R(1).LE.xe(1)) &	! beyond simulation range
   ! 			      .OR.(R(2).GE.ye(2)).OR.(R(2).LE.ye(1)) &
!			      .OR.(R(3).GE.ze(2)).OR.(R(3).LE.ze(1)))) THEN
!    print *, 'box extent exit'
!    IF (JTo4) write(49,*), 'B'
!    DO I = 1,3
!      RSTART(I)=R(I)
!    ENDDO
!    T2 = T
!    VPARSTART=VPAR
!    RETURN
!   ENDIF

   
   !exit the calculation if |B| decreases beyond some threshold - remember rg->infinity as b->0
   IF (modb.le.lowbthresh) THEN 
    print *, 'modb is too small'
    IF (JTo4) write(49,*), 'B'  
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    VPARSTART=VPAR
    RETURN
   ENDIF   
   
   !exit the calculation if the gyroradius approaches the lengthscale of the simulations
   IF (gyrorad.ge.lscl) THEN 
    print *, 'gyroradius is as large as imposed lengthscale!'
    IF (JTo4) write(49,*), 'G' 
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    VPARSTART=VPAR
    RETURN
   ENDIF
   
   IF ((H.lt.EPS).AND.(HNEXT.lt.EPS)) THEN ! both this and the next step are unbelievably small so quit before we get stuck!
    print *, 'timestep shrink'
    IF (JTo4) write(49,*), 'H'
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2 = T
    VPARSTART=VPAR
    RETURN
   ENDIF
   
   IF (UNDERFLOW.EQ.1) THEN	! JT fix to array overallocation in l.143
    print*, 'ERROR: timestep UNDERFLOW in RKQS, EXITING..'
    IF (JTo4) write(49,*), 'U'
    IF (JTo3) THEN
      WRITE(tempf,"(A,'O',I8.8,'.ufl')"),dlocN,pn    !
      open(55,file=tempf,recl=1024,status='unknown')
   
   
      CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
      bb=B/sqrt(dot(B,B))
      UE=Vscl*cross(E,B)/dot(B,B)
      Epar=dot(E,bb)

      e1=efct*0.5_num*M*Vscl*Vpar*Vscl*Vpar
      e2=efct*M*Vscl*Vscl*MU*sqrt(dot(B,B))
      e3=efct*0.5_num*M*dot(UE,UE)*Vscl*Vscl

      vtot_non = sqrt(Vpar*Vpar + 2.0_num*MU*sqrt(dot(B,B)))

      gyrofreq=AQoM*Bscl*sqrt(dot(B,B))
      gyroperiod=1.0_num/gyrofreq
      gyrorad=MoAQ*Vscl/Bscl*sqrt(2.0_num*MU/sqrt(dot(B,B)))
 
      write(55,*)Tscl*(T-T1),		&  !1
      Lscl*R,							&  !2,3,4
      Vscl*VPAR,						&  !5
      MU*sqrt(dot(B,B)),					&  !6
      vscl*sqrt(sum((DRDT-VPAR*bb)**2)),			&  !7
      Escl*E,							&  !8,9,10
      Bscl*B,							&  !11,12,13
      vscl*vtot_non,						&  !14   vtot_non
      e2,							&  !15
      e3,							&  !16
      0.5_num*vtot_non*vtot_non+dot(ue,ue),			&  !17  
      efct*Q*Escl*Vscl*dot(DRDT,E),				&  !18
      efct*M*Vscl*Lscl*MU*dot(B,DBDT)/sqrt(dot(B,B)),		&  !19 
      Vscl*Vf,							&  !20,21,22
      Escl*Epar,						&  !23
      vscl*UE,							&  !24,25,26
      gyrofreq,gyroperiod,gyrorad				   !27,28,29
      CLOSE(55)
    ENDIF
   
    DO I = 1,3
      RSTART(I)=R(I)
    ENDDO
    T2=T
    VPARSTART=VPAR
    RETURN 
   ENDIF

   
  !timestep limiter
   IF (HNEXT.lt.0.001*(T2-T1)) THEN	
    H=HNEXT
   ENDIF

  ENDDO                      !if we get to nstpmax...
  PRINT *, 'NSTP Loop ended'
  IF (JTo4) 	write(49,*), 'N'
  IF (JTo3)  	CLOSE(57)
  IF (JTo2)  	CLOSE(56)
  IF (writervs)	CLOSE(29)
  RETURN
   
 RETURN

END SUBROUTINE RKDRIVE

END MODULE M_driverN
