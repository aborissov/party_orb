MODULE M_rkqsR 

   USE global
   USE M_rkckR
   USE M_rkcolR
   IMPLICIT NONE

   CONTAINS

! SUBROUTINE RKQS (R,DRDT,U,DUDT,GAMMA,DGAMMADT, T,HTRY,MU,EPS,RSCAL,HDID,HNEXT,T1,T2)
 SUBROUTINE RKQS (R,DRDT,U,DUDT,GAMMA,DGAMMADT, T,HTRY,MU,EPS,RSCAL,HDID,HNEXT,T1,T2,B,uflag,RERR)

!##################################################################### 
 !this subroutine is the stepper and basically calls rkck to take one
 !step while monitoring truncation error to ensure accuracy. It adjusts
 !stepsize. Inputs are the what you'd expect. Outputs are new values of R,
 ! DRDT, T, hdid and hnext.
!###################################################################

 IMPLICIT NONE
 REAL(num), INTENT(INOUT), DIMENSION(3)	:: R, DRDT
 REAL(num), INTENT(INOUT)		:: U, DUDT, GAMMA, DGAMMADT
 REAL(num), INTENT(INOUT) 		:: T
 REAL(num), INTENT(IN)			:: EPS, HTRY, T1, T2
 REAL(num), INTENT(IN), DIMENSION(5)	:: RSCAL
 REAL(num), INTENT(OUT)			:: HDID, HNEXT
 REAL(num), DIMENSION(5)		:: RERR 
 REAL(num), DIMENSION(3)		:: RTEMP
 REAL(num)				:: UTEMP,GAMMATEMP,gamma_old,mu,mu_old,nu,nu0,alpha,alpha_old,modB,vtherm
 REAL(num)				:: ERRMAX, H, HTEMP, TNEW, HNEW, H2
 REAL(num), PARAMETER			:: SAFETY=0.9_num, PGROW=-0.2_num
 REAL(num), PARAMETER			:: PSHRINK=-0.25_num, ERRCON=1.89e-4
 !JT adds these for debugging:
 REAL(num), DIMENSION(3)		:: E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf
 REAL(num), DIMENSION(3)		:: bb
 REAL(num), DIMENSION(3)		:: ENERGY
 REAL(num), DIMENSION(3)		:: UE
 REAL(num) 				:: efct,e1,e2,e3, vtot_non, ek, vpar,rho,temperature,eta
 REAL(num) 				:: gyrofreq, gyroperiod, gyrorad, Epar
 CHARACTER(LEN=35)			:: tempfile
 INTEGER, INTENT(OUT)			:: uflag
 uflag=0
  
 call FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
 modB = sqrt(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
 UE = CROSS(E,B)/(modB*modB)
 
 H=HTRY                   !Initial value for stepsize
 mu_old = mu
 gamma_old = gamma
 alpha = atan(sqrt(2.0_num*mu*modB)/U)
 if (alpha .lt. 0) alpha = pi + alpha
 !print *, 'RKQS alpha = ', alpha
 alpha_old = alpha

 if (adjust_timestep .eq. 1) then
 DO 
 !PRINT *,'rkqs_R = ',R
  mu = mu_old
  gamma = gamma_old
  alpha = alpha_old


  CALL RKCK(R,DRDT,U,DUDT,GAMMA,DGAMMADT,T,H,MU,RTEMP,UTEMP,GAMMATEMP, RERR,T1,T2)
  !print 669, R, DRDT
  ! 669 format ('R3=[',ES9.2,',',ES9.2,',',ES9.2,'], DRDT=[',ES9.2,',',ES9.2,',',ES9.2,']')
  ERRMAX=maxval(abs(RERR(:)/RSCAL(:)))/EPS 
  !print *,'RKQS DUDT = ', DUDT
  !WRITE(111,*) RERR,ERRMAX,H,T
  !WRITE(*,1001)"time:",Tscl*(T-t1)
  !WRITE(*,1004)"RERR:",RERR
  !WRITE(*,1004)"RSCAL:",RSCAL
  !WRITE(*,1004)"abs(RE/RS):",abs(RERR(:)/RSCAL(:))
  !WRITE(*,1001)"ERRMAX:",ERRMAX
  !WRITE(*,*)'****'
  !1004 FORMAT (a,5E13.4)
  !1001 FORMAT (a,E13.4)
  !print *,NSTP,H,ERRMAX
  
  IF (ERRMAX <= 1.0) EXIT
   HTEMP=SAFETY*H*(ERRMAX**PSHRINK)
   
  HNEW=sign(max(abs(HTEMP),0.1_num*abs(H)),H)  
  !PRINT 1013,  abs(HTEMP), 0.1_num*abs(H), H, HNEW
  !1013 FORMAT ("abs(HTEMP)=",E11.4," 0.1|H|=",E11.4," H=",E11.4, " HNEW=",E11.4) 
  
  H=sign(max(abs(HTEMP),0.1_num*abs(H)),H)
  TNEW=T+H
  IF (TNEW == T) THEN
   PRINT *, 'STEPSIZE UNDERFLOW IN RKQS'
   PRINT *, 'Particle No.',pn           !particleno in old code
   PRINT *, 'T & TNEW',T,TNEW
   uflag=1
   EXIT
  END IF
 END DO

 IF (ERRMAX > ERRCON) THEN
  HNEXT=SAFETY*H*(ERRMAX**PGROW)
 ELSE
  HNEXT=5.0_num*H
 END IF
 HDID=H
 T=T+H
 R(:)=RTEMP(:) 
 U = UTEMP
 !print *, 'RKQS U = ', U
 GAMMA = GAMMATEMP

!print 667, R
!   667 format ('R=[',ES9.2,',',ES9.2,',',ES9.2,']')
 
 else if (adjust_timestep .eq. 0) then
  H2 = H
  !if (eta .ne. 0) H2 = H*1.0E-2_num
  CALL RKCK(R,DRDT,U,DUDT,GAMMA,DGAMMADT,T,H2,MU,RTEMP,UTEMP,GAMMATEMP, RERR,T1,T2)
  T = T+H2
  HDID = H2
  HNEXT = H
  R(:)=RTEMP(:) 
  U = UTEMP
  GAMMA = GAMMATEMP
  !print *, 'rkqs H = ', H2
 end if

 END SUBROUTINE RKQS

END MODULE M_rkqsR 
