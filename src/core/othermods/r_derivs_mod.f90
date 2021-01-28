MODULE M_derivsR
 USE global
 USE M_fields
 USE M_products 
 use bifrost_fields_mod, only: bifrost_grid
  
 IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: DERIVS
  PRIVATE :: JACOB
    
  CONTAINS
  !------------------------------------------------------------------------------!   
SUBROUTINE JACOB(T,R,U,GAMMA,MU,T1,T2,PD,BF_grid)
!+ requires fields module for dot and cross product functions to work.
!  calculates Jacobian for RHS's of particle equations
!
! Attempting to modify to be relativistic
 IMPLICIT NONE
 REAL(num), INTENT(IN)			:: T,MU, T1,T2, U, GAMMA
 !REAL(num), INTENT(OUT) 		:: DUDT, DGAMMADT
 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 !REAL(num), DIMENSION(3), INTENT(OUT)	:: DRDT
 REAL(num), DIMENSION(5,5), INTENT(OUT)	:: PD
 
 REAL(num), DIMENSION(3)		:: B,E,Vf
 REAL(num), DIMENSION(3)		:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3)		:: GRADB,DlittleBDT,DlittleBDX,DlittleBDY,DlittleBDZ
 REAL(num), DIMENSION(3)		:: UE, DUEDX,DUEDY,DUEDZ,DUEDT
 REAL(num), DIMENSION(3)		:: GRADDRIFT,DlittleBDS,UEGRADlittleB,UEGRADUE,DUEDS
 REAL(num), DIMENSION(3)		:: RELDRIFT1, RELDRIFT2, DRperpDT, d1, d2, d3, d4, d5
 REAL(num), DIMENSION(3)		:: DRPERPDT_DU, DRPERPDT_DGAMMA
 REAL(num)				:: MODB, oMODB, DMODBDS, EPAR,GRADBT, VPAR, DmodBDT
 REAL(num)				:: fac, facsq, ofac, ofacsq
 !LOGICAL, INTENT(OUT)			:: ERR 
 type(bifrost_grid)     :: BF_grid  ! bifrost grid

 CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2, BF_grid)
 MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
 
 oMODB=1.0_num/MODB					! 1/|B| - computationally less expensive to multiply..

 VPAR = U/GAMMA
 EPAR = DOT(E,B)*oMODB 
 
 fac=sqrt(1.0_num-vscl*vscl/c/c*EPAR*EPAR*oMODB*oMODB)
 facsq=fac*fac
 ofac=1.0_num/fac
 ofacsq=1.0_num/fac/fac


 GRADB(1) = DOT(B,DBDX)*oMODB
 GRADB(2) = DOT(B,DBDY)*oMODB
 GRADB(3) = DOT(B,DBDZ)*oMODB
 GRADBT   = DOT(B,DBDT)*oMODB

 DlittleBDX = DBDX*oMODB - B*GRADB(1)*oMODB*oMODB
 DlittleBDY = DBDY*oMODB - B*GRADB(2)*oMODB*oMODB
 DlittleBDZ = DBDZ*oMODB - B*GRADB(3)*oMODB*oMODB
 DlittleBDT = DBDT*oMODB - B*GRADBT*oMODB*oMODB

 UE = CROSS(E,B)*oMODB*oMODB	! UE = ExB/|B|^2

 !Need to get the derivatives of components of the UE:
 DUEDX = (CROSS(DEDX,B) + CROSS(E,DBDX) - 2.0_num*UE*DOT(B,DBDX))*oMODB*oMODB
 DUEDY = (CROSS(DEDY,B) + CROSS(E,DBDY) - 2.0_num*UE*DOT(B,DBDY))*oMODB*oMODB
 DUEDZ = (CROSS(DEDZ,B) + CROSS(E,DBDZ) - 2.0_num*UE*DOT(B,DBDZ))*oMODB*oMODB
 DUEDT = (CROSS(DEDT,B) + CROSS(E,DBDT) - 2.0_num*UE*DOT(B,DBDT))*oMODB*oMODB

 DMODBDS=dot(B,B(1)*DBDX+B(2)*DBDY+B(3)*DBDZ)*oMODB*oMODB

 GRADDRIFT = CROSS(B,GRADB)*oMODB*oMODB*ofacsq

 DlittleBDS = (B(1)*DlittleBDX + B(2)*DlittleBDY + B(3)*DlittleBDZ)*oMODB
 
 UEGRADlittleB = (UE(1)*DlittleBDX + UE(2)*DlittleBDY + UE(3)*DlittleBDZ)

 DUEDS = (B(1)*DUEDX + B(2)*DUEDY + B(3)*DUEDZ)*oMODB*oMODB

 UEGRADUE = UE(1)*DUEDX + UE(2)*DUEDY + UE(3)*DUEDZ

 
 ! DUDT = Omscl*tscl*EPAR - MU/gamma*DMODBDS*fac + gamma*DOT(UE,DlittleBDT)

  DmodBDT=GRADBT-dot(B,DlittleBDT)
  RELDRIFT1=U*DlittleBDT+gamma*DUEDT+vscl*vscl/c/c*MU/gamma*UE*DmodBDT*fac
  !RELDRIFT1=MU/gamma*GRADDRIFT+U*DlittleBDT+gamma*DUEDT+vscl*vscl/c/c*mu/gamma*UE*DmodBDT				!(assuming eperp<<B)
  RELDRIFT2=U/gamma*EPAR*UE
 
  d1=MU/gamma*GRADDRIFT
  d2=MU/gamma*vscl*vscl/c/c*CROSS(B,UE)*omodB*omodB*ofacsq*DmodBDT
  d3=U*CROSS(B,DlittleBDT)*omodB*omodB*ofacsq
  d4=gamma*CROSS(B,DUEDT)*omodB*omodB*ofacsq
  d5=Epar*U/gamma*vscl*vscl/c/c*CROSS(B,UE)*omodB*omodB*ofacsq
 
 !DRperpDT= UE + oneuponOmscl/tscl*(d1+d2+d3+d4)+d5
   
!  DRDT =  DRperpDT + (U/gamma)*B*oMODB
    
!  DgammaDT=Omscl*tscl*vscl*vscl/c/c*(dot(DRperpDT,E)+U/gamma*dot(B,E)*oMODB)+vscl*vscl/c/c*MU/gamma*DmodBDT*fac

!each element is du, dr(x3), dgamma
 PD(1,1)=0.0_num
 PD(1,2)=0.0_num
 PD(1,3)=0.0_num
 PD(1,4)=0.0_num
 PD(1,5)=DOT(UE,DlittleBDT)+MU/gamma/gamma*DMODBDS*fac
 
 !note extra bits on ends to reflect DRDT and not DRPERPDT
 DRPERPDT_DU(1)=oneuponOmscl/tscl*(d3(1)/U)+d5(1)/U
 DRPERPDT_DU(2)=oneuponOmscl/tscl*(d3(2)/U)+d5(2)/U
 DRPERPDT_DU(3)=oneuponOmscl/tscl*(d3(3)/U)+d5(3)/U
 DRPERPDT_DGAMMA(1)=oneuponOmscl/tscl*(-d1(1)/gamma-d2(1)/gamma+d4(1)/gamma)-d5(1)/gamma
 DRPERPDT_DGAMMA(2)=oneuponOmscl/tscl*(-d1(2)/gamma-d2(2)/gamma+d4(2)/gamma)-d5(2)/gamma
 DRPERPDT_DGAMMA(3)=oneuponOmscl/tscl*(-d1(3)/gamma-d2(3)/gamma+d4(3)/gamma)-d5(3)/gamma
 
 
 PD(2,1)=DRPERPDT_DU(1)+B(1)*oMODB/gamma
 PD(2,2)=0.0_num
 PD(2,3)=0.0_num
 PD(2,4)=0.0_num
 PD(2,5)=DRPERPDT_DGAMMA(1)-B(1)*oMODB/gamma/gamma*U
 PD(3,1)=DRPERPDT_DU(2)+B(2)*oMODB/gamma
 PD(3,2)=0.0_num
 PD(3,3)=0.0_num
 PD(3,4)=0.0_num
 PD(3,5)=DRPERPDT_DGAMMA(2)-B(2)*oMODB/gamma/gamma*U
 PD(4,1)=DRPERPDT_DU(3)+B(3)*oMODB/gamma
 PD(4,2)=0.0_num
 PD(4,3)=0.0_num
 PD(4,4)=0.0_num
 PD(4,5)=DRPERPDT_DGAMMA(3)-B(3)*oMODB/gamma/gamma*U

 
 PD(5,1)=vscl*vscl/c/c*(Omscl*tscl*DOT(DRPERPDT_DU,E)+Omscl*tscl*DOT(B,E)*oMODB/gamma)
 PD(5,2)=0.0_num
 PD(5,3)=0.0_num
 PD(5,4)=0.0_num
 PD(5,5)=vscl*vscl/c/c*(Omscl*tscl*DOT(DRPERPDT_DGAMMA,E)-Omscl*tscl/gamma/gamma*U*DOT(B,E)*oMODB-MU/gamma/gamma*DmodBDT*fac)

END SUBROUTINE JACOB
!------------------------------------------------------------------------------!   
SUBROUTINE DERIVS (T, R, DRDT, U, DUDT,GAMMA, DGAMMADT,MU, T1, T2, BF_grid)
!+ requires fields module for dot and cross product functions to work.
! Calculates individual components of equations (44) and (45) of Guiliani et al (2005)
!
! Attempting to modify to be relativistic
 IMPLICIT NONE
 REAL(num), INTENT(IN)			:: T,MU, T1,T2, U, GAMMA
 REAL(num), INTENT(OUT) 		:: DUDT, DGAMMADT
 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), DIMENSION(3), INTENT(OUT)	:: DRDT
 !REAL(num), DIMENSION(3), INTENT(OUT)	:: UGB, UC
 
 REAL(num), DIMENSION(3)		:: B,E,Vf
 REAL(num), DIMENSION(3)		:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3)		:: GRADB,DlittleBDT,DlittleBDX,DlittleBDY,DlittleBDZ
 REAL(num), DIMENSION(3)		:: VE, DVEDX,DVEDY,DVEDZ,DVEDT, UE
 REAL(num), DIMENSION(3)		:: GRADDRIFT,DlittleBDS,VEGRADlittleB,VEGRADVE,DVEDS
 REAL(num), DIMENSION(3)		:: OTHERS
 REAL(num), DIMENSION(3)		:: RELDRIFT1, RELDRIFT2, DRperpDT, d1, d2, d3, d4, d5
 !REAL(num), DIMENSION(3)		:: SCRAE
 REAL(num)				:: MODB, oMODB, DMODBDS, EPAR,GRADBT, VPAR, DmodBDT
 REAL(num)				:: fac, facsq, ofac, ofacsq
 !LOGICAL, INTENT(OUT)			:: ERR 

 type(bifrost_grid)     :: BF_grid  ! bifrost grid 

 ! ALEXEI: debugging
 integer  :: i

 ! ALEXEI: debugging
 call check_nan(T)
 call check_nan(U)
 call check_nan(GAMMA)
 call check_nan(DGAMMADT)
 call check_nan(MU)
 call check_nan(T1)
 call check_nan(T2)
 do i = 1, 3
   call check_nan(R(i))
 enddo

 CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2, BF_grid)
 MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
 
 oMODB=1.0_num/MODB					! 1/|B| - computationally less expensive to multiply..

 VPAR = U/GAMMA
 EPAR = DOT(E,B)*oMODB 
 
 fac=sqrt(1.0_num-vscl*vscl/c/c*EPAR*EPAR*oMODB*oMODB)
 facsq=fac*fac
 ofac=1.0_num/fac
 ofacsq=1.0_num/fac/fac


 GRADB(1) = DOT(B,DBDX)*oMODB
 GRADB(2) = DOT(B,DBDY)*oMODB
 GRADB(3) = DOT(B,DBDZ)*oMODB
 GRADBT   = DOT(B,DBDT)*oMODB

 DlittleBDX = DBDX*oMODB - B*GRADB(1)*oMODB*oMODB
 DlittleBDY = DBDY*oMODB - B*GRADB(2)*oMODB*oMODB
 DlittleBDZ = DBDZ*oMODB - B*GRADB(3)*oMODB*oMODB
 DlittleBDT = DBDT*oMODB - B*GRADBT*oMODB*oMODB

 VE = CROSS(E,B)*oMODB*oMODB	! VE = ExB/|B|^2

 !Need to get the derivatives of components of the VE:
 DVEDX = (CROSS(DEDX,B) + CROSS(E,DBDX) - 2.0_num*VE*DOT(B,DBDX))*oMODB*oMODB
 DVEDY = (CROSS(DEDY,B) + CROSS(E,DBDY) - 2.0_num*VE*DOT(B,DBDY))*oMODB*oMODB
 DVEDZ = (CROSS(DEDZ,B) + CROSS(E,DBDZ) - 2.0_num*VE*DOT(B,DBDZ))*oMODB*oMODB
 DVEDT = (CROSS(DEDT,B) + CROSS(E,DBDT) - 2.0_num*VE*DOT(B,DBDT))*oMODB*oMODB

 DMODBDS=dot(B,B(1)*DBDX+B(2)*DBDY+B(3)*DBDZ)*oMODB*oMODB

 GRADDRIFT = CROSS(B,GRADB)*oMODB*oMODB*ofacsq
 
! PRINT*, 'DER: B', B
! PRINT*, 'DER: GRADB', GRADB  
! PRINT*, 'DER: GRADDRIFT', GRADDRIFT 

 DlittleBDS = (B(1)*DlittleBDX + B(2)*DlittleBDY + B(3)*DlittleBDZ)*oMODB
 
 VEGRADlittleB = (VE(1)*DlittleBDX + VE(2)*DlittleBDY + VE(3)*DlittleBDZ)

 DVEDS = (B(1)*DVEDX + B(2)*DVEDY + B(3)*DVEDZ)*oMODB*oMODB

 VEGRADVE = VE(1)*DVEDX + VE(2)*DVEDY + VE(3)*DVEDZ
 
!All the acceleration drift terms, to be crossed with B
! ACCDRIFT = VPAR*DlittleBDT + VPAR*VPAR*DlittleBDS + VPAR*VEGRADlittleB + DVEDT + VPAR*DVEDS + VEGRADVE


!all the terms that make up the last bit of the parallel equation
 OTHERS = DlittleBDT + VPAR*DlittleBDS + VEGRADlittleB

!-------------------------------------------------------------------------------!
! The final equations of motion (eqs. 44 and 45 from Guiliani et al. (2005)	!
!-------------------------------------------------------------------------------!
!  DVPARDT = Omscl*tscl*EPAR - MU*DMODBDS + DOT(UE,OTHERS)		! My normalisation!!
!  DRDT = UE + oneuponOmscl/tscl*(MU*GRADDRIFT + CROSS(B,ACCDRIFT)*oMODB*oMODB) + VPAR*B*oMODB		! extra term due to Rperp->R!
!
!  UGB=oneuponOmscl/tscl*MU*GRADDRIFT
!  UC =oneuponOmscl/tscl*VPAR*VPAR*cross(B,DlittleBDS)*oMODB

!-----------------------------------------------------------------------!
! Relativistic versions (eqs. 1.76-1.79 from Northrop (1963)		!
!-----------------------------------------------------------------------!
!now 3 variables - U(=vpar*gamma), R and gamma
! will also attempt to include full 1-epar^2/b^2 thing
! - have included eperp<<B version so you can see where the non-dimensional version might have messed up things..
 
  DUDT = Omscl*tscl*EPAR - MU/gamma*DMODBDS*fac + gamma*DOT(VE,OTHERS)                          ! THRELFALLETAL 2015
  
  DmodBDT=GRADBT-dot(B,DlittleBDT)
  RELDRIFT1=U*DlittleBDT+gamma*DVEDT+vscl*vscl/c/c*MU/gamma*VE*DmodBDT*fac
  RELDRIFT2=U/gamma*EPAR*VE
  
  !PRINT*, 'DER: mu/gamma', MU/gamma
  d1=MU/gamma*GRADDRIFT
  d2=MU/gamma*vscl*vscl/c/c*CROSS(B,VE)*omodB*omodB*ofacsq*DmodBDT
  d3=U*CROSS(B,DlittleBDT)*omodB*omodB*ofacsq
  d4=gamma*CROSS(B,DVEDT)*omodB*omodB*ofacsq
  d5=Epar*U/gamma*vscl*vscl/c/c*CROSS(B,VE)*omodB*omodB*ofacsq
  
  DRperpDT= VE + oneuponOmscl/tscl*(d1+d2+d3+d4)+d5      
  
  DRDT =  DRperpDT + (U/gamma)*B*oMODB
  
  DgammaDT=Omscl*tscl*vscl*vscl/c/c*(dot(DRperpDT,E)+U/gamma*dot(B,E)*oMODB)+vscl*vscl/c/c*MU/gamma*DmodBDT*fac   !THRELFALL ET AL 2015
  return

!-----------------------------------------------------------------------!
! Relativistic versions WITH Birn et al (2004) Normalisation		!
!-----------------------------------------------------------------------!
!  UE=gamma*VE
!
!  !DUDT = sigma*DOT(UE,DlittleBDT)-oneosigma/myepsil*EPAR - oneosigma*MU/gamma*DMODBDS*fac			! full eqn
!  DUDT = -oneosigma/myepsil*EPAR - oneosigma*MU/gamma*DMODBDS*fac						 ! ignore sigmas
!  
!  GRADDRIFT = CROSS(B,GRADB)*oMODB*oMODB*ofacsq
!  
! 
!  d1=0.5_num*MU/gamma*GRADDRIFT
!  !d3=sigma*U*CROSS(B,DlittleBDT)*omodB*omodB*ofacsq
!  !d4=sigma*sigma*gamma*CROSS(B,DVEDT)*omodB*omodB*ofacsq
!  d2=vscl*vscl/c/c*MU/gamma*CROSS(B,VE)*omodB*omodB*ofacsq*DmodBDT
!  
!  d5=vscl*vscl/c/c*Epar*U/gamma*CROSS(B,VE)*omodB*omodB*ofacsq*oneosigma
!  
!  !DRperpDT= VE-myepsil*(d1+d2+d3+d4)+d5								! full equation
!  DRperpDT= VE-myepsil*(d1+d2)+d5									! ignore sigma, sigma^2
!  
!  DRDT =  DRperpDT + (U/gamma)*B*oMODB
!  
!  DgammaDT=vscl*vscl/c/c*oneosigma*oneosigma*(-oneomyepsil*(dot(DRperpDT,E)+oneosigma*U/gamma*dot(B,E)*oMODB)+MU/gamma*DmodBDT*fac)

END SUBROUTINE DERIVS
!-----------------------------------------------------------------------------!   
END MODULE M_derivsR
