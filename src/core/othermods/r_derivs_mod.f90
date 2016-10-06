MODULE M_derivsR
 USE global
 USE M_fields
 USE M_products 
  
 IMPLICIT NONE
  
  CONTAINS
!------------------------------------------------------------------------------!   
SUBROUTINE DERIVS (T, R, DRDT, U, DUDT,GAMMA, DGAMMADT,MU, T1, T2)
!+ requires fields module for dot and cross product functions to work.
! Calculates individual components of equations (44) and (45) of Guiliani et al (2005)
!
! Attempting to modify to be relativistic
 IMPLICIT NONE
 REAL(num), INTENT(IN)			:: T,MU, T1,T2, U
 REAL(num), INTENT(OUT) 		:: DUDT, DGAMMADT, GAMMA
 REAL(num)                              :: rho,temperature,eta
 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), DIMENSION(3), INTENT(OUT)	:: DRDT
 !REAL(num), DIMENSION(3), INTENT(OUT)	:: UGB, UC
 
 REAL(num), DIMENSION(3)		:: B,E,Vf
 REAL(num), DIMENSION(3)		:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3)		:: GRADB,DlittleBDT,DlittleBDX,DlittleBDY,DlittleBDZ
 REAL(num), DIMENSION(3)		:: VE, DVEDX,DVEDY,DVEDZ,DVEDT, UE
 REAL(num), DIMENSION(3)		:: GRADDRIFT,DlittleBDS,VEGRADlittleB,VEGRADVE,DVEDS
 REAL(num), DIMENSION(3)		:: ACCDRIFT, OTHERS
 REAL(num), DIMENSION(3)		:: RELDRIFT1, RELDRIFT2, DRperpDT, d1, d2, d3, d4, d5
 !REAL(num), DIMENSION(3)		:: SCRAE
 REAL(num)				:: MODB, oMODB, DMODBDS, EPAR,GRADBT, VPAR, DmodBDT, GAMMAOLD
 REAL(num)				:: fac, facsq, ofac, ofacsq
 !LOGICAL, INTENT(OUT)			:: ERR 

 CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
 MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
 
 oMODB=1.0_num/MODB					! 1/|B| - computationally less expensive to multiply..
 VE = CROSS(E,B)*oMODB*oMODB	! VE = ExB/|B|^2
 UE = gamma*VE                                          ! using gamma implicitly?
 CALL gammacalc(gamma,U,UE,mu,modB)

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

 !PRINT *,'E = ',E,'B = ',B,'OMODB = ',oMODB

 !Need to get the derivatives of components of the VE:
 DVEDX = (CROSS(DEDX,B) + CROSS(E,DBDX) - 2.0_num*VE*DOT(B,DBDX))*oMODB*oMODB
 DVEDY = (CROSS(DEDY,B) + CROSS(E,DBDY) - 2.0_num*VE*DOT(B,DBDY))*oMODB*oMODB
 DVEDZ = (CROSS(DEDZ,B) + CROSS(E,DBDZ) - 2.0_num*VE*DOT(B,DBDZ))*oMODB*oMODB
 DVEDT = (CROSS(DEDT,B) + CROSS(E,DBDT) - 2.0_num*VE*DOT(B,DBDT))*oMODB*oMODB

 DMODBDS=dot(B,B(1)*DBDX+B(2)*DBDY+B(3)*DBDZ)*oMODB*oMODB

 GRADDRIFT = CROSS(B,GRADB)*oMODB*oMODB*ofacsq

 DlittleBDS = (B(1)*DlittleBDX + B(2)*DlittleBDY + B(3)*DlittleBDZ)*oMODB
 
 VEGRADlittleB = (VE(1)*DlittleBDX + VE(2)*DlittleBDY + VE(3)*DlittleBDZ)

 DVEDS = (B(1)*DVEDX + B(2)*DVEDY + B(3)*DVEDZ)*oMODB*oMODB

 VEGRADVE = VE(1)*DVEDX + VE(2)*DVEDY + VE(3)*DVEDZ
 
!All the acceleration drift terms, to be crossed with B
! ACCDRIFT = VPAR*DlittleBDT + VPAR*VPAR*DlittleBDS + VPAR*VEGRADlittleB + DVEDT + VPAR*DVEDS + VEGRADVE


!all the terms that make up the last bit of the parallel equation
 OTHERS = DlittleBDT + VPAR*DlittleBDS + VEGRADlittleB
 !OTHERS = DlittleBDT + U*DlittleBDS + VEGRADlittleB

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
 
  DUDT = Omscl*tscl*EPAR - MU/gamma*DMODBDS*fac + gamma*DOT(VE,OTHERS)						! THRELFALLETAL 2015

  DmodBDT=GRADBT-dot(B,DlittleBDT)
  RELDRIFT1=U*DlittleBDT+gamma*DVEDT+vscl*vscl/c/c*MU/gamma*VE*DmodBDT*fac
  !RELDRIFT1=MU/gamma*GRADDRIFT+U*DlittleBDT+gamma*DVEDT+vscl*vscl/c/c*mu/gamma*VE*DmodBDT				!(assuming eperp<<B)
  RELDRIFT2=U/gamma*EPAR*VE
 
  d1=MU/gamma*GRADDRIFT
  d2=MU/gamma*vscl*vscl/c/c*CROSS(B,VE)*omodB*omodB*ofacsq*DmodBDT
  d3=U*CROSS(B,OTHERS)*omodB*omodB*ofacsq
  d4=gamma*CROSS(B,DVEDT)*omodB*omodB*ofacsq
  d5=Epar*U/gamma*vscl*vscl/c/c*CROSS(B,VE)*omodB*omodB*ofacsq
  
  
   DRperpDT= VE + oneuponOmscl/tscl*(d1+d2+d3+d4)+d5									! THRELFALL ET AL 2015
   !print *, 'derivs ve = ', ve*Vscl, ' drperpdt = ', DRperpDT*Vscl

 
  !DRperpDT=VE + oneuponOmscl/tscl*(CROSS(B,RELDRIFT1)*oMODB*oMODB)+vscl*vscl/c/c*(CROSS(B,RELDRIFT2)*oMODB*oMODB)	!(assuming eperp<<B)
   DRDT =  DRperpDT + (U/gamma)*B*oMODB
   !print *, 'derivs DRDT = ', DRDT*Vscl, ' drperpdt = ', DRperpDT*Vscl
  
  DgammaDT=Omscl*tscl*vscl*vscl/c/c*(dot(DRperpDT,E)+U/gamma*dot(B,E)*oMODB)+vscl*vscl/c/c*MU/gamma*DmodBDT*fac	!THRELFALL ET AL 2015
  

!-----------------------------------------------------------------------!
! Relativistic versions WITH Birn et al (2004) Normalisation		!
!-----------------------------------------------------------------------!
!  UE=gamma*VE
!
!  !DUDT = sigma*DOT(UE,DlittleBDT)-oneosigma/epsilon*EPAR - oneosigma*MU/gamma*DMODBDS*fac			! full eqn
!  DUDT = -oneosigma/epsilon*EPAR - oneosigma*MU/gamma*DMODBDS*fac						 ! ignore sigmas
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
!  !DRperpDT= VE-epsilon*(d1+d2+d3+d4)+d5								! full equation
!  DRperpDT= VE-epsilon*(d1+d2)+d5									! ignore sigma, sigma^2
!  
!  DRDT =  DRperpDT + (U/gamma)*B*oMODB
!  
!  DgammaDT=vscl*vscl/c/c*oneosigma*oneosigma*(-oneoepsilon*(dot(DRperpDT,E)+oneosigma*U/gamma*dot(B,E)*oMODB)+MU/gamma*DmodBDT*fac)
 
END SUBROUTINE DERIVS
!-------------------------------------------------------------------------------------------!   
subroutine gammacalc(gamma,U,UE,mu,modB)

implicit none
real(num), intent(in)                   :: U,mu,modB
real(num), dimension(3)                 :: UE,VE
real(num)                               :: gamma

VE = UE/gamma

gamma = sqrt(1.0_num + vscl*vscl*U*U/(c*c) + 2.0_num*m*vscl*vscl*mu*modB/(m*c*c))/sqrt(1.0_num - vscl*vscl*dot(VE,VE)/(c*c))

end subroutine gammacalc
!-------------------------------------------------------------------------------------------!   
END MODULE M_derivsR
