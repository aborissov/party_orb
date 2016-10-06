MODULE M_derivsN

 USE GLOBAL
 USE M_fields
 USE M_products
  
 IMPLICIT NONE
 
  PRIVATE
  PUBLIC :: DERIVS
  
  CONTAINS
!------------------------------------------------------------------------------!   
SUBROUTINE DERIVS(T, R, DRDT, VPAR, DVPARDT,MU, T1, T2, UGB, UC)
! Calculates individual bits of equations (44) and (45) of Guiliani et al (2005)
!
 IMPLICIT NONE
 REAL(num), INTENT(IN)			:: T,MU, T1,T2, VPAR
 REAL(num), INTENT(OUT) 		:: DVPARDT
 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), DIMENSION(3), INTENT(OUT)	:: DRDT
 REAL(num), DIMENSION(3), INTENT(OUT)	:: UGB, UC
 
 REAL(num), DIMENSION(3)		:: B,E,Vf
 REAL(num), DIMENSION(3)		:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3)		:: GRADB,DlittleBDT,DlittleBDX,DlittleBDY,DlittleBDZ
 REAL(num), DIMENSION(3)		:: UE, DUEDX,DUEDY,DUEDZ,DUEDT
 REAL(num), DIMENSION(3)		:: GRADDRIFT,DlittleBDS,UEGRADlittleB,UEGRADUE,DUEDS
 REAL(num), DIMENSION(3)		:: ACCDRIFT, OTHERS
 !REAL(num), DIMENSION(3)		:: SCRAE
 REAL(num)				:: MODB, oMODB, DMODBDS, EPAR,GRADBT, rho, temperature,eta
 

 !CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2)
 CALL FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
 MODB = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))		! |B|
 
 oMODB=1.0_num/MODB					! 1/|B| - computationally less expensive to multiply..

 EPAR = DOT(E,B)*oMODB 

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

 GRADDRIFT = CROSS(B,GRADB)*oMODB*oMODB

 DlittleBDS = (B(1)*DlittleBDX + B(2)*DlittleBDY + B(3)*DlittleBDZ)*oMODB
 
 UEGRADlittleB = (UE(1)*DlittleBDX + UE(2)*DlittleBDY + UE(3)*DlittleBDZ)

 DUEDS = (B(1)*DUEDX + B(2)*DUEDY + B(3)*DUEDZ)*oMODB*oMODB

 UEGRADUE = UE(1)*DUEDX + UE(2)*DUEDY + UE(3)*DUEDZ
 
!All the acceleration drift terms, to be crossed with B
 ACCDRIFT = VPAR*DlittleBDT + VPAR*VPAR*DlittleBDS + VPAR*UEGRADlittleB + DUEDT + VPAR*DUEDS + UEGRADUE


!all the terms that make up the last bit of the parallel equation
 OTHERS = DlittleBDT + VPAR*DlittleBDS + UEGRADlittleB

!-------------------------------------------------------------------------------!
! The final equations of motion (eqs. 44 and 45 from Guiliani et al. (2005)	!
!-------------------------------------------------------------------------------!
!problem seems to be with parallel electric field - QoM is BIG.

  DVPARDT = Omscl*tscl*EPAR - MU*DMODBDS + DOT(UE,OTHERS)		! My normalisation
! print *, "Omega_i=",Omega_i*tscl
 !*tscl*EPAR
! print *, MU*DMODBDS
 !print *, DOT(UE,OTHERS)
 !STOP
 !DVPARDT = QoM*Lscl*Bscl*Bscl/Escl*EPAR - MU*DMODBDS + DOT(UE,OTHERS)	! Original normalisation
 !DVPARDT = QoM*EPAR - MU*DMODBDS + DOT(UE,OTHERS)			! No normalisation

  !DRDT = UE + MoQ*Escl/Lscl/Bscl/Bscl*(MU*GRADDRIFT + CROSS(B,ACCDRIFT)*oMODB*oMODB) + VPAR*B*oMODB	! original normalisation
  DRDT = UE + oneuponOmscl/tscl*(MU*GRADDRIFT + CROSS(B,ACCDRIFT)*oMODB*oMODB) + VPAR*B*oMODB		! extra term due to Rperp->R

 !print*, "dRperp/dt=", vscl*(DRDT-VPAR*B*oMODB)*1e6, "micrometers/sec"
  !print*, "mu:",mu

  UGB=oneuponOmscl/tscl*MU*GRADDRIFT
  UC =oneuponOmscl/tscl*VPAR*VPAR*cross(B,DlittleBDS)*oMODB

! What are these temp values?

!  temp(1)=MU*dot(B,DBDT)*oMODB
!  temp(2)=MU*VPAR*DMODBDS
!  temp(3)=MU*dot(UE,GRADB)

!  temp(1)=dot(B,DBDT)*oMODB
!  temp(2)=VPAR*DMODBDS
!  temp(3)=dot(UE,GRADB)

!  temp2(1)=-MU*DMODBDS
!  temp2(2)=dot(UE,DlittleBDT)
!  temp2(3)=dot(UE,VPAR*DlittleBDS)
!  temp2(4)=dot(UE,UEGRADlittleB)
!  sum_temp2=temp2(1)+temp2(2)+temp2(3)+temp2(4)

! drift(1,:) = UE(:)
! drift(2,:) = MoQ*Escl/Lscl/Bscl/Bscl*MU*GRADDRIFT(:)
! drift(3,:) = MoQ*Escl/Lscl/Bscl/Bscl*CROSS(B,VPAR*DlittleBDT)*oMODB*oMODB
! drift(4,:) = MoQ*Escl/Lscl/Bscl/Bscl*CROSS(B,VPAR*VPAR*DlittleBDS)*oMODB*oMODB
! drift(5,:) = MoQ*Escl/Lscl/Bscl/Bscl*CROSS(B,VPAR*UEGRADlittleB)*oMODB*oMODB
! drift(6,:) = MoQ*Escl/Lscl/Bscl/Bscl*CROSS(B,DUEDT)*oMODB*oMODB
! drift(7,:) = MoQ*Escl/Lscl/Bscl/Bscl*CROSS(B,VPAR*DUEDS)*oMODB*oMODB
! drift(8,:) = MoQ*Escl/Lscl/Bscl/Bscl*CROSS(B,UEGRADUE)*oMODB*oMODB

END SUBROUTINE DERIVS
!----------------------------------------------------------------------------! 
END MODULE M_derivsN
