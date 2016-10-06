MODULE lare_fields
! lare interface which is dynamically called, depending on whether you are interested in 2d or 3d lare data
  
  USE GLOBAL
  USE lare_functions, ONLY: str_cmp, T2d, T3d

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: LAREFIELDS

  CONTAINS 

!----------------------------------------------------------------------;   
!SUBROUTINE LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
SUBROUTINE LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,rho,temperature,eta)
! wrapper module which calculates variables at sub grid particle locations, 
!  and outputs values in specific order

 REAL(num), INTENT(OUT)                        :: rho,temperature,eta
 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf, j
 REAL(num), INTENT(IN) 			:: T
 REAL(num), DIMENSION(39)		:: iquants
 !INTEGER				:: ij, jjx, jjy, jjz
 !LOGICAL				:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE.
 !LOGICAL, INTENT(OUT)			:: ERR=.FALSE.

  IF (((R(1).lt.myx(4)).and.(R(1).ge.myx(nx-3))).OR.((R(2).lt.myy(4)).and.(R(2).ge.myy(ny-3)))) THEN
   print*, 'returning as out of bounds'	! this should not get triggered but hey.
   RETURN
  ENDIF
  
  IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d"))) THEN
   IF ((R(3).lt.myz(4)).and.(R(3).ge.myz(nz-3))) THEN
    PRINT*, 'outside z bounds, returning'
    RETURN
   ENDIF
   !iquants=T3d(R,T)
   print *,'WARNING: need to change T3d to be compatible with more outputs. Stopped in lare_fields.f90'
   stop
  ELSE
   iquants=T2d(R,T)
  ENDIF
  ! compare array against itself for NANs         
  IF (ALL(iquants.eq.iquants)) THEN
  !everything is fine
  ELSE
   PRINT*, 'NAN DETECTED IN LARE OUTPUT'
   print *, 'lare_fields = ', iquants
   !ERR=.TRUE.
   STOP
  ENDIF
   
   B	=iquants(1:3)
  ! print*, iquants(3)
   
   Vf	=iquants(4:6)
   E	=iquants(7:9)
   j	=iquants(10:12)
   DBDX	=iquants(13:15)
   DBDY	=iquants(16:18)
   DBDZ	=iquants(19:21)
   DEDX	=iquants(22:24)
   DEDY	=iquants(25:27)
   DEDZ	=iquants(28:30)
   DBDT	=iquants(31:33)
   DEDT	=iquants(34:36)  
   rho = iquants(37)
   temperature = iquants(38)
   eta = iquants(39)

   !print *, 'lare_fields epar 0= ', (E(1)*B(1) + E(2)*B(2) + E(3)*B(3))/sqrt(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
   !print *, 'lare fields E = ', E, ' B = ', B
   !if (eta .ne. 0) print *, 'lare fields eta = ', eta, ' R = ', R, ' T = ', T

   

END SUBROUTINE LAREFIELDS

END MODULE lare_fields

