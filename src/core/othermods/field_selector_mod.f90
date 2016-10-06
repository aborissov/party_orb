MODULE M_fields
  
  USE global
  USE sep_fields, ONLY: SEPFIELDS
  USE lare_fields, ONLY: LAREFIELDS
  USE CMT_fields, ONLY: CMTFIELDS
  USE test_fields, ONLY: TESTFIELDS
  USE lare_functions, ONLY: str_cmp
  USE bourdin_fields, ONLY: BOURDINFIELDS
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: FIELDS

  CONTAINS 
!---------------------------------------------------------------------------!
SUBROUTINE FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,rho,temperature,eta)
!selects which fields module we actually use.

 INTEGER, SAVE :: route = -1
 REAL(num)                              :: eta
 REAL(num), INTENT(OUT)                 :: rho,temperature
 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf
 REAL(num), INTENT(IN) 			:: T, T1, T2

100 SELECT CASE(route)
  CASE(-1) ! setup on first use of subroutine
    IF ((str_cmp(FMOD, "L3D")).OR.(str_cmp(FMOD, "l3d"))) THEN
      route=1
    ELSE IF ((str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN
      route=2  
    ELSE IF ((str_cmp(FMOD, "SEP")).OR.(str_cmp(FMOD, "sep"))) THEN
      route=3   
    ELSE IF ((str_cmp(FMOD, "CMT")).OR.(str_cmp(FMOD, "cmt"))) THEN
      route=4
    ELSE IF ((str_cmp(FMOD, "TEST")).OR.(str_cmp(FMOD, "test"))) THEN
      route=5
    ELSE IF ((str_cmp(FMOD, "BOR")).OR.(str_cmp(FMOD, "bor"))) THEN
      route=6
    ELSE
      PRINT*, "incorrect module selection, choose from:"
      PRINT*, "['l3d','l2d','sep','CMT','test','bour']"
      STOP
    END IF
    GO TO 100 ! now actually head back and select case we want!
  CASE(1)
    !CALL LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
    CALL LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,rho,temperature,eta)
  CASE(2)
    !CALL LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
    CALL LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,rho,temperature,eta)
  CASE(3)
    CALL SEPFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(4)
    CALL CMTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(5)
    CALL TESTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(6)
    CALL BOURDINFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
END SELECT


END SUBROUTINE FIELDS

END MODULE M_fields

