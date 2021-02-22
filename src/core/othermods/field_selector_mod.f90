MODULE M_fields
  
  USE global
  USE sep_fields, ONLY: SEPFIELDS
  USE lare_fields, ONLY: LAREFIELDS
  USE CMT_fields, ONLY: CMTFIELDS
  USE test_fields, ONLY: TESTFIELDS
  USE bourdin_fields, ONLY: BOURDINFIELDS  
  USE NLFF_fields, ONLY: NLFFFIELDS
  USE MHDp_fields, ONLY: MHDpFIELDS
  USE FR_fields, ONLY: FRFIELDS
  USE bifrost_fields_mod, ONLY: bifrost_fields, bifrost_grid
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: FIELDS

  CONTAINS 
!---------------------------------------------------------------------------!
SUBROUTINE FIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,T1,T2,BF_grid)
!selects which fields module we actually use.

 INTEGER, SAVE :: route = -1
 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf
 REAL(num), INTENT(IN) 			:: T, T1, T2 ! ALEXEI: check if we actually need T1 and T2 here...

 ! ALEXEI: somehow make this work only when using bifrost fields
 type(bifrost_grid), intent(in) :: BF_grid

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
    ELSE IF ((str_cmp(FMOD, "FRE")).OR.(str_cmp(FMOD, "fre"))) THEN
      route=7
    ELSE IF ((str_cmp(FMOD, "NLFF")).OR.(str_cmp(FMOD, "nlff"))) THEN
      route=8    
    ELSE IF ((str_cmp(FMOD, "MHDp")).OR.(str_cmp(FMOD, "mhdp"))) THEN
      route=9
    ELSE IF ((str_cmp(FMOD, "BF"))) THEN
      route=10
    ELSE
      PRINT*, "FIELD SELECTOR MODULE: incorrect field choice, choose from:"
      PRINT*, "['l3d','l2d','sep','CMT','test','bour','FRE', 'NLFF', 'MHDp', 'BF']"
      STOP
    END IF
    GO TO 100 ! now actually head back and select case we want!
  CASE(1)
    CALL LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(2)
    CALL LAREFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(3)
    CALL SEPFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(4)
    CALL CMTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(5)
    CALL TESTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(6)
    CALL BOURDINFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(7)
    CALL FRFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(8)
    CALL NLFFFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(9)
    CALL MHDpFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
  CASE(10)
    CALL bifrost_fields(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,BF_grid)
END SELECT


END SUBROUTINE FIELDS

END MODULE M_fields

