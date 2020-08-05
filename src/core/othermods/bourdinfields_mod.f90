MODULE bourdin_fields
  
  USE global
  USE M_products
  USE lare_functions

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bour_ini, bour_fini, BOURDINFIELDS

  CONTAINS 
!------------------------------------------------------------------------------!
SUBROUTINE bour_ini
 
 INTEGER				:: i, j, k, rl, c, g, t1x,t1y,t1z
 INTEGER(kind=8) 			:: reclen
 REAL(num)				:: tmp1, tmp2, tmp3
 CHARACTER(len=6)			:: nom='VAR378'
 CHARACTER(len=6)			:: nomBabs='_B_ABS',nomEabs='_E_ABS',nomUabs='_U_ABS'
 CHARACTER(len=2)			:: nomB='_B',nomE='_E',nomU='_U' 
 CHARACTER(len=6)			:: nomDBDX='_DB_DX',nomDBDY='_DB_DY',nomDBDZ='_DB_DZ'
 !REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: data
 !CHARACTER(len=7) 			:: fmt1='(I4.4)', istring 	! format descriptor
 CHARACTER(LEN = 20+data_dir_max_length):: locB, locE, locU, locDBDX, locDBDY, locDBDZ 
 LOGICAL				:: ONLY=.FALSE.

 locB=trim(adjustl(sloc))//trim(nom)//trim(nomB)//filetypeb
 locE=trim(adjustl(sloc))//trim(nom)//trim(nomE)//filetypeb
 locU=trim(adjustl(sloc))//trim(nom)//trim(nomU)//filetypeb
 locDBDX=trim(adjustl(sloc))//trim(nom)//trim(nomDBDX)//filetypeb
 locDBDY=trim(adjustl(sloc))//trim(nom)//trim(nomDBDY)//filetypeb
 locDBDZ=trim(adjustl(sloc))//trim(nom)//trim(nomDBDZ)//filetypeb
 
 nx=1024
 ny=1024
 nz=256
 rl=nx*ny*nz*3			! typical size of B arrays: 3 = 3 components of field
 !reclen=long(8)*long(nx)*long(ny)	
 ! ALEXEI: apparently ifort doesn't like long()
 reclen=8*nx*ny
 
 ALLOCATE(vx(1:nx,1:ny,1:nz,0))
 ALLOCATE(vy(1:nx,1:ny,1:nz,0))
 ALLOCATE(vz(1:nx,1:ny,1:nz,0))
 ALLOCATE(bx(1:nx,1:ny,1:nz,0))
 ALLOCATE(by(1:nx,1:ny,1:nz,0))
 ALLOCATE(bz(1:nx,1:ny,1:nz,0))
 ALLOCATE(dbxdx(1:nx,1:ny,1:nz))
 ALLOCATE(dbydx(1:nx,1:ny,1:nz))
 ALLOCATE(dbzdx(1:nx,1:ny,1:nz))
 ALLOCATE(dbxdy(1:nx,1:ny,1:nz))
 ALLOCATE(dbydy(1:nx,1:ny,1:nz))
 ALLOCATE(dbzdy(1:nx,1:ny,1:nz))
 ALLOCATE(dbxdz(1:nx,1:ny,1:nz))
 ALLOCATE(dbydz(1:nx,1:ny,1:nz))
 ALLOCATE(dbzdz(1:nx,1:ny,1:nz))
 ALLOCATE(Ex(1:nx,1:ny,1:nz))
 ALLOCATE(Ey(1:nx,1:ny,1:nz))
 ALLOCATE(Ez(1:nx,1:ny,1:nz))
 ALLOCATE(myx(1:nx))
 ALLOCATE(myy(1:ny))
 ALLOCATE(myz(1:nz))

! read grid:
 PRINT*, 'reading grid'
 OPEN(66,file=locB, form="unformatted", recl=8,ACCESS='DIRECT')
 READ(66, rec=rl+1) tmp1
 DO i=1,nx
  READ(66, rec=rl+1+i) myx(i)
 ENDDO
 DO i=1,ny
  READ(66, rec=rl+1+nx+i) myy(i)
 ENDDO
 DO i=1,nz
  READ(66, rec=rl+1+nx+ny+i) myz(i)
 ENDDO
 CLOSE(66)
 !NORMALISE NOW:
 myx=myx/lscl
 myy=myy/lscl
 myz=myz/lscl
 WRITE(*,'(a,f7.2,a,f7.2,a,e9.2,a)') "xr=[", minval(myx),"->", maxval(myx),"]*", lscl, "m"
 WRITE(*,'(a,f7.2,a,f7.2,a,e9.2,a)') "yr=[", minval(myy),"->", maxval(myy),"]*", lscl, "m"
 WRITE(*,'(a,f7.2,a,f7.2,a,e9.2,a)') "zr=[", minval(myz),"->", maxval(myz),"]*", lscl, "m"

! read fields
 OPEN(66,file=locB, form="unformatted", recl=reclen,ACCESS='DIRECT')
 OPEN(67,file=locE, form="unformatted", recl=reclen,ACCESS='DIRECT')
 OPEN(68,file=locU, form="unformatted", recl=reclen,ACCESS='DIRECT')
 OPEN(76,file=locDBDX, form="unformatted", recl=reclen,ACCESS='DIRECT')
 OPEN(76,file=locDBDX, form="unformatted", recl=reclen,ACCESS='DIRECT')
 OPEN(77,file=locDBDY, form="unformatted", recl=reclen,ACCESS='DIRECT')
 OPEN(78,file=locDBDZ, form="unformatted", recl=reclen,ACCESS='DIRECT') 
 write (*,"(a)",advance="no") 'reading Bourdin field variables: '
 DO j=1,nz
   READ(66, rec=j) 	bx(1:nx,1:ny,j,0)
   READ(66, rec=nz+j) 	by(1:nx,1:ny,j,0)
   READ(66, rec=2*nz+j) bz(1:nx,1:ny,j,0)
   READ(67, rec=j) 	Ex(1:nx,1:ny,j)
   READ(67, rec=nz+j) 	Ey(1:nx,1:ny,j)
   READ(67, rec=2*nz+j) Ez(1:nx,1:ny,j)
   READ(68, rec=j) 	Vx(1:nx,1:ny,j,0)
   READ(68, rec=nz+j) 	Vy(1:nx,1:ny,j,0)
   READ(68, rec=2*nz+j) Vz(1:nx,1:ny,j,0)
   READ(76, rec=j) 	dbxdx(1:nx,1:ny,j)
   READ(76, rec=nz+j) 	dbydx(1:nx,1:ny,j)
   READ(76, rec=2*nz+j) dbzdx(1:nx,1:ny,j)
   READ(77, rec=j) 	dbxdy(1:nx,1:ny,j)
   READ(77, rec=nz+j) 	dbydy(1:nx,1:ny,j)
   READ(77, rec=2*nz+j) dbzdy(1:nx,1:ny,j)
   READ(78, rec=j) 	dbxdz(1:nx,1:ny,j)
   READ(78, rec=nz+j) 	dbydz(1:nx,1:ny,j)
   READ(78, rec=2*nz+j) dbzdz(1:nx,1:ny,j)
   IF (j.eq.nz/4) write (*,"(a)",advance="no") '.'
   IF (j.eq.nz/2) write (*,"(a)",advance="no") '.'
   IF (j.eq.3*nz/4) write (*,"(a)",advance="no") '.'
  ENDDO
  write (*,"(a)") 'done!' 
  CLOSE(66)
  CLOSE(67)
  CLOSE(68)
  CLOSE(76)
  CLOSE(77)
  CLOSE(78)
  
  !currently do not have E field derivs
 ! can easily include this, but for now calculate E-field derivs on the fly like for Lare3d..
 
END SUBROUTINE bour_ini
!------------------------------------------------------------------------------!
SUBROUTINE bour_fini

    DEALLOCATE(Ex)
    DEALLOCATE(Ey)
    DEALLOCATE(Ez)
    DEALLOCATE(vx)
    DEALLOCATE(vy)
    DEALLOCATE(vz)
    DEALLOCATE(bx)
    DEALLOCATE(by)
    DEALLOCATE(bz)
    DEALLOCATE(dbxdx)
    DEALLOCATE(dbydx)
    DEALLOCATE(dbzdx)
    DEALLOCATE(dbxdy)
    DEALLOCATE(dbydy)
    DEALLOCATE(dbzdy)
    DEALLOCATE(dbxdz)
    DEALLOCATE(dbydz)
    DEALLOCATE(dbzdz)
    DEALLOCATE(myx)
    DEALLOCATE(myy)
    DEALLOCATE(myz)

END SUBROUTINE bour_fini
!------------------------------------------------------------------------------!
FUNCTION b3d(R)
! function to calculate local quantities at the particle position
! (slimmed down version of f3d with only local interpolation!).
! OCT2014: updated to include higher order derivs to increase accuracy

   REAL(num), DIMENSION(3), INTENT(IN)	:: R		!actual position
   REAL(num), DIMENSION(27)		:: b3d
   REAL(num), DIMENSION(3)		:: dg, odg, coffset
   REAL(num)				:: bxt, byt, bzt
   REAL(num)				:: temp
   REAL(num)				:: vxt, vyt, vzt, oneotwelve=1.0_num/12.0_num
   REAL(num)				:: Ext, Eyt, Ezt
   REAL(num)				:: dbxdxt,dbxdyt,dbxdzt,dbydxt,dbydyt,dbydzt,dbzdxt,dbzdyt,dbzdzt
   REAL(num)				:: dExdxt,dExdyt,dExdzt,dEydxt,dEydyt,dEydzt,dEzdxt,dEzdyt,dEzdzt
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: dmExdx,dmExdy,dmExdz,dmEydx,dmEydy,dmEydz,dmEzdx,dmEzdy,dmEzdz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: mEx, mEy, mEz
   INTEGER, DIMENSION(3) 		:: l		! set at value it could never reach!
   INTEGER				:: jjx, jjy, jjz
   LOGICAL					:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE., tempflag
   !INTEGER, SAVE :: out=1
    temp=0.0_num

    l=(/-nx,-ny,-nz/)					! initial value of l, set to silly value (as >=0 triggers flag)   
    tempflag=.FALSE.  
  ! --STEP ONE-- !
  ! first, locate the (x,y,z) position on the grid
  ! particle found between (jjx,jjy,jjz) and (jjx+1,jjy+1,jjz+1)   
     
   DO jjx=5,nx-4
    IF ((R(1).ge.myx(jjx)).and.(R(1).lt.myx(jjx+1))) THEN 
      l(1)=jjx
      EXIT
    ENDIF
   ENDDO
   IF (l(1).eq.(-nx)) fxflag=.TRUE.
   DO jjy=5,ny-4
    IF ((R(2).ge.myy(jjy)).and.(R(2).lt.myy(jjy+1))) THEN
      l(2)=jjy
      EXIT
    ENDIF
   ENDDO
   IF (l(2).eq.(-ny)) fyflag=.TRUE.
   DO jjz=5,nz-4
    IF ((R(3).ge.myz(jjz)).and.(R(3).lt.myz(jjz+1))) THEN
      l(3)=jjz
      EXIT
    ENDIF
   ENDDO
   IF (l(3).eq.(-nz)) fzflag=.TRUE.
 

 ! I don't think this is an actual error since going beyond the grid will be detected and these results ignored
  IF ((fxflag).OR.(fyflag).OR.(fzflag)) THEN
    !print*, 'triggered'		
    b3d(:)=0.0_num
   RETURN
  ENDIF
  !call this cell "l"
  ! need some sort of condition here to guarantee that jj is far enough away from boundaries
  
  ! grid spacing may be VARIABLE. CHECK: use local grid spacing!!
    dg=(/abs(myx(l(1)+1)-myx(l(1))),abs(myy(l(2)+1)-myy(l(2))),abs(myz(l(3)+1)-myz(l(3)))/)	! grid spacing
    odg=(/1.0_num/dg(1),1.0_num/dg(2),1.0_num/dg(3)/)	! one over grid spacing
   coffset=(/(R(1)-myx(l(1)))*odg(1),(R(2)-myy(l(2)))*odg(2),(R(3)-myz(l(3)))*odg(3)/)
   
   !PRINT*, COFFSET
   
   ! --STEP TWO-- !
   ! begin interpolation of basic variables we already have, e.g. bx, by, bz, vx, vy, vz
       
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!1
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&	!1
   bx(l(1),l(2),l(3),0),bx(l(1)+1,l(2),l(3),0),bx(l(1),l(2)+1,l(3),0),&
   bx(l(1)+1,l(2)+1,l(3),0),bx(l(1),l(2),l(3)+1,0),bx(l(1)+1,l(2),l(3)+1,0), &
   bx(l(1),l(2)+1,l(3)+1,0),bx(l(1)+1,l(2)+1,l(3)+1,0))
   bxt=temp/bscl

   !print*, 'Bx:', temp
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!2
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   by(l(1),l(2),l(3),0),by(l(1)+1,l(2),l(3),0),by(l(1),l(2)+1,l(3),0),&
   by(l(1)+1,l(2)+1,l(3),0),by(l(1),l(2),l(3)+1,0),by(l(1)+1,l(2),l(3)+1,0), &
   by(l(1),l(2)+1,l(3)+1,0),by(l(1)+1,l(2)+1,l(3)+1,0))
   byt=temp/bscl

   !print*, 'By:', temp

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!3
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   bz(l(1),l(2),l(3),0),bz(l(1)+1,l(2),l(3),0),bz(l(1),l(2)+1,l(3),0),&
   bz(l(1)+1,l(2)+1,l(3),0),bz(l(1),l(2),l(3)+1,0),bz(l(1)+1,l(2),l(3)+1,0), &
   bz(l(1),l(2)+1,l(3)+1,0),bz(l(1)+1,l(2)+1,l(3)+1,0))
   bzt=temp/bscl
   
   !print*, 'Bz:', temp
   !STOP  
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!4
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   vx(l(1),l(2),l(3),0),vx(l(1)+1,l(2),l(3),0),vx(l(1),l(2)+1,l(3),0),&
   vx(l(1)+1,l(2)+1,l(3),0),vx(l(1),l(2),l(3)+1,0),vx(l(1)+1,l(2),l(3)+1,0), &
   vx(l(1),l(2)+1,l(3)+1,0),vx(l(1)+1,l(2)+1,l(3)+1,0))
   vxt=temp/vscl

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!5
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   vy(l(1),l(2),l(3),0),vy(l(1)+1,l(2),l(3),0),vy(l(1),l(2)+1,l(3),0),&
   vy(l(1)+1,l(2)+1,l(3),0),vy(l(1),l(2),l(3)+1,0),vy(l(1)+1,l(2),l(3)+1,0), &
   vy(l(1),l(2)+1,l(3)+1,0),vy(l(1)+1,l(2)+1,l(3)+1,0))
   vyt=temp/vscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!6
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   vz(l(1),l(2),l(3),0),vz(l(1)+1,l(2),l(3),0),vz(l(1),l(2)+1,l(3),0),&
   vz(l(1)+1,l(2)+1,l(3),0),vz(l(1),l(2),l(3)+1,0),vz(l(1)+1,l(2),l(3)+1,0), &
   vz(l(1),l(2)+1,l(3)+1,0),vz(l(1)+1,l(2)+1,l(3)+1,0))
   vzt=temp/vscl
     
   
  !IF (tempflag) THEN
   ALLOCATE(mEx(-4:5,-4:5,-4:5))
   ALLOCATE(mEy(-4:5,-4:5,-4:5))
   ALLOCATE(mEz(-4:5,-4:5,-4:5))

   mEx=Ex(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5)
   mEy=Ey(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5)
   mEz=Ez(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5)
   
   ALLOCATE(dmExdx(0:1,-2:3,-2:3))
   ALLOCATE(dmEydx(0:1,-2:3,-2:3))
   ALLOCATE(dmEzdx(0:1,-2:3,-2:3))
   DO iz=-2,3   
    DO iy=-2,3
     DO ix=0,1
      dmExdx(ix,iy,iz)=oneotwelve*(-mEx(ix+2,iy,iz)+8.0_num*mEx(ix+1,iy,iz)-8.0_num*mEx(ix-1,iy,iz)+mEx(ix-2,iy,iz))*odg(1)
      dmEydx(ix,iy,iz)=oneotwelve*(-mEy(ix+2,iy,iz)+8.0_num*mEy(ix+1,iy,iz)-8.0_num*mEy(ix-1,iy,iz)+mEy(ix-2,iy,iz))*odg(1)
      dmEzdx(ix,iy,iz)=oneotwelve*(-mEz(ix+2,iy,iz)+8.0_num*mEz(ix+1,iy,iz)-8.0_num*mEz(ix-1,iy,iz)+mEz(ix-2,iy,iz))*odg(1)
     ENDDO
    ENDDO
   ENDDO
   
   ALLOCATE(dmExdy(-2:3,0:1,-2:3))
   ALLOCATE(dmEydy(-2:3,0:1,-2:3))
   ALLOCATE(dmEzdy(-2:3,0:1,-2:3))  
   DO iz=-2,3   
    DO iy=0,1
     DO ix=-2,3  
      dmExdy(ix,iy,iz)=oneotwelve*(-mEx(ix,iy+2,iz)+8.0_num*mEx(ix,iy+1,iz)-8.0_num*mEx(ix,iy-1,iz)+mEx(ix,iy-2,iz))*odg(2)
      dmEydy(ix,iy,iz)=oneotwelve*(-mEy(ix,iy+2,iz)+8.0_num*mEy(ix,iy+1,iz)-8.0_num*mEy(ix,iy-1,iz)+mEy(ix,iy-2,iz))*odg(2)
      dmEzdy(ix,iy,iz)=oneotwelve*(-mEz(ix,iy+2,iz)+8.0_num*mEz(ix,iy+1,iz)-8.0_num*mEz(ix,iy-1,iz)+mEz(ix,iy-2,iz))*odg(2)
     ENDDO
    ENDDO
   ENDDO   
   
   ALLOCATE(dmExdz(-2:3,-2:3,0:1))
   ALLOCATE(dmEydz(-2:3,-2:3,0:1))
   ALLOCATE(dmEzdz(-2:3,-2:3,0:1))   
   DO iz=0,1   
    DO iy=-2,3
     DO ix=-2,3
      dmExdz(ix,iy,iz)=oneotwelve*(-mEx(ix,iy,iz+2)+8.0_num*mEx(ix,iy,iz+1)-8.0_num*mEx(ix,iy,iz-1)+mEx(ix,iy,iz-2))*odg(3)
      dmEydz(ix,iy,iz)=oneotwelve*(-mEy(ix,iy,iz+2)+8.0_num*mEy(ix,iy,iz+1)-8.0_num*mEy(ix,iy,iz-1)+mEy(ix,iy,iz-2))*odg(3)
      dmEzdz(ix,iy,iz)=oneotwelve*(-mEz(ix,iy,iz+2)+8.0_num*mEz(ix,iy,iz+1)-8.0_num*mEz(ix,iy,iz-1)+mEz(ix,iy,iz-2))*odg(3)
     ENDDO
    ENDDO
   ENDDO
     
   DEALLOCATE(mEx)
   DEALLOCATE(mEy)
   DEALLOCATE(mEz)
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!10
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmExdx(0,0,0), dmExdx(1,0,0), dmExdx(0,1,0), dmExdx(1,1,0),&
    dmExdx(0,0,1), dmExdx(1,0,1), dmExdx(0,1,1), dmExdx(1,1,1))
   dExdxt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!11
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmEydx(0,0,0), dmEydx(1,0,0), dmEydx(0,1,0), dmEydx(1,1,0),&
    dmEydx(0,0,1), dmEydx(1,0,1), dmEydx(0,1,1), dmEydx(1,1,1))
   dEydxt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!12
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmEzdx(0,0,0), dmEzdx(1,0,0), dmEzdx(0,1,0), dmEzdx(1,1,0),&
    dmEzdx(0,0,1), dmEzdx(1,0,1), dmEzdx(0,1,1), dmEzdx(1,1,1))
   dEzdxt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!13
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmExdy(0,0,0), dmExdy(1,0,0), dmExdy(0,1,0), dmExdy(1,1,0),&
    dmExdy(0,0,1), dmExdy(1,0,1), dmExdy(0,1,1), dmExdy(1,1,1))
   dExdyt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!14
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmEydy(0,0,0), dmEydy(1,0,0), dmEydy(0,1,0), dmEydy(1,1,0),&
    dmEydy(0,0,1), dmEydy(1,0,1), dmEydy(0,1,1), dmEydy(1,1,1))
   dEydyt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!15
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmEzdy(0,0,0), dmEzdy(1,0,0), dmEzdy(0,1,0), dmEzdy(1,1,0),&
    dmEzdy(0,0,1), dmEzdy(1,0,1), dmEzdy(0,1,1), dmEzdy(1,1,1))
   dEzdyt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!16
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmExdz(0,0,0), dmExdz(1,0,0), dmExdz(0,1,0), dmExdz(1,1,0),&
    dmExdz(0,0,1), dmExdz(1,0,1), dmExdz(0,1,1), dmExdz(1,1,1))
   dExdzt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!17
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmEydz(0,0,0), dmEydz(1,0,0), dmEydz(0,1,0), dmEydz(1,1,0),&
    dmEydz(0,0,1), dmEydz(1,0,1), dmEydz(0,1,1), dmEydz(1,1,1))
   dEydzt=temp/Escl
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!18
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
    dmEzdz(0,0,0), dmEzdz(1,0,0), dmEzdz(0,1,0), dmEzdz(1,1,0),&
    dmEzdz(0,0,1), dmEzdz(1,0,1), dmEzdz(0,1,1), dmEzdz(1,1,1))
   dEzdzt=temp/Escl
   
   DEALLOCATE(dmExdx)
   DEALLOCATE(dmEydx)
   DEALLOCATE(dmEzdx)
   DEALLOCATE(dmExdy)
   DEALLOCATE(dmEydy)
   DEALLOCATE(dmEzdy)
   DEALLOCATE(dmExdz)
   DEALLOCATE(dmEydz)
   DEALLOCATE(dmEzdz)
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!7
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   Ex(l(1),l(2),l(3)),Ex(l(1)+1,l(2),l(3)),Ex(l(1),l(2)+1,l(3)),&
   Ex(l(1)+1,l(2)+1,l(3)),Ex(l(1),l(2),l(3)+1),Ex(l(1)+1,l(2),l(3)+1), &
   Ex(l(1),l(2)+1,l(3)+1),Ex(l(1)+1,l(2)+1,l(3)+1))
   Ext=temp/Escl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!8
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   Ey(l(1),l(2),l(3)),Ey(l(1)+1,l(2),l(3)),Ey(l(1),l(2)+1,l(3)),&
   Ey(l(1)+1,l(2)+1,l(3)),Ey(l(1),l(2),l(3)+1),Ey(l(1)+1,l(2),l(3)+1), &
   Ey(l(1),l(2)+1,l(3)+1),Ey(l(1)+1,l(2)+1,l(3)+1))
   Eyt=temp/Escl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!9
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   Ez(l(1),l(2),l(3)),Ez(l(1)+1,l(2),l(3)),Ez(l(1),l(2)+1,l(3)),&
   Ez(l(1)+1,l(2)+1,l(3)),Ez(l(1),l(2),l(3)+1),Ez(l(1)+1,l(2),l(3)+1), &
   Ez(l(1),l(2)+1,l(3)+1),Ez(l(1)+1,l(2)+1,l(3)+1))
   Ezt=temp/Escl

  
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!10
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbxdx(l(1),l(2),l(3)),dbxdx(l(1)+1,l(2),l(3)),dbxdx(l(1),l(2)+1,l(3)),&
   dbxdx(l(1)+1,l(2)+1,l(3)),dbxdx(l(1),l(2),l(3)+1),dbxdx(l(1)+1,l(2),l(3)+1), &
   dbxdx(l(1),l(2)+1,l(3)+1),dbxdx(l(1)+1,l(2)+1,l(3)+1))
   dbxdxt=temp/bscl*lscl

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!11
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbydx(l(1),l(2),l(3)),dbydx(l(1)+1,l(2),l(3)),dbydx(l(1),l(2)+1,l(3)),&
   dbydx(l(1)+1,l(2)+1,l(3)),dbydx(l(1),l(2),l(3)+1),dbydx(l(1)+1,l(2),l(3)+1), &
   dbydx(l(1),l(2)+1,l(3)+1),dbydx(l(1)+1,l(2)+1,l(3)+1))
   dbydxt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!12
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbzdx(l(1),l(2),l(3)),dbzdx(l(1)+1,l(2),l(3)),dbzdx(l(1),l(2)+1,l(3)),&
   dbzdx(l(1)+1,l(2)+1,l(3)),dbzdx(l(1),l(2),l(3)+1),dbzdx(l(1)+1,l(2),l(3)+1), &
   dbzdx(l(1),l(2)+1,l(3)+1),dbzdx(l(1)+1,l(2)+1,l(3)+1))
   dbzdxt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!13
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbxdy(l(1),l(2),l(3)),dbxdy(l(1)+1,l(2),l(3)),dbxdy(l(1),l(2)+1,l(3)),&
   dbxdy(l(1)+1,l(2)+1,l(3)),dbxdy(l(1),l(2),l(3)+1),dbxdy(l(1)+1,l(2),l(3)+1), &
   dbxdy(l(1),l(2)+1,l(3)+1),dbxdy(l(1)+1,l(2)+1,l(3)+1))
   dbxdyt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!14
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbydy(l(1),l(2),l(3)),dbydy(l(1)+1,l(2),l(3)),dbydy(l(1),l(2)+1,l(3)),&
   dbydy(l(1)+1,l(2)+1,l(3)),dbydy(l(1),l(2),l(3)+1),dbydy(l(1)+1,l(2),l(3)+1), &
   dbydy(l(1),l(2)+1,l(3)+1),dbydy(l(1)+1,l(2)+1,l(3)+1))
   dbydyt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!15
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbzdy(l(1),l(2),l(3)),dbzdy(l(1)+1,l(2),l(3)),dbzdy(l(1),l(2)+1,l(3)),&
   dbzdy(l(1)+1,l(2)+1,l(3)),dbzdy(l(1),l(2),l(3)+1),dbzdy(l(1)+1,l(2),l(3)+1), &
   dbzdy(l(1),l(2)+1,l(3)+1),dbzdy(l(1)+1,l(2)+1,l(3)+1))
   dbzdyt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!16
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbxdz(l(1),l(2),l(3)),dbxdz(l(1)+1,l(2),l(3)),dbxdz(l(1),l(2)+1,l(3)),&
   dbxdz(l(1)+1,l(2)+1,l(3)),dbxdz(l(1),l(2),l(3)+1),dbxdz(l(1)+1,l(2),l(3)+1), &
   dbxdz(l(1),l(2)+1,l(3)+1),dbxdz(l(1)+1,l(2)+1,l(3)+1))
   dbxdzt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!17
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbydz(l(1),l(2),l(3)),dbydz(l(1)+1,l(2),l(3)),dbydz(l(1),l(2)+1,l(3)),&
   dbydz(l(1)+1,l(2)+1,l(3)),dbydz(l(1),l(2),l(3)+1),dbydz(l(1)+1,l(2),l(3)+1), &
   dbydz(l(1),l(2)+1,l(3)+1),dbydz(l(1)+1,l(2)+1,l(3)+1))
   dbydzt=temp/bscl*lscl
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!18
   temp=linterp3d(coffset(1),coffset(2),coffset(3),&
   dbzdz(l(1),l(2),l(3)),dbzdz(l(1)+1,l(2),l(3)),dbzdz(l(1),l(2)+1,l(3)),&
   dbzdz(l(1)+1,l(2)+1,l(3)),dbzdz(l(1),l(2),l(3)+1),dbzdz(l(1)+1,l(2),l(3)+1), &
   dbzdz(l(1),l(2)+1,l(3)+1),dbzdz(l(1)+1,l(2)+1,l(3)+1))
   dbzdzt=temp/bscl*lscl        

   b3d=(/bxt,byt,bzt,vxt,vyt,vzt,Ext,Eyt,Ezt, &
         dbxdxt,dbydxt,dbzdxt,dbxdyt,dbydyt,dbzdyt,dbxdzt,dbydzt,dbzdzt, &
	 dExdxt,dEydxt,dEzdxt,dExdyt,dEydyt,dEzdyt,dExdzt,dEydzt,dEzdzt/)
   RETURN


END FUNCTION b3d
!------------------------------------------------------------------------------!    
SUBROUTINE BOURDINFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
! Bourdin fields
 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf
 REAL(num), INTENT(IN) 			:: T
 REAL(num), DIMENSION(27)		:: iquants

! j has been removed from here since we don't ever use it!

   iquants=b3d(R)

   B	=iquants(1:3)
   Vf	=iquants(4:6)
   E	=iquants(7:9)
   !E	=0.0_num
   DBDX	=iquants(10:12)
   DBDY	=iquants(13:15)
   DBDZ	=iquants(16:18)
   DBDT	=0.0_num
   !DEDX	=0.0_num
   !DEDY	=0.0_num
   !DEDZ	=0.0_num
   DEDX	=iquants(19:21)
   DEDY	=iquants(22:24)
   DEDZ	=iquants(25:27)
   DEDT	=0.0_num

END SUBROUTINE BOURDINFIELDS

END MODULE bourdin_fields

