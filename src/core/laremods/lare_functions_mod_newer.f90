MODULE lare_functions
!+ contains destagger, interpolation, and sub grid calculation of variables

USE GLOBAL

 IMPLICIT NONE

  PRIVATE
  PUBLIC :: stagger_bx, stagger_by, stagger_bz, stagger_bx_2d, stagger_by_2d
  PUBLIC :: linterp3d, linterp2d, linterp1d, T2d, t3d, str_cmp, check_dims

CONTAINS
!--------------------------------------
FUNCTION stagger_bx(var)
!+ function for de-staggering bx
    REAL(num), DIMENSION(:,:), INTENT(IN)  		:: var
    REAL(num), DIMENSION(SIZE(var,1), SIZE(var,2))	:: stagger_bx 
    INTEGER	:: mny, mnz
    
    mny=size(var,1)
    mnz=size(var,2)
    
    !staggering right
    stagger_bx(1:mny-1,1:mnz)=1.5_num*(var(1:mny-1,1:mnz)+var(2:mny,1:mnz))
    stagger_bx(mny,1:mnz)=var(mny,1:mnz)
    !staggering up
    stagger_bx(1:mny,1:mnz-1)=1.5_num*(stagger_bx(1:mny,1:mnz-1)+stagger_bx(1:mny,2:mnz))
    stagger_bx(1:mny,mnz)=stagger_bx(1:mny,mnz)
       
    RETURN

END FUNCTION stagger_bx
!-------------------------------!
FUNCTION stagger_by(var)
!+ function for de-staggering by for arbitrary portion of grid

    REAL(num), DIMENSION(:,:), INTENT(IN)  		:: var
    REAL(num), DIMENSION(SIZE(var,1), SIZE(var,2))	:: stagger_by  
    INTEGER	:: mnx, mnz

    mnx=size(var,1)
    mnz=size(var,2)
    
    !staggering right
    stagger_by(1:mnx-1,1:mnz)=0.5_num*(var(1:mnx-1,1:mnz)+var(2:mnx,1:mnz))
    stagger_by(mnx,1:mnz)=var(mnx,1:mnz)
    !staggering up
    stagger_by(1:mnx,1:mnz-1)=0.5_num*(stagger_by(1:mnx,1:mnz-1)+stagger_by(1:mnx,2:mnz))
    stagger_by(1:mnx,mnz)=stagger_by(1:mnx,mnz)
       
    RETURN

END FUNCTION stagger_by
!-------------------------------!
FUNCTION stagger_bz(var)
!+ function for de-staggering bz (for arbitrary portion of grid)
    REAL(num), DIMENSION(:,:), INTENT(IN)  		:: var
    REAL(num), DIMENSION(SIZE(var,1), SIZE(var,2))	:: stagger_bz    
    INTEGER	:: mnx, mny
    
    mnx=size(var,1)
    mny=size(var,2)
    
!+ function for de-staggering bz
!  (separate function necessary since nx, ny, nz may be different values)    
    !staggering right
    stagger_bz(1:mnx-1,1:mny)=0.5_num*(var(1:mnx-1,1:mny)+var(2:mnx,1:mny))
    stagger_bz(mnx,1:mny)=var(mnx,1:mny)
    !staggering up
    stagger_bz(1:mnx,1:mny-1)=0.5_num*(stagger_bz(1:mnx,1:mny-1)+stagger_bz(1:mnx,2:mny))
    stagger_bz(1:mnx,mny)=stagger_bz(1:mnx,mny)
       
    RETURN

END FUNCTION stagger_bz
!--------------------------------------!
FUNCTION stagger_bx_2d(var)
!+ function for de-staggering bx
    REAL(num), DIMENSION(:), INTENT(IN)  	:: var
    REAL(num), DIMENSION(SIZE(var,1))		:: stagger_bx_2d 
    INTEGER	:: mny  
    mny=size(var,1)
    
    !staggering right
    stagger_bx_2d(1:mny-1)=0.5_num*(var(1:mny-1)+var(2:mny))
    stagger_bx_2d(mny)=var(mny)
       
    RETURN

END FUNCTION stagger_bx_2d
!-------------------------------!
FUNCTION stagger_by_2d(var)
!+ function for de-staggering in a 2D grid
    REAL(num), DIMENSION(:), INTENT(IN)  	:: var
    REAL(num), DIMENSION(SIZE(var,1))		:: stagger_by_2d 
    INTEGER	:: mnx
 
    mnx=size(var,1)
    !staggering right
    stagger_by_2d(1:mnx-1)=0.5_num*(var(1:mnx-1)+var(2:mnx))
    stagger_by_2d(mnx)=var(mnx)
       
    RETURN

END FUNCTION stagger_by_2d
!--------------------------------------!
FUNCTION linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
! peforms sub-grid interpolation of a given quantity
! April 2014: attempted update to evaluate derivs of variable passed in too, into 4d array
! - ultimately decided this add on was unecessary

    REAL(num)			:: linterp3d
    REAL(num), INTENT(IN)	:: dx, dy, dz
    REAL(num), INTENT(IN)	:: f000,f100,f010,f110,f001,f101,f011,f111
    REAL(num)			:: a,b,c,d,e,f,g,h

    IF((abs(dx).gt.1.0_num).OR.(abs(dy).gt.1.0_num).OR.(abs(dz).gt.1.0_num)) THEN
      PRINT*, 'CRITICAL ERROR: dx, dy or dz is TOO BIG'
      STOP
    ENDIF

    a=f000
    b=f100-f000
    c=f010-f000
    d=f110-f100-f010+f000
    e=f001-f000
    f=f101-f100-f001+f000
    g=f011-f010-f001+f000
    h=f111-f110-f101-f011+f100+f010+f001-f000
  
  linterp3d=a+b*dx+c*dy+d*dx*dy+e*dz+f*dx*dz+g*dy*dz+h*dx*dy*dz
  
END FUNCTION linterp3d
!-------------------------------!
FUNCTION linterp2d(dx,dy,f00,f10,f01,f11)
!2D (bilinear) interpolation, see Haynes & Parnell (2007)


    REAL(num)			:: linterp2d
    REAL(num), INTENT(IN)	:: dx, dy
    REAL(num), INTENT(IN)	:: f00,f10,f01,f11
    REAL(num)			:: a,b,c,d

    IF((abs(dx).gt.1.0_num).OR.(abs(dy).gt.1.0_num)) THEN
      PRINT*, 'CRITICAL ERROR: dx, dy is TOO BIG'
      STOP
    ENDIF

    a=f00
    b=f10-f00
    c=f01-f00
    d=f11-f10-f01+f00
    
  linterp2d=a+b*dx+c*dy+d*dx*dy

END FUNCTION linterp2d
!-------------------------------;
FUNCTION linterp1d(dx,f0,f1)
! performs 1d sub-grid interpolation of a given quantity

    REAL(num)			:: linterp1d
    REAL(num), INTENT(IN)	:: dx, f0, f1
 
    IF(abs(dx).gt.1.0_num) THEN
      PRINT*, 'CRITICAL ERROR: dx is TOO BIG, dx = ', dx
      !STOP
      linterp1d = -12345678.9_num
      return
    ENDIF
  
  linterp1d=(1.0_num-dx)*f0+dx*f1

END FUNCTION linterp1d
!-------------------------------------------------------------------------------------
FUNCTION T2d(R,T)
! function to calculate local quantities at the particle position on lare2d grid
! Dec 2015 - 2d version of T3d

   REAL(num), DIMENSION(3), INTENT(IN)		:: R		!actual position
   REAL(num), INTENT(IN)			:: T
   REAL(num), DIMENSION(39)			:: T2d
   REAL(num)					:: temp,  modj, dgt, odgt
   REAL(num), DIMENSION(2)			:: dg, odg, coffset
   REAL(num), DIMENSION(:), ALLOCATABLE		:: bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt, rhot, temperaturet, etat
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dbxdxt,dbxdyt,dbydxt,dbydyt,dbzdxt,dbzdyt
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dExdxt,dExdyt,dEydxt,dEydyt,dEzdxt,dEzdyt
   REAL(num), DIMENSION(:, :), ALLOCATABLE	:: mbx, mby, mbz, mEx, mEy, mEz, mjx, mjy, mjz, mvx, mvy, mvz
   REAL(num), DIMENSION(:, :), ALLOCATABLE	:: dmbxdx,dmbxdy,dmbydx,dmbydy,dmbzdx,dmbzdy
   REAL(num), DIMENSION(:, :), ALLOCATABLE	:: dmExdx,dmExdy,dmEydx,dmEydy,dmEzdx,dmEzdy
   REAL(num), DIMENSION(:, :), ALLOCATABLE	:: meta
   INTEGER, DIMENSION(3) 			:: l		! set at value it could never reach!
   INTEGER					:: jjx, jjy, jjz,jjt, rpt
   LOGICAL					:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE.
 
    temp=0.0_num
    dg=(/myx(2)-myx(1),myy(2)-myy(1)/)	! grid spacing
    odg=(/1.0_num/dg(1),1.0_num/dg(2)/)	! one over grid spacing
    l=(/-nx,-ny,-nframes/)					! initial value of l, set to silly value (as <=0 triggers flag)   

       
  ! --STEP ONE-- !
  ! first, locate the (x,y,z) position on the grid
  ! particle found between (jjx,jjy,jjz) and (jjx+1,jjy+1,jjz+1)
  ! need four gridpoints either side to guarantee we can create derivs successfully   
    
   !print *,'R = ',R
   if (isnan(R(1)) .or. isnan(R(2))) then
        print *, 'lare_function position nan, stopping.'
        stop
   endif
   DO jjx=4,nx-4
    IF ((R(1).ge.myx(jjx)).and.(R(1).lt.myx(jjx+1))) THEN 
      l(1)=jjx
      !print *,'FOUND X'
      EXIT
    !ELSE
    !  PRINT *, 'CANNOT FIND R(1) IN LARE x RANGE'
    ENDIF
   ENDDO
   DO jjy=4,ny-4
    IF ((R(2).ge.myy(jjy)).and.(R(2).lt.myy(jjy+1))) THEN
      l(2)=jjy
      !print *,'FOUND y'
      EXIT
    !ELSE
    !  PRINT *, 'CANNOT FIND R(2) IN LARE y RANGE'
    ENDIF
   ENDDO


! No guarantee we have more than one frame. IF we have one, this routine doesn't bother interpolating in time
  IF (nframes.gt.1) THEN
   dgt=ltimes(2)-ltimes(1)
   odgt=1.0_num/dgt
   !print*, ltimes(1), ltimes(0)
   DO jjt=1,nframes
    IF ((T.GE.ltimes(jjt)).AND.((T.LT.ltimes(jjt+1)))) THEN
      l(3)=jjt
      dgt=ltimes(jjt+1)-ltimes(jjt)
      odgt=1.0_num/dgt
      EXIT
    ELSE
      !PRINT *, 'CANNOT FIND TIME IN LARE TIME RANGE'
      !print *, 'T = ',T
      !print *, 'ltimes = ',ltimes
      !stop
    ENDIF
   ENDDO
   rpt=1
   ALLOCATE(bxt(2), byt(2), bzt(2),vxt(2), vyt(2), vzt(2), Ext(2), Eyt(2), Ezt(2), jxt(2), jyt(2), jzt(2), rhot(2), temperaturet(2),etat(2))
   ALLOCATE(dbxdxt(2),dbxdyt(2),dbydxt(2),dbydyt(2),dbzdxt(2),dbzdyt(2))
   ALLOCATE(dExdxt(2),dExdyt(2),dEydxt(2),dEydyt(2),dEzdxt(2),dEzdyt(2))
  ELSE
   l(3)=1
   rpt=0
   ALLOCATE(bxt(1), byt(1), bzt(1),vxt(1), vyt(1), vzt(1), Ext(1), Eyt(1), Ezt(1), jxt(1), jyt(1), jzt(1), rhot(1), temperaturet(1), etat(1))
   ALLOCATE(dbxdxt(1),dbxdyt(1),dbydxt(1),dbydyt(1),dbzdxt(1),dbzdyt(1))
   ALLOCATE(dExdxt(1),dExdyt(1),dEydxt(1),dEydyt(1),dEzdxt(1),dEzdyt(1))
  !what happens if there is one frame to work with?
  ENDIF 

  !print *,'l = ',l
  IF (minval(l).le.0) THEN
  ! CHECK: have we located the particle on the grid?
   PRINT*, 'particle identified at', L
   PRINT*, 'x/y grid bounds-> 1:', nx, ny
   PRINT*, 't grid bounds->', ltimes(1),':',ltimes(nframes)
   print*, 'particle position R = ', R
   STOP
  ENDIF

!  print*, l(2)
   coffset=(/(R(1)-myx(l(1)))*odg(1),(R(2)-myy(l(2)))*odg(2)/)
   
   ! --STEP TWO-- !
   ! begin interpolation of basic variables we already have, e.g. bx, by, bz, vx, vy, vz
       
  DO it=0,rpt	! NEED TO REPEAT FOR INTERPOLATION BETWEEN FRAMES, JT DEC 2015
  !(dx,dy,f00,f10,f01,f11)
   !print *,'l(1) = ',l(1), size(rho)
   !rhot(it+1) = rho(l(1),l(2),1,l(3)+it)
   !temperaturet(it+1) = temperature(l(1),l(2),1,l(3)+it)
   !etat(it+1) = eta(l(1),l(2),1,l(3)+it)

   temp=linterp2d(coffset(1),coffset(2), &
   bx(l(1),l(2),1,l(3)+it),bx(l(1)+1,l(2),1,l(3)+it),bx(l(1),l(2)+1,1,l(3)+it),bx(l(1)+1,l(2)+1,1,l(3)+it))
   bxt(it+1)=temp		!1
   
   temp=linterp2d(coffset(1),coffset(2), &
   by(l(1),l(2),1,l(3)+it),by(l(1)+1,l(2),1,l(3)+it),by(l(1),l(2)+1,1,l(3)+it),by(l(1)+1,l(2)+1,1,l(3)+it))
   byt(it+1)=temp		!2
   
   temp=linterp2d(coffset(1),coffset(2), &
   bz(l(1),l(2),1,l(3)+it),bz(l(1)+1,l(2),1,l(3)+it),bz(l(1),l(2)+1,1,l(3)+it),bz(l(1)+1,l(2)+1,1,l(3)+it))
   bzt(it+1)=temp		!3
 !  print*, 'it', it!, 'bzt(it)=', linterp2d(coffset(1),coffset(2), &
   !bz(l(1),l(2),0,l(3)+it),bz(l(1)+1,l(2),0,l(3)+it),bz(l(1),l(2)+1,0,l(3)+it),bz(l(1)+1,l(2)+1,0,l(3)+it))
 !  print*, bz(l(1),l(2),0,l(3)+it), bz(l(1)+1,l(2),0,l(3)+it)
 !  print*, bz(l(1),l(2)+1,0,l(3)+it),bz(l(1)+1,l(2)+1,0,l(3)+it)
   
   temp=linterp2d(coffset(1),coffset(2), &
   vx(l(1),l(2),1,l(3)+it),vx(l(1)+1,l(2),1,l(3)+it),vx(l(1),l(2)+1,1,l(3)+it),vx(l(1)+1,l(2)+1,1,l(3)+it))
   vxt(it+1)=temp		!4
   
   temp=linterp2d(coffset(1),coffset(2), &
   vy(l(1),l(2),1,l(3)+it),vy(l(1)+1,l(2),1,l(3)+it),vy(l(1),l(2)+1,1,l(3)+it),vy(l(1)+1,l(2)+1,1,l(3)+it))
   vyt(it+1)=temp		!5
   
   temp=linterp2d(coffset(1),coffset(2), &
   vz(l(1),l(2),1,l(3)+it),vz(l(1)+1,l(2),1,l(3)+it),vz(l(1),l(2)+1,1,l(3)+it),vz(l(1)+1,l(2)+1,1,l(3)+it))
   vzt(it+1)=temp		!6

   temp=linterp2d(coffset(1),coffset(2), &
   rho(l(1),l(2),1,l(3)+it),rho(l(1)+1,l(2),1,l(3)+it),rho(l(1),l(2)+1,1,l(3)+it),rho(l(1)+1,l(2)+1,1,l(3)+it))
   rhot(it+1)=temp		!6

   temp=linterp2d(coffset(1),coffset(2), &
   temperature(l(1),l(2),1,l(3)+it),temperature(l(1)+1,l(2),1,l(3)+it),temperature(l(1),l(2)+1,1,l(3)+it),temperature(l(1)+1,l(2)+1,1,l(3)+it))
   temperaturet(it+1)=temp		!6

   !print *, 'lare_func eta', eta(l(1),l(2),1,l(3)+it),eta(l(1)+1,l(2),1,l(3)+it),eta(l(1),l(2)+1,1,l(3)+it),eta(l(1)+1,l(2)+1,1,l(3)+it), 'pos', myx(l(1)), myy(l(2)), R(1), R(2)
   !print *, 'lare_func eta', eta(l(1),l(2),1,l(3)+it),eta(l(1)+1,l(2),1,l(3)+it),eta(l(1),l(2)+1,1,l(3)+it),eta(l(1)+1,l(2)+1,1,l(3)+it), l(3) + it
   temp=linterp2d(coffset(1),coffset(2), &
   eta(l(1),l(2),1,l(3)+it),eta(l(1)+1,l(2),1,l(3)+it),eta(l(1),l(2)+1,1,l(3)+it),eta(l(1)+1,l(2)+1,1,l(3)+it))
   etat(it+1)=temp		!6
   !etat = 0.005_num
   !
   !! --STEP THREE--!
   ! create temporary arrays around target cell and calculate derivs using extra cells
   ! - how many cells depend on which finite difference routine used
   ! - need db/dr->j and etadj/dr for E-derivs
   ! JT OCT2014: commented out O(h) derivs, and insert O(h^2)
   
   ! begin with temporary B arrays:
   ALLOCATE(mbx(-4:5,-4:5))
   ALLOCATE(mby(-4:5,-4:5))
   ALLOCATE(mbz(-4:5,-4:5))

   mbx=bx(l(1)-4:l(1)+5,l(2)-4:l(2)+5,1,l(3)+it)
   mby=by(l(1)-4:l(1)+5,l(2)-4:l(2)+5,1,l(3)+it)
   mbz=bz(l(1)-4:l(1)+5,l(2)-4:l(2)+5,1,l(3)+it)


   ! calculate x-derivs first..
   ALLOCATE(dmbxdx(-2:3,-4:5))
   ALLOCATE(dmbydx(-2:3,-4:5))
   ALLOCATE(dmbzdx(-2:3,-4:5))
   
    DO iy=-4,5
     DO ix=-2,3
      dmbxdx(ix,iy)=(-mbx(ix+2,iy)+8.0_num*mbx(ix+1,iy)-8.0_num*mbx(ix-1,iy)+mbx(ix-2,iy))*odg(1)*oneotwelve
      dmbydx(ix,iy)=(-mby(ix+2,iy)+8.0_num*mby(ix+1,iy)-8.0_num*mby(ix-1,iy)+mby(ix-2,iy))*odg(1)*oneotwelve
      dmbzdx(ix,iy)=(-mbz(ix+2,iy)+8.0_num*mbz(ix+1,iy)-8.0_num*mbz(ix-1,iy)+mbz(ix-2,iy))*odg(1)*oneotwelve
     END DO
    END DO

   
   ! ..followed by y-derivs..
   ALLOCATE(dmbxdy(-4:5,-2:3))
   ALLOCATE(dmbydy(-4:5,-2:3))
   ALLOCATE(dmbzdy(-4:5,-2:3))
 
    DO iy=-2,3
     DO ix=-4,5
      dmbxdy(ix,iy)=(-mbx(ix,iy+2)+8.0_num*mbx(ix,iy+1)-8.0_num*mbx(ix,iy-1)+mbx(ix,iy-2))*odg(2)*oneotwelve
      dmbydy(ix,iy)=(-mby(ix,iy+2)+8.0_num*mby(ix,iy+1)-8.0_num*mby(ix,iy-1)+mby(ix,iy-2))*odg(2)*oneotwelve
      dmbzdy(ix,iy)=(-mbz(ix,iy+2)+8.0_num*mbz(ix,iy+1)-8.0_num*mbz(ix,iy-1)+mbz(ix,iy-2))*odg(2)*oneotwelve
     END DO
    END DO

   ! .. and finally z derivs (EXCEPT WE DON'T HAVE ANY!)

   ! NB the arrays are odd shapes depending on deriv direction but the target cell remains at [0:1,0:1,0:1]
   ! also need temporary velocity arrays, in order to calc E=-etaJ+vxB

   ALLOCATE(mvx(-2:3,-2:3))
   ALLOCATE(mvy(-2:3,-2:3))
   ALLOCATE(mvz(-2:3,-2:3))
   mvx=vx(l(1)-2:l(1)+3,l(2)-2:l(2)+3,1,l(3)+it)
   mvy=vy(l(1)-2:l(1)+3,l(2)-2:l(2)+3,1,l(3)+it)
   mvz=vz(l(1)-2:l(1)+3,l(2)-2:l(2)+3,1,l(3)+it)
   
   
   ! now calculate temporary currents and also local eta values

   ALLOCATE(mjx(-2:3,-2:3))
   ALLOCATE(mjy(-2:3,-2:3))
   ALLOCATE(mjz(-2:3,-2:3))
   ALLOCATE(mEx(-2:3,-2:3))
   ALLOCATE(mEy(-2:3,-2:3))
   ALLOCATE(mEz(-2:3,-2:3))
   ALLOCATE(meta(-2:3,-2:3))
   meta=eta(l(1)-2:l(1)+3,l(2)-2:l(2)+3,1,l(3)+it)
   
   
    DO iy=-2,3
     DO ix=-2,3
      mjx(ix,iy)=dmbzdy(ix,iy)!-dmbydz(ix,iy,iz)
      mjy(ix,iy)=-dmbzdx(ix,iy)!+dmbxdz(ix,iy,iz)
      mjz(ix,iy)=dmbydx(ix,iy)-dmbxdy(ix,iy)
      modj=(mjx(ix,iy)*mjx(ix,iy)+mjy(ix,iy)*mjy(ix,iy)+mjz(ix,iy)*mjz(ix,iy))**0.5_num   
      !meta(ix,iy)=0.5_num*(tanh((modj-jcrit)/rwidth)+1.0_num)*eta          ! smooth transition between eta = 0 and non-zero
      !meta(ix,iy)=eta
      !if (eta .ne. 0.0_num) print *, 'lare_functions eta = ', eta
     END DO
    END DO


   ! now interpolate derivs to correct point:  
   dbxdxt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmbxdx(0,0), dmbxdx(1,0), dmbxdx(0,1), dmbxdx(1,1))	! 7
   dbydxt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmbydx(0,0), dmbydx(1,0), dmbydx(0,1), dmbydx(1,1))	! 8
   dbzdxt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmbzdx(0,0), dmbzdx(1,0), dmbzdx(0,1), dmbzdx(1,1))	! 9
   dbxdyt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmbxdy(0,0), dmbxdy(1,0), dmbxdy(0,1), dmbxdy(1,1)) ! 10
   dbydyt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmbydy(0,0), dmbydy(1,0), dmbydy(0,1), dmbydy(1,1)) ! 11
   dbzdyt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmbzdy(0,0), dmbzdy(1,0), dmbzdy(0,1), dmbzdy(1,1)) ! 12            
   
   DEALLOCATE(dmbxdx)
   DEALLOCATE(dmbydx)
   DEALLOCATE(dmbzdx)
   DEALLOCATE(dmbxdy)
   DEALLOCATE(dmbydy)
   DEALLOCATE(dmbzdy)


   ! interpolate j too at this stage:
   !mjx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   jxt(it+1)=linterp2d(coffset(1),coffset(2), &
    mjx(0,0), mjx(1,0), mjx(0,1), mjx(1,1))
   jyt(it+1)=linterp2d(coffset(1),coffset(2), &
    mjy(0,0), mjy(1,0), mjy(0,1), mjy(1,1))
   jzt(it+1)=linterp2d(coffset(1),coffset(2), &
    mjz(0,0), mjz(1,0), mjz(0,1), mjz(1,1))
   
   ! --STEP FOUR--
   ! Ohm's Law: E=vxB-etaJ OR E=-etaJ
    !DO iy=-2,3
    ! DO ix=-2,3
    !  mEx(ix,iy)=meta(ix,iy)*mjx(ix,iy)-(mvy(ix,iy)*mbz(ix,iy)-mvz(ix,iy)*mby(ix,iy))
    !  mEy(ix,iy)=meta(ix,iy)*mjy(ix,iy)-(mvz(ix,iy)*mbx(ix,iy)-mvx(ix,iy)*mbz(ix,iy))
    !  mEz(ix,iy)=meta(ix,iy)*mjz(ix,iy)-(mvx(ix,iy)*mby(ix,iy)-mvy(ix,iy)*mbx(ix,iy))
    ! ENDDO
    !ENDDO
    Ext(it+1) = etat(it+1)*jxt(it+1) - (vyt(it+1)*bzt(it+1) - vzt(it+1)*byt(it+1))
    Eyt(it+1) = etat(it+1)*jyt(it+1) - (vzt(it+1)*bxt(it+1) - vxt(it+1)*bzt(it+1))
    Ezt(it+1) = etat(it+1)*jzt(it+1) - (vxt(it+1)*byt(it+1) - vyt(it+1)*bxt(it+1))

      
   
   DEALLOCATE(mvx)
   DEALLOCATE(mvy)
   DEALLOCATE(mvz)
   DEALLOCATE(mbx)
   DEALLOCATE(mby)
   DEALLOCATE(mbz)
   DEALLOCATE(mjx)
   DEALLOCATE(mjy)
   DEALLOCATE(mjz)    
   DEALLOCATE(meta)
   
   ! --STEP FIVE-- !
   ! finally, calculate derivatives of electric field:
   
   ALLOCATE(dmExdx(0:1,-2:3))
   ALLOCATE(dmEydx(0:1,-2:3))
   ALLOCATE(dmEzdx(0:1,-2:3))  
    DO iy=-2,3
     DO ix=0,1
      dmExdx(ix,iy)=oneotwelve*(-mEx(ix+2,iy)+8.0_num*mEx(ix+1,iy)-8.0_num*mEx(ix-1,iy)+mEx(ix-2,iy))*odg(1)
      dmEydx(ix,iy)=oneotwelve*(-mEy(ix+2,iy)+8.0_num*mEy(ix+1,iy)-8.0_num*mEy(ix-1,iy)+mEy(ix-2,iy))*odg(1)
      dmEzdx(ix,iy)=oneotwelve*(-mEz(ix+2,iy)+8.0_num*mEz(ix+1,iy)-8.0_num*mEz(ix-1,iy)+mEz(ix-2,iy))*odg(1)
     ENDDO
    ENDDO

   

   ALLOCATE(dmExdy(-2:3,0:1))
   ALLOCATE(dmEydy(-2:3,0:1))
   ALLOCATE(dmEzdy(-2:3,0:1))   
    DO iy=0,1
     DO ix=-2,3  
      dmExdy(ix,iy)=oneotwelve*(-mEx(ix,iy+2)+8.0_num*mEx(ix,iy+1)-8.0_num*mEx(ix,iy-1)+mEx(ix,iy-2))*odg(2)
      dmEydy(ix,iy)=oneotwelve*(-mEy(ix,iy+2)+8.0_num*mEy(ix,iy+1)-8.0_num*mEy(ix,iy-1)+mEy(ix,iy-2))*odg(2)
      dmEzdy(ix,iy)=oneotwelve*(-mEz(ix,iy+2)+8.0_num*mEz(ix,iy+1)-8.0_num*mEz(ix,iy-1)+mEz(ix,iy-2))*odg(2)
     ENDDO
    ENDDO  
   
   
 
   dExdxt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmExdx(0,0), dmExdx(1,0), dmExdx(0,1), dmExdx(1,1))
   dEydxt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmEydx(0,0), dmEydx(1,0), dmEydx(0,1), dmEydx(1,1))
   dEzdxt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmEzdx(0,0), dmEzdx(1,0), dmEzdx(0,1), dmEzdx(1,1))
   dExdyt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmExdy(0,0), dmExdy(1,0), dmExdy(0,1), dmExdy(1,1))
   dEydyt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmEydy(0,0), dmEydy(1,0), dmEydy(0,1), dmEydy(1,1))
   dEzdyt(it+1)=linterp2d(coffset(1),coffset(2), &
    dmEzdy(0,0), dmEzdy(1,0), dmEzdy(0,1), dmEzdy(1,1))

   
   DEALLOCATE(dmExdx)
   DEALLOCATE(dmEydx)
   DEALLOCATE(dmEzdx)
   DEALLOCATE(dmExdy)
   DEALLOCATE(dmEydy)
   DEALLOCATE(dmEzdy)


   !Ext(it+1)=linterp2d(coffset(1),coffset(2), &
   ! mEx(0,0), mEx(1,0), mEx(0,1), mEx(1,1))
   !Eyt(it+1)=linterp2d(coffset(1),coffset(2), &
   ! mEy(0,0), mEy(1,0), mEy(0,1), mEy(1,1))
   !Ezt(it+1)=linterp2d(coffset(1),coffset(2), &
   ! mEz(0,0), mEz(1,0), mEz(0,1), mEz(1,1))
      
   !if (Ext(it+1)*bxt(it+1) + Eyt(it+1)*byt(it+1) + Ezt(it+1)*bzt(it+1) .ne. 0) print *, 'lare_functions E parallel = ', &
   !  (Ext(it+1)*bxt(it+1) + Eyt(it+1)*byt(it+1) + Ezt(it+1)*bzt(it+1))/sqrt(bxt(it+1)*bxt(it+1) + byt(it+1)*byt(it+1) + bzt(it+1)*bzt(it+1))*Escl
    
   DEALLOCATE(mEx)
   DEALLOCATE(mEy)
   DEALLOCATE(mEz)
   
  END DO


   ! its possible to run lare3d in real units - lets make sure things are normalised if this happens
   IF (.not. lare_norm) THEN
    bxt=bxt/bscl
    byt=byt/bscl
    bzt=bzt/bscl
    vxt=vxt/vscl
    vyt=vyt/vscl
    vzt=vzt/vscl
    Ext=Ext/Escl
    Eyt=Eyt/Escl
    Ezt=Ezt/Escl
    jxt=jxt/bscl*lscl
    jyt=jyt/bscl*lscl
    jzt=jzt/bscl*lscl	
    dbxdxt=dbxdxt/bscl*lscl
    dbydxt=dbydxt/bscl*lscl
    dbzdxt=dbzdxt/bscl*lscl
    dbxdyt=dbxdyt/bscl*lscl
    dbydyt=dbydyt/bscl*lscl
    dbzdyt=dbzdyt/bscl*lscl
    dExdxt=dExdxt/Escl*lscl
    dEydxt=dEydxt/Escl*lscl
    dEzdxt=dEzdxt/Escl*lscl
    dExdyt=dExdyt/Escl*lscl
    dEydyt=dEydyt/Escl*lscl
    dEzdyt=dEzdyt/Escl*lscl
   ENDIF


   IF (nframes.gt.1) THEN

    T2d(1)=linterp1d((T-ltimes(l(3)))*odgt,bxt(1),bxt(2))
    if(T2d(1) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 1, T = ', T, ' ltimes = ', ltimes(l(3)), ' dgt = ', dgt
        stop
    endif
    T2d(2)=linterp1d((T-ltimes(l(3)))*odgt,byt(1),byt(2))
    if(T2d(2) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 2'
        stop
    endif
    T2d(3)=linterp1d((T-ltimes(l(3)))*odgt,bzt(1),bzt(2))
    if(T2d(3) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 3'
        stop
    endif
    T2d(4)=linterp1d((T-ltimes(l(3)))*odgt,vxt(1),vxt(2))
    if(T2d(4) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 4'
        stop
    endif
    T2d(5)=linterp1d((T-ltimes(l(3)))*odgt,vyt(1),vyt(2))
    if(T2d(5) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 5'
        stop
    endif
    T2d(6)=linterp1d((T-ltimes(l(3)))*odgt,vzt(1),vzt(2))
    if(T2d(6) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 6'
        stop
    endif
    T2d(7)=linterp1d((T-ltimes(l(3)))*odgt,Ext(1),Ext(2))
    if(T2d(7) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 7'
        stop
    endif
    T2d(8)=linterp1d((T-ltimes(l(3)))*odgt,Eyt(1),Eyt(2))
    if(T2d(8) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 8'
        stop
    endif
    T2d(9)=linterp1d((T-ltimes(l(3)))*odgt,Ezt(1),Ezt(2))   
    if(T2d(9) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 9'
        stop
    endif
    T2d(10)=linterp1d((T-ltimes(l(3)))*odgt,jxt(1),jxt(2))
    if(T2d(10) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 10'
        stop
    endif
    T2d(11)=linterp1d((T-ltimes(l(3)))*odgt,jyt(1),jyt(2))
    if(T2d(11) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 11'
        stop
    endif
    T2d(12)=linterp1d((T-ltimes(l(3)))*odgt,jzt(1),jzt(2))
    if(T2d(12) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 12'
        stop
    endif
    T2d(13)=linterp1d((T-ltimes(l(3)))*odgt,dbxdxt(1),dbxdxt(2))
    if(T2d(13) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 13'
        stop
    endif
    T2d(14)=linterp1d((T-ltimes(l(3)))*odgt,dbydxt(1),dbydxt(2))
    if(T2d(14) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 14'
        stop
    endif
    T2d(15)=linterp1d((T-ltimes(l(3)))*odgt,dbzdxt(1),dbzdxt(2))
    if(T2d(15) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 15'
        stop
    endif
    T2d(16)=linterp1d((T-ltimes(l(3)))*odgt,dbxdyt(1),dbxdyt(2))
    if(T2d(16) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 16'
        stop
    endif
    T2d(17)=linterp1d((T-ltimes(l(3)))*odgt,dbydyt(1),dbydyt(2))
    if(T2d(17) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 17'
        stop
    endif
    T2d(18)=linterp1d((T-ltimes(l(3)))*odgt,dbzdyt(1),dbzdyt(2))
    if(T2d(18) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 18'
        stop
    endif
    T2d(19)=0.0_num		! Z derivs explicitly set to zero in 2d case
    T2d(20)=0.0_num
    T2d(21)=0.0_num
    T2d(22)=linterp1d((T-ltimes(l(3)))*odgt,dExdxt(1),dExdxt(2))
    if(T2d(22) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 22'
        stop
    endif
    T2d(23)=linterp1d((T-ltimes(l(3)))*odgt,dEydxt(1),dEydxt(2))
    if(T2d(23) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 23'
        stop
    endif
    T2d(24)=linterp1d((T-ltimes(l(3)))*odgt,dEzdxt(1),dEzdxt(2))
    if(T2d(24) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 24'
        stop
    endif
    T2d(25)=linterp1d((T-ltimes(l(3)))*odgt,dExdyt(1),dExdyt(2))
    if(T2d(25) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 25'
        stop
    endif
    T2d(26)=linterp1d((T-ltimes(l(3)))*odgt,dEydyt(1),dEydyt(2))
    if(T2d(26) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 26'
        stop
    endif
    T2d(27)=linterp1d((T-ltimes(l(3)))*odgt,dEzdyt(1),dEzdyt(2))
    if(T2d(27) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 27'
        stop
    endif
    T2d(28)=0.0_num
    T2d(29)=0.0_num
    T2d(30)=0.0_num
    T2d(31)=(bxt(2)-bxt(1))*odgt
    T2d(32)=(byt(2)-byt(1))*odgt
    T2d(33)=(bzt(2)-bzt(1))*odgt
    T2d(34)=(Ext(2)-Ext(1))*odgt	! as we have multiple snapshots, we can calculate time derivs finally!
    T2d(35)=(Eyt(2)-Eyt(1))*odgt
    T2d(36)=(Ezt(2)-Ezt(1))*odgt
    T2d(37)=linterp1d((T-ltimes(l(3)))*odgt,rhot(1),rhot(2))
    if(T2d(37) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 37'
        stop
    endif
    T2d(38)=linterp1d((T-ltimes(l(3)))*odgt,temperaturet(1),temperaturet(2))
    if(T2d(38) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 38'
        stop
    endif
    T2d(39)=linterp1d((T-ltimes(l(3)))*odgt,etat(1),etat(2))
    if(T2d(39) .eq. -12345678.9_num) then 
        print *, 'lare_functions failed at 39'
        stop
    endif

    print *, 'larefields epar 0= ', (Ext(1)*Bxt(2) + Eyt(1)*Byt(2) + Ezt(1)*Bzt(2))/sqrt(Bxt(2)*Bxt(2) + Byt(2)*Byt(2) + Bzt(2)*Bzt(2))
    print *, 'larefields epar 1= ', (Ext(2)*Bxt(1) + Eyt(2)*Byt(1) + Ezt(2)*Bzt(1))/sqrt(Bxt(1)*Bxt(1) + Byt(1)*Byt(1) + Bzt(1)*Bzt(1))
    print *, 'larefields epar 2= ', (T2d(7)*T2d(1) + T2d(8)*T2d(2) + T2d(9)*T2d(3))/sqrt(T2d(1)*T2d(1) + T2d(2)*T2d(2) + T2d(3)*T2d(3))

    DEALLOCATE(bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt)
    DEALLOCATE(dbxdxt,dbxdyt,dbydxt,dbydyt,dbzdxt,dbzdyt)
    DEALLOCATE(dExdxt,dExdyt,dEydxt,dEydyt,dEzdxt,dEzdyt)

   ELSE
   
    T2d(1)=bxt(1)
    T2d(2)=byt(1)
    T2d(3)=bzt(1)
    T2d(4)=vxt(1)
    T2d(5)=vyt(1)
    T2d(6)=vzt(1)
    T2d(7)=Ext(1)
    T2d(8)=Eyt(1)
    T2d(9)=Ezt(1)   
    T2d(10)=jxt(1)
    T2d(11)=jyt(1)
    T2d(12)=jzt(1)
    T2d(13)=dbxdxt(1)
    T2d(14)=dbydxt(1)
    T2d(15)=dbzdxt(1)
    T2d(16)=dbxdyt(1)
    T2d(17)=dbydyt(1)
    T2d(18)=dbzdyt(1)
    T2d(19)=0.0_num		! Z derivs explicitly set to zero in 2d case
    T2d(20)=0.0_num
    T2d(21)=0.0_num
    T2d(22)=dExdxt(1)
    T2d(23)=dEydxt(1)
    T2d(24)=dEzdxt(1)
    T2d(25)=dExdyt(1)
    T2d(26)=dEydyt(1)
    T2d(27)=dEzdyt(1)
    T2d(28)=0.0_num
    T2d(29)=0.0_num
    T2d(30)=0.0_num
    T2d(31:36)=0.0_num
    T2d(37)=rhot(1)
    T2d(38)=temperaturet(1)
    T2d(39)=etat(1)

   ENDIF
   

   RETURN


END FUNCTION T2d
!-------------------------------------------------------------------------------------
FUNCTION T3d(R,T)
! function to calculate local quantities at the particle position
! (slimmed down version of f3d with only local interpolation!).
! OCT2014: updated to include higher order derivs to increase accuracy

   REAL(num), DIMENSION(3), INTENT(IN)		:: R		!actual position
   REAL(num), INTENT(IN)			:: T
   REAL(num), DIMENSION(36)			:: T3d
   REAL(num)					:: temp,  modj, dgt, odgt
   REAL(num), DIMENSION(3)			:: dg, odg, coffset
   REAL(num), DIMENSION(:), ALLOCATABLE		:: bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dbxdxt,dbxdyt,dbxdzt,dbydxt,dbydyt,dbydzt,dbzdxt,dbzdyt,dbzdzt
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dExdxt,dExdyt,dExdzt,dEydxt,dEydyt,dEydzt,dEzdxt,dEzdyt,dEzdzt
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: mbx, mby, mbz, mEx, mEy, mEz, mjx, mjy, mjz, mvx, mvy, mvz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: dmbxdx,dmbxdy,dmbxdz,dmbydx,dmbydy,dmbydz,dmbzdx,dmbzdy,dmbzdz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: dmExdx,dmExdy,dmExdz,dmEydx,dmEydy,dmEydz,dmEzdx,dmEzdy,dmEzdz
   REAL(num), DIMENSION(:, :, :), ALLOCATABLE	:: meta
   INTEGER, DIMENSION(4) 			:: l		! set at value it could never reach!
   INTEGER					:: jjx, jjy, jjz,jjt, rpt
   LOGICAL					:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE.
 
      temp=0.0_num
    dg=(/myx(2)-myx(1),myy(2)-myy(1),myz(2)-myz(1)/)	! grid spacing
    odg=(/1.0_num/dg(1),1.0_num/dg(2),1.0_num/dg(3)/)	! one over grid spacing
    l=(/-nx,-ny,-nz,-nframes/)					! initial value of l, set to silly value (as >=0 triggers flag)   
    
    
       
  ! --STEP ONE-- !
  ! first, locate the (x,y,z) position on the grid
  ! particle found between (jjx,jjy,jjz) and (jjx+1,jjy+1,jjz+1)   
     
   DO jjx=4,nx-4
    IF ((R(1).ge.myx(jjx)).and.(R(1).lt.myx(jjx+1))) THEN 
      l(1)=jjx
      EXIT
    ENDIF
   ENDDO
   !IF (l(1).eq.(-nx)) fxflag=.TRUE.
   DO jjy=4,ny-4
    IF ((R(2).ge.myy(jjy)).and.(R(2).lt.myy(jjy+1))) THEN
      l(2)=jjy
      EXIT
    ENDIF
   ENDDO
   !IF (l(2).eq.(-ny)) fyflag=.TRUE.
   DO jjz=4,nz-4
    IF ((R(3).ge.myz(jjz)).and.(R(3).lt.myz(jjz+1))) THEN
      l(3)=jjz
      EXIT
    ENDIF
   ENDDO

! No guarantee we have more than one frame. IF we have one, this routine doesn't bother interpolating in time
  IF (nframes.gt.1) THEN
   dgt=ltimes(2)-ltimes(1)
   odgt=1.0_num/dgt
   DO jjt=1,nframes
    IF ((T.GE.ltimes(jjt)).AND.((T.LT.ltimes(jjt+1)))) THEN
      l(4)=jjt
      EXIT
    ELSE
      PRINT *, 'CANNOT FIND TIME IN LARE TIME RANGE'
    ENDIF
   ENDDO
   rpt=1
   ALLOCATE(bxt(2), byt(2), bzt(2),vxt(2), vyt(2), vzt(2), Ext(2), Eyt(2), Ezt(2), jxt(2), jyt(2), jzt(2))
   ALLOCATE(dbxdxt(2),dbxdyt(2),dbxdzt(2),dbydxt(2),dbydyt(2),dbydzt(2),dbzdxt(2),dbzdyt(2),dbzdzt(2))
   ALLOCATE(dExdxt(2),dExdyt(2),dExdzt(2),dEydxt(2),dEydyt(2),dEydzt(2),dEzdxt(2),dEzdyt(2),dEzdzt(2))
  ELSE
   l(4)=1
   rpt=0
   ALLOCATE(bxt(1), byt(1), bzt(1),vxt(1), vyt(1), vzt(1), Ext(1), Eyt(1), Ezt(1), jxt(1), jyt(1), jzt(1))
   ALLOCATE(dbxdxt(1),dbxdyt(1),dbxdzt(1),dbydxt(1),dbydyt(1),dbydzt(1),dbzdxt(1),dbzdyt(1),dbzdzt(1))
   ALLOCATE(dExdxt(1),dExdyt(1),dExdzt(1),dEydxt(1),dEydyt(1),dEydzt(1),dEzdxt(1),dEzdyt(1),dEzdzt(1))
  !what happens if there is one frame to work with?
  ENDIF 


   coffset=(/(R(1)-myx(l(1)))*odg(1),(R(2)-myy(l(2)))*odg(2),(R(3)-myz(l(3)))*odg(3)/)

   
   ! --STEP TWO-- !
   ! begin interpolation of basic variables we already have, e.g. bx, by, bz, vx, vy, vz
       
  DO it=0,rpt	! NEED TO REPEAT FOR INTERPOLATION BETWEEN FRAMES, JT DEC 2015
  
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!1   
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   bx(l(1),l(2),l(3),l(4)+it),bx(l(1)+1,l(2),l(3),l(4)+it),bx(l(1),l(2)+1,l(3),l(4)+it),&
   bx(l(1)+1,l(2)+1,l(3),l(4)+it),bx(l(1),l(2),l(3)+1,l(4)+it),bx(l(1)+1,l(2),l(3)+1,l(4)+it), &
   bx(l(1),l(2)+1,l(3)+1,l(4)+it),bx(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   bxt(it+1)=temp
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!2
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   by(l(1),l(2),l(3),l(4)+it),by(l(1)+1,l(2),l(3),l(4)+it),by(l(1),l(2)+1,l(3),l(4)+it),&
   by(l(1)+1,l(2)+1,l(3),l(4)+it),by(l(1),l(2),l(3)+1,l(4)+it),by(l(1)+1,l(2),l(3)+1,l(4)+it), &
   by(l(1),l(2)+1,l(3)+1,l(4)+it),by(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   byt(it+1)=temp

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!3
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   bz(l(1),l(2),l(3),l(4)+it),bz(l(1)+1,l(2),l(3),l(4)+it),bz(l(1),l(2)+1,l(3),l(4)+it),&
   bz(l(1)+1,l(2)+1,l(3),l(4)+it),bz(l(1),l(2),l(3)+1,l(4)+it),bz(l(1)+1,l(2),l(3)+1,l(4)+it), &
   bz(l(1),l(2)+1,l(3)+1,l(4)+it),bz(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   bzt(it+1)=temp

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!4
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   vx(l(1),l(2),l(3),l(4)+it),vx(l(1)+1,l(2),l(3),l(4)+it),vx(l(1),l(2)+1,l(3),l(4)+it),&
   vx(l(1)+1,l(2)+1,l(3),l(4)+it),vx(l(1),l(2),l(3)+1,l(4)+it),vx(l(1)+1,l(2),l(3)+1,l(4)+it), &
   vx(l(1),l(2)+1,l(3)+1,l(4)+it),vx(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   vxt(it+1)=temp

   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!5
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   vy(l(1),l(2),l(3),l(4)+it),vy(l(1)+1,l(2),l(3),l(4)+it),vy(l(1),l(2)+1,l(3),l(4)+it),&
   vy(l(1)+1,l(2)+1,l(3),l(4)+it),vy(l(1),l(2),l(3)+1,l(4)+it),vy(l(1)+1,l(2),l(3)+1,l(4)+it), &
   vy(l(1),l(2)+1,l(3)+1,l(4)+it),vy(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   vyt(it+1)=temp
   
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!6
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   vz(l(1),l(2),l(3),l(4)+it),vz(l(1)+1,l(2),l(3),l(4)+it),vz(l(1),l(2)+1,l(3),l(4)+it),&
   vz(l(1)+1,l(2)+1,l(3),l(4)+it),vz(l(1),l(2),l(3)+1,l(4)+it),vz(l(1)+1,l(2),l(3)+1,l(4)+it), &
   vz(l(1),l(2)+1,l(3)+1,l(4)+it),vz(l(1)+1,l(2)+1,l(3)+1,l(4)+it))
   vzt(it+1)=temp
   
   ! --STEP THREE--!
   ! create temporary arrays around target cell and calculate derivs using extra cells
   ! - how many cells depend on which finite difference routine used
   ! - need db/dr->j and etadj/dr for E-derivs
   ! JT OCT2014: commented out O(h) derivs, and insert O(h^2)
   
   ! begin with temporary B arrays:
   ALLOCATE(mbx(-4:5,-4:5,-4:5))
   ALLOCATE(mby(-4:5,-4:5,-4:5))
   ALLOCATE(mbz(-4:5,-4:5,-4:5))
   !ALLOCATE(mbx(-2:3,-2:3,-2:3))
   !ALLOCATE(mby(-2:3,-2:3,-2:3))
   !ALLOCATE(mbz(-2:3,-2:3,-2:3))
     
   !mbx=bx(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3)
   !mby=by(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3)
   !mbz=bz(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3)
   mbx=bx(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5,l(4)+it)
   mby=by(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5,l(4)+it)
   mbz=bz(l(1)-4:l(1)+5,l(2)-4:l(2)+5,l(3)-4:l(3)+5,l(4)+it)

   !ALLOCATE(dmbxdx(-1:2,-2:3,-2:3))
   !ALLOCATE(dmbydx(-1:2,-2:3,-2:3))
   !ALLOCATE(dmbzdx(-1:2,-2:3,-2:3))
   !DO iz=-2,3   
   ! DO iy=-2,3
   !  DO ix=-1,2
   !   dmbxdx(ix,iy,iz)=0.5_num*(mbx(ix+1,iy,iz)-mbx(ix-1,iy,iz))*odg(1)
   !   dmbydx(ix,iy,iz)=0.5_num*(mby(ix+1,iy,iz)-mby(ix-1,iy,iz))*odg(1)
   !   dmbzdx(ix,iy,iz)=0.5_num*(mbz(ix+1,iy,iz)-mbz(ix-1,iy,iz))*odg(1)
   !  ENDDO
   ! ENDDO
   !ENDDO
   ! calculate x-derivs first..
   ALLOCATE(dmbxdx(-2:3,-4:5,-4:5))
   ALLOCATE(dmbydx(-2:3,-4:5,-4:5))
   ALLOCATE(dmbzdx(-2:3,-4:5,-4:5))

   DO iz=-4,5   
    DO iy=-4,5
     DO ix=-2,3
      dmbxdx(ix,iy,iz)=(-mbx(ix+2,iy,iz)+8.0_num*mbx(ix+1,iy,iz)-8.0_num*mbx(ix-1,iy,iz)+mbx(ix-2,iy,iz))*odg(1)*oneotwelve
      dmbydx(ix,iy,iz)=(-mby(ix+2,iy,iz)+8.0_num*mby(ix+1,iy,iz)-8.0_num*mby(ix-1,iy,iz)+mby(ix-2,iy,iz))*odg(1)*oneotwelve
      dmbzdx(ix,iy,iz)=(-mbz(ix+2,iy,iz)+8.0_num*mbz(ix+1,iy,iz)-8.0_num*mbz(ix-1,iy,iz)+mbz(ix-2,iy,iz))*odg(1)*oneotwelve
     END DO
    END DO
   END DO

   
   !ALLOCATE(dmbxdy(-2:3,-1:2,-2:3))
   !ALLOCATE(dmbydy(-2:3,-1:2,-2:3))
   !ALLOCATE(dmbzdy(-2:3,-1:2,-2:3))
   !DO iz=-2,3   
   ! DO iy=-1,2
   !  DO ix=-2,3
   !   dmbxdy(ix,iy,iz)=0.5_num*(mbx(ix,iy+1,iz)-mbx(ix,iy-1,iz))*odg(2)
   !   dmbydy(ix,iy,iz)=0.5_num*(mby(ix,iy+1,iz)-mby(ix,iy-1,iz))*odg(2)
   !   dmbzdy(ix,iy,iz)=0.5_num*(mbz(ix,iy+1,iz)-mbz(ix,iy-1,iz))*odg(2)
   !  ENDDO
   ! ENDDO
   !ENDDO
   ! ..followed by y-derivs..
   ALLOCATE(dmbxdy(-4:5,-2:3,-4:5))
   ALLOCATE(dmbydy(-4:5,-2:3,-4:5))
   ALLOCATE(dmbzdy(-4:5,-2:3,-4:5))

   DO iz=-4,5   
    DO iy=-2,3
     DO ix=-4,5
      dmbxdy(ix,iy,iz)=(-mbx(ix,iy+2,iz)+8.0_num*mbx(ix,iy+1,iz)-8.0_num*mbx(ix,iy-1,iz)+mbx(ix,iy-2,iz))*odg(2)*oneotwelve
      dmbydy(ix,iy,iz)=(-mby(ix,iy+2,iz)+8.0_num*mby(ix,iy+1,iz)-8.0_num*mby(ix,iy-1,iz)+mby(ix,iy-2,iz))*odg(2)*oneotwelve
      dmbzdy(ix,iy,iz)=(-mbz(ix,iy+2,iz)+8.0_num*mbz(ix,iy+1,iz)-8.0_num*mbz(ix,iy-1,iz)+mbz(ix,iy-2,iz))*odg(2)*oneotwelve
     END DO
    END DO
   END DO

   !ALLOCATE(dmbxdz(-2:3,-2:3,-1:2))
   !ALLOCATE(dmbydz(-2:3,-2:3,-1:2))
   !ALLOCATE(dmbzdz(-2:3,-2:3,-1:2))
   !DO iz=-1,2
   ! DO iy=-2,3   
   !  DO ix=-2,3
   !   dmbxdz(ix,iy,iz)=0.5_num*(mbx(ix,iy,iz+1)-mbx(ix,iy,iz-1))*odg(3)
   !   dmbydz(ix,iy,iz)=0.5_num*(mby(ix,iy,iz+1)-mby(ix,iy,iz-1))*odg(3)
   !   dmbzdz(ix,iy,iz)=0.5_num*(mbz(ix,iy,iz+1)-mbz(ix,iy,iz-1))*odg(3)
   !  ENDDO
   ! ENDDO
   !ENDDO
   ! .. and finally z derivs.
   ALLOCATE(dmbxdz(-4:5,-4:5,-2:3))
   ALLOCATE(dmbydz(-4:5,-4:5,-2:3))
   ALLOCATE(dmbzdz(-4:5,-4:5,-2:3))

   DO iz=-2,3
    DO iy=-4,5   
     DO ix=-4,5
      dmbxdz(ix,iy,iz)=oneotwelve*(-mbx(ix,iy,iz+2)+8.0_num*mbx(ix,iy,iz+1)-8.0_num*mbx(ix,iy,iz-1)+mbx(ix,iy,iz-2))*odg(3)
      dmbydz(ix,iy,iz)=oneotwelve*(-mby(ix,iy,iz+2)+8.0_num*mby(ix,iy,iz+1)-8.0_num*mby(ix,iy,iz-1)+mby(ix,iy,iz-2))*odg(3)
      dmbzdz(ix,iy,iz)=oneotwelve*(-mbz(ix,iy,iz+2)+8.0_num*mbz(ix,iy,iz+1)-8.0_num*mbz(ix,iy,iz-1)+mbz(ix,iy,iz-2))*odg(3)
     ENDDO
    ENDDO
   ENDDO

   
   ! NB the arrays are odd shapes depending on deriv direction but the target cell remains at [0:1,0:1,0:1]
   ! also need temporary velocity arrays, in order to calc E=-etaJ+vxB
   !ALLOCATE(mvx(-1:2,-1:2,-1:2))
   !ALLOCATE(mvy(-1:2,-1:2,-1:2))
   !ALLOCATE(mvz(-1:2,-1:2,-1:2))
   !mvx=vx(l(1)-1:l(1)+2,l(2)-1:l(2)+2,l(3)-1:l(3)+2)
   !mvy=vy(l(1)-1:l(1)+2,l(2)-1:l(2)+2,l(3)-1:l(3)+2)
   !mvz=vz(l(1)-1:l(1)+2,l(2)-1:l(2)+2,l(3)-1:l(3)+2)
   ALLOCATE(mvx(-2:3,-2:3,-2:3))
   ALLOCATE(mvy(-2:3,-2:3,-2:3))
   ALLOCATE(mvz(-2:3,-2:3,-2:3))
   mvx=vx(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   mvy=vy(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   mvz=vz(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   
   
   ! now calculate temporary currents and also local eta values
   !ALLOCATE(mjx(-1:2,-1:2,-1:2))
   !ALLOCATE(mjy(-1:2,-1:2,-1:2))
   !ALLOCATE(mjz(-1:2,-1:2,-1:2))
   !ALLOCATE(mEx(-1:2,-1:2,-1:2))
   !ALLOCATE(mEy(-1:2,-1:2,-1:2))
   !ALLOCATE(mEz(-1:2,-1:2,-1:2))
   !ALLOCATE(meta(-1:2,-1:2,-1:2))
   !DO iz=-1,2   
   ! DO iy=-1,2
   !  DO ix=-1,2
   ALLOCATE(mjx(-2:3,-2:3,-2:3))
   ALLOCATE(mjy(-2:3,-2:3,-2:3))
   ALLOCATE(mjz(-2:3,-2:3,-2:3))
   ALLOCATE(mEx(-2:3,-2:3,-2:3))
   ALLOCATE(mEy(-2:3,-2:3,-2:3))
   ALLOCATE(mEz(-2:3,-2:3,-2:3))
   ALLOCATE(meta(-2:3,-2:3,-2:3))
   meta=eta(l(1)-2:l(1)+3,l(2)-2:l(2)+3,l(3)-2:l(3)+3,l(4)+it)
   

   DO iz=-2,3  
    DO iy=-2,3
     DO ix=-2,3
      mjx(ix,iy,iz)=dmbzdy(ix,iy,iz)-dmbydz(ix,iy,iz)
      mjy(ix,iy,iz)=dmbxdz(ix,iy,iz)-dmbzdx(ix,iy,iz)
      mjz(ix,iy,iz)=dmbydx(ix,iy,iz)-dmbxdy(ix,iy,iz)
      modj=(mjx(ix,iy,iz)*mjx(ix,iy,iz)+mjy(ix,iy,iz)*mjy(ix,iy,iz)+mjz(ix,iy,iz)*mjz(ix,iy,iz))**0.5_num   
      !!OCT 2014: replace step function with smooth eta ramp using tanh profile:
      !IF (modj.gt.jcrit) THEN
      ! meta(ix,iy,iz)=eta
      !ELSE
      ! meta(ix,iy,iz)=0.0_num
      !ENDIF
      !meta(ix,iy,iz)=0.5_num*(tanh((modj-jcrit)/rwidth)+1.0_num)*eta          
     END DO
    END DO
   END DO


   ! now interpolate derivs to correct point:
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)    
   !dbxdxt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!22
   dbxdxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdx(0,0,0), dmbxdx(1,0,0), dmbxdx(0,1,0), dmbxdx(1,1,0),&
    dmbxdx(0,0,1), dmbxdx(1,0,1), dmbxdx(0,1,1), dmbxdx(1,1,1))
   !dbydxt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!23
   dbydxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydx(0,0,0), dmbydx(1,0,0), dmbydx(0,1,0), dmbydx(1,1,0),&
    dmbydx(0,0,1), dmbydx(1,0,1), dmbydx(0,1,1), dmbydx(1,1,1))
   !dbzdxt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!24
   dbzdxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdx(0,0,0), dmbzdx(1,0,0), dmbzdx(0,1,0), dmbzdx(1,1,0),&
    dmbzdx(0,0,1), dmbzdx(1,0,1), dmbzdx(0,1,1), dmbzdx(1,1,1))
   !dbxdyt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!25
   dbxdyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdy(0,0,0), dmbxdy(1,0,0), dmbxdy(0,1,0), dmbxdy(1,1,0),&
    dmbxdy(0,0,1), dmbxdy(1,0,1), dmbxdy(0,1,1), dmbxdy(1,1,1))
   !dbydyt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!26
   dbydyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydy(0,0,0), dmbydy(1,0,0), dmbydy(0,1,0), dmbydy(1,1,0),&
    dmbydy(0,0,1), dmbydy(1,0,1), dmbydy(0,1,1), dmbydy(1,1,1))
   !dbzdyt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!27
   dbzdyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdy(0,0,0), dmbzdy(1,0,0), dmbzdy(0,1,0), dmbzdy(1,1,0),&
    dmbzdy(0,0,1), dmbzdy(1,0,1), dmbzdy(0,1,1), dmbzdy(1,1,1))            
   !dbxdzt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!28
   dbxdzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdz(0,0,0), dmbxdz(1,0,0), dmbxdz(0,1,0), dmbxdz(1,1,0),&
    dmbxdz(0,0,1), dmbxdz(1,0,1), dmbxdz(0,1,1), dmbxdz(1,1,1))
   !dbydzt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!29
   dbydzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydz(0,0,0), dmbydz(1,0,0), dmbydz(0,1,0), dmbydz(1,1,0),&
    dmbydz(0,0,1), dmbydz(1,0,1), dmbydz(0,1,1), dmbydz(1,1,1))
   !dbzdzt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!30
   dbzdzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdz(0,0,0), dmbzdz(1,0,0), dmbzdz(0,1,0), dmbzdz(1,1,0),&
    dmbzdz(0,0,1), dmbzdz(1,0,1), dmbzdz(0,1,1), dmbzdz(1,1,1))
   
   DEALLOCATE(dmbxdx)
   DEALLOCATE(dmbydx)
   DEALLOCATE(dmbzdx)
   DEALLOCATE(dmbxdy)
   DEALLOCATE(dmbydy)
   DEALLOCATE(dmbzdy)
   DEALLOCATE(dmbxdz)
   DEALLOCATE(dmbydz)
   DEALLOCATE(dmbzdz)

   ! interpolate j too at this stage:
   !mjx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!7
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mjx(0,0,0), mjx(1,0,0), mjx(0,1,0), mjx(1,1,0), &
    mjx(0,0,1), mjx(1,0,1), mjx(0,1,1), mjx(1,1,1))
   jxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!8
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mjy(0,0,0), mjy(1,0,0), mjy(0,1,0), mjy(1,1,0), &
    mjy(0,0,1), mjy(1,0,1), mjy(0,1,1), mjy(1,1,1))
   jyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!9
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mjz(0,0,0), mjz(1,0,0), mjz(0,1,0), mjz(1,1,0), &
    mjz(0,0,1), mjz(1,0,1), mjz(0,1,1), mjz(1,1,1))
   jzt(it+1)=temp
   
   ! --STEP FOUR--
   ! Ohm's Law: E=vxB-etaJ OR E=-etaJ
   !DO iz=-1,2  
   ! DO iy=-1,2
   !  DO ix=-1,2
   !   mEx(ix,iy,iz)=eta*(dmbzdy(ix,iy,iz)-dmbydz(ix,iy,iz))-(mvy(ix,iy,iz)*mbz(ix,iy,iz)-mvz(ix,iy,iz)*mby(ix,iy,iz))
   !   mEy(ix,iy,iz)=eta*(dmbxdz(ix,iy,iz)-dmbzdx(ix,iy,iz))-(mvz(ix,iy,iz)*mbx(ix,iy,iz)-mvx(ix,iy,iz)*mbz(ix,iy,iz))
   !   mEz(ix,iy,iz)=eta*(dmbydx(ix,iy,iz)-dmbxdy(ix,iy,iz))-(mvx(ix,iy,iz)*mby(ix,iy,iz)-mvy(ix,iy,iz)*mbx(ix,iy,iz))
   DO iz=-2,3 
    DO iy=-2,3
     DO ix=-2,3
      mEx(ix,iy,iz)=meta(ix,iy,iz)*mjx(ix,iy,iz)-(mvy(ix,iy,iz)*mbz(ix,iy,iz)-mvz(ix,iy,iz)*mby(ix,iy,iz))
      mEy(ix,iy,iz)=meta(ix,iy,iz)*mjy(ix,iy,iz)-(mvz(ix,iy,iz)*mbx(ix,iy,iz)-mvx(ix,iy,iz)*mbz(ix,iy,iz))
      mEz(ix,iy,iz)=meta(ix,iy,iz)*mjz(ix,iy,iz)-(mvx(ix,iy,iz)*mby(ix,iy,iz)-mvy(ix,iy,iz)*mbx(ix,iy,iz))
     ENDDO
    ENDDO
   ENDDO
      
   
   DEALLOCATE(mvx)
   DEALLOCATE(mvy)
   DEALLOCATE(mvz)
   DEALLOCATE(mbx)
   DEALLOCATE(mby)
   DEALLOCATE(mbz)
   DEALLOCATE(mjx)
   DEALLOCATE(mjy)
   DEALLOCATE(mjz)    
   DEALLOCATE(meta)
   
   ! --STEP FIVE-- !
   ! finally, calculate derivatives of electric field:
   
   !ALLOCATE(dmExdx(0:1,-1:2,-1:2))
   !ALLOCATE(dmEydx(0:1,-1:2,-1:2))
   !ALLOCATE(dmEzdx(0:1,-1:2,-1:2))
   !DO iz=-1,2   
   ! DO iy=-1,2
   !  DO ix=0,1
      !dmExdx(ix,iy,iz)=0.5_num*(mEx(ix+1,iy,iz)-mEx(ix-1,iy,iz))*odg(1)
      !dmEydx(ix,iy,iz)=0.5_num*(mEy(ix+1,iy,iz)-mEy(ix-1,iy,iz))*odg(1)
      !dmEzdx(ix,iy,iz)=0.5_num*(mEz(ix+1,iy,iz)-mEz(ix-1,iy,iz))*odg(1)
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
   
 !  ALLOCATE(dmExdy(-1:2,0:1,-1:2))
 !  ALLOCATE(dmEydy(-1:2,0:1,-1:2))
 !  ALLOCATE(dmEzdy(-1:2,0:1,-1:2))  
 !  DO iz=-1,2   
 !   DO iy=0,1
 !    DO ix=-1,2
 !     dmExdy(ix,iy,iz)=0.5_num*(mEx(ix,iy+1,iz)-mEx(ix,iy-1,iz))*odg(2)
 !     dmEydy(ix,iy,iz)=0.5_num*(mEy(ix,iy+1,iz)-mEy(ix,iy-1,iz))*odg(2)
 !     dmEzdy(ix,iy,iz)=0.5_num*(mEz(ix,iy+1,iz)-mEz(ix,iy-1,iz))*odg(2)
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
   
   !ALLOCATE(dmExdz(-1:2,-1:2,0:1))
   !ALLOCATE(dmEydz(-1:2,-1:2,0:1))
   !ALLOCATE(dmEzdz(-1:2,-1:2,0:1))   
   !DO iz=0,1   
   ! DO iy=-1,2
   !  DO ix=-1,2
   !   dmExdz(ix,iy,iz)=0.5_num*(mEx(ix,iy,iz+1)-mEx(ix,iy,iz-1))*odg(3)
   !   dmEydz(ix,iy,iz)=0.5_num*(mEy(ix,iy,iz+1)-mEy(ix,iy,iz-1))*odg(3)
   !   dmEzdz(ix,iy,iz)=0.5_num*(mEz(ix,iy,iz+1)-mEz(ix,iy,iz-1))*odg(3)
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
   
   !dmExdx(1:2,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
         
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!10
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdx(0,0,0), dmExdx(1,0,0), dmExdx(0,1,0), dmExdx(1,1,0),&
    dmExdx(0,0,1), dmExdx(1,0,1), dmExdx(0,1,1), dmExdx(1,1,1))
   dExdxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!11
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydx(0,0,0), dmEydx(1,0,0), dmEydx(0,1,0), dmEydx(1,1,0),&
    dmEydx(0,0,1), dmEydx(1,0,1), dmEydx(0,1,1), dmEydx(1,1,1))
   dEydxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!12
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdx(0,0,0), dmEzdx(1,0,0), dmEzdx(0,1,0), dmEzdx(1,1,0),&
    dmEzdx(0,0,1), dmEzdx(1,0,1), dmEzdx(0,1,1), dmEzdx(1,1,1))
   dEzdxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!13
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdy(0,0,0), dmExdy(1,0,0), dmExdy(0,1,0), dmExdy(1,1,0),&
    dmExdy(0,0,1), dmExdy(1,0,1), dmExdy(0,1,1), dmExdy(1,1,1))
   dExdyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!14
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydy(0,0,0), dmEydy(1,0,0), dmEydy(0,1,0), dmEydy(1,1,0),&
    dmEydy(0,0,1), dmEydy(1,0,1), dmEydy(0,1,1), dmEydy(1,1,1))
   dEydyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!15
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdy(0,0,0), dmEzdy(1,0,0), dmEzdy(0,1,0), dmEzdy(1,1,0),&
    dmEzdy(0,0,1), dmEzdy(1,0,1), dmEzdy(0,1,1), dmEzdy(1,1,1))
   dEzdyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!16
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdz(0,0,0), dmExdz(1,0,0), dmExdz(0,1,0), dmExdz(1,1,0),&
    dmExdz(0,0,1), dmExdz(1,0,1), dmExdz(0,1,1), dmExdz(1,1,1))
   dExdzt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!17
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydz(0,0,0), dmEydz(1,0,0), dmEydz(0,1,0), dmEydz(1,1,0),&
    dmEydz(0,0,1), dmEydz(1,0,1), dmEydz(0,1,1), dmEydz(1,1,1))
   dEydzt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!18
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdz(0,0,0), dmEzdz(1,0,0), dmEzdz(0,1,0), dmEzdz(1,1,0),&
    dmEzdz(0,0,1), dmEzdz(1,0,1), dmEzdz(0,1,1), dmEzdz(1,1,1))
   dEzdzt(it+1)=temp
   
   DEALLOCATE(dmExdx)
   DEALLOCATE(dmEydx)
   DEALLOCATE(dmEzdx)
   DEALLOCATE(dmExdy)
   DEALLOCATE(dmEydy)
   DEALLOCATE(dmEzdy)
   DEALLOCATE(dmExdz)
   DEALLOCATE(dmEydz)
   DEALLOCATE(dmEzdz)
   
   !mEx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
         
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!19
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEx(0,0,0), mEx(1,0,0), mEx(0,1,0), mEx(1,1,0),&
    mEx(0,0,1), mEx(1,0,1), mEx(0,1,1), mEx(1,1,1))
   Ext(it+1)=temp 
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!20
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEy(0,0,0), mEy(1,0,0), mEy(0,1,0), mEy(1,1,0),&
    mEy(0,0,1), mEy(1,0,1), mEy(0,1,1), mEy(1,1,1))
   Eyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!21
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEz(0,0,0), mEz(1,0,0), mEz(0,1,0), mEz(1,1,0),&
    mEz(0,0,1), mEz(1,0,1), mEz(0,1,1), mEz(1,1,1))
   Ezt(it+1)=temp 
    
   DEALLOCATE(mEx)
   DEALLOCATE(mEy)
   DEALLOCATE(mEz)
   
  END DO
   
   !we can try and sort other variables later (e.g. density/temp) when we're happy
             
   !rhot=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!31
   !rho(l(1),l(2),l(3)),rho(l(1)+1,l(2),l(3)),rho(l(1),l(2)+1,l(3)),&
   !rho(l(1)+1,l(2)+1,l(3)),rho(l(1),l(2),l(3)+1),rho(l(1)+1,l(2),l(3)+1), &
   !rho(l(1),l(2)+1,l(3)+1),rho(l(1)+1,l(2)+1,l(3)+1))
   
   !tempt=linterp3d(R(1)-myx(l(1)),R(2)-myy(l(2)),R(3)-myz(l(3)),&	!31
   !temperature(l(1),l(2),l(3)),temperature(l(1)+1,l(2),l(3)),temperature(l(1),l(2)+1,l(3)),&
   !temperature(l(1)+1,l(2)+1,l(3)),temperature(l(1),l(2),l(3)+1),temperature(l(1)+1,l(2),l(3)+1), &
   !temperature(l(1),l(2)+1,l(3)+1),temperature(l(1)+1,l(2)+1,l(3)+1))

   ! its possible to run lare3d in real units - lets make sure things are normalised if this happens
   IF (.not. lare_norm) THEN
    bxt=bxt/bscl
    byt=byt/bscl
    bzt=bzt/bscl
    vxt=vxt/vscl
    vyt=vyt/vscl
    vzt=vzt/vscl
    Ext=Ext/Escl
    Eyt=Eyt/Escl
    Ezt=Ezt/Escl
    jxt=jxt/bscl*lscl
    jyt=jyt/bscl*lscl
    jzt=jzt/bscl*lscl	
    dbxdxt=dbxdxt/bscl*lscl
    dbydxt=dbydxt/bscl*lscl
    dbzdxt=dbzdxt/bscl*lscl
    dbxdyt=dbxdyt/bscl*lscl
    dbydyt=dbydyt/bscl*lscl
    dbzdyt=dbzdyt/bscl*lscl
    dbxdzt=dbxdzt/bscl*lscl
    dbydzt=dbydzt/bscl*lscl
    dbzdzt=dbzdzt/bscl*lscl
    dExdxt=dExdxt/Escl*lscl
    dEydxt=dEydxt/Escl*lscl
    dEzdxt=dEzdxt/Escl*lscl
    dExdyt=dExdyt/Escl*lscl
    dEydyt=dEydyt/Escl*lscl
    dEzdyt=dEzdyt/Escl*lscl
    dExdzt=dExdzt/Escl*lscl
    dEydzt=dEydzt/Escl*lscl
    dEzdzt=dEzdzt/Escl*lscl
   ENDIF


   IF (nframes.gt.1) THEN

    T3d(1)=linterp1d((T-ltimes(l(4)))*odgt,bxt(1),bxt(2))
    T3d(2)=linterp1d((T-ltimes(l(4)))*odgt,byt(1),byt(2))
    T3d(3)=linterp1d((T-ltimes(l(4)))*odgt,bzt(1),bzt(2))
    T3d(4)=linterp1d((T-ltimes(l(4)))*odgt,vxt(1),vxt(2))
    T3d(5)=linterp1d((T-ltimes(l(4)))*odgt,vyt(1),vyt(2))
    T3d(6)=linterp1d((T-ltimes(l(4)))*odgt,vzt(1),vzt(2))
    T3d(7)=linterp1d((T-ltimes(l(4)))*odgt,Ext(1),Ext(2))
    T3d(8)=linterp1d((T-ltimes(l(4)))*odgt,Eyt(1),Eyt(2))
    T3d(9)=linterp1d((T-ltimes(l(4)))*odgt,Ezt(1),Ezt(2))   
    T3d(10)=linterp1d((T-ltimes(l(4)))*odgt,jxt(1),jxt(2))
    T3d(11)=linterp1d((T-ltimes(l(4)))*odgt,jyt(1),jyt(2))
    T3d(12)=linterp1d((T-ltimes(l(4)))*odgt,jzt(1),jzt(2))
    T3d(13)=linterp1d((T-ltimes(l(4)))*odgt,dbxdxt(1),dbxdxt(2))
    T3d(14)=linterp1d((T-ltimes(l(4)))*odgt,dbydxt(1),dbydxt(2))
    T3d(15)=linterp1d((T-ltimes(l(4)))*odgt,dbzdxt(1),dbzdxt(2))
    T3d(16)=linterp1d((T-ltimes(l(4)))*odgt,dbxdyt(1),dbxdyt(2))
    T3d(17)=linterp1d((T-ltimes(l(4)))*odgt,dbydyt(1),dbydyt(2))
    T3d(18)=linterp1d((T-ltimes(l(4)))*odgt,dbzdyt(1),dbzdyt(2))
    T3d(19)=linterp1d((T-ltimes(l(4)))*odgt,dbxdzt(1),dbxdzt(2))
    T3d(20)=linterp1d((T-ltimes(l(4)))*odgt,dbydzt(1),dbydzt(2))
    T3d(21)=linterp1d((T-ltimes(l(4)))*odgt,dbzdzt(1),dbzdzt(2))
    T3d(22)=linterp1d((T-ltimes(l(4)))*odgt,dExdxt(1),dExdxt(2))
    T3d(23)=linterp1d((T-ltimes(l(4)))*odgt,dEydxt(1),dEydxt(2))
    T3d(24)=linterp1d((T-ltimes(l(4)))*odgt,dEzdxt(1),dEzdxt(2))
    T3d(25)=linterp1d((T-ltimes(l(4)))*odgt,dExdyt(1),dExdyt(2))
    T3d(26)=linterp1d((T-ltimes(l(4)))*odgt,dEydyt(1),dEydyt(2))
    T3d(27)=linterp1d((T-ltimes(l(4)))*odgt,dEzdyt(1),dEzdyt(2))
    T3d(28)=linterp1d((T-ltimes(l(4)))*odgt,dExdzt(1),dExdzt(2))
    T3d(29)=linterp1d((T-ltimes(l(4)))*odgt,dEydzt(1),dEydzt(2))
    T3d(30)=linterp1d((T-ltimes(l(4)))*odgt,dEzdzt(1),dEzdzt(2))
    T3d(31)=(bxt(2)-bxt(1))*odgt
    T3d(32)=(byt(2)-byt(1))*odgt
    T3d(33)=(bzt(2)-bzt(1))*odgt
    T3d(34)=(Ext(2)-Ext(1))*odgt	! as we have multiple snapshots, we can calculate time derivs finally!
    T3d(35)=(Eyt(2)-Eyt(1))*odgt
    T3d(36)=(Ezt(2)-Ezt(1))*odgt


    DEALLOCATE(bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt)
    DEALLOCATE(dbxdxt,dbxdyt,dbxdzt,dbydxt,dbydyt,dbydzt,dbzdxt,dbzdyt,dbzdzt)
    DEALLOCATE(dExdxt,dExdyt,dExdzt,dEydxt,dEydyt,dEydzt,dEzdxt,dEzdyt,dEzdzt)

   ELSE
   
    T3d(1)=bxt(1)
    T3d(2)=byt(1)
    T3d(3)=bzt(1)
    T3d(4)=vxt(1)
    T3d(5)=vyt(1)
    T3d(6)=vzt(1)
    T3d(7)=Ext(1)
    T3d(8)=Eyt(1)
    T3d(9)=Ezt(1)   
    T3d(10)=jxt(1)
    T3d(11)=jyt(1)
    T3d(12)=jzt(1)
    T3d(13)=dbxdxt(1)
    T3d(14)=dbydxt(1)
    T3d(15)=dbzdxt(1)
    T3d(16)=dbxdyt(1)
    T3d(17)=dbydyt(1)
    T3d(18)=dbzdyt(1)
    T3d(19)=dbxdzt(1)
    T3d(20)=dbydzt(1)
    T3d(21)=dbzdzt(1)
    T3d(22)=dExdxt(1)
    T3d(23)=dEydxt(1)
    T3d(24)=dEzdxt(1)
    T3d(25)=dExdyt(1)
    T3d(26)=dEydyt(1)
    T3d(27)=dEzdyt(1)
    T3d(28)=dExdzt(1)
    T3d(29)=dEydzt(1)
    T3d(30)=dEzdzt(1)
    T3d(31:36)=0.0_num
   ENDIF
   
   RETURN

END FUNCTION T3d
!-----------------------------------------------;
  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) :: str_in, str_test
    CHARACTER(30) :: str_trim
    LOGICAL :: str_cmp
    INTEGER :: l

    str_trim = TRIM(ADJUSTL(str_in))
    l = LEN(str_test)

    IF (l > LEN(str_in)) THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    IF (str_trim(l+1:l+1) /= ' ') THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    str_cmp = str_trim(1:l) == str_test

  END FUNCTION str_cmp
!---------------------------------------------------  
    SUBROUTINE check_dims(dims)

    INTEGER, INTENT(IN) :: dims(:)
    INTEGER, DIMENSION(c_ndims) :: global_dims
    INTEGER :: ierr

    global_dims = (/ nx_global, ny_global, nz_global /)

    IF (ALL(dims(1:c_ndims) == global_dims(1:c_ndims))) RETURN

    IF (rank == 0) THEN
      PRINT*, '*** ERROR ***'
      PRINT*, 'Number of gridpoints in restart dump does not match', &
          ' the control parameters.'
      PRINT*, 'Control grid: ', nx_global, ',', ny_global, ',', nz_global
      PRINT*, 'Restart dump grid: ', dims(1), ',', dims(2), ',', dims(3)
    END IF

    CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
    STOP

  END SUBROUTINE check_dims
  
END MODULE lare_functions
