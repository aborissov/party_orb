MODULE lare_functions
!+ contains destagger, interpolation, and sub grid calculation of variables

USE GLOBAL

 IMPLICIT NONE

  PRIVATE :: smoothstep, clamp ! Alternative to tanh, see wikipedia article
  PUBLIC :: stagger_bx, stagger_by, stagger_bz, stagger_bx_2d, stagger_by_2d
  PUBLIC :: stagger_j
  PUBLIC :: linterp3d, linterp2d, linterp1d, T2d, t3d, check_dims!,str_cmp

CONTAINS

function smoothstep(x, centre, half_width) 
  real(num), intent(in) :: x, centre, half_width
  real(num) :: scaled_x, edge0, edge1, smoothstep

  edge0 = centre - half_width
  edge1 = centre + half_width
  scaled_x = clamp((x - edge0) / (edge1 - edge0), 0.0_num, 1.0_num)
  smoothstep = scaled_x * scaled_x * (3 - 2 * scaled_x)
end function smoothstep

function clamp(x, lower_limit, upper_limit)
  real(num), intent(in) :: x, lower_limit, upper_limit
  real(num) :: clamp
  if (x < lower_limit) then
    clamp = lower_limit
    return
  else if (x > upper_limit) then
    clamp = upper_limit
    return
  else
    clamp = x
    return
  endif
end function clamp

!--------------------------------------
FUNCTION stagger_j(var)
!+ function for de-staggering any j
    REAL(num), DIMENSION(:,:), INTENT(IN)  		:: var
    REAL(num), DIMENSION(SIZE(var,1), SIZE(var,2))	:: stagger_j 
    INTEGER	:: mny, mnz
    
    mny=size(var,1)
    mnz=size(var,2)
    
    !staggering right
    stagger_j(1:mny-1,1:mnz)=0.5_num*(var(1:mny-1,1:mnz)+var(2:mny,1:mnz))
    stagger_j(mny,1:mnz)=var(mny,1:mnz)
    
    RETURN

END FUNCTION stagger_j
!-------------------------------!
FUNCTION stagger_bx(var)
!+ function for de-staggering bx
    REAL(num), DIMENSION(:,:), INTENT(IN)  		:: var
    REAL(num), DIMENSION(SIZE(var,1), SIZE(var,2))	:: stagger_bx 
    INTEGER	:: mny, mnz
    
    mny=size(var,1)
    mnz=size(var,2)
    
    !staggering right
    stagger_bx(1:mny-1,1:mnz)=0.5_num*(var(1:mny-1,1:mnz)+var(2:mny,1:mnz))
    stagger_bx(mny,1:mnz)=var(mny,1:mnz)
    !staggering up
    stagger_bx(1:mny,1:mnz-1)=0.5_num*(stagger_bx(1:mny,1:mnz-1)+stagger_bx(1:mny,2:mnz))
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
      PRINT*, 'CRITICAL ERROR: dx, dy or dz is TOO BIG (in linterp3d: lare_functions_mod)'
      IF (abs(dx).gt.1.0_num) PRINT*, 'dx:', abs(dx)
      IF (abs(dy).gt.1.0_num) PRINT*, 'dy:', abs(dy)      
      IF (abs(dz).gt.1.0_num) PRINT*, 'dz:', abs(dz)
      Rlost=.TRUE.
      linterp3d=0.0_num
      RETURN
      !STOP
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
      PRINT*, 'CRITICAL ERROR: dx, dy is TOO BIG (in linterp2d: lare_functions_mod)'      
      IF (abs(dx).gt.1.0_num) PRINT*, 'dx:', abs(dx)
      IF (abs(dy).gt.1.0_num) PRINT*, 'dy:', abs(dy)   
      Rlost=.TRUE.
      linterp2d=0.0_num
      RETURN
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
      PRINT*, 'CRITICAL ERROR: dx is TOO BIG'
      print*, abs(dx)
      Rlost=.TRUE.
      linterp1d=0.0_num
      RETURN
    ENDIF
  
  linterp1d=(1.0_num-dx)*f0+dx*f1

END FUNCTION linterp1d
!-------------------------------------------------------------------------------------
FUNCTION T2d(R,T)
! function to calculate local quantities at the particle position on lare2d grid
! Dec 2015 - 2d version of T3d

   REAL(num), DIMENSION(3), INTENT(IN)		:: R		!actual position
   REAL(num), INTENT(IN)			:: T
   REAL(num), DIMENSION(36)			:: T2d
   REAL(num)					:: temp,  modj
   REAL(num), DIMENSION(2)			:: grid_spacing, grid_spacing_inverse, coffset
   REAL(num), DIMENSION(:), ALLOCATABLE		:: dgt, odgt
   REAL(num), DIMENSION(2)		:: bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt
   REAL(num), DIMENSION(2)		:: dbxdxt,dbxdyt,dbydxt,dbydyt,dbzdxt,dbzdyt
   REAL(num), DIMENSION(2)		:: dExdxt,dExdyt,dEydxt,dEydyt,dEzdxt,dEzdyt
   REAL(num), DIMENSION(-4:5, -4:5)	:: mbx, mby, mbz 
   REAL(num), DIMENSION(-2:3, -2:3)	:: mEx, mEy, mEz 
   REAL(num), DIMENSION(-2:3, -2:3)	:: mjx, mjy, mjz 
   REAL(num), DIMENSION(-2:3, -2:3)	:: mvx, mvy, mvz
   REAL(num), DIMENSION(-2:3, -4:5)	:: dmbxdx,dmbydx,dmbzdx
   REAL(num), DIMENSION(-4:5, -2:3)	:: dmbxdy,dmbydy,dmbzdy
   REAL(num), DIMENSION(0:1, -2:3)	:: dmExdx,dmEydx,dmEzdx
   REAL(num), DIMENSION(-2:3, 0:1)	:: dmExdy,dmEydy,dmEzdy
   REAL(num), DIMENSION(-2:3, -2:3)	:: meta
   INTEGER, DIMENSION(3) 			:: part_grid_index		! set at value it could never reach!
   INTEGER, DIMENSION(3) 			:: part_grid_index_check		! ALEXEI: debugging!
   INTEGER					:: jjx, jjy, jjz,jjt, rpt
   LOGICAL					:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE.

   if (isnan(R(1)) .or. isnan(R(2)) .or. isnan(R(3))) then
     print *, "R nan ", R
     call abort
   endif
 
    temp=0.0_num
    grid_spacing=(/myx(2)-myx(1),myy(2)-myy(1)/)	! grid spacing
    grid_spacing_inverse=(/1.0_num/grid_spacing(1),1.0_num/grid_spacing(2)/)	! one over grid spacing
    part_grid_index=(/-nx,-ny,-nframes/)					! initial value of part_grid_index, set to silly value (as <=0 triggers flag)   
    part_grid_index_check=(/-nx,-ny,-nframes/)					! ALEXEI: debugging!   
    
       
  ! --STEP ONE-- !
  ! first, locate the (x,y,z) position on the grid
  ! particle found between (jjx,jjy,jjz) and (jjx+1,jjy+1,jjz+1)
  ! need four gridpoints either side to guarantee we can create derivs successfully   
     
  ! ALEXEI: why are we starting with 4 here? Do more testing to make sure this
  ! is correct
  ! ALEXEI: This assumes uniform grid! make sure to add check for this.
  part_grid_index(1) = floor((R(1) - myx(4))*grid_spacing_inverse(1)) + 4
  part_grid_index(2) = floor((R(2) - myy(4))*grid_spacing_inverse(2)) + 4
  
   !DO jjx=4,nx-4
   ! IF ((R(1).ge.myx(jjx)).and.(R(1).lt.myx(jjx+1))) THEN 
   !   part_grid_index(1)=jjx
   !   !part_grid_index_check(1)=jjx
   !   EXIT
   ! !ELSE
   ! !  PRINT *, 'CANNOT FIND R(1) IN LARE x RANGE'
   ! ENDIF
   !ENDDO
   !DO jjy=4,ny-4
   ! IF ((R(2).ge.myy(jjy)).and.(R(2).lt.myy(jjy+1))) THEN
   !   part_grid_index(2)=jjy
   !   !part_grid_index_check(2)=jjy
   !   EXIT
   ! !ELSE
   ! !  PRINT *, 'CANNOT FIND R(2) IN LARE y RANGE'
   ! ENDIF
   !ENDDO

  !if (any(part_grid_index(:) .ne. part_grid_index_check(:))) then
  !  print *, "calculated and brute force indices don't match."
  !  print *, "R ", R(:)
  !  print *, "min max x ", myx(4), myx(nx-4)
  !  print *, "min max y ", myy(4), myy(ny-4)
  !  print *, "calculated  ", part_grid_index(:)
  !  print *, "brute force ", part_grid_index_check(:)
  !  call abort
  !endif


! No guarantee we have more than one frame. IF we have one, this routine doesn't bother interpolating in time
  IF (nframes.gt.1) THEN
   ALLOCATE(dgt(nframes-1),odgt(nframes-1))
   dgt=ltimes(2:nframes)-ltimes(1:nframes-1)
   odgt=1.0_num/dgt
   !print*, ltimes(1), ltimes(0)
   DO jjt=1,nframes
    IF ((T.GE.ltimes(jjt)).AND.((T.LT.ltimes(jjt+1)))) THEN
      part_grid_index(3)=jjt
      EXIT
    !ELSE
    !  PRINT *, 'CANNOT FIND TIME IN LARE TIME RANGE'
    ENDIF
   ENDDO
   rpt=1
  ELSE
   part_grid_index(3)=1
   rpt=0
  ENDIF 

  IF (minval(part_grid_index).le.0) THEN
  ! CHECK: have we located the particle on the grid?
   PRINT*, 'particle identified at', L
   PRINT*, 'x/y grid bounds-> 1:', nx, ny
   PRINT*, 't grid bounds->', ltimes(1),':',ltimes(nframes)
   print*, 'position ', R
   call abort
   STOP
  ENDIF

   coffset=(/(R(1)-myx(part_grid_index(1)))*grid_spacing_inverse(1),(R(2)-myy(part_grid_index(2)))*grid_spacing_inverse(2)/)
   
   ! --STEP TWO-- !
   ! begin interpolation of basic variables we already have, e.g. bx, by, bz, vx, vy, vz
       
  DO it=0,rpt	! NEED TO REPEAT FOR INTERPOLATION BETWEEN FRAMES, JT DEC 2015
   temp=linterp2d(coffset(1),coffset(2), &
   bx(part_grid_index(1),part_grid_index(2),1,part_grid_index(3)+it),bx(part_grid_index(1)+1,part_grid_index(2),1,part_grid_index(3)+it),bx(part_grid_index(1),part_grid_index(2)+1,1,part_grid_index(3)+it),bx(part_grid_index(1)+1,part_grid_index(2)+1,1,part_grid_index(3)+it))
   bxt(it+1)=temp		!1
   
   temp=linterp2d(coffset(1),coffset(2), &
   by(part_grid_index(1),part_grid_index(2),1,part_grid_index(3)+it),by(part_grid_index(1)+1,part_grid_index(2),1,part_grid_index(3)+it),by(part_grid_index(1),part_grid_index(2)+1,1,part_grid_index(3)+it),by(part_grid_index(1)+1,part_grid_index(2)+1,1,part_grid_index(3)+it))
   byt(it+1)=temp		!2
   
   temp=linterp2d(coffset(1),coffset(2), &
   bz(part_grid_index(1),part_grid_index(2),1,part_grid_index(3)+it),bz(part_grid_index(1)+1,part_grid_index(2),1,part_grid_index(3)+it),bz(part_grid_index(1),part_grid_index(2)+1,1,part_grid_index(3)+it),bz(part_grid_index(1)+1,part_grid_index(2)+1,1,part_grid_index(3)+it))
   bzt(it+1)=temp		!3
 !  print*, 'it', it!, 'bzt(it)=', linterp2d(coffset(1),coffset(2), &
   !bz(part_grid_index(1),part_grid_index(2),0,part_grid_index(3)+it),bz(part_grid_index(1)+1,part_grid_index(2),0,part_grid_index(3)+it),bz(part_grid_index(1),part_grid_index(2)+1,0,part_grid_index(3)+it),bz(part_grid_index(1)+1,part_grid_index(2)+1,0,part_grid_index(3)+it))
 !  print*, bz(part_grid_index(1),part_grid_index(2),0,part_grid_index(3)+it), bz(part_grid_index(1)+1,part_grid_index(2),0,part_grid_index(3)+it)
 !  print*, bz(part_grid_index(1),part_grid_index(2)+1,0,part_grid_index(3)+it),bz(part_grid_index(1)+1,part_grid_index(2)+1,0,part_grid_index(3)+it)
   
   temp=linterp2d(coffset(1),coffset(2), &
   vx(part_grid_index(1),part_grid_index(2),1,part_grid_index(3)+it),vx(part_grid_index(1)+1,part_grid_index(2),1,part_grid_index(3)+it),vx(part_grid_index(1),part_grid_index(2)+1,1,part_grid_index(3)+it),vx(part_grid_index(1)+1,part_grid_index(2)+1,1,part_grid_index(3)+it))
   vxt(it+1)=temp		!4
   
   temp=linterp2d(coffset(1),coffset(2), &
   vy(part_grid_index(1),part_grid_index(2),1,part_grid_index(3)+it),vy(part_grid_index(1)+1,part_grid_index(2),1,part_grid_index(3)+it),vy(part_grid_index(1),part_grid_index(2)+1,1,part_grid_index(3)+it),vy(part_grid_index(1)+1,part_grid_index(2)+1,1,part_grid_index(3)+it))
   vyt(it+1)=temp		!5
   
   temp=linterp2d(coffset(1),coffset(2), &
   vz(part_grid_index(1),part_grid_index(2),1,part_grid_index(3)+it),vz(part_grid_index(1)+1,part_grid_index(2),1,part_grid_index(3)+it),vz(part_grid_index(1),part_grid_index(2)+1,1,part_grid_index(3)+it),vz(part_grid_index(1)+1,part_grid_index(2)+1,1,part_grid_index(3)+it))
   vzt(it+1)=temp		!6

   
   ! --STEP THREE--!
   ! create temporary arrays around target cell and calculate derivs using extra cells
   ! - how many cells depend on which finite difference routine used
   ! - need db/dr->j and etadj/dr for E-derivs
   ! JT OCT2014: commented out O(h) derivs, and insert O(h^2)
   
   ! begin with temporary B arrays:
   mbx=bx(part_grid_index(1)-4:part_grid_index(1)+5,part_grid_index(2)-4:part_grid_index(2)+5,1,part_grid_index(3)+it)
   mby=by(part_grid_index(1)-4:part_grid_index(1)+5,part_grid_index(2)-4:part_grid_index(2)+5,1,part_grid_index(3)+it)
   mbz=bz(part_grid_index(1)-4:part_grid_index(1)+5,part_grid_index(2)-4:part_grid_index(2)+5,1,part_grid_index(3)+it)


   ! calculate x-derivs first..
   DO iy=-4,5
    DO ix=-2,3
     dmbxdx(ix,iy)=(-mbx(ix+2,iy)+8.0_num*mbx(ix+1,iy)-8.0_num*mbx(ix-1,iy)+mbx(ix-2,iy))*grid_spacing_inverse(1)*oneotwelve
     dmbydx(ix,iy)=(-mby(ix+2,iy)+8.0_num*mby(ix+1,iy)-8.0_num*mby(ix-1,iy)+mby(ix-2,iy))*grid_spacing_inverse(1)*oneotwelve
     dmbzdx(ix,iy)=(-mbz(ix+2,iy)+8.0_num*mbz(ix+1,iy)-8.0_num*mbz(ix-1,iy)+mbz(ix-2,iy))*grid_spacing_inverse(1)*oneotwelve
    END DO
   END DO

   
   ! ..followed by y-derivs..
   DO iy=-2,3
    DO ix=-4,5
     dmbxdy(ix,iy)=(-mbx(ix,iy+2)+8.0_num*mbx(ix,iy+1)-8.0_num*mbx(ix,iy-1)+mbx(ix,iy-2))*grid_spacing_inverse(2)*oneotwelve
     dmbydy(ix,iy)=(-mby(ix,iy+2)+8.0_num*mby(ix,iy+1)-8.0_num*mby(ix,iy-1)+mby(ix,iy-2))*grid_spacing_inverse(2)*oneotwelve
     dmbzdy(ix,iy)=(-mbz(ix,iy+2)+8.0_num*mbz(ix,iy+1)-8.0_num*mbz(ix,iy-1)+mbz(ix,iy-2))*grid_spacing_inverse(2)*oneotwelve
    END DO
   END DO

   ! NB the arrays are odd shapes depending on deriv direction but the target cell remains at [0:1,0:1,0:1]
   ! also need temporary velocity arrays, in order to calc E=-etaJ+vxB

   mvx=vx(part_grid_index(1)-2:part_grid_index(1)+3,part_grid_index(2)-2:part_grid_index(2)+3,1,part_grid_index(3)+it)
   mvy=vy(part_grid_index(1)-2:part_grid_index(1)+3,part_grid_index(2)-2:part_grid_index(2)+3,1,part_grid_index(3)+it)
   mvz=vz(part_grid_index(1)-2:part_grid_index(1)+3,part_grid_index(2)-2:part_grid_index(2)+3,1,part_grid_index(3)+it)
   
   
   ! now calculate temporary currents and also local eta values
   DO iy=-2,3
    DO ix=-2,3
     mjx(ix,iy)=dmbzdy(ix,iy)!-dmbydz(ix,iy,iz)
     mjy(ix,iy)=-dmbzdx(ix,iy)!+dmbxdz(ix,iy,iz)
     mjz(ix,iy)=dmbydx(ix,iy)-dmbxdy(ix,iy)
     modj=sqrt(mjx(ix,iy)*mjx(ix,iy)+mjy(ix,iy)*mjy(ix,iy)+mjz(ix,iy)*mjz(ix,iy))
     meta(ix,iy)=smoothstep(modj, jcrit, rwidth)*eta+etabkg          
     ! ALEXEI: original implementation. remove when happy with smoothstep
     !meta(ix,iy)=0.5_num*(tanh((modj-jcrit)/rwidth)+1.0_num)*eta+etabkg          
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
    DO iy=-2,3
     DO ix=-2,3
      mEx(ix,iy)=meta(ix,iy)*mjx(ix,iy)-(mvy(ix,iy)*mbz(ix,iy)-mvz(ix,iy)*mby(ix,iy))
      mEy(ix,iy)=meta(ix,iy)*mjy(ix,iy)-(mvz(ix,iy)*mbx(ix,iy)-mvx(ix,iy)*mbz(ix,iy))
      mEz(ix,iy)=meta(ix,iy)*mjz(ix,iy)-(mvx(ix,iy)*mby(ix,iy)-mvy(ix,iy)*mbx(ix,iy))
     ENDDO
    ENDDO

   ! --STEP FIVE-- !
   ! finally, calculate derivatives of electric field:
   
   DO iy=-2,3
    DO ix=0,1
     dmExdx(ix,iy)=oneotwelve*(-mEx(ix+2,iy)+8.0_num*mEx(ix+1,iy)-8.0_num*mEx(ix-1,iy)+mEx(ix-2,iy))*grid_spacing_inverse(1)
     dmEydx(ix,iy)=oneotwelve*(-mEy(ix+2,iy)+8.0_num*mEy(ix+1,iy)-8.0_num*mEy(ix-1,iy)+mEy(ix-2,iy))*grid_spacing_inverse(1)
     dmEzdx(ix,iy)=oneotwelve*(-mEz(ix+2,iy)+8.0_num*mEz(ix+1,iy)-8.0_num*mEz(ix-1,iy)+mEz(ix-2,iy))*grid_spacing_inverse(1)
    ENDDO
   ENDDO

   DO iy=0,1
    DO ix=-2,3  
     dmExdy(ix,iy)=oneotwelve*(-mEx(ix,iy+2)+8.0_num*mEx(ix,iy+1)-8.0_num*mEx(ix,iy-1)+mEx(ix,iy-2))*grid_spacing_inverse(2)
     dmEydy(ix,iy)=oneotwelve*(-mEy(ix,iy+2)+8.0_num*mEy(ix,iy+1)-8.0_num*mEy(ix,iy-1)+mEy(ix,iy-2))*grid_spacing_inverse(2)
     dmEzdy(ix,iy)=oneotwelve*(-mEz(ix,iy+2)+8.0_num*mEz(ix,iy+1)-8.0_num*mEz(ix,iy-1)+mEz(ix,iy-2))*grid_spacing_inverse(2)
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

   Ext(it+1)=linterp2d(coffset(1),coffset(2), &
    mEx(0,0), mEx(1,0), mEx(0,1), mEx(1,1))
   Eyt(it+1)=linterp2d(coffset(1),coffset(2), &
    mEy(0,0), mEy(1,0), mEy(0,1), mEy(1,1))
   Ezt(it+1)=linterp2d(coffset(1),coffset(2), &
    mEz(0,0), mEz(1,0), mEz(0,1), mEz(1,1))
    
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

    T2d(1)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),bxt(1),bxt(2))
    T2d(2)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),byt(1),byt(2))
    T2d(3)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),bzt(1),bzt(2))
    T2d(4)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),vxt(1),vxt(2))
    T2d(5)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),vyt(1),vyt(2))
    T2d(6)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),vzt(1),vzt(2))
    T2d(7)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),Ext(1),Ext(2))
    T2d(8)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),Eyt(1),Eyt(2))
    T2d(9)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),Ezt(1),Ezt(2))   
    T2d(10)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),jxt(1),jxt(2))
    T2d(11)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),jyt(1),jyt(2))
    T2d(12)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),jzt(1),jzt(2))
    T2d(13)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dbxdxt(1),dbxdxt(2))
    T2d(14)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dbydxt(1),dbydxt(2))
    T2d(15)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dbzdxt(1),dbzdxt(2))
    T2d(16)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dbxdyt(1),dbxdyt(2))
    T2d(17)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dbydyt(1),dbydyt(2))
    T2d(18)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dbzdyt(1),dbzdyt(2))
    T2d(19)=0.0_num		! Z derivs explicitly set to zero in 2d case
    T2d(20)=0.0_num
    T2d(21)=0.0_num
    T2d(22)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dExdxt(1),dExdxt(2))
    T2d(23)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dEydxt(1),dEydxt(2))
    T2d(24)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dEzdxt(1),dEzdxt(2))
    T2d(25)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dExdyt(1),dExdyt(2))
    T2d(26)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dEydyt(1),dEydyt(2))
    T2d(27)=linterp1d((T-ltimes(part_grid_index(3)))*odgt(part_grid_index(3)),dEzdyt(1),dEzdyt(2))
    T2d(28)=0.0_num
    T2d(29)=0.0_num
    T2d(30)=0.0_num
    T2d(31)=(bxt(2)-bxt(1))*odgt(part_grid_index(3))
    T2d(32)=(byt(2)-byt(1))*odgt(part_grid_index(3))
    T2d(33)=(bzt(2)-bzt(1))*odgt(part_grid_index(3))
    T2d(34)=(Ext(2)-Ext(1))*odgt(part_grid_index(3))	! as we have multiple snapshots, we can calculate time derivs finally!
    T2d(35)=(Eyt(2)-Eyt(1))*odgt(part_grid_index(3))
    T2d(36)=(Ezt(2)-Ezt(1))*odgt(part_grid_index(3))
    
    DEALLOCATE(dgt,odgt)
  
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
  REAL(num)					:: temp,  modj
  REAL(num), DIMENSION(3)			:: grid_spacing, grid_spacing_inverse, coffset
  REAL(num), DIMENSION(:), ALLOCATABLE		:: dgt, odgt
  REAL(num), DIMENSION(2)		:: bxt, byt, bzt,vxt, vyt, vzt, Ext, Eyt, Ezt, jxt, jyt, jzt
  REAL(num), DIMENSION(2)		:: dbxdxt,dbxdyt,dbxdzt,dbydxt,dbydyt,dbydzt,dbzdxt,dbzdyt,dbzdzt
  REAL(num), DIMENSION(2)		:: dExdxt,dExdyt,dExdzt,dEydxt,dEydyt,dEydzt,dEzdxt,dEzdyt,dEzdzt
  REAL(num), DIMENSION(-4:5, -4:5, -4:5)	:: mbx, mby, mbz
  REAL(num), DIMENSION(-2:3, -2:3, -2:3)	:: mEx, mEy, mEz
  REAL(num), DIMENSION(-2:3, -2:3, -2:3)	:: mjx, mjy, mjz
  REAL(num), DIMENSION(-2:3, -2:3, -2:3)	:: mvx, mvy, mvz
  REAL(num), DIMENSION(-2:3, -4:5, -4:5)	:: dmbxdx,dmbydx,dmbzdx
  REAL(num), DIMENSION(-4:5, -2:3, -4:5)	:: dmbxdy,dmbydy,dmbzdy
  REAL(num), DIMENSION(-4:5, -4:5, -2:3)	:: dmbxdz,dmbydz,dmbzdz
  REAL(num), DIMENSION(0:1, -2:3, -2:3)	    :: dmExdx,dmEydx,dmEzdx
  REAL(num), DIMENSION(-2:3, 0:1, -2:3)     :: dmExdy,dmEydy,dmEzdy
  REAL(num), DIMENSION(-2:3, -2:3, 0:1)     :: dmExdz,dmEydz,dmEzdz
  REAL(num), DIMENSION(-2:3, -2:3, -2:3)	:: meta
  INTEGER, DIMENSION(4) 			:: part_grid_index		! set at value it could never reach!
  INTEGER, DIMENSION(4) 			:: part_grid_index_check		! ALEXEI: debugging!
  INTEGER					:: jjx, jjy, jjz,jjt, rpt
  LOGICAL					:: fxflag=.FALSE., fyflag=.FALSE., fzflag=.FALSE.
 
  temp=0.0_num
  grid_spacing=(/myx(2)-myx(1),myy(2)-myy(1),myz(2)-myz(1)/)	! grid spacing
  grid_spacing_inverse=(/1.0_num/grid_spacing(1),1.0_num/grid_spacing(2),1.0_num/grid_spacing(3)/)	! one over grid spacing
  part_grid_index=(/-nx,-ny,-nz,-nframes/)					! initial value of part_grid_index, set to silly value (as >=0 triggers flag)   
  part_grid_index_check=(/-nx,-ny,-nz,-nframes/)					! ALEXEI: debugging!   
    
          
  ! --STEP ONE-- !
  ! first, locate the (x,y,z) position on the grid
  ! particle found between (jjx,jjy,jjz) and (jjx+1,jjy+1,jjz+1)   
     
  ! ALEXEI: why are we starting with 4 here? Do more testing to make sure this
  ! is correct
  ! ALEXEI: This assumes uniform grid! make sure to add check for this.
  part_grid_index(1) = floor((R(1) - myx(4))*grid_spacing_inverse(1)) + 4
  part_grid_index(2) = floor((R(2) - myy(4))*grid_spacing_inverse(2)) + 4
  part_grid_index(3) = floor((R(3) - myz(4))*grid_spacing_inverse(3)) + 4
  
  !DO jjx=4,nx-4
  ! IF ((R(1).ge.myx(jjx)).and.(R(1).lt.myx(jjx+1))) THEN 
  !   part_grid_index_check(1)=jjx
  !   EXIT
  ! ENDIF
  !ENDDO
  !!IF (part_grid_index(1).eq.(-nx)) fxflag=.TRUE.
  !DO jjy=4,ny-4
  ! IF ((R(2).ge.myy(jjy)).and.(R(2).lt.myy(jjy+1))) THEN
  !   part_grid_index_check(2)=jjy
  !   EXIT
  ! ENDIF
  !ENDDO
  !!IF (part_grid_index(2).eq.(-ny)) fyflag=.TRUE.
  !DO jjz=4,nz-4
  ! IF ((R(3).ge.myz(jjz)).and.(R(3).lt.myz(jjz+1))) THEN
  !   part_grid_index_check(3)=jjz
  !   EXIT
  ! ENDIF
  !ENDDO

  !if (any(part_grid_index(:) .ne. part_grid_index_check(:))) then
  !  print *, "calculated and brute force indices don't match."
  !  print *, "R ", R(:)
  !  print *, "min max x ", myx(4), myx(nx-4)
  !  print *, "min max y ", myy(4), myy(ny-4)
  !  print *, "min max z ", myz(4), myz(nz-4)
  !  print *, "calculated  ", part_grid_index(:)
  !  print *, "brute force ", part_grid_index_check(:)
  !  call abort
  !endif

  IF ((part_grid_index(1).eq.(-nx)).or.(part_grid_index(2).eq.(-ny)).or.(part_grid_index(3).eq.(-nz))) THEN 
   !LETS TRY TO CATCH THIS BEFORE NOW
   !print*, R
   T3d(1:36)=-999.0_num
  ENDIF

! No guarantee we have more than one frame. IF we have one, this routine doesn't bother interpolating in time
  IF (nframes.gt.1) THEN
   ALLOCATE(dgt(nframes-1),odgt(nframes-1))
   dgt=ltimes(2:nframes)-ltimes(1:nframes-1)
   odgt=1.0_num/dgt
   DO jjt=1,nframes
    IF ((T.GE.ltimes(jjt)).AND.((T.LT.ltimes(jjt+1)))) THEN
      part_grid_index(4)=jjt 
      !print*, (T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4))
      EXIT
    !ELSE
    !  PRINT *, 'CANNOT FIND TIME IN LARE TIME RANGE'
    ENDIF
   ENDDO
   rpt=1
  ELSE
   part_grid_index(4)=1
   rpt=0
  !what happens if there is one frame to work with?
  ENDIF 


   coffset=(/(R(1)-myx(part_grid_index(1)))*grid_spacing_inverse(1),(R(2)-myy(part_grid_index(2)))*grid_spacing_inverse(2),(R(3)-myz(part_grid_index(3)))*grid_spacing_inverse(3)/)

   
   ! --STEP TWO-- !
   ! begin interpolation of basic variables we already have, e.g. bx, by, bz, vx, vy, vz
       
  DO it=0,rpt	! NEED TO REPEAT FOR INTERPOLATION BETWEEN FRAMES, JT DEC 2015
  
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!1   
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   bx(part_grid_index(1),part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),bx(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),bx(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),&
   bx(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),bx(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it),bx(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it), &
   bx(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it),bx(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it))
   bxt(it+1)=temp
   
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!2
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   by(part_grid_index(1),part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),by(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),by(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),&
   by(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),by(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it),by(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it), &
   by(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it),by(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it))
   byt(it+1)=temp

   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!3
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   bz(part_grid_index(1),part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),bz(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),bz(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),&
   bz(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),bz(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it),bz(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it), &
   bz(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it),bz(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it))
   bzt(it+1)=temp

   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!4
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   vx(part_grid_index(1),part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),vx(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),vx(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),&
   vx(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),vx(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it),vx(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it), &
   vx(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it),vx(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it))
   vxt(it+1)=temp

   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!5
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   vy(part_grid_index(1),part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),vy(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),vy(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),&
   vy(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),vy(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it),vy(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it), &
   vy(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it),vy(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it))
   vyt(it+1)=temp
   
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!6
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
   vz(part_grid_index(1),part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),vz(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3),part_grid_index(4)+it),vz(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),&
   vz(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3),part_grid_index(4)+it),vz(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it),vz(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1,part_grid_index(4)+it), &
   vz(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it),vz(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1,part_grid_index(4)+it))
   vzt(it+1)=temp
   
   ! --STEP THREE--!
   ! create temporary arrays around target cell and calculate derivs using extra cells
   ! - how many cells depend on which finite difference routine used
   ! - need db/dr->j and etadj/dr for E-derivs
   ! JT OCT2014: commented out O(h) derivs, and insert O(h^2)
   
   ! begin with temporary B arrays:
   mbx=bx(part_grid_index(1)-4:part_grid_index(1)+5,part_grid_index(2)-4:part_grid_index(2)+5,part_grid_index(3)-4:part_grid_index(3)+5,part_grid_index(4)+it)
   mby=by(part_grid_index(1)-4:part_grid_index(1)+5,part_grid_index(2)-4:part_grid_index(2)+5,part_grid_index(3)-4:part_grid_index(3)+5,part_grid_index(4)+it)
   mbz=bz(part_grid_index(1)-4:part_grid_index(1)+5,part_grid_index(2)-4:part_grid_index(2)+5,part_grid_index(3)-4:part_grid_index(3)+5,part_grid_index(4)+it)

   ! calculate x-derivs first..
   DO iz=-4,5   
    DO iy=-4,5
     DO ix=-2,3
      dmbxdx(ix,iy,iz)=(-mbx(ix+2,iy,iz)+8.0_num*mbx(ix+1,iy,iz)-8.0_num*mbx(ix-1,iy,iz)+mbx(ix-2,iy,iz))*grid_spacing_inverse(1)*oneotwelve
      dmbydx(ix,iy,iz)=(-mby(ix+2,iy,iz)+8.0_num*mby(ix+1,iy,iz)-8.0_num*mby(ix-1,iy,iz)+mby(ix-2,iy,iz))*grid_spacing_inverse(1)*oneotwelve
      dmbzdx(ix,iy,iz)=(-mbz(ix+2,iy,iz)+8.0_num*mbz(ix+1,iy,iz)-8.0_num*mbz(ix-1,iy,iz)+mbz(ix-2,iy,iz))*grid_spacing_inverse(1)*oneotwelve
     END DO
    END DO
   END DO
   
   ! ..followed by y-derivs..
   DO iz=-4,5   
    DO iy=-2,3
     DO ix=-4,5
      dmbxdy(ix,iy,iz)=(-mbx(ix,iy+2,iz)+8.0_num*mbx(ix,iy+1,iz)-8.0_num*mbx(ix,iy-1,iz)+mbx(ix,iy-2,iz))*grid_spacing_inverse(2)*oneotwelve
      dmbydy(ix,iy,iz)=(-mby(ix,iy+2,iz)+8.0_num*mby(ix,iy+1,iz)-8.0_num*mby(ix,iy-1,iz)+mby(ix,iy-2,iz))*grid_spacing_inverse(2)*oneotwelve
      dmbzdy(ix,iy,iz)=(-mbz(ix,iy+2,iz)+8.0_num*mbz(ix,iy+1,iz)-8.0_num*mbz(ix,iy-1,iz)+mbz(ix,iy-2,iz))*grid_spacing_inverse(2)*oneotwelve
     END DO
    END DO
   END DO

   ! .. and finally z derivs.
   DO iz=-2,3
    DO iy=-4,5   
     DO ix=-4,5
      dmbxdz(ix,iy,iz)=oneotwelve*(-mbx(ix,iy,iz+2)+8.0_num*mbx(ix,iy,iz+1)-8.0_num*mbx(ix,iy,iz-1)+mbx(ix,iy,iz-2))*grid_spacing_inverse(3)
      dmbydz(ix,iy,iz)=oneotwelve*(-mby(ix,iy,iz+2)+8.0_num*mby(ix,iy,iz+1)-8.0_num*mby(ix,iy,iz-1)+mby(ix,iy,iz-2))*grid_spacing_inverse(3)
      dmbzdz(ix,iy,iz)=oneotwelve*(-mbz(ix,iy,iz+2)+8.0_num*mbz(ix,iy,iz+1)-8.0_num*mbz(ix,iy,iz-1)+mbz(ix,iy,iz-2))*grid_spacing_inverse(3)
     ENDDO
    ENDDO
   ENDDO

   
   ! NB the arrays are odd shapes depending on deriv direction but the target cell remains at [0:1,0:1,0:1]
   ! also need temporary velocity arrays, in order to calc E=-etaJ+vxB
   mvx=vx(part_grid_index(1)-2:part_grid_index(1)+3,part_grid_index(2)-2:part_grid_index(2)+3,part_grid_index(3)-2:part_grid_index(3)+3,part_grid_index(4)+it)
   mvy=vy(part_grid_index(1)-2:part_grid_index(1)+3,part_grid_index(2)-2:part_grid_index(2)+3,part_grid_index(3)-2:part_grid_index(3)+3,part_grid_index(4)+it)
   mvz=vz(part_grid_index(1)-2:part_grid_index(1)+3,part_grid_index(2)-2:part_grid_index(2)+3,part_grid_index(3)-2:part_grid_index(3)+3,part_grid_index(4)+it)
   
   
   ! now calculate temporary currents and also local eta values
   DO iz=-2,3  
    DO iy=-2,3
     DO ix=-2,3
      mjx(ix,iy,iz)=dmbzdy(ix,iy,iz)-dmbydz(ix,iy,iz)
      mjy(ix,iy,iz)=dmbxdz(ix,iy,iz)-dmbzdx(ix,iy,iz)
      mjz(ix,iy,iz)=dmbydx(ix,iy,iz)-dmbxdy(ix,iy,iz)
      modj=sqrt(mjx(ix,iy,iz)*mjx(ix,iy,iz)+mjy(ix,iy,iz)*mjy(ix,iy,iz)+mjz(ix,iy,iz)*mjz(ix,iy,iz))
      meta(ix,iy,iz)=smoothstep(modj, jcrit, rwidth)*eta+etabkg          
      ! ALEXEI: original implementation, remove when happy with smoothstep
      !meta(ix,iy,iz)=0.5_num*(tanh((modj-jcrit)/rwidth)+1.0_num)*eta+etabkg          
     END DO
    END DO
   END DO


   ! now interpolate derivs to correct point:
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)    
   !dbxdxt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!22
   dbxdxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdx(0,0,0), dmbxdx(1,0,0), dmbxdx(0,1,0), dmbxdx(1,1,0),&
    dmbxdx(0,0,1), dmbxdx(1,0,1), dmbxdx(0,1,1), dmbxdx(1,1,1))
   !dbydxt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!23
   dbydxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydx(0,0,0), dmbydx(1,0,0), dmbydx(0,1,0), dmbydx(1,1,0),&
    dmbydx(0,0,1), dmbydx(1,0,1), dmbydx(0,1,1), dmbydx(1,1,1))
   !dbzdxt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!24
   dbzdxt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdx(0,0,0), dmbzdx(1,0,0), dmbzdx(0,1,0), dmbzdx(1,1,0),&
    dmbzdx(0,0,1), dmbzdx(1,0,1), dmbzdx(0,1,1), dmbzdx(1,1,1))
   !dbxdyt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!25
   dbxdyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdy(0,0,0), dmbxdy(1,0,0), dmbxdy(0,1,0), dmbxdy(1,1,0),&
    dmbxdy(0,0,1), dmbxdy(1,0,1), dmbxdy(0,1,1), dmbxdy(1,1,1))
   !dbydyt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!26
   dbydyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydy(0,0,0), dmbydy(1,0,0), dmbydy(0,1,0), dmbydy(1,1,0),&
    dmbydy(0,0,1), dmbydy(1,0,1), dmbydy(0,1,1), dmbydy(1,1,1))
   !dbzdyt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!27
   dbzdyt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdy(0,0,0), dmbzdy(1,0,0), dmbzdy(0,1,0), dmbzdy(1,1,0),&
    dmbzdy(0,0,1), dmbzdy(1,0,1), dmbzdy(0,1,1), dmbzdy(1,1,1))            
   !dbxdzt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!28
   dbxdzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbxdz(0,0,0), dmbxdz(1,0,0), dmbxdz(0,1,0), dmbxdz(1,1,0),&
    dmbxdz(0,0,1), dmbxdz(1,0,1), dmbxdz(0,1,1), dmbxdz(1,1,1))
   !dbydzt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!29
   dbydzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbydz(0,0,0), dmbydz(1,0,0), dmbydz(0,1,0), dmbydz(1,1,0),&
    dmbydz(0,0,1), dmbydz(1,0,1), dmbydz(0,1,1), dmbydz(1,1,1))
   !dbzdzt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!30
   dbzdzt(it+1)=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmbzdz(0,0,0), dmbzdz(1,0,0), dmbzdz(0,1,0), dmbzdz(1,1,0),&
    dmbzdz(0,0,1), dmbzdz(1,0,1), dmbzdz(0,1,1), dmbzdz(1,1,1))
   
   ! interpolate j too at this stage:
   !mjx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!7
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mjx(0,0,0), mjx(1,0,0), mjx(0,1,0), mjx(1,1,0), &
    mjx(0,0,1), mjx(1,0,1), mjx(0,1,1), mjx(1,1,1))
   jxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!8
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mjy(0,0,0), mjy(1,0,0), mjy(0,1,0), mjy(1,1,0), &
    mjy(0,0,1), mjy(1,0,1), mjy(0,1,1), mjy(1,1,1))
   jyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!9
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
      
   
   ! --STEP FIVE-- !
   ! finally, calculate derivatives of electric field:
   
   DO iz=-2,3   
    DO iy=-2,3
     DO ix=0,1
      dmExdx(ix,iy,iz)=oneotwelve*(-mEx(ix+2,iy,iz)+8.0_num*mEx(ix+1,iy,iz)-8.0_num*mEx(ix-1,iy,iz)+mEx(ix-2,iy,iz))*grid_spacing_inverse(1)
      dmEydx(ix,iy,iz)=oneotwelve*(-mEy(ix+2,iy,iz)+8.0_num*mEy(ix+1,iy,iz)-8.0_num*mEy(ix-1,iy,iz)+mEy(ix-2,iy,iz))*grid_spacing_inverse(1)
      dmEzdx(ix,iy,iz)=oneotwelve*(-mEz(ix+2,iy,iz)+8.0_num*mEz(ix+1,iy,iz)-8.0_num*mEz(ix-1,iy,iz)+mEz(ix-2,iy,iz))*grid_spacing_inverse(1)
     ENDDO
    ENDDO
   ENDDO
   
   DO iz=-2,3   
    DO iy=0,1
     DO ix=-2,3  
      dmExdy(ix,iy,iz)=oneotwelve*(-mEx(ix,iy+2,iz)+8.0_num*mEx(ix,iy+1,iz)-8.0_num*mEx(ix,iy-1,iz)+mEx(ix,iy-2,iz))*grid_spacing_inverse(2)
      dmEydy(ix,iy,iz)=oneotwelve*(-mEy(ix,iy+2,iz)+8.0_num*mEy(ix,iy+1,iz)-8.0_num*mEy(ix,iy-1,iz)+mEy(ix,iy-2,iz))*grid_spacing_inverse(2)
      dmEzdy(ix,iy,iz)=oneotwelve*(-mEz(ix,iy+2,iz)+8.0_num*mEz(ix,iy+1,iz)-8.0_num*mEz(ix,iy-1,iz)+mEz(ix,iy-2,iz))*grid_spacing_inverse(2)
     ENDDO
    ENDDO
   ENDDO   
   
   DO iz=0,1   
    DO iy=-2,3
     DO ix=-2,3
      dmExdz(ix,iy,iz)=oneotwelve*(-mEx(ix,iy,iz+2)+8.0_num*mEx(ix,iy,iz+1)-8.0_num*mEx(ix,iy,iz-1)+mEx(ix,iy,iz-2))*grid_spacing_inverse(3)
      dmEydz(ix,iy,iz)=oneotwelve*(-mEy(ix,iy,iz+2)+8.0_num*mEy(ix,iy,iz+1)-8.0_num*mEy(ix,iy,iz-1)+mEy(ix,iy,iz-2))*grid_spacing_inverse(3)
      dmEzdz(ix,iy,iz)=oneotwelve*(-mEz(ix,iy,iz+2)+8.0_num*mEz(ix,iy,iz+1)-8.0_num*mEz(ix,iy,iz-1)+mEz(ix,iy,iz-2))*grid_spacing_inverse(3)
     ENDDO
    ENDDO
   ENDDO
   
   !dmExdx(1:2,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
         
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!10
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdx(0,0,0), dmExdx(1,0,0), dmExdx(0,1,0), dmExdx(1,1,0),&
    dmExdx(0,0,1), dmExdx(1,0,1), dmExdx(0,1,1), dmExdx(1,1,1))
   dExdxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!11
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydx(0,0,0), dmEydx(1,0,0), dmEydx(0,1,0), dmEydx(1,1,0),&
    dmEydx(0,0,1), dmEydx(1,0,1), dmEydx(0,1,1), dmEydx(1,1,1))
   dEydxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!12
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdx(0,0,0), dmEzdx(1,0,0), dmEzdx(0,1,0), dmEzdx(1,1,0),&
    dmEzdx(0,0,1), dmEzdx(1,0,1), dmEzdx(0,1,1), dmEzdx(1,1,1))
   dEzdxt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!13
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdy(0,0,0), dmExdy(1,0,0), dmExdy(0,1,0), dmExdy(1,1,0),&
    dmExdy(0,0,1), dmExdy(1,0,1), dmExdy(0,1,1), dmExdy(1,1,1))
   dExdyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!14
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydy(0,0,0), dmEydy(1,0,0), dmEydy(0,1,0), dmEydy(1,1,0),&
    dmEydy(0,0,1), dmEydy(1,0,1), dmEydy(0,1,1), dmEydy(1,1,1))
   dEydyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!15
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdy(0,0,0), dmEzdy(1,0,0), dmEzdy(0,1,0), dmEzdy(1,1,0),&
    dmEzdy(0,0,1), dmEzdy(1,0,1), dmEzdy(0,1,1), dmEzdy(1,1,1))
   dEzdyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!16
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmExdz(0,0,0), dmExdz(1,0,0), dmExdz(0,1,0), dmExdz(1,1,0),&
    dmExdz(0,0,1), dmExdz(1,0,1), dmExdz(0,1,1), dmExdz(1,1,1))
   dExdzt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!17
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEydz(0,0,0), dmEydz(1,0,0), dmEydz(0,1,0), dmEydz(1,1,0),&
    dmEydz(0,0,1), dmEydz(1,0,1), dmEydz(0,1,1), dmEydz(1,1,1))
   dEydzt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!18
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    dmEzdz(0,0,0), dmEzdz(1,0,0), dmEzdz(0,1,0), dmEzdz(1,1,0),&
    dmEzdz(0,0,1), dmEzdz(1,0,1), dmEzdz(0,1,1), dmEzdz(1,1,1))
   dEzdzt(it+1)=temp
   
   !mEx(0:3,0:3,0:3)	<-x=1-2,y=1-2,z=1-2
   !linterp3d(dx,dy,dz,f000,f100,f010,f110,f001,f101,f011,f111)
         
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!19
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEx(0,0,0), mEx(1,0,0), mEx(0,1,0), mEx(1,1,0),&
    mEx(0,0,1), mEx(1,0,1), mEx(0,1,1), mEx(1,1,1))
   Ext(it+1)=temp 
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!20
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEy(0,0,0), mEy(1,0,0), mEy(0,1,0), mEy(1,1,0),&
    mEy(0,0,1), mEy(1,0,1), mEy(0,1,1), mEy(1,1,1))
   Eyt(it+1)=temp
   !temp=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!21
   temp=linterp3d(coffset(1),coffset(2),coffset(3), &
    mEz(0,0,0), mEz(1,0,0), mEz(0,1,0), mEz(1,1,0),&
    mEz(0,0,1), mEz(1,0,1), mEz(0,1,1), mEz(1,1,1))
   Ezt(it+1)=temp 
    
  END DO
   
   !we can try and sort other variables later (e.g. density/temp) when we're happy
             
   !rhot=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!31
   !rho(part_grid_index(1),part_grid_index(2),part_grid_index(3)),rho(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)),rho(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)),&
   !rho(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)),rho(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1),rho(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1), &
   !rho(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1),rho(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1))
   
   !tempt=linterp3d(R(1)-myx(part_grid_index(1)),R(2)-myy(part_grid_index(2)),R(3)-myz(part_grid_index(3)),&	!31
   !temperature(part_grid_index(1),part_grid_index(2),part_grid_index(3)),temperature(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)),temperature(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)),&
   !temperature(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)),temperature(part_grid_index(1),part_grid_index(2),part_grid_index(3)+1),temperature(part_grid_index(1)+1,part_grid_index(2),part_grid_index(3)+1), &
   !temperature(part_grid_index(1),part_grid_index(2)+1,part_grid_index(3)+1),temperature(part_grid_index(1)+1,part_grid_index(2)+1,part_grid_index(3)+1))

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

    T3d(1)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),bxt(1),bxt(2))
    T3d(2)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),byt(1),byt(2))
    T3d(3)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),bzt(1),bzt(2))
    T3d(4)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),vxt(1),vxt(2))
    T3d(5)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),vyt(1),vyt(2))
    T3d(6)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),vzt(1),vzt(2))
    T3d(7)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),Ext(1),Ext(2))
    T3d(8)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),Eyt(1),Eyt(2))
    T3d(9)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),Ezt(1),Ezt(2))   
    T3d(10)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),jxt(1),jxt(2))
    T3d(11)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),jyt(1),jyt(2))
    T3d(12)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),jzt(1),jzt(2))
    T3d(13)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbxdxt(1),dbxdxt(2))
    T3d(14)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbydxt(1),dbydxt(2))
    T3d(15)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbzdxt(1),dbzdxt(2))
    T3d(16)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbxdyt(1),dbxdyt(2))
    T3d(17)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbydyt(1),dbydyt(2))
    T3d(18)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbzdyt(1),dbzdyt(2))
    T3d(19)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbxdzt(1),dbxdzt(2))
    T3d(20)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbydzt(1),dbydzt(2))
    T3d(21)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dbzdzt(1),dbzdzt(2))
    T3d(22)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dExdxt(1),dExdxt(2))
    T3d(23)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dEydxt(1),dEydxt(2))
    T3d(24)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dEzdxt(1),dEzdxt(2))
    T3d(25)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dExdyt(1),dExdyt(2))
    T3d(26)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dEydyt(1),dEydyt(2))
    T3d(27)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dEzdyt(1),dEzdyt(2))
    T3d(28)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dExdzt(1),dExdzt(2))
    T3d(29)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dEydzt(1),dEydzt(2))
    T3d(30)=linterp1d((T-ltimes(part_grid_index(4)))*odgt(part_grid_index(4)),dEzdzt(1),dEzdzt(2))
    T3d(31)=(bxt(2)-bxt(1))*odgt(part_grid_index(4))
    T3d(32)=(byt(2)-byt(1))*odgt(part_grid_index(4))
    T3d(33)=(bzt(2)-bzt(1))*odgt(part_grid_index(4))
    T3d(34)=(Ext(2)-Ext(1))*odgt(part_grid_index(4))	! as we have multiple snapshots, we can calculate time derivs finally!
    T3d(35)=(Eyt(2)-Eyt(1))*odgt(part_grid_index(4))
    T3d(36)=(Ezt(2)-Ezt(1))*odgt(part_grid_index(4))
    DEALLOCATE(dgt,odgt)
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
!  FUNCTION str_cmp(str_in, str_test)!
!
!    CHARACTER(*), INTENT(IN) :: str_in, str_test
!    CHARACTER(30) :: str_trim
!    LOGICAL :: str_cmp
!    INTEGER :: l!
!
!    str_trim = TRIM(ADJUSTL(str_in))
!    l = LEN(str_test)!
!
!    IF (l > LEN(str_in)) THEN
!      str_cmp = .FALSE.
!      RETURN
!    END IF
!
!    IF (str_trim(l+1:l+1) /= ' ') THEN
!      str_cmp = .FALSE.
!      RETURN
!    END IF
!
!    str_cmp = str_trim(1:l) == str_test
!
!  END FUNCTION str_cmp
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
