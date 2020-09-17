module cell_struct

use global
use m_products, only: dot, cross
use m_driverr, only: rkdrive
use m_fields, only: fields
implicit none

type, public :: cell
  ! Particle quantities
  real, dimension(:,:), allocatable :: particle_coords
  real, dimension(:),  allocatable :: particle_v_par, particle_v_perp, particle_t

  ! Grid quantities
  real, dimension(3) :: grid_upper_bound, grid_lower_bound
contains
  procedure :: process => process_cell
end type cell

contains

subroutine process_cell(self, n_part_per_cell, cell_index, nok, nbad)
  implicit none

  class(cell) :: self
  integer :: i, n_part_per_cell, particle_index, cell_index
  integer :: nok, nbad ! ALEXEI: check out what these are for.
  logical :: reset_flag
  real(num), dimension(3) :: rstart

  !$omp parallel do private(mu,ustart,gammastart,Ekin,alpha,rstart,particle_index,reset_flag)
  do i = 1, n_part_per_cell
    particle_index = i + n_part_per_cell*cell_index
    rstart = self % particle_coords(:,i)
    Ekin = sqrt(self % particle_v_par(i) * self % particle_v_par(i) + &
                self % particle_v_perp(i) * self % particle_v_perp(i))
    ! Calculate magnetic moment, initial 4-velocity, gamma
    call JTmucalc(mu,ustart,gammastart,Ekin,alpha,rstart,t1,t2, reset_flag)

    ! ALEXEI: In original version there is a check that B is not too small. Need
    ! to replicate this!!! (use reset_flag)

    ! Integrate the particle orbit
    call rkdrive(rstart,ustart,gammastart,mu,t1,t2,eps,h1, particle_index,nok,nbad)
  end do
  !$omp end parallel do
end subroutine process_cell

SUBROUTINE JTMUcalc(mu,USTART,GAMMASTART, Ekin,alpha,RSTART,T1,T2, resetflag)

  REAL(num), DIMENSION(3),INTENT(IN) :: RSTART
  REAL(num), INTENT(IN) :: T1,T2, Ekin, Alpha				! 
  REAL(num), INTENT(OUT) :: mu, gammastart, Ustart
  REAL(num), DIMENSION(3) :: B,El,a2,a3,a4,a5,a6,a7,a8,a9,a10,ue, RT
  !REAL(num) :: magB,vtot,vperp, vparstart!,Erest
  REAL(num) :: modB,vtot, gamma, Etemp
  LOGICAL, INTENT(OUT)   :: resetflag
 
 resetflag=.FALSE.
 
 !calculate B, E, V at this point/time:
 CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2)

 !calculate magnitude of B
 modB=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
 modB=sqrt(modB)
 RT=RSTART

 IF (modb.le.lowbthresh) THEN 
  resetflag=.TRUE.
 ENDIF

 !print*, 'RT=', RT
 !PRINT*, 'B=', B
 !STOP
 !print*, modB

 !Erest = (M*c*c)*oneuponAQ
 
 ! E X B drift
 ue=cross(El,B)/dot(B,B)  !*0.5
 Etemp=Ekin
 
! print*, 'Etemp=', Ekin
 
 !need to check if 1/2mUE^2 is covered by the initial KE
 IF (Ekin.lt.0.5d0*m*vscl*vscl*dot(ue,ue)*6.242e18) THEN
  PRINT*, 'WARNING: raising Initial KE to account for local UE drift'
  Etemp=0.5*m*vscl*vscl*dot(ue,ue)*6.242e18
 ENDIF
 
 gamma=Etemp/m/c/c*AQ+1.0d0
 
 ! vtot^2=vpar^2+vperp^2+UE^2!!, UE should NOT come out of Gamma!
 !vtot=c/gamma*sqrt(gamma*gamma-1)			! I *think* vtot is dimensional
 !vtot=sqrt(c*c/gamma/gamma*(gamma*gamma-1)-vscl*vscl*dot(ue,ue))
 !vtot=sqrt(c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue))			! THRELFALL ET AL 2015
 vtot=sqrt(c*c*(1.0d0-1.0d0/gamma/gamma)-vscl*vscl*dot(ue,ue))
 ! vtot=sqrt(c*c-c/gamma*c/gamma-vce*vce*dot(ue,ue))			! THRELFALL ET AL 2015
 
 !print*, 'gamma ', gamma
 !print*, 'uE^2:', vscl*vscl*dot(ue,ue)
 !print*, 'c*c-c/gamma*c/gamma', c*c-c/gamma*c/gamma
 !print*, 'c*c*(1.0d0-1.0d0/gamma/gamma)', c*c*(1.0d0-1.0d0/gamma/gamma)
 !PRINT*, c*c-c/gamma*c/gamma-vscl*vscl*dot(ue,ue)
 !print*, 'vtot ', vtot 
 
 IF (vtot/=vtot) THEN
  PRINT*, VTOT
  PRINT*, 'Vtot NaN ENCOUNTERED'
  STOP
 ENDIF
 
 !Ustart=(c*sqrt(gamma*gamma-1))*cos(alpha)
 Ustart=gamma*vtot*cos(alpha)
 
 
 !print*, vtot*cos(alpha)
 
 mu=0.5_num*m*vtot*vtot*sin(alpha)*sin(alpha)*gamma*gamma/modB		
 
 USTART=Ustart/vscl					! hence, non-dimensionalising..
 GAMMASTART=gamma
 mu=mu/m/vscl/vscl					! no bscl normalising factor - using normalised B's already!
 !PRINT*, 'vtot=', vtot
 !PRINT*, 'gamma=', gamma
 !PRINT*, 'mu=', mu
 !PRINT*, 'modB=', modB
  !calculate mu
  !mu=0.5_num*vperp*vperp/magB*gammastart*gammastart
 !STOP
 ! WRITE (19,*) RStart,T1,Ekin,Alpha, mu*magB, 0.5_num*vparstart*vparstart
  !WRITE (19,*) vtot,vperp,vparstart,El,B,magB,mu

END SUBROUTINE

end module cell_struct
