module cell_struct

use global
use m_products, only: dot, cross
use m_driverr, only: rkdrive
use m_fields, only: fields
use bifrost_fields_mod, only: bifrost_grid
implicit none

type, public :: cell
  ! Particle quantities
  ! ALEXEI: sort out which ones we need to minimise computation and memory
  ! access times
  real(num), dimension(:,:), allocatable :: particle_coords  ! Location of the guiding center
  real(num), dimension(:),  allocatable  :: gamma            ! particle lorentz factor
  real(num), dimension(:),  allocatable  :: alpha            ! particle pitch angle
  real(num), dimension(:),  allocatable  :: particle_t       ! particle time
  real(num)                              :: particle_mass    ! particle mass
  type(bifrost_grid)                     :: bifrost_grid     ! bifrost grid structure

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
  real(num) ::  v_tot_2 ! square of total velocity
  real(num), dimension(3) :: rstart

  !$omp parallel do private(mu,ustart,gammastart,Ekin,alpha,rstart,particle_index,reset_flag)
  do i = 1, n_part_per_cell
    particle_index = i + n_part_per_cell*cell_index
    rstart = self % particle_coords(:,i)
    
    ! Calculate magnetic moment, initial 4-velocity, gamma
    ! ALEXEI: do we really need to call JTmucalc, or can we put this
    ! functionality into init_cells?
    call JTmucalc(mu,ustart,self % gamma(i),alpha,rstart,t1,t2,reset_flag,self % bifrost_grid)

    ! ALEXEI: In original version there is a check that B is not too small. Need
    ! to replicate this!!! (use reset_flag)

    ! Integrate the particle orbit
    call rkdrive(rstart,ustart,self % gamma(i),mu,t1,t2,eps,h1,particle_index,nok,nbad,self % bifrost_grid)
  end do
  !$omp end parallel do
end subroutine process_cell

SUBROUTINE JTMUcalc(mu,USTART,GAMMASTART,alpha,RSTART,T1,T2,resetflag,BF_grid)

  REAL(num), DIMENSION(3),INTENT(IN) :: RSTART
  REAL(num), INTENT(IN) :: T1,T2, Alpha, gammastart
  REAL(num), INTENT(OUT) :: mu, Ustart
  REAL(num), DIMENSION(3) :: B,El,a2,a3,a4,a5,a6,a7,a8,a9,a10,ue, RT
  REAL(num) :: modB,vtot, gamma, Etemp
  LOGICAL, INTENT(OUT)   :: resetflag
  type(bifrost_grid)    :: BF_grid      ! bifrost grid structure
 
 resetflag=.FALSE.
 
 !calculate B, E, V at this point/time:
 CALL FIELDS(RSTART,T1,El,B,a2,a3,a4,a5,a6,a7,a8,a9,a10,T1,T2,BF_grid)

 !calculate magnitude of B
 modB=sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
 RT=RSTART

 ! ALEXEI: make sure this is properly handled
 IF (modb.le.lowbthresh) THEN 
  resetflag=.TRUE.
 ENDIF
 
 ! E X B drift
 ue=cross(El,B)/dot(B,B)  !*0.5
 !Etemp=Ekin
 !
 !!need to check if 1/2mUE^2 is covered by the initial KE
 !IF (Ekin.lt.0.5d0*m*vscl*vscl*dot(ue,ue)*6.242e18) THEN
 ! PRINT*, 'WARNING: raising Initial KE to account for local UE drift'
 ! Etemp=0.5*m*vscl*vscl*dot(ue,ue)*6.242e18
 !ENDIF
 
 gamma=gammastart
 
 vtot=sqrt(c*c*(1.0d0-1.0d0/gamma/gamma)-vscl*vscl*dot(ue,ue))
 
 IF (vtot/=vtot) THEN
  PRINT*, VTOT
  PRINT*, 'Vtot NaN ENCOUNTERED'
  STOP
 ENDIF
 
 !Ustart=(c*sqrt(gamma*gamma-1))*cos(alpha)
 Ustart=gamma*vtot*cos(alpha)
 
 mu=0.5_num*m*vtot*vtot*sin(alpha)*sin(alpha)*gamma*gamma/modB		
 
 USTART=Ustart/vscl					! hence, non-dimensionalising..
 mu=mu/m/vscl/vscl					! no bscl normalising factor - using normalised B's already!

 if (isnan(mu)) then
   print *, m, vtot, sin(alpha), gamma, modB, vscl
   call abort
 endif

END SUBROUTINE

end module cell_struct
