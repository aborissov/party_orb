PROGRAM reltest		! Relativistic Particle Code: Main Program (JT Dec 2015)
!+ sets up initial field conditions, particle grid and hands over to rkdriver

USE GLOBAL
!USE lare_functions
USE mpi_routines
USE M_DRIVERR, ONLY: RKDRIVE
USE bourdin_fields, ONLY: bour_ini, bour_fini
USE NLFF_fields, ONLY: NLFFF_ini, NLFFF_fini
USE MHDp_fields, ONLY: MHDp_ini, MHDp_fini
USE M_products, ONLY: DOT, CROSS
USE M_fields, ONLY: FIELDS
USE gammadist_mod, ONLY: random_gamma
USE omp_lib
USE cell_struct
!USE io

IMPLICIT NONE
 
  INTEGER :: NOK, NBAD ! ALEXEI: Check what these are
  REAL(num), DIMENSION(3) :: spatial_grid_step, box_size
  integer :: n_cells = 1, n_part_per_cell
  type(cell), dimension(:), allocatable :: cells

  !-----------------------------------------------------------------------------
  
  !welcome screen
  WRITE(*,*) "====RELATIVISTIC particle code====="

  ! read in input parameters
  ! This provides the initial particle grid bounds R, AlphaMax, AlphaMin
  CALL read_param
     
  ! Specify box limits
  IF ((str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN
   l2dflag=.TRUE.
   c_ndims=2
   spatial_extent_lower_bound = (/-1., -1., -1./)
   spatial_extent_upper_bound = (/1., 1., 1./)
   CALL MPI_INIT(errcode)
   CALL mpi_initialise_2d
   IF ((R1(3).NE.R2(3)).OR.(RSTEPS(3).GT.1)) THEN
    PRINT*, '-FATAL ERROR-'
    PRINT*, 'Lare2D data requires single z position value in initial grid!!'
    STOP
   ENDIF
   PRINT*, '..evaluating particle array against lare grid..'
  ELSE
    ! ALEXEI: temporary hardcoded limits added to make things run. Make this
    ! whole limit selector work better somehow.
    spatial_extent_lower_bound = (/-1., -1., -1./)
    spatial_extent_upper_bound = (/1., 1., 1./)
  ENDIF
   
  ! Check that particles do not start outside the box extent!
  DO i = 1,3
    IF (((R1(i)/lscl).le.spatial_extent_lower_bound(i)).OR.((R2(i)/lscl).ge.spatial_extent_upper_bound(i)))  THEN
      WRITE(*,*) '..particles not within extent of dimension ', i
      print*, "particle grid lower bound", R1(i)/lscl
      print*, "box lower bound", spatial_extent_lower_bound(i) 
      print*, "particle grid upper bound", R2(i)/lscl
      print*, "box upper bound", spatial_extent_upper_bound(i) 
    ENDIF
  END DO

  ! ALEXEI: debugging
  print *, "particle grid within field grid extent"
   
  ! Calculate total number of particles and step sizes for spatial and alpha grid
  nparticles=RSTEPS(1)*RSTEPS(2)*RSTEPS(3)*(AlphaSteps-1)*EkinSteps
  dAlpha = (AlphaMax-AlphaMin)/(AlphaSteps - 1.0d0) 
  box_size = (/R2(1)-R1(1),R2(2)-R1(2),R2(3)-R1(3)/)
  DO i=1,3
    IF (RSTEPS(i).EQ.1) THEN 	;! if rsteps=1, 1/(rsteps-1)=1/0!!
      spatial_grid_step(i)=box_size(i)
    ELSE 
      spatial_grid_step(i)=box_size(i)/REAL(RSTEPS(i)-1)
    ENDIF
  ENDDO

  ! ALEXEI: debugging
  print *, "nparticles calculated"
   
  ! Init cells
  ! ALEXEI: for testing just using one cell for now, so n_part_per_cell = nparticles
  n_part_per_cell = nparticles
  cells = init_cells(n_cells, n_part_per_cell, RSTEPS(1)*RSTEPS(2)*RSTEPS(3), AlphaSteps-1, EkinSteps)

  ! ALEXEI: debugging
  print *, "cells initialised"

  ! Integrate particle orbits in the cells.
  print *, "n_cells n_part_per_cell ", n_cells, n_part_per_cell
  do i = 1, n_cells
    call cells(i) % process(n_part_per_cell, i-1, nok, nbad)
  end do

  ! ALEXEI: debugging
  print *, "cells processed"

 IF ((str_cmp(FMOD, "L2D")).OR.(str_cmp(FMOD, "l2d"))) THEN   !forget arrays at end
  CALL mpi_close                     ! mpi_routines.f90
  CALL MPI_FINALIZE(errcode)
 ENDIF

!ALEXEI: remember to deallocate all the cells

!------------------------------------------------------------------------------!
 contains
!------------------------------------------------------------------------------!
function init_cells(n_cells, n_part_per_cell, n_pos_per_cell, n_gamma_per_cell, n_alpha_per_cell) result(cells)
  integer, intent(in)   :: n_cells              ! number of cells
  integer, intent(in)   :: n_pos_per_cell       ! number of particle positions per cell
  integer, intent(in)   :: n_gamma_per_cell     ! number of gammas per cell
  integer, intent(in)   :: n_alpha_per_cell     ! number of pitch angles per cell
  integer, intent(in)   :: n_part_per_cell      ! total number of particles per cell
  type(cell), dimension(:), allocatable :: cells    ! cell structs array
  integer               :: i,j,k,l,part_index   ! counters

  allocate(cells(n_cells))

  ! ALEXEI: need to implement randomised ICs so crash for now
  if (randomise_r .or. randomise_a .or. randomise_e) then
    print *, "We can't do random ICs yet..."
    call abort
  end if

  do i = 1,n_cells
    ! ALEXEI: crash if using more than one cell, need to check that we're
    ! initialising correctly
    if (n_cells > 1) then
      print *, "make sure we initialise correctly."
      call abort
    end if

    ! Allocate the coordinate arrays, particle velocities and times
    allocate(cells(i) % particle_coords(3, n_part_per_cell))
    allocate(cells(i) % gamma(n_part_per_cell))
    allocate(cells(i) % alpha(n_part_per_cell))
    allocate(cells(i) % particle_t(n_part_per_cell))

    ! Set the allocated quantities to their initial values
    do j = 1,n_pos_per_cell
      do k = 1,n_gamma_per_cell
        do l = 1,n_alpha_per_cell
          part_index = j*k*l
          cells(i) % particle_mass = M ! ALEXEI: I'm sure we don't need this here.
          cells(i) % particle_coords(:, part_index) = R1(:) + spatial_grid_step(:) * ((i-1)*n_part_per_cell + (j-1))
          ! ALEXEI: make sure gamma is high enough to cover the ExB drift !!!
          ! (see old JTmucalc
          cells(i) % gamma(part_index) = (EkinLow + k * (EkinHigh - EkinLow)) &
                                         / (cells(i) % particle_mass *c*c) + 1.0_num
          cells(i) % alpha(part_index) = AlphaMin + l * (AlphaMax - AlphaMin) 
          cells(i) % particle_t(part_index) = t1
        end do
      end do
    end do

    ! Set the grid bounds
    cells(i) % grid_upper_bound = spatial_extent_upper_bound
    cells(i) % grid_lower_bound = spatial_extent_lower_bound
  end do
end function init_cells

END PROGRAM reltest
