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
  ! ALEXEI: temporary hardcoded limits added to make things run. Make this
  ! whole limit selector work better somehow.
  spatial_extent_lower_bound = (/-1., -1., -1./)
  spatial_extent_upper_bound = (/1., 1., 1./)
   
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
   
  ! Calculate total number of particles and step sizes for spatial and alpha grid
  nparticles=RSTEPS(1)*RSTEPS(2)*RSTEPS(3)*(AlphaSteps-1)*EkinSteps
  dAlpha = (AlphaMax-AlphaMin)/(AlphaSteps - 1.0d0) 
  DO i=1,3
    IF (RSTEPS(i).EQ.1) THEN 	;! if rsteps=1, 1/(rsteps-1)=1/0!!
      spatial_grid_step(i)=1.0_num 
    ELSE 
      spatial_grid_step(i)=1.0_num/REAL(RSTEPS(i)-1)
    ENDIF
  ENDDO
  box_size = (/R2(1)-R1(1),R2(2)-R1(2),R2(3)-R1(3)/)

  ! Init cells
  ! ALEXEI: for testing just using one cell for now, so n_part_per_cell = nparticles
  n_part_per_cell = nparticles
  cells = init_cells(n_cells, n_part_per_cell)

  do i = 1, n_cells
    !call process_cell(cells(i))
    call cells(i) % process(n_part_per_cell, nok, nbad)
  end do

!ALEXEI: remember to deallocate all the cells

!------------------------------------------------------------------------------!
 contains
!------------------------------------------------------------------------------!
function init_cells(n_cells, n_part_per_cell) result(cells)
  integer, intent(in) :: n_cells, n_part_per_cell
  type(cell), dimension(:), allocatable :: cells
  allocate(cells(n_cells))

  do i = 1,n_cells
    ! Allocate the coordinate arrays, particle velocities and times
    allocate(cells(i) % particle_coords(3, n_part_per_cell))
    allocate(cells(i) % particle_v_par(n_part_per_cell))
    allocate(cells(i) % particle_v_perp(n_part_per_cell))
    allocate(cells(i) % particle_t(n_part_per_cell))

    ! Set the allocated quantities to their initial values
    do j = 1,n_part_per_cell
      cells(i) % particle_coords(:, j) = R1(:) + spatial_grid_step(:) * (i*n_part_per_cell + j)
      cells(i) % particle_v_par = 0 ! ALEXEI: place holder put some actual math in here!!!
      cells(i) % particle_v_perp = 0 ! ALEXEI: place holder put some actual math in here!!!
      cells(i) % particle_t = t1
    end do

    ! Set the grid bounds
    cells(i) % grid_upper_bound = spatial_extent_upper_bound
    cells(i) % grid_lower_bound = spatial_extent_lower_bound
  end do
end function init_cells

END PROGRAM reltest
