module bifrost_fields_mod

  use global

  implicit none

  public    :: bifrost_fields
  private   :: read_bifrost_snapshots!, read_bifrost_grid, find_in_grid

  ! Data structure for holding grid associated with Bifrost simulation
  type, public :: bifrost_grid
    character(len=100)        :: filename   ! location of mesh file
    integer      :: nx, ny, nz   ! number of grid points in each direction
    real(num)    :: dx, dy, dz   ! grid spacing in each direction
    real(num), dimension(:), allocatable :: x    ! grid points in x direction
    real(num), dimension(:), allocatable :: y    ! grid points in y direction
    real(num), dimension(:), allocatable :: z    ! grid points in z direction
  end type bifrost_grid

  ! Data structure for holding data associated with a Bifrost snapshot
  ! ALEXEI: some of these variables should probably be in a common structure for
  ! parameters related to the reading
  type, public :: bifrost_snap
    character(len=100)  :: filename     ! location of snapshot ! ALEXEI: figure out 
                                        ! something to do with the constant length string
    integer             :: n_fields     ! number of fields in snapshot
    integer             :: grid_size    ! number of elements per field
    integer             :: dtype_size   ! size of datatype written in snapshot
    real(num)           :: T            ! time of snapshot
    real(num), dimension(:,:,:), allocatable :: Bx, By, Bz   ! Magnetic field components
    real(num), dimension(:,:,:), allocatable :: vx, vy, vz   ! Velocity field components
  end type bifrost_snap

  contains

  ! Returns the field components and derivatives in each direction at the
  ! location of the guiding centre from a set of Bifrost snapshots
  subroutine bifrost_fields(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf,BF_grid)
    REAL(num), DIMENSION(3),INTENT(OUT) :: B,E  ! Magnetic and electric fields 
                                                ! at the location of the particle 
                                                ! guiding centre
    REAL(num), DIMENSION(3),INTENT(OUT) :: DBDX,DBDY,DBDZ,DBDT   ! Derivatives of B 
    REAL(num), DIMENSION(3),INTENT(OUT) :: DEDX,DEDY,DEDZ,DEDT   ! Derivatives of E
    REAL(num), DIMENSION(3),INTENT(IN)  :: R    ! Location of the particle guiding centre 
    REAL(num), DIMENSION(3)             :: Vf   ! No clue...    
    REAL(num)                           :: T    ! current time

    type(bifrost_grid), intent(in)      :: BF_grid  ! grid associated with the bifrost simulation
    type(bifrost_snap)      :: BF_snap1, BF_snap2   ! Bifrost snapshots before and after 
                                                    ! current particle time
    integer, dimension(3)   :: index            ! Array containing the bifrost grid index 
                                                ! immediately below the location of the 
                                                ! guiding center
    real(num), dimension(3) :: offset           ! Array containg the offset of
                                                ! the particle within the bifrost grid cell
    
    ! Check that we have the correct two snapshots read in
    if (BF_snap1 % T > T .or. BF_snap2 % T <= T) then
      call read_bifrost_snapshots(T, BF_snap1, BF_snap2)
    endif

    ! Find position within Bifrost grid
    call find_in_grid(R, BF_grid, index, offset)

    ! Interpolate nearest gridpoints to get B, E

    ! Compute derivatives

  end subroutine bifrost_fields

  ! Reads in Bifrost snapshots before and after the specified time
  subroutine read_bifrost_snapshots(time, snap1, snap2)
    real(num), intent(in)   :: time         ! particle time
    type(bifrost_snap)      :: snap1, snap2 ! bifrost snapshots before and after
                                            ! particle time
    integer     :: snap1_fh, snap2_fh       ! filehandles for bifrost snapshots
    integer(kind=MPI_OFFSET_KIND)   :: offset   ! offset for data in bifrost file
    integer     :: mpierr

    ! Open the files
    call mpi_file_open(MPI_COMM_WORLD, snap1 % filename, MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, snap1_fh, mpierr)
    call check_mpi_error(mpierr)
    call mpi_file_open(MPI_COMM_WORLD, snap2 % filename, MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, snap2_fh, mpierr)
    call check_mpi_error(mpierr)

    ! Read in each of the fields
    ! ALEXEI: this is done explicitly for now, see if there is a better way to
    ! do this.

    ! Bx
    offset = 0
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % Bx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % Bx)
    ! By
    offset = snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % By)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % By)
    ! Bz
    offset = 2 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % Bz)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % Bz)

    ! vx
    offset = 3 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % vx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % vx)
    ! vy
    offset = 4 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % vy)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % vy)
    ! vz
    offset = 5 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % vz)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % vz)

    ! Close the files
    call mpi_file_close(snap1_fh, mpierr)
    call check_mpi_error(mpierr)
    call mpi_file_close(snap2_fh, mpierr)
    call check_mpi_error(mpierr)

  end subroutine read_bifrost_snapshots

  subroutine read_data_field(fh, offset, data_size, data_array)
    integer, intent(in)     :: fh           ! filehandle associated with the file being read
    integer(kind=MPI_OFFSET_KIND), intent(in)     :: offset       ! offset within the file
    integer, intent(in)     :: data_size    ! size of data being read
    real(num), dimension(:,:,:)  :: data_array     ! data being read

    integer     :: mpierr

    ! Set the offsets
    call mpi_file_set_view(fh, offset, MPI_REAL, MPI_REAL, "native", &
                           MPI_INFO_NULL, mpierr)
    call check_mpi_error(mpierr)

    ! Read the data
    call mpi_file_read(fh, data_array, data_size, MPI_REAL, status, mpierr)
    call check_mpi_error(mpierr)

  end subroutine read_data_field

  ! Locates the particle guiding center within the Bifrost grid and returns the
  ! largest index below the guiding center location (in each direction) and the
  ! corresponding offset within the grid cell. This assumes a uniform grid
  subroutine find_in_grid(R, grid, index, offset)
    real(num), dimension(3), intent(in)  :: R      ! Location of guiding center
    type(bifrost_grid), intent(in)            :: grid   ! The bifrost grid we're using
    integer, dimension(3), intent(out)   :: index  ! the indices of the grid cell
                                                   ! containing the test particle
    real(num), dimension(3), intent(out) :: offset ! the particle's offset within
                                                   !the cell
    index(1) = ceiling(R(1) / grid % dx)
    index(2) = ceiling(R(2) / grid % dx)
    index(3) = ceiling(R(3) / grid % dx)
    offset(1) = R(1) - grid % x(index(1))
    offset(2) = R(2) - grid % y(index(2))
    offset(3) = R(3) - grid % z(index(3))
  end subroutine find_in_grid

  ! Reads the mesh file to define the Bifrost grid we're using and assigns
  ! associated quantities
  subroutine init_bifrost_grid(grid)
    type(bifrost_grid)  :: grid     ! the bifrost grid we're initialising
    integer             :: fh       ! the filehandle to the grid file
    integer             :: grid_precision = 4 ! size in bytes of grid data
    integer(kind=MPI_OFFSET_KIND)   :: file_size    ! size of grid file
    integer(kind=MPI_OFFSET_KIND)   :: offset ! offset within file
    integer             :: mpierr

    call mpi_file_open(MPI_COMM_WORLD, grid % filename, MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, fh, mpierr)
    call check_mpi_error(mpierr)

    call mpi_file_get_size(fh, file_size, mpierr)
    call check_mpi_error(mpierr)

    ! If grid size is not specified in parameter file assume uniform grid
    grid % nx = file_size / (3 * grid_precision)
    grid % ny = file_size / (3 * grid_precision)
    grid % nz = file_size / (3 * grid_precision)
    ! ALEXEI: add ability to read grid size from parameter file

    allocate(grid % x(grid % nx))
    allocate(grid % y(grid % ny))
    allocate(grid % z(grid % nz))

    ! Set the offsets
    call mpi_file_set_view(fh, offset, MPI_REAL, MPI_REAL, "native", &
                           MPI_INFO_NULL, mpierr)
    call check_mpi_error(mpierr)

    ! Read the data
    call mpi_file_read(fh, grid % x, grid % nx, MPI_REAL, status, mpierr)
    call check_mpi_error(mpierr)
    call mpi_file_read(fh, grid % x, grid % ny, MPI_REAL, status, mpierr)
    call check_mpi_error(mpierr)
    call mpi_file_read(fh, grid % x, grid % nz, MPI_REAL, status, mpierr)
    call check_mpi_error(mpierr)

  end subroutine init_bifrost_grid

end module bifrost_fields_mod
