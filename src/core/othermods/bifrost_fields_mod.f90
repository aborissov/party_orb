module bifrost_fields_mod

  use global

  implicit none

  public    :: bifrost_fields
  private   :: read_bifrost_snapshots!, read_bifrost_grid, find_in_grid

  ! Data structure for holding grid associated with Bifrost simulation
  type, public :: bifrost_grid
    integer     :: nx, ny, nz   ! number of grid points in each direction
    real(num)   :: dx, dy, dz   ! grid spacing in each direction
    real(num), dimension(:), allocatable :: x    ! grid points in x
    real(num), dimension(:), allocatable :: y    ! grid points in y
    real(num), dimension(:), allocatable :: z    ! grid points in z
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
  subroutine bifrost_fields(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
    REAL(num), DIMENSION(3),INTENT(OUT) :: B,E  ! Magnetic and electric fields 
                                                ! at the location of the particle 
                                                ! guiding centre
    REAL(num), DIMENSION(3),INTENT(OUT) :: DBDX,DBDY,DBDZ,DBDT   ! Derivatives of B 
    REAL(num), DIMENSION(3),INTENT(OUT) :: DEDX,DEDY,DEDZ,DEDT   ! Derivatives of E
    REAL(num), DIMENSION(3),INTENT(IN)  :: R    ! Location of the particle guiding centre 
    REAL(num), DIMENSION(3)             :: Vf   ! No clue...    
    REAL(num)                           :: T    ! current time

    type(bifrost_grid)      :: BF_grid          ! grid associated with the bifrost simulation
    type(bifrost_snap)      :: BF_snap1, BF_snap2   ! Bifrost snapshots before and after 
                                                    ! current particle time
    integer, dimension(3)   :: index, offset    ! Arrays holding the bifrost grid index 
                                                ! immediately below the location of the 
                                                ! guiding center, and the offset
                                                ! from it.
    
    ! Check that we have the correct two snapshots read in
    if (BF_snap1 % T > T .or. BF_snap2 % T <= T) then
      call read_bifrost_snapshots(T, BF_snap1, BF_snap2)
    endif

    ! Find position within Bifrost grid
    !call find_in_grid(R, BF_grid, index, offset)

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
    call mpi_file_read(fh, data_size, data_size, MPI_REAL, status, mpierr)
    call check_mpi_error(mpierr)

  end subroutine read_data_field

  ! Locates the particle guiding center within the Bifrost grid and returns the
  ! largest index below the guiding center location (in each direction) and the
  ! corresponding offset within the grid cell. This assumes a uniform grid
  !subroutine find_in_grid(R, grid, index, offset)
  !  real(num), dimension(3), intent(in) :: R      ! Location of guiding center
  !  type(BF_grid), intent(in)           :: grid   ! The bifrost grid we're using
  !  integer, dimension(3)               :: index  ! the indices of the grid cell
  !                                                ! containing the test particle
  !  real(num), dimension(3)             :: offset ! the particle's offset within
  !                                                !the cell



  !end subroutine find_in_grid

end module bifrost_fields_mod
