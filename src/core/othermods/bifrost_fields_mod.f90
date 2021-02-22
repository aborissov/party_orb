module bifrost_fields_mod

  use global
  use lare_functions, only: linterp3d

  implicit none

  public    :: bifrost_fields
  private   :: read_bifrost_snapshots, find_in_grid

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

  ! Data structure for holding grid associated with Bifrost simulation
  type, public :: bifrost_grid
    character(len=100)        :: filename   ! location of mesh file
    integer      :: nx, ny, nz   ! number of grid points in each direction
    real(num)    :: dx, dy, dz   ! grid spacing in each direction
    real(num), dimension(:), allocatable :: x    ! grid points in x direction
    real(num), dimension(:), allocatable :: y    ! grid points in y direction
    real(num), dimension(:), allocatable :: z    ! grid points in z direction

    type(bifrost_snap)  :: snap1, snap2          ! bifrost snapshots before and
                                                 ! after current time
  end type bifrost_grid

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
    real(num)               :: time_offset      ! the time offset of the current
                                                ! time between the bounding snapshots
    real(num), dimension(2) :: snapshot_time_list   ! list of snapshot times

    ! ALEXEI: change this to be read in as a parameter
    snapshot_time_list = (/ 0.0_num, 1000.0_num/)
    
    ! Check that we have the correct two snapshots read in
    if (BF_grid % snap1 % T > T .or. BF_grid % snap2 % T <= T) then
      print *, "Need to read snapshot time T snap1 snap2 ", T, BF_grid % snap1 % T, BF_grid % snap2 % T, BF_grid % snap1 % T > T, BF_grid % snap2 % T <= T
      call find_bounding_snapshots(T, BF_grid % snap1, BF_grid % snap2, snapshot_time_list,2)
      call read_bifrost_snapshots(T, BF_grid % snap1, BF_grid % snap2, BF_grid)
    endif
    print *, "finished reading snapshots, proceeding with calculation"

    ! Find position within Bifrost grid
    call find_in_grid(R, BF_grid, index, offset)

    ! Compute the time offset
    time_offset = (T - BF_grid % snap1 % T)/(BF_grid % snap2 % T - BF_grid % snap1 % T)

    ! Interpolate nearest gridpoints to get B, E
    call interpolate_fields(E, B, index, offset, time_offset, BF_grid % snap1, BF_grid % snap2)

    ! Compute derivatives

  end subroutine bifrost_fields

  ! Reads in Bifrost snapshots before and after the specified time
  subroutine read_bifrost_snapshots(time, snap1, snap2, grid)
    real(num), intent(in)   :: time         ! particle time
    type(bifrost_snap)      :: snap1, snap2 ! bifrost snapshots before and after
                                            ! particle time
    type(bifrost_grid)      :: grid         ! bifrost grid being used
    integer     :: snap1_fh, snap2_fh       ! filehandles for bifrost snapshots
    integer(kind=MPI_OFFSET_KIND)   :: offset   ! offset for data in bifrost file
    integer     :: mpierr

    ! ALEXEI: timing
    real    :: t0, t1

    ! ALEXEI: timing 
    call cpu_time(t0)

    ! ALEXEI: temporarily hardcode the snapshot location
    snap1 % grid_size = grid % nx * grid % ny * grid % nz
    snap2 % grid_size = grid % nx * grid % ny * grid % nz

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
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % Bx, grid % nx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % Bx, grid % nx)
    ! By
    offset = snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % By, grid % nx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % By, grid % nx)
    ! Bz
    offset = 2 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % Bz, grid % nx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % Bz, grid % nx)

    ! vx
    offset = 3 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % vx, grid % nx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % vx, grid % nx)
    ! vy
    offset = 4 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % vy, grid % nx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % vy, grid % nx)
    ! vz
    offset = 5 * snap1 % grid_size * snap1 % dtype_size
    call read_data_field(snap1_fh, offset, snap1 % grid_size, snap1 % vz, grid % nx)
    call read_data_field(snap2_fh, offset, snap2 % grid_size, snap2 % vz, grid % nx)

    ! Close the files
    call mpi_file_close(snap1_fh, mpierr)
    call check_mpi_error(mpierr)
    call mpi_file_close(snap2_fh, mpierr)
    call check_mpi_error(mpierr)

    ! ALEXEI: timing
    call cpu_time(t1)
    print *, "finished reading snapshot data in ", t1 - t0
  end subroutine read_bifrost_snapshots

  subroutine read_data_field(fh, offset, data_size, data_array, n_grid)
    integer, intent(in)     :: fh           ! filehandle associated with the file being read
    integer(kind=MPI_OFFSET_KIND), intent(in)     :: offset       ! offset within the file
    integer, intent(in)     :: data_size    ! size of data being read
    real(num), dimension(:,:,:)  :: data_array     ! data being read
    integer                      :: n_grid  ! ALEXEI: size of debugging array
    real, dimension(n_grid,n_grid,n_grid) :: temp_data    !  ALEXEI: debugging temporary array

    integer     :: mpierr
    real    :: t0, t1   ! timing

    call cpu_time(t0)
    ! Set the offsets
    call mpi_file_set_view(fh, offset, MPI_REAL, MPI_REAL, "native", &
                           MPI_INFO_NULL, mpierr)
    call check_mpi_error(mpierr)

    ! Read the data
    call mpi_file_read(fh, temp_data, data_size, MPI_REAL, status, mpierr)
    call check_mpi_error(mpierr)
    ! ALEXEI: this is stupid, need to figure out how to cast on the fly
    ! (although this could be useful for normalisation, so might need to keep it
    ! like this)
    data_array = temp_data
    
    call cpu_time(t1)
    print *, "finished reading one data field in ", t1 - t0
  end subroutine read_data_field

  ! Locates the particle guiding center within the Bifrost grid and returns the
  ! largest index below the guiding center location (in each direction) and the
  ! corresponding offset within the grid cell. This assumes a uniform grid
  subroutine find_in_grid(R, grid, index, offset)
    real(num), dimension(3), intent(in)  :: R      ! Location of guiding center
    type(bifrost_grid), intent(in)       :: grid   ! The bifrost grid we're using
    integer, dimension(3), intent(out)   :: index  ! the indices of the grid cell
                                                   ! containing the test particle
    real(num), dimension(3), intent(out) :: offset ! the particle's offset within
                                                   !the cell
    index(1) = ceiling((R(1) - grid % x(1)) / grid % dx)
    index(2) = ceiling((R(2) - grid % y(1)) / grid % dy)
    index(3) = ceiling((R(3) - grid % z(1)) / grid % dz)
    offset(1) = (R(1) - grid % x(index(1))) / grid % dx
    offset(2) = (R(2) - grid % y(index(2))) / grid % dy
    offset(3) = (R(3) - grid % z(index(3))) / grid % dz
  end subroutine find_in_grid

  ! Computes E and B based on the particle's position and time offsets
  subroutine interpolate_fields(E, B, i, offset, time_offset, snap1, snap2)
    real(num), dimension(3)     :: E, B, V  ! The E, B, and V fields
                                            ! at the location of the particle
    integer, dimension(3)       :: i        ! the indices in each spatial
                                            ! direction of the grid cell
                                            ! containing the particle.
    real(num), dimension(3)     :: offset   ! the offset of the particle in each 
                                            ! spatial direction within the grid
                                            ! cell
    real(num)                   :: time_offset  ! the temporal offset of the
                                                ! particle between the two snapshots
    type(bifrost_snap)          :: snap1, snap2 ! snapshots before and after
                                                ! current time
    
    B(1) = interpolate_single_field(snap1 % Bx, snap2 % Bx, i, offset, time_offset)
    B(2) = interpolate_single_field(snap1 % By, snap2 % By, i, offset, time_offset)
    B(3) = interpolate_single_field(snap1 % Bz, snap2 % Bz, i, offset, time_offset)
    V(1) = interpolate_single_field(snap1 % Vx, snap2 % Vx, i, offset, time_offset)
    V(2) = interpolate_single_field(snap1 % Vy, snap2 % Vy, i, offset, time_offset)
    V(3) = interpolate_single_field(snap1 % Vz, snap2 % Vz, i, offset, time_offset)
    
    ! ALEXEI: Note that this is the ideal Ohm's law for now
    E(1) = V(3)*B(2) - V(2)*B(3)
    E(2) = V(1)*B(3) - V(3)*B(1) 
    E(3) = V(2)*B(1) - V(1)*B(2) 
  end subroutine interpolate_fields

  ! Interpolates an individual field in space and time
  function interpolate_single_field(field1, field2, i, dx, dt)
    real(num), dimension(:,:,:)     :: field1   ! field at time 1 being interpolated
    real(num), dimension(:,:,:)     :: field2   ! field at time 2 being interpolated
    integer, dimension(3)           :: i        ! indices of cell to interpolate in
    real(num), dimension(3)         :: dx       ! offsets of cell to interpolate
    real(num)                       :: dt       ! time offset to interpolate
    real(num)   :: interpolate_single_field     ! result of interpolation
    real(num)   :: temp1, temp2     ! accumulators for spatially interpolated
                                    ! quantities at snapshot times

    temp1 = linterp3d(dx(1), dx(2), dx(3), &
                      field1(i(1),i(2),i(3)), field1(i(1)+1,i(2),i(3)), &
                      field1(i(1),i(2)+1,i(3)), field1(i(1)+1,i(2)+1,i(3)), &
                      field1(i(1),i(2),i(3)+1), field1(i(1)+1,i(2),i(3)+1), &
                      field1(i(1),i(2)+1,i(3)+1), field1(i(1)+1,i(2)+1,i(3)+1))    
    temp2 = linterp3d(dx(1), dx(2), dx(3), &
                      field2(i(1),i(2),i(3)), field2(i(1)+1,i(2),i(3)), &
                      field2(i(1),i(2)+1,i(3)), field2(i(1)+1,i(2)+1,i(3)), &
                      field2(i(1),i(2),i(3)+1), field2(i(1)+1,i(2),i(3)+1), &
                      field2(i(1),i(2)+1,i(3)+1), field2(i(1)+1,i(2)+1,i(3)+1))    
    if (temp1*(1.0_num - dt) + temp2*dt == 0) print *, "interopolation zero ", temp1, (1.0_num - dt), temp2, dt
    interpolate_single_field = temp1*(1.0_num - dt) + temp2*dt
  end function interpolate_single_field

  ! Finds the snapshots immediately before and after the current particle time
  ! based on a list of snapshot times. It then sets the appropriate filenames
  ! and times for the snapshots in preparation for reading them in
  subroutine find_bounding_snapshots(T, snap1, snap2, time_list, n_snap)
    real(num)           :: T                ! Particle time
    type(bifrost_snap)  :: snap1, snap2     ! bifrost snapshot structures to read into
    real(num), dimension(:)     :: time_list    ! list of snapshot times
    integer             :: n_snap           ! number of snapshots
    real(num)           :: dt
    character(len=22)   :: root
    integer             :: snap_index

    ! directory and filename root for where to find the snapshots
    ! ALEXEI: make this an input of some sort
    root = "../bifrost_atmosphere_"

    ! Calculate
    dt = time_list(2) - time_list(1)
    snap_index = floor((T - time_list(1))/dt)
    if (snap_index > n_snap) then
      print *, "calculated snapshot index greater than number of snapshots"
      call abort
    end if

    ! Set the filenames and snapshot times
    write(snap1 % filename, "(A22,I4.4,A4)") root, snap_index, ".bin"
    write(snap2 % filename, "(A22,I4.4,A4)") root, snap_index+1, ".bin"
    snap1 % T = time_list(snap_index+1)
    snap2 % T = time_list(snap_index+2)
  end subroutine find_bounding_snapshots

  ! Reads the mesh file to define the Bifrost grid we're using and assigns
  ! associated quantities
  subroutine init_bifrost_grid(grid)
    type(bifrost_grid)  :: grid     ! the bifrost grid we're initialising
    integer             :: fh       ! the filehandle to the grid file
    integer             :: grid_precision = 4 ! size in bytes of grid data
    integer(kind=MPI_OFFSET_KIND)   :: file_size    ! size of grid file
    integer(kind=MPI_OFFSET_KIND)   :: offset ! offset within file
    integer             :: mpierr
    real(num), parameter    :: Mm_to_m = 1.0e6_num  ! conversion factor 

    ! ALEXEI: temporarily just hardcode the filename and number of grid points in.
    grid % filename = "../bifrost_atmosphere.mesh"
    grid % nx = 768
    grid % ny = 768
    grid % nz = 768

    ! Allocate dimensions of grid
    allocate(grid % x(grid % nx))
    allocate(grid % y(grid % ny))
    allocate(grid % z(grid % nz))

    ! Allocate space for snapshots
    allocate(grid % snap1 % Bx(grid % nx, grid % ny, grid % nz))
    allocate(grid % snap1 % By(grid % ny, grid % ny, grid % nz))
    allocate(grid % snap1 % Bz(grid % nz, grid % ny, grid % nz))
    allocate(grid % snap1 % vx(grid % nx, grid % ny, grid % nz))
    allocate(grid % snap1 % vy(grid % ny, grid % ny, grid % nz))
    allocate(grid % snap1 % vz(grid % nz, grid % ny, grid % nz))
    allocate(grid % snap2 % Bx(grid % nx, grid % ny, grid % nz))
    allocate(grid % snap2 % By(grid % ny, grid % ny, grid % nz))
    allocate(grid % snap2 % Bz(grid % nz, grid % ny, grid % nz))
    allocate(grid % snap2 % vx(grid % nx, grid % ny, grid % nz))
    allocate(grid % snap2 % vy(grid % ny, grid % ny, grid % nz))
    allocate(grid % snap2 % vz(grid % nz, grid % ny, grid % nz))

    ! Use one rank to open and read the file. Have to do it like this because
    ! the mesh is a formatted file.
    if (rank == 0) then
      open(unit=15,file=grid % filename, form="FORMATTED")
      read(15,*) grid % x, grid % y, grid % z
      close(15)
      grid % x = grid % x * Mm_to_m / Lscl
      grid % y = grid % y * Mm_to_m / Lscl
      grid % z = grid % z * Mm_to_m / Lscl
      print *, "min max x ", grid % x(1), grid % x(grid % nx)
      print *, "min max y ", grid % y(1), grid % y(grid % ny)
      print *, "min max z ", grid % z(1), grid % z(grid % nz)
    end if

    call mpi_bcast(grid % x, grid % nx, mpireal, 0, mpi_comm_world, mpierr)
    call check_mpi_error(mpierr)

    ! Compute grid spacing
    grid % dx = grid % x(2) - grid % x(1)
    grid % dy = grid % y(2) - grid % y(1)
    grid % dz = grid % z(2) - grid % z(1)

    print *, "initialised bifrost grid"
  end subroutine init_bifrost_grid

  subroutine deallocate_bifrost_grid(grid)
    type(bifrost_grid)  :: grid

    ! deallocate dimensions of grid
    deallocate(grid % x)
    deallocate(grid % y)
    deallocate(grid % z)

    ! deallocate space for snapshots
    deallocate(grid % snap1 % Bx)
    deallocate(grid % snap1 % By)
    deallocate(grid % snap1 % Bz)
    deallocate(grid % snap1 % vx)
    deallocate(grid % snap1 % vy)
    deallocate(grid % snap1 % vz)
    deallocate(grid % snap2 % Bx)
    deallocate(grid % snap2 % By)
    deallocate(grid % snap2 % Bz)
    deallocate(grid % snap2 % vx)
    deallocate(grid % snap2 % vy)
    deallocate(grid % snap2 % vz)
  end subroutine deallocate_bifrost_grid

end module bifrost_fields_mod
