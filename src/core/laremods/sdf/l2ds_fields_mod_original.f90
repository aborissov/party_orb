MODULE l2ds_fields

  USE GLOBAL
  USE sdf
  USE sdf_job_info
  USE sdf_common
  USE lare_functions

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: L2DSGRID, L2DSINIFIELDS

CONTAINS
!****************************************************************************
  SUBROUTINE L2DSGRID
  ! subroutine to read in ONLY the Lare grid of a *.sdf file
! features several horrible hacks of the lare3d code. Sorry Tony Arber!

    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1
    INTEGER :: blocktype, datatype, code_io_version, string_len
    INTEGER :: ierr, iblock, nblocks, geometry
    INTEGER, DIMENSION(4) :: dims
    LOGICAL :: restart_flag    
    INTEGER :: mpireal = MPI_DOUBLE_PRECISION
    INTEGER :: sdf_num = c_datatype_real8
    REAL(num) :: dx, dy
    !TYPE(jobid_type) :: jobid
    TYPE(sdf_file_handle) :: sdf_handle

    INTEGER :: step

    CALL get_job_id(jobid)
    step = -1
    !c_ndims=ndims
    
    !ALLOCATE(coordinates(c_ndims), n_global_min(c_ndims), n_global_max(c_ndims))
    !ALLOCATE(extents(2*c_ndims))
    !ALLOCATE(global_dims(c_ndims))

    IF (rank == 0) THEN
      PRINT*,'Attempting to restart from file: ', sdfloc
    END IF

    CALL sdf_open(sdf_handle, sdfloc, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
    string_len, restart_flag)

IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
      STOP
    END IF

    IF (.NOT.str_cmp(code_name, TRIM(c_code_name2))) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by ', &
            TRIM(c_code_name2) // '.', 'Unable to ', 'continue.'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
      STOP
    END IF

    IF (string_len > c_max_string_length) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file string lengths are too large to read.'
        PRINT*, 'Please increase the size of "c_max_string_length" in ', &
            'shared_data.F90 to ','be at least ', string_len
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
      STOP
    END IF

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      !CALL create_ascii_header
    END IF

    !IF (rank == 0) PRINT*, 'Input file contains', nblocks, 'blocks';

    CALL sdf_read_blocklist(sdf_handle)	!!!UP TO READING THE BLOCKS. 
    CALL sdf_seek_start(sdf_handle)

    global_dims = (/ nx_global+1, ny_global+1/)

    DO iblock = 1, nblocks
      !print*, iblock
      
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_constant)
        CYCLE

      CASE(c_blocktype_plain_mesh)
        IF (ndims /= c_ndims .OR. datatype /= sdf_num &
            .OR. .NOT.str_cmp(block_id, 'grid')) CYCLE
        !print*, 'GRID FOUND'
        
	CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims, extents)

        IF (geometry /= c_geometry_cartesian &
            .OR. ALL(dims(1:c_ndims) /= global_dims(1:c_ndims))) CYCLE

        ! Should read the grid from file at this point?
        x_min = extents(1)
        x_max = extents(c_ndims+1)
        y_min = extents(2)
        y_max = extents(c_ndims+2)


	!print*, x_min, x_max, y_min, y_max

        length_x = x_max - x_min
        length_y = y_max - y_min

    ! Initially assume uniform grid
        dx = length_x / REAL(nx_global, num)
        dy = length_y / REAL(ny_global, num)

    ! Grid cell boundary for x coordinates
    DO ix = 0, nx_global
      myx(ix) = x_min + REAL(ix, num) * dx
    END DO
    DO iy = 0, ny_global
      myy(iy) = y_min + REAL(iy, num) * dy
    END DO
    
      CASE(c_blocktype_plain_variable)
       CYCLE

      END SELECT
    END DO


    CALL sdf_close(sdf_handle)
    CALL MPI_BARRIER(comm, errcode)
  ! PRINT*, 'SUCCESSFULLY OPENED AND CLOSED SDF FILE!!'
  !print*, vx(0:10,0,0)

  END SUBROUTINE L2DSGRID
!****************************************************************************
  SUBROUTINE L2DSINIFIELDS
  ! subroutine to read in ONLY the Lare grid of a *.sdf file
! features several horrible hacks of the lare3d code. Sorry Tony Arber!

    INTEGER				:: ii
    INTEGER 				:: blocktype, datatype, code_io_version, string_len
    INTEGER 				:: ierr, iblock, nblocks, geometry
    LOGICAL				:: restart_flag 
    INTEGER				:: mpireal = MPI_DOUBLE_PRECISION
    INTEGER				:: sdf_num = c_datatype_real8
    INTEGER, DIMENSION(4) 		:: dims
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: data
    CHARACTER(LEN=c_id_length) 		:: code_name, block_id, mesh_id, str1
    

    !TYPE(jobid_type) :: jobid
    TYPE(sdf_file_handle) :: sdf_handle

    INTEGER :: step

    CALL get_job_id(jobid)

    step = -1

    IF (rank == 0) THEN
      PRINT*,'reading from file: ', sdfloc
    END IF

    CALL sdf_open(sdf_handle, sdfloc, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
    string_len, restart_flag)

IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
      STOP
    END IF

    IF (.NOT.str_cmp(code_name, TRIM(c_code_name2))) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by ', &
            TRIM(c_code_name2) // '.', 'Unable to ', 'continue.'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
      STOP
    END IF

    IF (string_len > c_max_string_length) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file string lengths are too large to read.'
        PRINT*, 'Please increase the size of "c_max_string_length" in ', &
            'shared_data.F90 to ','be at least ', string_len
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
      STOP
    END IF

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      ltimes(frame)=time
      !CALL create_ascii_header
    END IF

    !IF (rank == 0) PRINT*, 'Input file contains', nblocks, 'blocks';

    CALL sdf_read_blocklist(sdf_handle)	!!!UP TO READING THE BLOCKS. 
    CALL sdf_seek_start(sdf_handle)

    global_dims = (/ nx_global+1, ny_global+1, nz_global+1 /)

    DO iblock = 1, nblocks
      !print*, iblock
      
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_constant)
        !print*, 'constant:', block_id
        CYCLE
        !IF (str_cmp(block_id, 'dt')) THEN
        !  CALL sdf_read_srl(sdf_handle, dt_from_restart)
        !ELSE IF (str_cmp(block_id, 'time_prev')) THEN
        !  CALL sdf_read_srl(sdf_handle, time_prev)
        !ELSE IF (str_cmp(block_id, 'visc_heating')) THEN
        !  CALL sdf_read_srl(sdf_handle, total_visc_heating)
        !  IF (rank /= 0) total_visc_heating = 0
        !END IF
      CASE(c_blocktype_plain_mesh)
        IF (ndims /= c_ndims .OR. datatype /= sdf_num &
            .OR. .NOT.str_cmp(block_id, 'grid')) CYCLE
        !print*, 'grid:', block_id
        
	CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims, extents)

        IF (geometry /= c_geometry_cartesian &
            .OR. ALL(dims(1:c_ndims) /= global_dims(1:c_ndims))) CYCLE

        ! Should read the grid from file at this point?
        x_min = extents(1)
        x_max = extents(c_ndims+1)
        y_min = extents(2)
        y_max = extents(c_ndims+2)

      CASE(c_blocktype_plain_variable)
        IF (ndims /= c_ndims .OR. datatype /= sdf_num) CYCLE
        CALL sdf_read_plain_variable_info(sdf_handle, dims, str1, mesh_id)
        IF (.NOT.str_cmp(mesh_id, 'grid')) CYCLE
        IF (str_cmp(block_id, 'Rho')) THEN
	   CYCLE 
!          CALL check_dims(dims)
!          CALL sdf_read_plain_variable(sdf_handle, rho, &
!              cell_distribution, cell_subarray)
        ELSE IF (str_cmp(block_id, 'Energy')) THEN
 	   CYCLE
 !         CALL check_dims(dims)
 !         CALL sdf_read_plain_variable(sdf_handle, energy, &
 !             cell_distribution, cell_subarray)
        ELSE IF (str_cmp(block_id, 'Vx')) THEN
          dims = dims - 1
!	  print*, 'variable:', block_id
          CALL check_dims(dims)
	!  print*, 'dims:', dims
	!  print*, nx
          ALLOCATE(data(-2:nx+2, -2:ny+2))
          CALL sdf_read_plain_variable(sdf_handle, data, &
              node_distribution, node_subarray)
          vx(1:nx,1:ny,1,frame)=data(1:nx,1:ny)
          DEALLOCATE(data)
        ELSE IF (str_cmp(block_id, 'Vy')) THEN
!   	  print*, 'variable:', block_id
          dims = dims - 1
          CALL check_dims(dims)
	  ALLOCATE(data(-2:nx+2, -2:ny+2))
          CALL sdf_read_plain_variable(sdf_handle, data, &
              node_distribution, node_subarray)
	  vy(1:nx,1:ny,1,frame)=data(1:nx, 1:ny)
          DEALLOCATE(data)
        ELSE IF (str_cmp(block_id, 'Vz')) THEN
!	  print*, 'variable:', block_id
          dims = dims - 1
          CALL check_dims(dims)
	  ALLOCATE(data(-2:nx+2, -2:ny+2))
          CALL sdf_read_plain_variable(sdf_handle, data, &
              node_distribution, node_subarray)
          vz(1:nx,1:ny,1,frame)=data(1:nx, 1:ny)
	  DEALLOCATE(data)
        ELSE IF (str_cmp(block_id, 'Bx')) THEN
 !  	  print*, 'variable:', block_id
          dims(1) = dims(1) - 1
          CALL check_dims(dims)
	  ALLOCATE(data(-2:nx+2, -1:ny+2))
          CALL sdf_read_plain_variable(sdf_handle, data, &
              bx_distribution, bx_subarray)
	  ! destagger routines differ in 2d and 3d for B.    
	  DO ii=1,nx
	   bx(ii,1:ny,1,frame)=stagger_bx_2d(data(ii,1:ny))
	  ENDDO
	  DEALLOCATE(data)
        ELSE IF (str_cmp(block_id, 'By')) THEN
  !	  print*, 'variable:', block_id
          IF (c_ndims >= 2) dims(2) = dims(2) - 1
          CALL check_dims(dims)
          ALLOCATE(data(-1:nx+2, -2:ny+2))
          CALL sdf_read_plain_variable(sdf_handle, data, &
              by_distribution, by_subarray)
	  DO ii=1,ny
	   by(1:nx,ii,1,frame)=stagger_by_2d(data(1:nx,ii))
	  ENDDO 
	  DEALLOCATE(data)
        ELSE IF (str_cmp(block_id, 'Bz')) THEN
  !        print*, 'variable:', block_id
	  IF (c_ndims >= 3) dims(3) = dims(3) - 1
          CALL check_dims(dims)
	  ALLOCATE(data(-1:nx+2, -1:ny+2))
          CALL sdf_read_plain_variable(sdf_handle, data, &
              bz_distribution, bz_subarray)
	  bz(1:nx,1:ny,1,frame)=stagger_bz(data(1:nx,1:ny))
	  DEALLOCATE(data)
        END IF
      END SELECT
    END DO
   

    CALL sdf_close(sdf_handle)
    CALL MPI_BARRIER(comm, errcode)
  ! PRINT*, 'SUCCESSFULLY OPENED AND CLOSED SDF FILE!!'

  END SUBROUTINE L2DSINIFIELDS


END MODULE l2ds_fields
