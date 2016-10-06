MODULE l3dc_fields
! Module contains several routines which set up the global lare grid and variable arrays.
  
  USE GLOBAL
  USE iocommon
  USE iocontrol
  USE input
  USE input_cartesian
  USE lare_functions


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: L3DCGRID, L3DCINIFIELDS!, L3DFIELDS

  CONTAINS 
!-------------------------------------------------------------------  
SUBROUTINE L3DCGRID
! subroutine to read in ONLY the Lare grid of a *.cfd file
! features several horrible hacks of the lare3d code. Sorry Tony Arber!

 CHARACTER(LEN = 20) :: name, class, mesh_name, mesh_class
 INTEGER, DIMENSION(3) :: dims
 REAL(num) :: time_d
 REAL(num), DIMENSION(3) :: extent
 REAL(num), DIMENSION(3) :: stagger
 INTEGER:: nblocks, type, nd, sof, snap

  
  print*, 'reading grid from:', cfdloc
  
! BEGIN HACK OF LAREXD READ:  

   !CALL cfd_open(cfdloc, rank, comm, MPI_MODE_RDONLY)
   cfd_comm = comm
   cfd_rank = rank
   cfd_mode = MPI_MODE_RDONLY
   CALL cfd_open_read(cfdloc)
   nblocks = cfd_get_nblocks()

    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      !IF (rank == 0) PRINT *, ix,nblocks, name, class, type
      !print*, name
      
      IF (type == TYPE_SNAPSHOT) THEN	!if type=3
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
      END IF

      IF (type == TYPE_MESH) THEN	!if type=1

        CALL cfd_get_3d_cartesian_grid_all(myx, myy, myz)

        !CALL cfd_skip_block()
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN		!if type=2
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)

        IF (nd /= DIMENSION_3D) THEN
          IF (rank == 0) PRINT *, "Non 3D Dataset found in input file, ", &
              "ignoring and continuing."
          CALL cfd_skip_block()
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT *, "Non - Cartesian variable block found ", &
              "in file, ignoring and continuing"
          CALL cfd_skip_block()
          CYCLE
        END IF

        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block()
      ELSE
        ! Unknown block, just skip it
        CALL cfd_skip_block()
      END IF
    END DO

    CALL cfd_close()
    CALL MPI_BARRIER(comm, errcode)


END SUBROUTINE L3DCGRID
!------------------------------------------------------------------------------! 
SUBROUTINE L3DCINIFIELDS
! designed to read in ALL field variables from Lare3d cfd files. Apologies to T. Arber.

    INTEGER				:: ii
    CHARACTER(LEN = 20) :: name, class, mesh_name, mesh_class
    !INTEGER :: nblocks, type, nd, sof, snap
    INTEGER		:: nblocks, type, nd, sof, snap
    INTEGER, DIMENSION(3) :: dims
    REAL(num) :: time_d
    REAL(num), DIMENSION(3) :: extent
    REAL(num), DIMENSION(3) :: stagger
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: data
    
! BEGIN HACK OF LAREXD READ:
  PRINT*, 'reading variables from: ', cfdloc   

   ALLOCATE(data(0:nx, 0:ny, 0:nz))

   cfd_comm = comm
   cfd_rank = rank
   cfd_mode = MPI_MODE_RDONLY
   CALL cfd_open(cfdloc, rank, comm, MPI_MODE_RDONLY)
   
   nblocks = cfd_get_nblocks()

    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      IF (rank == 0) PRINT *, ix,nblocks, name, class, type
      !print*, name
      
      IF (type == TYPE_SNAPSHOT) THEN	!if type=3
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
	ltimes(frame)=time_d
      END IF

      IF (type == TYPE_MESH) THEN	!if type=1
        ! Strangely, LARE doesn't actually read in the grid from a file
        ! This can be fixed, but for the moment, just go with the flow and
        ! Replicate the old behaviour
 
 	!cfd_get_2d_cartesian_variable_parallel(data, subtype)
        !CALL cfd_get_3d_cartesian_grid_all(myx, myy, myz)
	
        CALL cfd_skip_block()
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN		!if type=2
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)
			
        IF (nd /= DIMENSION_3D) THEN
          IF (rank == 0) PRINT *, "Non 3D Dataset found in input file, ", &
              "ignoring and continuing."
          CALL cfd_skip_block()
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT *, "Non - Cartesian variable block found ", &
              "in file, ignoring and continuing"
          CALL cfd_skip_block()
          CYCLE
        END IF

        ! We now have a valid variable, let's load it up
        ! First error trapping
        CALL cfd_get_nd_cartesian_variable_metadata_all(nd, dims, extent, &
            stagger, mesh_name, mesh_class)

        IF (dims(1) /= nx_global+1 &
            .OR. dims(2) /= ny_global+1 .OR. dims(3) /= nz_global+1) THEN
          IF (rank == 0) PRINT *, "Size of grid represented by one more ", &
              "variables invalid. Continuing"
	  print*, name
          CALL cfd_skip_block
          CYCLE
        END IF

        IF (sof /= num) THEN
          IF (rank == 0) PRINT *, "Precision of data does not match ", &
              "precision of code. Continuing."
          CALL cfd_skip_block
        END IF

        ! We're not interested in the other parameters, so if we're here,
        ! load up the data
        CALL cfd_get_3d_cartesian_variable_parallel(data, subtype)
        ! Now have the data, just copy it to correct place
	!IF (str_cmp(name(1:3), "Rho")) THEN
        !  rho(0:nx, 0:ny, 0:nz) = data
        !END IF
        !IF (str_cmp(name(1:6), "Energy")) THEN
        !  energy(0:nx, 0:ny, 0:nz) = data
        !END IF

        IF (str_cmp(name(1:2), "Vx")) THEN
!  	 PRINT *, ix,nblocks, name, class
          vx(1:nx, 1:ny, 1:nz,frame) = data(1:nx,1:ny,1:nz)
        END IF
        IF (str_cmp(name(1:2), "Vy")) THEN
! 	 PRINT *, ix,nblocks, name, class
          vy(1:nx, 1:ny, 1:nz,frame) = data(1:nx,1:ny,1:nz)
        END IF
        IF (str_cmp(name(1:2), "Vz")) THEN
!	 PRINT *, ix,nblocks, name, class
          vz(1:nx, 1:ny, 1:nz,frame) = data(1:nx,1:ny,1:nz)
        END IF
        IF (str_cmp(name(1:2), "Bx")) THEN
	! destagger on the fly to same locations as vx, vy, vz
	  DO ii=1,nx
	   bx(ii,1:ny,1:nz,frame)=stagger_bx(data(ii,1:ny,1:nz))
	  ENDDO  		
	END IF
        IF (str_cmp(name(1:2), "By")) THEN
!	 PRINT *, ix,nblocks, name, class
        !  by(0:nx, 0:ny, 0:nz) = data
	 DO ii=0,ny
	  by(1:nx,ii,1:nz,frame)=stagger_by(data(1:nx,ii,1:nz))
	 ENDDO 	 
        END IF
        IF (str_cmp(name(1:2), "Bz")) THEN
!	 PRINT *, ix,nblocks, name, class

	 DO ii=1,nz
	  bz(1:nx,1:ny,ii,frame)=stagger_bz(data(1:nx,1:ny,ii))
	 ENDDO 
        END IF
	
	!IF (str_cmp(name(1:11), "Temperature")) THEN
        !  temperature(0:nx, 0:ny, 0:nz) = data
        !END IF
	!IF (str_cmp(name(1:8), "Pressure")) THEN
        !  pressure(0:nx, 0:ny, 0:nz) = data
        !END IF
	!IF (str_cmp(name(1:3), "eta")) THEN
        !  eta(0:nx, 0:ny, 0:nz) = data
        !END IF

        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block()
      ELSE
        ! Unknown block, just skip it
	
	!print*, name(1:8)
	!IF (str_cmp(name(1:8), "Pressure")) THEN
        !  pressure(0:nx, 0:ny, 0:nz) = data
        !END IF
	
        CALL cfd_skip_block()
      END IF
    END DO

    DEALLOCATE(data)
    CALL cfd_close()
    CALL MPI_BARRIER(comm, errcode)

END SUBROUTINE L3DCINIFIELDS


END MODULE l3dc_fields

