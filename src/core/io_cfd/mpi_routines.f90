MODULE mpi_routines

  USE global
  USE l3dc_fields, ONLY:  L3DCINIFIELDS, L3DCGRID
  USE l3ds_fields, ONLY:  L3DSINIFIELDS, L3DSGRID
  USE l2dc_fields, ONLY:  L2DCINIFIELDS, L2DCGRID
  USE l2ds_fields, ONLY:  L2DSINIFIELDS, L2DSGRID 
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_initialise_2d, mpi_close
  PRIVATE:: mpi_create_types_2d, mpi_create_types_3d,mpi_destroy_types

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_initialise
!  Subroutine which sets up the Lare3d grids from either cfd or sdf files, and related settings

    INTEGER :: icoord
    INTEGER :: x_coords, y_coords, z_coords
    LOGICAL :: periods(3), reorder, ce, se, ce2, se2
    INTEGER :: starts(3), sizes(3), subsizes(3)
    INTEGER :: nx0, ny0, nz0
    INTEGER :: nxp, nyp, nzp
    INTEGER :: cx, cy, cz

101 format(a) 

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)
    
    ALLOCATE(coordinates(c_ndims), n_global_min(c_ndims), n_global_max(c_ndims))
    ALLOCATE(global_dims(c_ndims),local_dims(c_ndims))
    ALLOCATE(extents(2*c_ndims))

    dims = (/ nprocy, nprocx,nprocz /)
    local_dims = (/nx, ny,nz/)
    global_dims = (/nx_global, ny_global,nz_global/)


    IF (MAX(dims(1), 1)*MAX(dims(2), 1)*MAX(dims(3), 1) .GT. nproc) THEN
      dims = 0
      IF (rank .EQ. 0) THEN
        PRINT *, "Too many processors requested in override."
        PRINT *, "Reverting to automatic decomposition."
        PRINT *, "******************************************"
        PRINT *, ""
      END IF
    END IF

    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    nprocx = dims(3)
    nprocy = dims(2)
    nprocz = dims(1)

    periods=.FALSE.
    reorder=.FALSE.


    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, &
        reorder, comm, errcode)

    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, 3, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_y_min, proc_y_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_z_min, proc_z_max, errcode)

    cx = coordinates(3)
    cy = coordinates(2)
    cz = coordinates(1)

    ! Create the subarray for this problem: subtype decribes where this
    ! process's data fits into the global picture.

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nz0 = nz_global / nprocz
    nx = nx0
    ny = ny0
    nz = nz0

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx .NE. nx_global) THEN
      nxp = (nx0 + 1) * nprocx - nx_global
      IF (cx .GE. nxp) nx = nx0 + 1
    ELSE
      nxp = nprocx
    ENDIF
    IF (ny0 * nprocy .NE. ny_global) THEN
      nyp = (ny0 + 1) * nprocy - ny_global
      IF (cy .GE. nyp) ny = ny0 + 1
    ELSE
      nyp = nprocy
    ENDIF
    IF (nz0 * nprocz .NE. nz_global) THEN
      nzp = (nz0 + 1) * nprocz - nz_global
      IF (cz .GE. nzp) nz = nz0 + 1
    ELSE
      nzp = nprocz
    ENDIF

    ALLOCATE(cell_nx_mins(0:nprocx-1), cell_nx_maxs(0:nprocx-1))
    ALLOCATE(cell_ny_mins(0:nprocy-1), cell_ny_maxs(0:nprocy-1))
    ALLOCATE(cell_nz_mins(0:nprocz-1), cell_nz_maxs(0:nprocz-1))
    ! Set up the starting point for my subgrid (assumes arrays start at 0)

    DO icoord = 0, nxp - 1
      cell_nx_mins(icoord) = icoord * nx0 + 1
      cell_nx_maxs(icoord) = (icoord + 1) * nx0
    END DO
    DO icoord = nxp, nprocx - 1
      cell_nx_mins(icoord) = nxp * nx0 + (icoord - nxp) * (nx0 + 1) + 1
      cell_nx_maxs(icoord) = nxp * nx0 + (icoord - nxp + 1) * (nx0 + 1)
    END DO

    DO icoord = 0, nyp - 1
      cell_ny_mins(icoord) = icoord * ny0 + 1
      cell_ny_maxs(icoord) = (icoord + 1) * ny0
    END DO
    DO icoord = nyp, nprocy - 1
      cell_ny_mins(icoord) = nyp * ny0 + (icoord - nyp) * (ny0 + 1) + 1
      cell_ny_maxs(icoord) = nyp * ny0 + (icoord - nyp + 1) * (ny0 + 1)
    END DO

    DO icoord = 0, nzp - 1
      cell_nz_mins(icoord) = icoord * nz0 + 1
      cell_nz_maxs(icoord) = (icoord + 1) * nz0
    END DO
    DO icoord = nzp, nprocz - 1
      cell_nz_mins(icoord) = nzp * nz0 + (icoord - nzp) * (nz0 + 1) + 1
      cell_nz_maxs(icoord) = nzp * nz0 + (icoord - nzp + 1) * (nz0 + 1)
    END DO

    n_global_min(1) = cell_nx_mins(cx) - 1
    n_global_max(1) = cell_nx_maxs(cx)

    n_global_min(2) = cell_ny_mins(cy) - 1
    n_global_max(2) = cell_ny_maxs(cy)

    n_global_min(3) = cell_nz_mins(cz) - 1
    n_global_max(3) = cell_nz_maxs(cz)

    nx = n_global_max(1) - n_global_min(1)
    ny = n_global_max(2) - n_global_min(2)
    nz = n_global_max(3) - n_global_min(3)


    IF (cx .LT. nxp) THEN
      starts(1) = cx * nx0
    ELSE
      starts(1) = nxp * nx0 + (cx - nxp) * (nx0 + 1)
    ENDIF
    IF (cy .LT. nyp) THEN
      starts(2) = cy * ny0
    ELSE
      starts(2) = nyp * ny0 + (cy - nyp) * (ny0 + 1)
    ENDIF
    IF (cz .LT. nzp) THEN
      starts(3) = cz * nz0
    ELSE
      starts(3) = nzp * nz0 + (cz - nzp) * (nz0 + 1)
    ENDIF

    ! the grid sizes
    subsizes = (/ nx+1, ny+1, nz+1 /)
    sizes = (/ nx_global+1, ny_global+1, nz_global+1 /)

    ! set up and commit the subarray type
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subtype, errcode)

    CALL MPI_TYPE_COMMIT(subtype, errcode)

    ! Calculate initial displacement value:
    ! nx, ny, nz, (xb, yb, zb, time) * size of float
    initialdisp = 12 + (nx_global + ny_global + 3) * num

   nx=nx_global
   ny=ny_global
   nz=nz_global
   
   !ALLOCATE(rho(0:nx,0:ny,0:nz))
   !ALLOCATE(energy(0:nx,0:ny,0:nz)) 
   ALLOCATE(vx(0:nx,0:ny,0:nz,0:nframes-1))
   ALLOCATE(vy(0:nx,0:ny,0:nz,0:nframes-1))
   ALLOCATE(vz(0:nx,0:ny,0:nz,0:nframes-1))
   ALLOCATE(bx(0:nx,0:ny,0:nz,0:nframes-1))
   ALLOCATE(by(0:nx,0:ny,0:nz,0:nframes-1))
   ALLOCATE(bz(0:nx,0:ny,0:nz,0:nframes-1))
   ALLOCATE(myx(0:nx))
   ALLOCATE(myy(0:ny))
   ALLOCATE(myz(0:nz))
   ALLOCATE(ltimes(0:nframes-1))

   !ALLOCATE(eta(0:nx,0:ny,0:nz))
   !ALLOCATE(temperature(0:nx,0:ny,0:nz))
   !ALLOCATE(pressure(0:nx,0:ny,0:nz))   gflag=.false.
  
  WRITE (istring,fmt1) mysnap 			! converting integer to string using an 'internal file'
  cfdloc=trim(adjustl(sloc))//trim(istring)//filetype1		! store new filename
  sdfloc=trim(adjustl(sloc))//trim(istring)//filetype2		! store new filename
  
  ce=.FALSE.
  se=.FALSE.
  
  INQUIRE(file=TRIM(cfdloc),exist=ce)
  INQUIRE(file=TRIM(sdfloc),exist=se)
  

  IF ((ce).NEQV.(se)) THEN	! do flags differ? (if so, one is ON!)
   IF (ce) THEN
    print*, 'found old lare3d file (cfd)'
    CALL L3DCGRID()		! lare grid reads in the grid from the output file
   END IF
   IF (se) THEN
    print*, 'found new lare3d file (sdf)'
    ALLOCATE(local_dims(c_ndims))
    local_dims = (/nx, ny, nz/)
    CALL mpi_create_types_3d
    CALL L3DSGRID()
   END IF
 ELSE
  IF (ce.AND.se) THEN
   print*, 'ERROR: BOTH cfd and sfd formats at same location, terminating'
   STOP
  ELSE
  !IF ((.NOT.ce).AND.(.NOT.se)) THEN
   print*, 'ERROR: NEITHER cfd and sfd formats found at location, terminating'
   print*, '-OR SOMETHING ELSE HAS GONE HORRIBLY WRONG-'
   STOP
  END IF
 ENDIF
   
   WRITE(*,101)  'checking grid extent: '
   
   IF (evenlarefield) THEN
    WRITE(*,101,advance='no') '..EVEN lare grid used..'
    IF (.not. ((maxval(myx).ge.x_end).OR.(minval(myx).le.(-1.0_num*x_end)))) THEN
      WRITE(*,*) '..lare x extent < x_end - check boundaries are consistent'
      gflag=.true.
    ENDIF
    IF (.not. ((maxval(myy).ge.y_end).OR.(minval(myy).le.(-1.0_num*y_end)))) THEN
      WRITE(*,*) '..lare y extent < y_end - check boundaries are consistent'
      gflag=.true.
    ENDIF
    IF (.not. ((maxval(myz).ge.z_end).OR.(minval(myz).le.(-1.0_num*z_end)))) THEN
      WRITE(*,*) '..lare z extent < z_end - check boundaries are consistent'
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating as fields module will not work beyond bounds set in Lare.'
      STOP
    ELSE 
      WRITE(*,101,advance='no') 'fine!'
    ENDIF 
   ELSE
    WRITE(*,101,advance='no') '..ODD lare grid used..'
    IF (.not. ((maxval(myx).ge.xe(2)).OR.(minval(myx).le.xe(1)))) THEN
      WRITE(*,*) '..particle start positions beyond specified lare grid range in x'
      gflag=.true.
    ENDIF
    IF (.not. ((maxval(myy).ge.ye(2)).OR.(minval(myy).le.ye(1)))) THEN
      WRITE(*,*) '..particle start positions beyond specified lare grid range in y'
      gflag=.true.
    ENDIF
    IF (.not. ((maxval(myz).ge.ze(2)).OR.(minval(myz).le.ze(1)))) THEN
      WRITE(*,*) '..particle start positions beyond specified lare grid range in z'
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particles cannot begin beyond bounds set by Lare fields.'
      STOP
    ELSE 
      WRITE(*,101,advance='no') 'fine!'
    ENDIF 
   ENDIF
    !WRITE(*,101,advance='no')  '...done!'
    WRITE(*,*)  '..now loading in Lare variables..'
   
   frame=0
   DO WHILE (frame.LT.(nframes))
    WRITE (istring,fmt1) mysnap+frame 			! converting integer to string using an 'internal file'
    cfdloc=trim(adjustl(sloc))//trim(istring)//filetype1		! store new filename
    sdfloc=trim(adjustl(sloc))//trim(istring)//filetype2		! store new filename
  
    IF (ce) THEN
      INQUIRE(file=TRIM(cfdloc),exist=ce2)
      IF (ce2) THEN
       CALL L3DCINIFIELDS()
      ELSE
       PRINT*, 'file "', TRIM(cfdloc), '" missing'
      ENDIF
    ENDIF
    IF (se) THEN
      INQUIRE(file=TRIM(sdfloc),exist=se2)
      IF (se2) THEN
       CALL L3DSINIFIELDS()
      ELSE
       PRINT*, 'file "', TRIM(sdfloc), '" missing'
      ENDIF
    ENDIF
    frame=frame+1
   END DO
   
    IF (rank == 0) start_time = MPI_WTIME()


  END SUBROUTINE mpi_initialise
!----------------------------------------------------!
  SUBROUTINE mpi_initialise_2d
!+ Initialisation for 2d Lare2d environment setup.  

    INTEGER :: icoord
    INTEGER :: x_coords, y_coords
    LOGICAL :: periods(2), reorder, ce, se, ce2, se2
    INTEGER :: starts(2), sizes(2), subsizes(2)
    INTEGER :: nx0, ny0
    INTEGER :: nxp, nyp
    INTEGER :: cx, cy

101 format(a) 

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)
    
    ALLOCATE(coordinates(c_ndims), n_global_min(c_ndims), n_global_max(c_ndims))
    ALLOCATE(global_dims(c_ndims))
    ALLOCATE(extents(2*c_ndims))

    dims = (/ nprocy, nprocx /)
    global_dims = (/nx_global, ny_global/)

    IF (MAX(dims(1), 1)*MAX(dims(2), 1) .GT. nproc) THEN
      dims = 0
      IF (rank .EQ. 0) THEN
        PRINT *, "Too many processors requested in override."
        PRINT *, "Reverting to automatic decomposition."
        PRINT *, "******************************************"
        PRINT *, ""
      END IF
    END IF

    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    nprocx = dims(2)
    nprocy = dims(1)

    periods=.FALSE.
    reorder=.FALSE.


    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, &
        reorder, comm, errcode)

    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, 2, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_y_min, proc_y_max, errcode)


    cx = coordinates(2)
    cy = coordinates(1)

    ! Create the subarray for this problem: subtype decribes where this
    ! process's data fits into the global picture.

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy

    nx = nx0
    ny = ny0

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx .NE. nx_global) THEN
      nxp = (nx0 + 1) * nprocx - nx_global
      IF (cx .GE. nxp) nx = nx0 + 1
    ELSE
      nxp = nprocx
    ENDIF
    IF (ny0 * nprocy .NE. ny_global) THEN
      nyp = (ny0 + 1) * nprocy - ny_global
      IF (cy .GE. nyp) ny = ny0 + 1
    ELSE
      nyp = nprocy
    ENDIF


    ALLOCATE(cell_nx_mins(0:nprocx-1), cell_nx_maxs(0:nprocx-1))
    ALLOCATE(cell_ny_mins(0:nprocy-1), cell_ny_maxs(0:nprocy-1))

    ! Set up the starting point for my subgrid (assumes arrays start at 0)

    DO icoord = 0, nxp - 1
      cell_nx_mins(icoord) = icoord * nx0 + 1
      cell_nx_maxs(icoord) = (icoord + 1) * nx0
    END DO
    DO icoord = nxp, nprocx - 1
      cell_nx_mins(icoord) = nxp * nx0 + (icoord - nxp) * (nx0 + 1) + 1
      cell_nx_maxs(icoord) = nxp * nx0 + (icoord - nxp + 1) * (nx0 + 1)
    END DO

    DO icoord = 0, nyp - 1
      cell_ny_mins(icoord) = icoord * ny0 + 1
      cell_ny_maxs(icoord) = (icoord + 1) * ny0
    END DO
    DO icoord = nyp, nprocy - 1
      cell_ny_mins(icoord) = nyp * ny0 + (icoord - nyp) * (ny0 + 1) + 1
      cell_ny_maxs(icoord) = nyp * ny0 + (icoord - nyp + 1) * (ny0 + 1)
    END DO

  

    n_global_min(1) = cell_nx_mins(cx) - 1
    n_global_max(1) = cell_nx_maxs(cx)

    n_global_min(2) = cell_ny_mins(cy) - 1
    n_global_max(2) = cell_ny_maxs(cy)


    nx = n_global_max(1) - n_global_min(1)
    ny = n_global_max(2) - n_global_min(2)



    IF (cx .LT. nxp) THEN
      starts(1) = cx * nx0
    ELSE
      starts(1) = nxp * nx0 + (cx - nxp) * (nx0 + 1)
    ENDIF
    IF (cy .LT. nyp) THEN
      starts(2) = cy * ny0
    ELSE
      starts(2) = nyp * ny0 + (cy - nyp) * (ny0 + 1)
    ENDIF


    ! the grid sizes
    subsizes = (/ nx+1, ny+1 /)
    sizes = (/ nx_global+1, ny_global+1 /)

    ! set up and commit the subarray type
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subtype, errcode)

    CALL MPI_TYPE_COMMIT(subtype, errcode)

    ! Calculate initial displacement value:
    ! nx, ny, nz, (xb, yb, zb, time) * size of float
    initialdisp = 12 + (nx_global + ny_global + 3) * num

   nx=nx_global
   ny=ny_global

   
   !ALLOCATE(rho(0:nx,0:ny,0:nz))
   !ALLOCATE(energy(0:nx,0:ny,0:nz)) 
   ALLOCATE(vx(0:nx,0:ny,0:0,0:nframes-1))
   ALLOCATE(vy(0:nx,0:ny,0:0,0:nframes-1))
   ALLOCATE(vz(0:nx,0:ny,0:0,0:nframes-1))
   ALLOCATE(bx(0:nx,0:ny,0:0,0:nframes-1))
   ALLOCATE(by(0:nx,0:ny,0:0,0:nframes-1))
   ALLOCATE(bz(0:nx,0:ny,0:0,0:nframes-1))
   ALLOCATE(myx(0:nx))
   ALLOCATE(myy(0:ny))
   ALLOCATE(ltimes(0:nframes-1))

   !ALLOCATE(eta(0:nx,0:ny,0:nz))
   !ALLOCATE(temperature(0:nx,0:ny,0:nz))
   !ALLOCATE(pressure(0:nx,0:ny,0:nz))   gflag=.false.
  
  WRITE (istring,fmt1) mysnap 			! converting integer to string using an 'internal file'
  cfdloc=trim(adjustl(sloc))//trim(istring)//filetype1		! store new filename
  sdfloc=trim(adjustl(sloc))//trim(istring)//filetype2		! store new filename
  
  ce=.FALSE.
  se=.FALSE.
  
  INQUIRE(file=TRIM(cfdloc),exist=ce)
  INQUIRE(file=TRIM(sdfloc),exist=se)
  

  IF ((ce).NEQV.(se)) THEN	! do flags differ? (if so, one is ON!)
   IF (ce) THEN
    print*, 'found old lare2d file (cfd)'
    CALL L2DCGRID()		! lare grid reads in the grid from the output file
   END IF
   IF (se) THEN
    print*, 'found new lare2d file (sdf)'
    ALLOCATE(local_dims(c_ndims))
    local_dims = (/nx, ny/)
    CALL mpi_create_types_2d
    !print*, 'create_types_finishes'
    CALL L2DSGRID()
   END IF
 ELSE
  IF (ce.AND.se) THEN
   print*, 'ERROR: BOTH cfd and sfd formats at same location, terminating'
   STOP
  ELSE
  !IF ((.NOT.ce).AND.(.NOT.se)) THEN
   print*, 'ERROR: NEITHER cfd and sfd formats found at location, terminating'
   print*, '-OR SOMETHING ELSE HAS GONE HORRIBLY WRONG-'
   STOP
  END IF
 ENDIF
   
   WRITE(*,101)  'checking grid extent: '
   
   IF (evenlarefield) THEN
    WRITE(*,101,advance='no') '..EVEN lare grid used..'
    IF (.not. ((maxval(myx).ge.x_end).OR.(minval(myx).le.(-1.0_num*x_end)))) THEN
      !print*, maxval(myx), x_end
      WRITE(*,*) '..lare x extent < x_end - check boundaries are consistent'
      gflag=.true.
    ENDIF
    IF (.not. ((maxval(myy).ge.y_end).OR.(minval(myy).le.(-1.0_num*y_end)))) THEN
      WRITE(*,*) '..lare y extent < y_end - check boundaries are consistent'
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating as fields module will not work beyond bounds set in Lare.'
      STOP
    ELSE 
      WRITE(*,101,advance='no') 'fine!'
    ENDIF 
   ELSE
    WRITE(*,101,advance='no') '..ODD lare grid used..'
    IF (.not. ((maxval(myx).ge.xe(2)).OR.(minval(myx).le.xe(1)))) THEN
      WRITE(*,*) '..particle start positions beyond specified lare grid range in x'
      gflag=.true.
    ENDIF
    IF (.not. ((maxval(myy).ge.ye(2)).OR.(minval(myy).le.ye(1)))) THEN
      WRITE(*,*) '..particle start positions beyond specified lare grid range in y'
      gflag=.true.
    ENDIF
    IF (gflag) THEN
      WRITE(*,*) 'terminating: particles cannot begin beyond bounds set by Lare fields.'
      STOP
    ELSE 
      WRITE(*,101,advance='no') 'fine!'
    ENDIF 
   ENDIF
    !WRITE(*,101,advance='no')  '...done!'
    WRITE(*,*)  '..now loading in Lare variables..'
   
   frame=0
   DO WHILE (frame.LT.(nframes))
    WRITE (istring,fmt1) mysnap+frame 			! converting integer to string using an 'internal file'
    cfdloc=trim(adjustl(sloc))//trim(istring)//filetype1		! store new filename
    sdfloc=trim(adjustl(sloc))//trim(istring)//filetype2		! store new filename
  
    IF (ce) THEN
      INQUIRE(file=TRIM(cfdloc),exist=ce2)
      IF (ce2) THEN
       CALL L2DCINIFIELDS()
      ELSE
       PRINT*, 'file "', TRIM(cfdloc), '" missing'
      ENDIF
    ENDIF
    IF (se) THEN
      INQUIRE(file=TRIM(sdfloc),exist=se2)
      IF (se2) THEN
       CALL L2DSINIFIELDS()
      ELSE
       PRINT*, 'file "', TRIM(sdfloc), '" missing'
      ENDIF
    ENDIF
    frame=frame+1
   END DO
   
   !STOP
    IF (rank == 0) start_time = MPI_WTIME()


  END SUBROUTINE mpi_initialise_2d

  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

    IF (rank == 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
      WRITE(20, *)
      WRITE(20, '("runtime = ", i4, "h ", i2, "m ", i2, &
          & "s on ", i4, " process elements.")') hours, minutes, seconds, nproc
    END IF

    CALL MPI_BARRIER(comm, errcode)

    DEALLOCATE(vx)
    DEALLOCATE(vy)
    DEALLOCATE(vz)
    DEALLOCATE(bx)
    DEALLOCATE(by)
    DEALLOCATE(bz)
    DEALLOCATE(myx)
    DEALLOCATE(myy)
    DEALLOCATE(myz)

  END SUBROUTINE mpi_close
  
  SUBROUTINE mpi_create_types_3d

    !INTEGER, DIMENSION(3) :: nx, ny, nz 
    INTEGER :: sizes(c_ndims), subsizes(c_ndims), starts(c_ndims)
    !INTEGER :: local_dims(c_ndims), global_dims(c_ndims)
    INTEGER :: idir, vdir, mpitype
    INTEGER, PARAMETER :: ng = 2 ! Number of ghost cells


    ! File view for cell-centred variables (excluding the ghost cells)
    sizes = global_dims
    subsizes = local_dims
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_distribution = mpitype

    ! Subarray for cell-centred variable which has no ghost cells
    sizes = subsizes
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cellng_subarray = mpitype

    ! Cell-centred array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for cell-centred variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_subarray = mpitype

    ! MPI subtypes for communication of cell-centred variables

    ! ng cells, 1d slice of cell-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_yface = mpitype
    
   IF (c_ndims.eq.3) THEN
    idir = 3
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_zface = mpitype
   ENDIF
   
    ! File view for node-centred variables (excluding the ghost cells)
    sizes = global_dims + 1
    subsizes = local_dims + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_distribution = mpitype

    ! Subarray for node-centred variable which has no ghost cells
    sizes = subsizes
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    nodeng_subarray = mpitype

    ! Node-centred array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for node-centred variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_subarray = mpitype


    ! Array sizes for Bx-sized variables
    vdir = 1
    sizes = global_dims
    sizes(vdir) = sizes(vdir) + 1

    ! File view for Bx-sized variables (excluding the ghost cells)
    subsizes = local_dims
    subsizes(vdir) = subsizes(vdir) + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_distribution = mpitype

    ! Bx-sized array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for Bx-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_subarray = mpitype


    ! Array sizes for By-sized variables
    vdir = 2
    sizes = global_dims
    sizes(vdir) = sizes(vdir) + 1

    ! File view for By-sized variables (excluding the ghost cells)
    subsizes = local_dims
    subsizes(vdir) = subsizes(vdir) + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_distribution = mpitype

    ! By-sized array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for By-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_subarray = mpitype


    ! Array sizes for Bz-sized variables
    vdir = 3
    sizes = global_dims
    sizes(vdir) = sizes(vdir) + 1

    ! File view for Bz-sized variables (excluding the ghost cells)
    subsizes = local_dims
    subsizes(vdir) = subsizes(vdir) + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_distribution = mpitype

    ! Bz-sized array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for Bz-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bz_subarray = mpitype



  END SUBROUTINE mpi_create_types_3d
  
   SUBROUTINE mpi_create_types_2d

    INTEGER :: sizes(c_ndims), subsizes(c_ndims), starts(c_ndims)
    !INTEGER :: local_dims(c_ndims), global_dims(c_ndims)
    INTEGER :: idir, vdir, mpitype
    INTEGER, PARAMETER :: ng = 2 ! Number of ghost cells

    !local_dims = (/nx, ny/)
    !global_dims = (/nx_global, ny_global/)

    ! File view for cell-centred variables (excluding the ghost cells)
    sizes = global_dims
    subsizes = local_dims
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_distribution = mpitype

    ! Subarray for cell-centred variable which has no ghost cells
    sizes = subsizes
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cellng_subarray = mpitype

    ! Cell-centred array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for cell-centred variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_subarray = mpitype

    ! MPI subtypes for communication of cell-centred variables

    ! ng cells, 1d slice of cell-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    cell_yface = mpitype

    ! File view for node-centred variables (excluding the ghost cells)
    sizes = global_dims + 1
    subsizes = local_dims + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_distribution = mpitype

    ! Subarray for node-centred variable which has no ghost cells
    sizes = subsizes
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    nodeng_subarray = mpitype

    ! Node-centred array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for node-centred variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_subarray = mpitype

    ! MPI subtypes for communication of node-centred variables

    ! ng cells, 1d slice of node-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_yface = mpitype

    ! ng+1 cells, 1d slice of node-centred variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_xface1 = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    node_yface1 = mpitype

    ! Array sizes for Bx-sized variables
    vdir = 1
    sizes = global_dims
    sizes(vdir) = sizes(vdir) + 1

    ! File view for Bx-sized variables (excluding the ghost cells)
    subsizes = local_dims
    subsizes(vdir) = subsizes(vdir) + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_distribution = mpitype

    ! Bx-sized array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for Bx-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_subarray = mpitype

    ! MPI subtypes for communication of Bx-sized variables

    ! ng cells, 1d slice of Bx-sized variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_yface = mpitype

    ! ng+1 cells, 1d slice of Bx-sized variable

    idir = vdir
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    bx_xface1 = mpitype

    ! Array sizes for By-sized variables
    vdir = 2
    sizes = global_dims
    sizes(vdir) = sizes(vdir) + 1

    ! File view for By-sized variables (excluding the ghost cells)
    subsizes = local_dims
    subsizes(vdir) = subsizes(vdir) + 1
    starts = n_global_min

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_distribution = mpitype

    ! By-sized array dimensions
    sizes = subsizes + 2 * ng

    ! Subarray for By-sized variable which excludes the ghost cells
    starts = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_subarray = mpitype

    ! MPI subtypes for communication of By-sized variables

    ! ng cells, 1d slice of By-sized variable

    idir = 1
    subsizes = sizes
    subsizes(idir) = ng
    starts = 0

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_xface = mpitype

    idir = 2
    subsizes = sizes
    subsizes(idir) = ng

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_yface = mpitype

    ! ng+1 cells, 1d slice of By-sized variable

    idir = vdir
    subsizes = sizes
    subsizes(idir) = ng + 1

    mpitype = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)

    by_yface1 = mpitype

    ! MPI types for Bz-sized variables - same as cell-centred variable in 2d

    bz_distribution = cell_distribution
    bz_subarray = cell_subarray
    bz_xface = cell_xface
    bz_yface = cell_yface

  END SUBROUTINE mpi_create_types_2d
  
  
  SUBROUTINE mpi_destroy_types

    !CALL MPI_TYPE_FREE(cell_xface, errcode)
    !CALL MPI_TYPE_FREE(cell_yface, errcode)
    !CALL MPI_TYPE_FREE(cell_zface, errcode)
    !CALL MPI_TYPE_FREE(node_xface, errcode)
    !CALL MPI_TYPE_FREE(node_yface, errcode)
    !CALL MPI_TYPE_FREE(node_zface, errcode)
    !CALL MPI_TYPE_FREE(node_xface1, errcode)
    !CALL MPI_TYPE_FREE(node_yface1, errcode)
    !CALL MPI_TYPE_FREE(node_zface1, errcode)
    !CALL MPI_TYPE_FREE(bx_xface, errcode)
    !CALL MPI_TYPE_FREE(bx_yface, errcode)
    !CALL MPI_TYPE_FREE(bx_zface, errcode)
    !CALL MPI_TYPE_FREE(by_xface, errcode)
    !CALL MPI_TYPE_FREE(by_yface, errcode)
    !CALL MPI_TYPE_FREE(by_zface, errcode)
    !CALL MPI_TYPE_FREE(bz_xface, errcode)
    !CALL MPI_TYPE_FREE(bz_yface, errcode)
    !CALL MPI_TYPE_FREE(bz_zface, errcode)
    !CALL MPI_TYPE_FREE(bx_xface1, errcode)
    !CALL MPI_TYPE_FREE(by_yface1, errcode)
    !CALL MPI_TYPE_FREE(bz_zface1, errcode)
    CALL MPI_TYPE_FREE(cell_subarray, errcode)
    CALL MPI_TYPE_FREE(node_subarray, errcode)
    !CALL MPI_TYPE_FREE(cellng_subarray, errcode)
    !CALL MPI_TYPE_FREE(nodeng_subarray, errcode)
    CALL MPI_TYPE_FREE(cell_distribution, errcode)
    CALL MPI_TYPE_FREE(node_distribution, errcode)
    CALL MPI_TYPE_FREE(bx_subarray, errcode)
    CALL MPI_TYPE_FREE(by_subarray, errcode)
    CALL MPI_TYPE_FREE(bz_subarray, errcode)
    CALL MPI_TYPE_FREE(bx_distribution, errcode)
    CALL MPI_TYPE_FREE(by_distribution, errcode)
    CALL MPI_TYPE_FREE(bz_distribution, errcode)

  END SUBROUTINE mpi_destroy_types
  

END MODULE mpi_routines
