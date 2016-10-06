MODULE global
!+ contains lots of constants and global variables
  
  USE sdf_job_info
  USE MPI
 
 IMPLICIT NONE
  ! INCLUDE 'mpif.h'
 
 SAVE 
  
!########################################################################## 
 CHARACTER(Len = 4), PARAMETER	:: FMOD='l2d' ! SWITCH BETWEEN FIELDS: "l3d","l2d", "SEP","CMT","test", or "bor"
 INTEGER, PARAMETER		:: mysnap=0001	!  no. of ****.cfd/****.sdf file (if "l3d")
 INTEGER, PARAMETER		:: nframes=1	! no. of frames
 CHARACTER(Len = 40)		:: sloc='lare2d_runs/Data/'

!##########################################################################
! now some stuff required to plug in lare data 
 
 INTEGER, PARAMETER		:: num = KIND(1.0d0), dbl = KIND(1.D0)
 INTEGER, PARAMETER		:: NKEEPMAX = 1000001		! max no of dumps
 INTEGER(KIND=8), PARAMETER	:: NSTPMAX  = 5E9	! max no of steps
 INTEGER(KIND=8)		:: NSTP				! step counter
 INTEGER, PARAMETER		:: NSTORE =50, MAXTIME=200	! how often does the RK solver output data, every NSTORE steps?
 INTEGER 			:: ix, iy, iz,it, iix, iiy, iiz,iit, i
 INTEGER 			:: frame 
 
!JT DEBUGGING SWITCHES:
 LOGICAL, PARAMETER		:: writervs=.TRUE.					! ARE WE WRITING? (ALWAYS TRUE!) 
 LOGICAL, PARAMETER		:: JTo=.TRUE., JTo2=.FALSE., JTo3=.FALSE., JTO4=.TRUE.	! various debugging switches (2&3 output every NSTP)
 LOGICAL, PARAMETER		:: FIELDDUMP=.FALSE.					! switch to dump the lare fields to unformatted data files.
 LOGICAL, PARAMETER		:: everystepswitch=.TRUE.				! dumps EVERY NSTP to each particle data file.
 
!PARTICLE quantities: 								
 REAL(num), DIMENSION(3) 	:: R1,R2, tempr
 REAL(num) 			:: time, tempa, tempe
 REAL(num) 			:: Ekin,Alpha,AlphaMax,AlphaMin,dalpha, Ekinlow,Ekinhigh,T1Keep,T2Keep
 REAL(num) 			:: T1,T2, H1,EPS, USTART, mu, USTARTKEEP, GAMMASTART,GAMMASTARTKEEP, VPARSTART
 INTEGER, DIMENSION(3) 		:: RSteps 
 INTEGER			:: pn, nparticles
 INTEGER			:: p_restart_no, p_stop_no
 INTEGER 			:: EKinSteps, AlphaSteps

 INTEGER                        :: adjust_timestep = 0                          ! AB 0 for fixed timestep, 1 for adjustable
 INTEGER                        :: collisions = 1                               ! AB 0 for no collisions, 1 for collisions
 INTEGER                        :: print_number = 1                             ! AB 0 for no printing particle number, 1 for printing particle number
 LOGICAL			:: p_restart=.FALSE., p_stop=.FALSE. 		! are we starting or stopping midway through the arrays?
 LOGICAL			:: RANDOMISE_R, RANDOMISE_A, RANDOMISE_E	! switches for randomising position, angle and energy (TWO WORK CURRENTLY!)

 LOGICAL		:: maxwellEfirst
 REAL(num), PARAMETER	:: maxwellpeaktemp= 1e6_num
 ! Local parameters
REAL (num), PARAMETER  	:: one = 1.0_num, zero = 0.0_num
 
! CONSTANTS
 REAL(num), PARAMETER	:: pi = 3.1415926535897932_num
 REAL(num), PARAMETER 	:: mu0_si =  4.0e-7_num * pi
 REAL(num), PARAMETER 	:: kb=1.380650500e-23_num
 REAL(num), PARAMETER	:: Mp = 1.6726E-27_num, Me=9.1093826E-31_num
 REAL(num), PARAMETER	:: Qp = 1.60217653E-19_num, Qe = -1.60217653E-19_num
 REAL(num), PARAMETER	:: Q = Qe, M = Me					!SELECT PROTON OR ELECTRON HERE
 REAL(num), PARAMETER	:: oneuponQ = 1.d0/Q, AQ = abs(Q), oneuponAQ= 1.d0/AQ
 REAL(num), PARAMETER	:: MoQ=M/Q, QoM=Q/M, oM=1.0d0/M, MoAQ=M/AQ, AQoM=AQ/M	!protons
 REAL(num), PARAMETER	:: c = 2.99792458E8_num, oneuponc=1.0_num/c   ! speed of light
 REAL(num), PARAMETER	:: oneotwelve=1.0_num/12.0_num 

! NORMALISING SCALES 
 REAL(num), PARAMETER	:: Lscl = 1.0E1_num     		! 10 Mega meters (1e7)
 REAL(num), PARAMETER	:: Bscl = 0.03_num		! 100 Gauss 	 (0.01)
 REAL(num), PARAMETER	:: Tscl = 4.82E-7_num		! 100s	
 REAL(num), PARAMETER	:: Escl = Lscl*Bscl/Tscl	! 10V/cm	 (1e3)
 REAL(num), PARAMETER	:: Vscl = Lscl/Tscl		! 10^7/10^2=10^5m/s=100km/s
 REAL(num), PARAMETER	:: Ekscl = M*Vscl*Vscl		! 10^10*9e-31=9e-21joules
 REAL(num), PARAMETER	:: Omscl=QoM*Bscl, oneuponOmscl=1.d0/Omscl		! qB/M
 REAL(num), PARAMETER	:: jscl = Bscl/Lscl/mu0_si
 REAL(num), PARAMETER	:: etascl = mu0_si*lscl*vscl
 REAL(num), PARAMETER	:: myT1=0.0_num/Tscl
 REAL(num), PARAMETER	:: pscl=bscl*bscl/mu0_si			! Lare normalisation of pressure
 REAL(num), PARAMETER	:: rhoscl=bscl*bscl/mu0_si/vscl/vscl		! Lare normalisation of mass density
 !REAL(num), PARAMETER	:: tempscl=pscl/rhoscl/Mp*kb			! Lare normalisation of temperature
 REAL(num), PARAMETER	:: tempscl=vscl*vscl*1.2*Mp/kb			! Lare normalisation of temperature AB

 LOGICAL, PARAMETER	:: lare_norm=.TRUE.	! TRUE if Lare output IS NORMALISED, else FALSE (so normalise varibles and derivs).
 						! (usually true, unless someone has messed with control.f90)

! BIRN ET AL NORMALISATION 
 REAL(num), PARAMETER	:: sigma=sqrt(Me/Mp), oneosigma=sqrt(Mp/Me) 				!sigma is small, oneosigma is big. 
 REAL(num), PARAMETER	:: myepsil=Mp/abs(Qe)/Bscl/tscl, oneomyepsil=abs(Qe)*Bscl*tscl/Mp	! epsilon small
 REAL(num), PARAMETER	:: vce=vscl/sigma


! parameters for Guiliani et al model
! (required by CMT module)
 REAL(num), PARAMETER	:: d=Lscl      ! depth of monopoles
 REAL(num), PARAMETER	:: Lv = Lscl   ! Height where stretching starts
 REAL(num), PARAMETER	:: Bzfinal =  0. !.1426754410 
 REAL(num), PARAMETER	:: esp=1.  ! esponent in log transformation
 REAL(num), PARAMETER	:: cc=0.4  ! coeff in log transformation
 REAL(num), PARAMETER	:: c1=-0.15E6   ! Magnetic Charge
 
 !Separator model parameters (Wilmot-Smith & Hornig, 2011)
 ! (required by SEP module)
 REAL(num), PARAMETER	:: rat1=20.0_num	!r1=ratio of B1 to B0, B1=r1*B0 (in paper B1=20B0)
 REAL(num), PARAMETER	:: rat2=1e-6_num	!r2=ratio of a to z0, a=r2*z0	(in paper a=1/2,z0=5)
 REAL(num), PARAMETER	:: rat3=0.2_num		!r3=ratio of l to z0, l=r3*z0	(in paper ll=1,z0=5)
 REAL(num), PARAMETER	:: b0=Bscl, E0=Escl		! field setup in WS&H (2011)
 REAL(num), PARAMETER	:: L=Lscl, oL=1.0_num/L		! field setup in WS&H (2011)
 REAL(num), PARAMETER	:: z0=5.0_num*L
 REAL(num), PARAMETER	:: b1=rat1*b0, a=rat2*z0, ll=rat3*z0, oa=1.0_num/a, oll=1.0_num/ll 
 REAL(num), PARAMETER	:: xc=0.0_num, yc=0.0_num, zc=0.0_num		! center of flux ring
 REAL(num), PARAMETER	:: tau_n=1.0_num, tau=1.0_num*Tscl
! REAL(num), PARAMETER	:: xymult=10.0_num, zmult=1.0_num

 ! Lare field parameters
 ! (required by LARE modules)

 !REAL(num), DIMENSION(2), PARAMETER	:: xe=(/0.1_num,99.9_num/),ye=(/-0.1_num,99.9_num/),ze=(/-19.5_num,79.5_num/)
 !REAL(num), DIMENSION(2), PARAMETER	:: xe=(/-0.9_num,0.9_num/),ye=(/-0.9_num,0.9_num/),ze=(/-100.00_num,100.0_num/)
 REAL(num), DIMENSION(2), PARAMETER	:: xe=(/-8.0_num,8.0_num/),ye=(/-6.1_num,6.1_num/),ze=(/-100.00_num,100.0_num/) ! working collision lare 2d grid
 !REAL(num), DIMENSION(2), PARAMETER	:: xe=(/-80.0_num,80.0_num/),ye=(/-60.1_num,60.1_num/),ze=(/-100.00_num,100.0_num/) ! for testing
 !REAL(num), PARAMETER			:: eta=0.001_num, jcrit=25.0_num
 !REAL(num), PARAMETER			:: eta=0.005_num, jcrit=0.9_num, rwidth=0.05_num
 REAL(num), PARAMETER			:: jcrit=0.9_num, rwidth=0.05_num
 !REAL(num), PARAMETER			:: eta=0.001_num, jcrit=20.0_num, rwidth=0.5_num

 CHARACTER(Len = 5), PARAMETER 	:: dloc='Data/'
 CHARACTER(Len = 11) 		:: dlocN=dloc//'DataN/',dlocR=dloc//'DataR/'

 INTEGER, PARAMETER				:: mpireal = MPI_DOUBLE_PRECISION
 !INTEGER, PARAMETER				:: nx_global=64, ny_global=64, nz_global=64		! julie's lare3d cfd config
 !INTEGER, PARAMETER				:: nx_global=120, ny_global=120, nz_global=480		! alan's lare3d sdf config
 !INTEGER, PARAMETER				:: nx_global=128, ny_global=128, nz_global=1		! James's lare2d cfd config
 !INTEGER, PARAMETER                             :: nx_global=256, ny_global=256, nz_global=1            ! James's lare2d cfd config
 INTEGER, PARAMETER                             :: nx_global=512, ny_global=512, nz_global=1            ! James's lare2d cfd config
 INTEGER, PARAMETER 				:: data_dir_max_length = 64
 INTEGER 					:: nx, ny, nz	
 INTEGER, DIMENSION(:), ALLOCATABLE		:: dims

 !REAL(num), DIMENSION(:, :, :,:), ALLOCATABLE 	:: bx, by, bz, vx, vy, vz, rho, temperature
 !REAL(num), DIMENSION(:, :, :), ALLOCATABLE 	:: jx, jy, jz, Ex, Ey, Ez
 REAL(num), DIMENSION(:, :, :,:), ALLOCATABLE 	:: bx, by, bz, vx, vy, vz, rho, temperature, eta, jx, jy, jz 
 REAL(num), DIMENSION(:, :, :), ALLOCATABLE 	:: Ex, Ey, Ez
 REAL(num), DIMENSION(:, :, :), ALLOCATABLE 	:: dbxdx,dbxdy,dbxdz,dbydx,dbydy,dbydz,dbzdx,dbzdy,dbzdz
 REAL(num), DIMENSION(:, :, :), ALLOCATABLE 	:: dExdx,dExdy,dExdz,dEydx,dEydy,dEydz,dEzdx,dEzdy,dEzdz
 REAL(num), DIMENSION(:), 	ALLOCATABLE 	:: myx, myy, myz, ltimes
  
 CHARACTER(LEN = data_dir_max_length) :: data_dir
 CHARACTER(Len = 14)	:: mydataloc
 CHARACTER(Len = 4)	:: filetype1='.cfd', filetype2='.sdf'
 CHARACTER(LEN = 20) 	:: name, class, mesh_name, mesh_class
 
 !---BOURDIN DATA DEFINITIONS---! 
 CHARACTER(Len = 8)	:: filetypeb='.bin_f77'
 LOGICAL		:: bourdinflag, l3dflag, analyticalflag, l2dflag
 REAL(num), DIMENSION(2)	:: xee, yee, zee
  
  ! MPI data
 INTEGER :: rank, proc_x_min, proc_x_max, proc_y_min, proc_y_max, proc_z_min, proc_z_max
 INTEGER :: errcode, comm, tag, nproc, nprocx, nprocy, nprocz
 INTEGER :: status(MPI_STATUS_SIZE)

  ! file handling
 INTEGER :: subtype, obstype
 INTEGER(KIND = MPI_OFFSET_KIND) :: initialdisp
 INTEGER :: initial
 INTEGER :: n_zeros = 4
 INTEGER :: output_file = 0
  
 LOGICAL :: gflag, pflag1=.true.
 
 CHARACTER(len=7) 				:: fmt1='(I4.4)', istring 	! format descriptor
 CHARACTER(LEN = 20+data_dir_max_length) 	:: cfdloc, sdfloc
 
 ! NEW LARE3D PLUGIN INFO
   !-- cut and pasted from "version info"--!
  CHARACTER(LEN=*), PARAMETER :: c_code_name = 'Lare3D'
  CHARACTER(LEN=*), PARAMETER :: c_code_name2 = 'Lare2D'
  TYPE(jobid_type) :: jobid
  INTEGER :: ndims, c_ndims
  INTEGER, DIMENSION(:), ALLOCATABLE :: coordinates, n_global_min, n_global_max
  INTEGER :: cell_subarray, cellng_subarray, cell_distribution
  INTEGER :: node_subarray, nodeng_subarray, node_distribution
  INTEGER :: bx_subarray, bx_distribution
  INTEGER :: by_subarray, by_distribution
  INTEGER :: bz_subarray, bz_distribution
  INTEGER :: cell_xface, node_xface, node_xface1
  INTEGER :: cell_yface, node_yface, node_yface1
  INTEGER :: cell_zface, node_zface, node_zface1
  INTEGER :: bx_xface, by_xface, bz_xface, bx_xface1
  INTEGER :: bx_yface, by_yface, bz_yface, by_yface1
  INTEGER :: bx_zface, by_zface, bz_zface, bz_zface1
  REAL(num) :: x_min, x_max, length_x
  REAL(num) :: y_min, y_max, length_y
  REAL(num) :: z_min, z_max, length_z
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_nx_mins, cell_nx_maxs
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_ny_mins, cell_ny_maxs
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_nz_mins, cell_nz_maxs
  
  INTEGER, DIMENSION(:), ALLOCATABLE	:: global_dims, local_dims
  REAL(num), DIMENSION(:), ALLOCATABLE :: extents, extent, stagger
 
!----------------------------------------------------
 contains
!----------------------------------------------------!
SUBROUTINE read_param
 Namelist/inputdata/T1,T2,H1,EPS,AlphaSteps,AlphaMin,AlphaMax,R1,R2,RSteps,EkinLow,&
 EKinHigh,EkinSteps,p_restart, p_restart_no, p_stop, p_stop_no, &
 RANDOMISE_R, RANDOMISE_A, RANDOMISE_E!

  OPEN(20,file='newinput.dat',status='unknown')
  ! JT: Alpha - angle between initial V & B
  !
  READ(20,nml=inputdata)
  CLOSE(20)

END SUBROUTINE read_param
!----------------------------------------------------!
SUBROUTINE init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
END SUBROUTINE
!----------------------------------------------------!

END MODULE global
