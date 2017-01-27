module M_rkcolR

use global
use M_products
use M_derivsR

implicit none

contains

        subroutine rkcol(mu,gamma,U,B,DUDT,DBDT,UE,dt,T,eta,temperature,rho)
        
        implicit none

        real(num), dimension(3)         :: B,DBDT,UE,flow,VE
        real(num)                       :: mu,gamma,U,DmodBDT,DUDT,dt,modB,vpar,T
        real(num)                       :: nu,beta,betadot,vdot,gammadot,vtot,A_g,B_g,A_beta,B_beta,dnudv,theta,lambda,dlambdadv
        real(num)                       :: Fg,Dg,dDgdv,Fbeta,Dbeta,dDbeta,vtherm,eta,temperature,old_gamma,old_beta,rho
        integer                         :: iter

        vtherm = sqrt(3.0_num*kb*temperature*tempscl/m)/vscl
        !print *, 'vtherm = ', vtherm*vscl, ' temperature = ', temperature*tempscl, ' eta = ', eta
        !vtherm = 1.57E7_num/vscl
        iter = 1

        VE = UE/gamma
        vpar = U/gamma
        modB = sqrt(dot(B,B))
        DmodBDT = 2.0_num*dot(DBDT,B)/modB
        !vtot = sqrt(vpar*vpar + 2.0_num*mu*modB)
        vtot = sqrt(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE))/vscl
        !beta = vpar/vtot
        theta = atan(sqrt(2.0_num*mu*modB)/U)
        if (theta .lt. 0.0_num) theta = pi + theta
	!print *, 'rkcol  theta = ', theta
        beta = cos(theta)
        if (isnan(vtot)) then
          U = vtot*beta*gamma
          mu = gamma*gamma*vtot*vtot*(1.0_num - beta*beta)/(2.0_num*modB)
          print *, 'rkcol vtot is nan. stopping particle orbit.'
          return
        endif

        !betadot = -beta/(2.0_num*mu)*(mu/modB*DmodBDT - 2.0_num*mu/U*DUDT)
        betadot = beta*(1.0_num-beta*beta)*(-1.0_num/(2.0_num*modB)*DmodBDT + 1.0_num/U*DUDT)
        !print *,'rkcol betadot error', (betadot-(beta*(1.0_num-beta*beta)*(-1.0_num/(2.0_num*modB)*DmodBDT + 1.0_num/U*DUDT)))/(beta*(1.0_num-beta*beta)*(-1.0_num/(2.0_num*modB)*DmodBDT + 1.0_num/U*DUDT))
        vdot = DUDT*vpar/vtot + DmodBDT*mu/vtot
        gammadot = (2.0_num*u/(c*c)*DUDT + 2.0_num*mu/(c*c)*DmodBDT)*vscl*vscl/(2.0_num*sqrt(1.0_num + vscl*vscl*U*U/(c*c) + 2.0_num*m*vscl*vscl*mu*modB/(c*c))*sqrt(1.0_num - vscl*vscl*dot(VE,VE)/(c*c)))
        !print *, 'rkcol',u/(c*c),DUDT
        ! note that the above expression for gammadot computes the derivative of gamma given its definition in r_derivs file, assuming that dVe/dt is negligible

        call collision_frequency(dt,nu,lambda,dnudv,dlambdadv,gamma,UE,vtherm,mu,B,modB,U,flow,temperature,eta,rho)

        if (nu .eq. -1.0_num) then
                U = 0.0_num/0.0_num
                return
        endif

        call collision_coefficients(vtot,gamma,beta,nu,lambda,dnudv,dlambdadv,mu,U,B,flow,Fg,Dg,dDgdv,Fbeta,Dbeta,dDbeta)

        A_g = gammadot + Fg + dDgdv
        B_g = sqrt(2.0_num*Dg)
        A_beta = betadot + Fbeta + dDbeta
        B_beta = sqrt(2.0_num*Dbeta)

        old_gamma = gamma
        old_beta = beta
        !if (eta .ne. 0) print *, 'rkcol 1 gamma = ', gamma, ' UE = ',dot(UE,UE)
        call scatter(beta,gamma,A_g,B_g,A_beta,B_beta,dt,nu,T,iter)
        !if (eta .ne. 0) print *, 'rkcol 2 gamma = ', gamma, ' UE = ',dot(UE,UE)
        
        vtot = sqrt(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE))/vscl
        !do while (isnan(vtot)) 
        !  gamma = old_gamma
        !  beta = old_beta
        !  call scatter(beta,gamma,A_g,B_g,A_beta,B_beta,dt,nu,T,iter)
        !  vtot = sqrt(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE))/vscl
        !  !print *, 'rkcol vtot was nan, called scatter again. vtot^2 = ',(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE)), gamma, iter
        !  print *, 'rkcol gamma = ', gamma, A_g*dt, B_g
        !  iter = iter + 1
        !enddo

        vpar = vtot*beta
        mu = gamma*gamma*vtot*vtot*(1.0_num - beta*beta)/(2.0_num*modB)
        U = vpar*gamma

        if (isnan(U)) then
                print *, 'rkcol u is nan. vtot = ', vtot, ' beta = ', beta, ' gamma  = ', gamma, ' A_g = ', A_g, ' B_g = ', B_g, 'dt = ', dt
                print *, 'rkcol vtot is nan. stopping particle orbit.'
                return
                !stop
                !print *, 'rkcol U is nan. stopping particle orbit'
        endif

        end subroutine rkcol

        !************************************************************

        subroutine scatter(beta,gamma,A_g,B_g,A_beta,B_beta,dt,nu,T,iter)

        implicit none

        real(num)               :: beta,beta1,gamma,modB,A_g,B_g,A_beta,B_beta,dt,d_beta,d_g,y,y1,nu,old_vtot
        real(num)               :: dw_g,dw_beta,r1,r2,T,TT,TTT,zeta1,zeta2
        integer                 :: j,rank,ierr,iter

        call random_number(TT)
        call random_number(TTT)
        !print *, 'rkcol check TT = ', TT, TTT
        zeta1 = 0.0_num
        zeta2 = 0.0_num
        do while (zeta1 .eq. zeta2)
                zeta1 = random_normal(TT*T*iter)
                zeta2 = random_normal(TTT*T*T*iter)
        enddo
        !if (iter .gt. 1) print *, 'rkcol random numbers: ', zeta1, zeta2, TT, TT*TT, TT*TT*TT
        dw_g = zeta1*sqrt(dt)*1.0E0_num
        dw_beta = zeta2*sqrt(dt)*1.0E0_num

        d_g = A_g*dt + B_g*dw_g
        d_beta = A_beta*dt + B_beta*dw_beta
        !print *, 'rkcol A_beta = ', A_beta*dt, ' B_ beta = ', B_beta*dw_beta, ' d_beta = ', d_beta

        y = beta + A_beta*dt + B_beta*dw_beta 
        if (y .gt. 1.0_num) then
                y = -y + floor(y) + 1.0_num
        elseif (y .lt. -1.0_num) then
                y = -y + ceiling(y) - 1.0_num
        endif

        gamma = gamma + d_g
        if (gamma .le. 1.0_num) gamma = gamma - 2.0_num*d_g
	if (isnan(gamma)) then
		print *, 'rkcol gamma is nan, stopping. dgamma = ', d_g
		stop
	endif
        
        beta = beta + d_beta
        !beta = beta + (sqrt((1.0_num - y*y)*nu) - B_beta)/2.0_num*(dw_beta*dw_beta - dt)/sqrt(dt)
        !if (abs(beta) .gt. 1000) then
        !        print *,'rkcol line 132'
        !        beta = sign(1.0_num - 1.0E-5_num,beta)
        !elseif (beta .gt. 1.0_num) then
        if (beta .gt. 1.0_num) then
                beta = -beta + floor(beta) + 1.0_num
        elseif (beta .lt. -1.0_num) then
                beta = -beta + ceiling(beta) - 1.0_num
        endif

        end subroutine scatter

        !************************************************************

        subroutine collision_coefficients(vtot,gamma,beta,nu,lambda,dnudv,dlambdadv,mu,U,B,flow,Fg,Dg,dDgdv,Fbeta,Dbeta,dDbeta)

        implicit none

        real(num),dimension(3)          :: flow,B,vel_par
        real(num)                       :: Fg,Dg,dDgdv,Fbeta,Dbeta,dDbeta,nu,mu,U,modB,ne,sg,sgv,sgb,dnudv,vtot,beta,vfr,vtr,b1,b2,b3,afv,adv,gamma
        real(num)                       :: Temp,biggamma,biggammaeff,xx,yy,ndens,epsilon0,lambda,dlambdadv
        real(num)                       :: vtherm

        modB = sqrt(dot(B,B))
        vel_par = U*B/modB
        
        vfr = abs(U/gamma - dot(flow,B)/modB)
        vtr = abs(vtot - vtherm)

        afv = 1.0E-6_num
        adv = 1.0E-6_num

        b1 = 1.0_num
        b2 = 1.0_num
        b3 = 1.0_num

        if (beta .gt. 0) then
                if ((dot(flow,B) .gt. 0) .and. (abs(U/gamma) .lt. sqrt(dot(flow,flow)))) then
                        sgv = 1.0_num
                        sgb = 1.0_num
                else
                        sgv = -1.0_num
                        sgb = -1.0_num
                endif
        else
                if ((dot(flow,B) .lt. 0) .and. (abs(U/gamma) .lt. sqrt(dot(flow,flow)))) then
                        sgv = 1.0_num
                        sgb = -1.0_num
                else
                        sgv = -1.0_num
                        sgb = 1.0_num
                endif
        endif

        ! variables for coulomb collisions in gamma
        Temp = 1.0E6_num
        ndens = 1.0E16_num
	epsilon0 = 8.85E-12_num
        xx = sqrt((gamma - 1.0_num)*m*c*c/(kb*Temp))
        biggamma = 4.0_num*pi*qe*qe*qe*qe*20.0_num*ndens/(me*me)/(m*m*c*c*c*c)/(4.0_num*pi*epsilon0)**2*vtot*vscl
        biggammaeff = biggamma*m*m*1.0E0_num
        Fg = 0.0_num
        Dg = 0.0_num
        !Fg = -biggammaeff/(2.0_num*(gamma - 1.0_num))*(erf(xx) - 2.0_num*xx*(2.0_num*exp(-xx*xx)/sqrt(pi)))
        !Dg = biggammaeff*(erf(xx) - xx*(2.0_num*exp(-xx*xx)/sqrt(pi)))/(2.0_num*xx*xx)
        dDgdv = 0.0_num

	!if (xx .lt. 0.01_num) then
	!	Fg = biggammaeff/(sqrt(pi)*(gamma - 1.0_num))*xx
	!	Dg = 2.0_num*biggammaeff*xx/(3.0_num*sqrt(pi))
	!	!print *, 'rkcol here'
	!endif

        yy = erf(xx) - (erf(xx) - xx*(2.0_num*exp(-xx*xx)/sqrt(pi)))/(2.0_num*xx*xx)
        !Fbeta = beta*biggammaeff/(4.0_num*(gamma - 1.0_num)*(gamma - 1.0_num)*m*m*c*c*c*c)*yy
        !Dbeta = (1.0_num - beta*beta)/(8.0_num*(gamma - 1.0_num)*(gamma - 1.0_num)*m*m*c*c*c*c)*biggammaeff*yy
        !print *, 'rkcol Fg = ', Fg, ' Dg = ', Dg, ' energy difference = ', (gamma - 1.0_num)*m*c*c - (kb*Temp)
        !print *, 'rkcol gamma', gamma, ' energy = ', (gamma-1.0_num)*m*c*c, ' xx = ', xx

        Fbeta = -nu*beta
        Dbeta = nu*(1.0_num - beta*beta)/2.0_num
        !Fbeta = -vtot/lambda*beta
        !Dbeta = vtot/lambda*(1.0_num - beta*beta)/2.0_num
        !Fbeta = 0.0_num
        !Dbeta = 0.0_num
        dDbeta = 0.0_num

        !print *, 'rkcol Fg = ', Fg, ' Dg = ', Dg, ' Fbeta = ', Fbeta, ' Dbeta = ', Dbeta, ' beta = ', beta 

        if (isnan(Fbeta) .or. isnan(Dbeta)) then
                print *, 'rkcol Fbeta or Dbeta is nan. stopping.', vfr,vtot,nu,beta
                stop
        endif


        end subroutine collision_coefficients

        !************************************************************

        subroutine collision_frequency(dt,nu,lambda,dnudv,dlambdadv,gamma,UE,vtherm,mu,B,modB,U,flow,temperature,eta,rho)

        implicit none

        real(num),dimension(3)  :: UE,flow,B
        real(num)               :: nu,gamma,vtot,vtherm,dnudv,mu,modB,U,a1,a2,vfr,vtr,temperature,dt
        real(num)               :: lambda_ei,eta_spitzer,lambda,dlambdadv,alpha,lambda0,eta,kappa,rho

        vtot = sqrt(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE))/vscl
        vfr = abs(U/gamma - dot(flow,B)/modB)
        vtr = abs(vtot - vtherm)
        
        eta_spitzer = 2.4E3_num*(temperature*tempscl)**(-1.5_num)/(Vscl*Lscl*mu0_si)
        !print *, 'eta ', eta, ' eta_spitzer ', eta_spitzer, ' eta scaling ', Vscl*Lscl*mu0_si, ' temperature ', (temperature*tempscl)**(-1.5_num), ' ratio = ', eta/eta_spitzer
        !stop
        !eta_spitzer = eta*1.0E-9_num
        !eta = eta*Vscl*Lscl*mu0_si
        !if (eta .ne. 0) print *, 'rkcol code eta_sp/eta = ', eta_spitzer/eta, eta_spitzer, eta, temperature*tempscl

        if (eta .eq. 0.0_num) then
          nu = 0.0_num
        else 
          alpha = 0.0_num
          kappa = 1.0E-6_num
          !kappa = eta_spitzer/eta
          !kappa = max(eta_spitzer/eta,1.0E-6_num)
          !print *, 'rkcol kappa = ', kappa, ' temperature = ', temperature*tempscl, ' eta ratio ', eta/eta_spitzer
          lambda_ei = 2.0E8_num/Lscl
          !lambda_ei = 2.0E5_num/Lscl
          lambda0 = lambda_ei*kappa
          !lambda0 = lambda_ei*sqrt(kappa)
          lambda = lambda0/(1.0_num + vtot/vtherm)**alpha
          nu = vtot/lambda
        end if
        nu = nu*1.0E0_num
        !nu = min(nu,30.0_num)
        dnudv = 0.0_num
        dlambdadv = 0.0_num
        dt = min(1.0_num/(10.0_num*nu),dt_min)
        !if (eta .ne. 0) print *, 'rkcol 1/nu = ', tscl/nu, ' dt = ', dt*tscl, ' 1/w_g = ', m/(q*sqrt(B(1)**2 + B(2)**2 + B(3)**2)*Bscl), ' ratio ', 1.0_num/kappa
        if (tscl/nu .lt. abs(m/(q*sqrt(B(1)**2 + B(2)**2 + B(3)**2)*Bscl))) then
                print *, 'scattering more frequent than gyration, stopping orbit,scattering timescale', tscl/nu, ' gyroperiod ', abs(m/(q*sqrt(B(1)**2 + B(2)**2 + B(3)**2)*Bscl))
                nu = -1.0_num
        endif

        end subroutine collision_frequency

        !************************************************************
        ! taken from http://jblevins.org/mirror/amiller/random.f90

        FUNCTION random_normal(time) RESULT(fn_val)
        
        ! Adapted from the following Fortran 77 code
        !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
        !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
        !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
        
        !  The function random_normal() returns a normally distributed pseudo-random
        !  number with zero mean and unit variance.
        
        !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
        !  and J.F. Monahan augmented with quadratic bounding curves.
        
        REAL(num) :: fn_val
        
        !     Local variables
        REAL(num)     :: s = 0.449871_num, t = -0.386595_num, a = 0.19600_num, b = 0.25472_num,    &
                    r1 = 0.27597_num, r2 = 0.27846_num, u, v, x, y, q, time
        integer       :: rank, ierr,n,iter
        integer, allocatable :: new_seed(:), old_seed(:)
        
        !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
        call random_seed(size = n)
        allocate(old_seed(n))
        allocate(new_seed(n))
        
        iter = 0
        DO
          !CALL RANDOM_SEED()
          call random_seed(get = old_seed)
          CALL init_random_seed2(time,iter)
          call random_seed(get = new_seed)
          !print *, 'rkcol check seeds, old seed = ', old_seed, ' new seed = ', new_seed
          CALL RANDOM_NUMBER(u)
          CALL RANDOM_NUMBER(v)
               !call mpi_comm_rank(mpi_comm_world, rank, ierr)
               !print *, 'rkcol rank = ',rank, ' random number = ', u,v
               !call mpi_barrier(mpi_comm_world,ierr)
               !!if (rank .eq. 0) stop
          v = 1.7156_num * (v - 0.5_num)
        
        !     Evaluate the quadratic form
          x = u - s
          y = ABS(v) - t
          q = x**2 + y*(a*y - b*x)
            iter = iter+1
        
        !     Accept P if inside inner ellipse
          IF (q < r1) EXIT
        !     Reject P if outside outer ellipse
          IF (q > r2) CYCLE
        !     Reject P if outside acceptance region
          IF (v**2 < -4.0*LOG(u)*u**2) EXIT
        END DO
        !print *,'rkcol number of iterations to get random number: ', iter
        
        !     Return ratio of P's coordinates as the normal deviate
        fn_val = v/u

        deallocate(old_seed)
        deallocate(new_seed)
        !print *, 'rkcol check gaussian random number = ', fn_val
               
        RETURN
        
        END FUNCTION random_normal

        !function random_normal() result(v)
        !
        !real(num)       :: u,v

        !call init_random_seed()
        !call random_number(u)
        !v = u*2.0_num - 1.0_num

        !return 

        !end function

        !************************************************************
        ! taken from
        ! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED, not really sure if this actually works
        subroutine init_random_seed2(time,iter)
            implicit none
            integer, allocatable :: seed(:), old_seed(:)
            integer             :: i, n, un, istat, dt(8), pid,iter
            integer(num)        :: t
            integer             :: ierr, rank
            real(num)           :: time
          
            call random_seed(size = n)
            allocate(seed(n))
            allocate(old_seed(n))
            call random_seed(get = old_seed)
            !print *,'rkcol rank', rank, 'seed = ', seed, ' old_seed = ', old_seed
            ! First try if the OS provides a random number generator
            !open(newunit=un, file="/dev/urandom", access="stream", &
            !     form="unformatted", action="read", status="old", iostat=istat)
            !if (istat == 0) then
            !   read(un) seed
            !   close(un)
            !else
            !   ! Fallback to XOR:ing the current time and pid. The PID is
            !   ! useful in case one launches multiple instances of the same
            !   ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_num * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_num * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_num * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
                        !call mpi_comm_rank(mpi_comm_world, rank, ierr)
                        !print *,'rkcol rank = ', rank,' t = ', t
                        !call mpi_barrier(mpi_comm_world,ierr)
               t = t + time + iter
               !print *,'rkcol t = ', t, ' time = ', time, ' iter = ', iter
               do i = 1, n
                  seed(i) = lcg(t)
                  if (any(old_seed .eq. seed(i))) then
                        seed(i) = 2.0_num*seed(i)
                        call mpi_comm_rank(mpi_comm_world, rank, ierr)
                        print *,'rkcol rank = ', rank, 'here'
                  endif
               end do
                        call mpi_comm_rank(mpi_comm_world, rank, ierr)
                        !print *,'rkcol rank = ', rank,' seed = ', seed
                        !call sleep(2)
            !end if
            call random_seed(put=seed)
            deallocate(seed)
            deallocate(old_seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(num) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_num)
              end if
              s = mod(s * 279470273_num, 4294967291_num)
              lcg = int(mod(s, int(huge(0), num)), kind(0))
            end function lcg
          end subroutine init_random_seed2


end module M_rkcolR
