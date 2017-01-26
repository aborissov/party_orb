module M_rkcolR

use global
use M_products

implicit none

contains

        subroutine rkcol(mu,gamma,U,B,DUDT,DBDT,UE,dt,rho,temperature,flow)
        
        implicit none

        real(num), dimension(3)         :: B,DBDT,UE,flow
        real(num)                       :: mu,gamma,U,DmodBDT,DUDT,dt,modB,vpar,rho,temperature
        real(num)                       :: nu,beta,betadot,vdot,vtot,A_v,B_v,A_beta,B_beta,dnudv
        real(num)                       :: Fv,Dv,dDvdv,Fbeta,Dbeta,dDbeta,vtherm

        !print *,'here'

        vtherm = sqrt(3.0_num*kb*temperature/mp)

        vpar = U/gamma
        modB = sqrt(dot(B,B))
        DmodBDT = 2.0_num*dot(DBDT,B)/modB
        vtot = sqrt(vpar*vpar + 2.0_num*mu*modB)
        !vtot = sqrt(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE))/vscl
        beta = vpar/sqrt(2.0_num*mu*modB + vpar*vpar)

        betadot = beta*beta*beta/2.0_num*(DmodBDT/((1.0_num - beta*beta)*modB) - 2.0_num*DUDT/((1.0_num - beta*beta)*vpar))
        vdot = DUDT*vpar/vtot + DmodBDT*mu/vtot
        
        call collision_frequency(nu,dnudv,gamma,UE,vtherm,mu,B,modB,U,flow)

        call collision_coefficients(vtot,gamma,beta,nu,dnudv,mu,U,B,flow,Fv,Dv,dDvdv,Fbeta,Dbeta,dDbeta,vtherm)

        A_v = vdot + Fv + dDvdv
        B_v = sqrt(2.0_num*Dv)
        A_beta = betadot + Fbeta + dDbeta
        B_beta = sqrt(2.0_num*Dbeta)
        if (isnan(B_beta)) then
                print *, 'rkcol B_beta is nan, stopping. Dbeta = ', Dbeta
                stop
        endif

        !A_beta = -beta*nu
        !B_beta = sqrt(nu*(1.0_num - beta*beta))

        call scatter(beta,vtot,A_v,B_v,A_beta,B_beta,dt)

        gamma = sqrt(1.0_num + vscl*vscl*U*U/(c*c) + 2.0_num*m*vscl*vscl*mu*modB/(m*c*c))/sqrt(1.0_num - vscl*vscl*dot(UE/gamma,UE/gamma)/(c*c))
        vpar = vtot*beta
        U = vpar*gamma
        mu = vtot*vtot*(1.0_num - beta*beta)/(2.0_num*modB)

        if (isnan(U)) then
                print *, 'rkcol U is nan, stopping.'
                stop
        endif

        end subroutine rkcol

        !************************************************************

        subroutine scatter(beta,vtot,A_v,B_v,A_beta,B_beta,dt)

        implicit none

        real(num)               :: beta,beta1,vtot,modB,A_v,B_v,A_beta,B_beta,dt,d_beta,d_v
        real(num)               :: dw_v,dw_beta,r1,r2
        integer                 :: j

        dw_v = random_normal()*sqrt(dt)*1.0E0_num
        dw_beta = random_normal()*sqrt(dt)*1.0E0_num

        !A_v = 0.0_num
        !B_v = 0.0_num
        !A_beta = 0.0_num
        !A_beta = 0.0_num

        d_v = A_v*dt + B_v*dw_v
        d_beta = A_beta*dt + B_beta*dw_beta

        vtot = abs(vtot + d_v)
        if (vtot .lt. 0.0_num) vtot = vtot - 2.0_num*d_v
        !beta = beta + d_beta
        !if (abs(beta) .gt. 1.0_num) beta = beta - 2.0_num*d_beta
        beta1 = modulo(beta + d_beta,1.0_num)
        if (beta1*beta .lt. 0.0_num) then
                beta = beta1 - 1.0_num
        else
                beta = beta1
        endif
        ! Note: this scheme is only order 1/2 (I think), to do order 1 see Eduard's code and wikipedia article (or higher orders)
        
        if (isnan(vtot)) then
                print *, 'rkcol vtot is nan, stopping.'
                stop
        endif

        end subroutine scatter

        !************************************************************

        subroutine collision_coefficients(vtot,gamma,beta,nu,dnudv,mu,U,B,flow,Fv,Dv,dDvdv,Fbeta,Dbeta,dDbeta,vtherm)

        implicit none

        real(num),dimension(3)          :: flow,B,vel_par
        real(num)                       :: Fv,Dv,dDvdv,Fbeta,Dbeta,dDbeta,nu,mu,U,modB,ne,sg,sgv,sgb,dnudv,vtot,beta,vfr,vtr,b1,b2,b3,afv,adv,gamma,vtherm

        if (abs(beta) .gt. 1.0_num) then
                print *, 'beta magnitude greater than 1 stopping soon'
        endif

        modB = sqrt(dot(B,B))
        vel_par = U*B/modB
        
        vfr = abs(U/gamma - dot(flow,B)/modB)
        vtr = abs(vtot - vtherm)

        afv = 1.0_num
        adv = 1.0_num

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

        !Fv = sgv*10.0_num
        !Dv = 1.0_num
        !dDvdv = 0.0_num
        !Fbeta = sgb*10.0_num
        !Dbeta = 1.0_num
        Fv = afv*sgv*nu*vfr**b1*vtr**b2
        Dv = adv*nu*vtr**b3
        dDvdv = adv*sgv*(dnudv*vtr**b3 + sign(b3*vtr**(b3 - 1.0_num),vtot - vtherm))
        Fbeta = 1.0E0_num*sign(abs(beta*nu),sgb)
        Dbeta = nu*(1.0_num - beta*beta)/2.0_num
        !dDbeta = sg*dnudv*(1.0_num - beta*beta)/2.0_num
        dDbeta = 0.0_num
        
        if (isnan(Dbeta)) then
                print *, 'rkcol Dbeta is nan, stopping. beta = ', beta, ' nu = ', nu
                stop
        endif

        if (Dbeta .lt. 0.0_num) then
                print *, 'rkcol Dbeta = ', Dbeta, ' beta = ', beta
                stop
        endif

        end subroutine collision_coefficients

        !************************************************************

        subroutine collision_frequency(nu,dnudv,gamma,UE,vtherm,mu,B,modB,U,flow)
        !subroutine collision_frequency(nu,dnudv,vtot,vtherm)

        implicit none

        real(num),dimension(3)  :: UE,flow,B
        real(num)               :: nu,nu0,gamma,vtot,vtherm,dnudv,mu,modB,U,a1,a2,vfr,vtr,rho,temperature,vpar


        !vtot = sqrt(c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE))/vscl
        vpar = U/gamma
        vtot = sqrt(vpar*vpar + 2.0_num*mu*modB)
        vfr = abs(U/gamma - dot(flow,B)/modB)
        vtr = abs(vtot - vtherm)

        a1 = 0.0_num
        a2 = 1.0_num
        
        nu = nu0*exp(-a1*vfr - a2*vtr)
        if (isnan(nu)) then
                print *, 'rkcol nu is nan, stopping.'
                print *,c*c - c*c/(gamma*gamma) - vscl*vscl*dot(UE,UE)
                print *, dot(UE,UE)
                stop
        endif
        if (vtot .gt. vtherm) then
                dnudv = -nu
        else
                dnudv = nu0*exp(vtot - vtherm)
        endif
        ! WARNING: dnudv calculated assuming a1 = 0. Need to change (calculate more derivatives) when a1 != 0.

        end subroutine collision_frequency

        !************************************************************
        ! taken from http://jblevins.org/mirror/amiller/random.f90

        FUNCTION random_normal() RESULT(fn_val)
        
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
                    r1 = 0.27597_num, r2 = 0.27846_num, u, v, x, y, q
        
        !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
        
        DO
          !CALL RANDOM_SEED()
          !CALL initialize_random_seed()
          CALL RANDOM_NUMBER(u)
          CALL RANDOM_NUMBER(v)
          v = 1.7156_num * (v - 0.5_num)
        
        !     Evaluate the quadratic form
          x = u - s
          y = ABS(v) - t
          q = x**2 + y*(a*y - b*x)
        
        !     Accept P if inside inner ellipse
          IF (q < r1) EXIT
        !     Reject P if outside outer ellipse
          IF (q > r2) CYCLE
        !     Reject P if outside acceptance region
          IF (v**2 < -4.0*LOG(u)*u**2) EXIT
        END DO
        
        !     Return ratio of P's coordinates as the normal deviate
        fn_val = v/u
        RETURN
        
        END FUNCTION random_normal

        !function random_normal() result(v)
        !
        !real(num)       :: u,v

        !call initialize_random_seed()
        !call random_number(u)
        !v = u*2.0_num - 1.0_num

        !return 

        !end function

        !************************************************************
        ! taken from
        ! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED, not really sure if this actually works
        !subroutine initialize_random_seed()
        !    implicit none
        !    integer, allocatable :: seed(:)
        !    integer :: i, n, un, istat, dt(8), pid
        !    integer(num) :: t
        !  
        !    call random_seed(size = n)
        !    allocate(seed(n))
        !    ! First try if the OS provides a random number generator
        !    !open(newunit=un, file="/dev/urandom", access="stream", &
        !    !     form="unformatted", action="read", status="old", iostat=istat)
        !    !if (istat == 0) then
        !    !   read(un) seed
        !    !   close(un)
        !    !else
        !       ! Fallback to XOR:ing the current time and pid. The PID is
        !       ! useful in case one launches multiple instances of the same
        !       ! program in parallel.
        !       call system_clock(t)
        !       if (t == 0) then
        !          call date_and_time(values=dt)
        !          t = (dt(1) - 1970) * 365_num * 24 * 60 * 60 * 1000 &
        !               + dt(2) * 31_num * 24 * 60 * 60 * 1000 &
        !               + dt(3) * 24_num * 60 * 60 * 1000 &
        !               + dt(5) * 60 * 60 * 1000 &
        !               + dt(6) * 60 * 1000 + dt(7) * 1000 &
        !               + dt(8)
        !       end if
        !       pid = getpid()
        !       t = ieor(t, int(pid, kind(t)))
        !       do i = 1, n
        !          seed(i) = lcg(t)
        !       end do
        !    !end if
        !    call random_seed(put=seed)
        !  contains
        !    ! This simple PRNG might not be good enough for real work, but is
        !    ! sufficient for seeding a better PRNG.
        !    function lcg(s)
        !      integer :: lcg
        !      integer(num) :: s
        !      if (s == 0) then
        !         s = 104729
        !      else
        !         s = mod(s, 4294967296_num)
        !      end if
        !      s = mod(s * 279470273_num, 4294967291_num)
        !      lcg = int(mod(s, int(huge(0), num)), kind(0))
        !    end function lcg
        !  end subroutine initalize_random_seed


end module M_rkcolR
