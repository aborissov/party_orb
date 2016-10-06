MODULE test_fields
  
  USE global
  USE M_products

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: TESTFIELDS

  CONTAINS 
!------------------------------------------------------------------------------!    
SUBROUTINE TESTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)
! testfields - routines to test how the relativistic and nonrelativistic versions conform.

 REAL(num), DIMENSION(3),INTENT(OUT)	:: B,E
 REAL(num), DIMENSION(3),INTENT(OUT)	:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3),INTENT(IN)	:: R
 REAL(num), DIMENSION(3)  		:: Vf, dBedx,dBedy,dBedz
 REAL(num), INTENT(IN) 			:: T
 REAL(num)				:: ntemp, ontemp,x,y,z,L,vyinf,Lvx,Lvy,Letax,Letay,eta0,lambda
  

 ntemp=64.0_num
 ontemp=1.0_num/ntemp 

 x = R(1)
 y = R(2)
 z = R(3)
 L = 10.0_num    !Since length is already normalised this is ratio between L (see Paul Wood's thesis) and Lscl. 
 Letax = 1.666E0_num*L
 Letay = 1.666E-1_num*L
 eta0 = 0.00029_num
 Lvx = L
 Lvy = 1.0_num
 vyinf = 0.6E-1_num
 lambda = 1.0_num
  
  !full B=Be+B1
 B(1)=lambda*sinh(y)/(cosh(y) + exp(-x*x/(L*L)))
 B(2)=2.0_num*lambda/(L*L)*(x*exp(-x*x/(L*L)))/(cosh(y) + exp(-x*x/(L*L)))
 B(3)=1.0_num

 !print*, T, T1, tau
 
 ! non-dimensionalised
 !B=B/Bscl

 ! no velocity?
 Vf(1)=vyinf * Lvx * tanh(x / Lvx) / Lvy / cosh(y / Lvy) ** 2
 Vf(2)=-vyinf * tanh(y / Lvy) / cosh(x / Lvx) ** 2
 Vf(3)=0.0_num

 !Vf = Vf*tscl

 ! Electric fields
 E(1)=vyinf*tanh(y/Lvy)/cosh(x/Lvx)**2
 E(2)=vyinf * Lvx * tanh(x / Lvx) / Lvy / cosh(y / Lvy) ** 2
 E(3)=eta0 * (0.2D1 * lambda * exp(-x ** 2 / L ** 2) / L ** 2 / (                       &
     cosh(y) + exp(-x ** 2 / L ** 2)) - 0.4D1 * lambda * x ** 2 * exp(-                 &
     x ** 2 / L ** 2) / L ** 4 / (cosh(y) + exp(-x ** 2 / L ** 2)) +                    &
     0.4D1 * lambda * x ** 2 * exp(-x ** 2 / L ** 2) ** 2 / L ** 4 / (                  &
     cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 - lambda * cosh(y) / (cosh(                  &
     y) + exp(-x ** 2 / L ** 2)) + lambda * sinh(y) ** 2 / (cosh(y) +                   &
      exp(-x ** 2 / L ** 2)) ** 2) / cosh(x / Letax) ** 2 / cosh(Y /                    &
     Letay) ** 2 - 0.2D1 * vyinf * Lvx * tanh(x / Lvx) * lambda * x * exp               &
     (-x ** 2 / L ** 2) / Lvy / cosh(y / Lvy) ** 2 / L ** 2 / (cosh(y                   &
     ) + exp(-x ** 2 / L ** 2)) - vyinf * tanh(y / Lvy) * lambda * sinh                 &
     (y) / cosh(x / Lvx) ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2))
 
 !print*, B1*a
 !print*, 'testfields E dot B = ', E(1)*B(1) + E(2)*B(2) + E(3)*B(3)
 !print*, Escl
 !print*, E/Escl
 !STOP
! non-dimensionalised 
 !E=E/Escl
 !E = 0.0_num

 !---------DERIVATIVES-------!
 ! calculate the spatial derivatives of individual components - in this case, only the x derivatives matter, but make sure only the z component is assigned..
 dBdx(1)=0.2D1 * lambda * sinh(y) * x * exp(-x ** 2 / L ** 2) / (cosh    &
     (y) + exp(-x ** 2 / L ** 2)) ** 2 / L ** 2
 dBdx(2)=0.2D1 * lambda * exp(-x ** 2 / L ** 2) / L ** 2 / (cosh(y) +                            &
      exp(-x ** 2 / L ** 2)) - 0.4D1 * lambda * x ** 2 * exp(-x ** 2 /                  &
     L ** 2) / L ** 4 / (cosh(y) + exp(-x ** 2 / L ** 2)) + 0.4D1 *                     &
     lambda * x ** 2 * exp(-x ** 2 / L ** 2) ** 2 / L ** 4 / (cosh(y) + exp             &
     (-x ** 2 / L ** 2)) ** 2
 dBdx(3)=0.0_num
 
 dBdy(1)=lambda * cosh(y) / (cosh(y) + exp(-x ** 2 / L ** 2)) - lambda   &
      * sinh(y) ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2
 dBdy(2)=0.2D1 * lambda * sinh(y) * x * exp(-x ** 2 / L ** 2) / (cosh    &
     (y) + exp(-x ** 2 / L ** 2)) ** 2 / L ** 2
 dBdy(3)=0.0_num

 dBdz(1)=0.0_num
 dBdz(2)=0.0_num
 dBdz(3)=0.0_num
 
 ! not sure quite how to deal with time dependence in here yet. Or normalisation.
 !DBDX = dBdx/Bscl		! lengths normalised, fields NOT 
 !DBDY = dBdy/Bscl
 !DBDZ = dBdz/Bscl

 
 !dBdt= Bf/tau
 dBdt= 0.0_num		
! non-dimensionalised
 DBDT= dBdt*Tscl/Bscl

! E is purely in z, therefore derivatives only vary in z too
 dEdx(1)=-0.2D1 * vyinf * tanh(y / Lvy) * sinh(x / Lvx) / cosh(x / Lvx   &
     ) ** 3 / Lvx
 dEdx(2)=vyinf * (0.1D1 - tanh(x / Lvx) ** 2) / Lvy / cosh(y / Lvy) ** 2
 dEdx(3)=-0.2D1 * eta0 * (0.2D1 * lambda * exp(-x ** 2 / L ** 2) / L                     &
      ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)) - 0.4D1 * lambda * x **          &
      2 * exp(-x ** 2 / L ** 2) / L ** 4 / (cosh(y) + exp(-x ** 2 / L           &
     ** 2)) + 0.4D1 * lambda * x ** 2 * exp(-x ** 2 / L ** 2) ** 2 / L          &
     ** 4 / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 - lambda * cosh(y)           &
      / (cosh(y) + exp(-x ** 2 / L ** 2)) + lambda * sinh(y) ** 2 / (           &
     cosh(y) + exp(-x ** 2 / L ** 2)) ** 2) * sinh(x / Letax) / cosh(x          &
      / Letax) ** 3 / cosh(Y / Letay) ** 2 / Letax + eta0 * (-0.12D2 *          &
     lambda * x * exp(-x ** 2 / L ** 2) / L ** 4 / (cosh(y) + exp(-x            &
     ** 2 / L ** 2)) + 0.12D2 * lambda * exp(-x ** 2 / L ** 2) ** 2 * x         &       
     / L ** 4 / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 + 0.8D1 * lambda         &
      * x ** 3 * exp(-x ** 2 / L ** 2) / L ** 6 / (cosh(y) + exp(-x             &
     ** 2 / L ** 2)) - 0.24D2 * lambda * x ** 3 * exp(-x ** 2 / L ** 2)         &       
     ** 2 / L ** 6 / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 + 0.16D2 *          &
      lambda * x ** 3 * exp(-x ** 2 / L ** 2) ** 3 / L ** 6 / (cosh(y)          &
      + exp(-x ** 2 / L ** 2)) ** 3 - 0.2D1 * lambda * cosh(y) * x *            &
     exp(-x ** 2 / L ** 2) / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 / L         &
      ** 2 + 0.4D1 * lambda * sinh(y) ** 2 * x * exp(-x ** 2 / L ** 2)          &
      / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 3 / L ** 2) / cosh(x /             &
     Letax) ** 2 / cosh(Y / Letay) ** 2 - 0.2D1 * vyinf * (0.1D1 - tanh(x       &       
      / Lvx) ** 2) * lambda * x * exp(-x ** 2 / L ** 2) / Lvy / cosh(y          &
      / Lvy) ** 2 / L ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)) - 0.2D1          &
      * vyinf * Lvx * tanh(x / Lvx) * lambda * exp(-x ** 2 / L ** 2) /          &
     Lvy / cosh(y / Lvy) ** 2 / L ** 2 / (cosh(y) + exp(-x ** 2 / L             &
     ** 2)) + 0.4D1 * vyinf * Lvx * tanh(x / Lvx) * lambda * x ** 2 * exp       &       
     (-x ** 2 / L ** 2) / Lvy / cosh(y / Lvy) ** 2 / L ** 4 / (cosh(            &
     y) + exp(-x ** 2 / L ** 2)) - 0.4D1 * vyinf * Lvx * tanh(x / Lvx)          &
     * lambda * x ** 2 * exp(-x ** 2 / L ** 2) ** 2 / Lvy / cosh(y /            &
     Lvy) ** 2 / L ** 4 / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 + 0.2D1        &
      * vyinf * tanh(y / Lvy) * lambda * sinh(y) * sinh(x / Lvx) /              &
     cosh(x / Lvx) ** 3 / (cosh(y) + exp(-x ** 2 / L ** 2)) / Lvx - 0.2D1       &
      * vyinf * tanh(y / Lvy) * lambda * sinh(y) * x * exp(-x ** 2              &
     / L ** 2) / cosh(x / Lvx) ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)          &
     ) ** 2 / L ** 2

 !DEDX=DEDX/Escl			! again, all lengths in problem ALREADY normalised, so 1/Escl
 
 dEdy(1)=vyinf * (0.1D1 - tanh(y / Lvy) ** 2) / Lvy / cosh(x / Lvx) ** 2
 dEdy(2)=0.2D1 * vyinf * Lvx * tanh(x / Lvx) * sinh(y / Lvy) / Lvy **    &
      2 / cosh(y / Lvy) ** 3
 dEdy(3)=eta0 * (-0.2D1 * lambda * exp(-x ** 2 / L ** 2) * sinh(y)                               &
     / L ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)) ** 2 + 0.4D1 * lambda                 &
      * x ** 2 * exp(-x ** 2 / L ** 2) * sinh(y) / L ** 4 / (cosh(y)                    &
      + exp(-x ** 2 / L ** 2)) ** 2 - 0.8D1 * lambda * x ** 2 * exp(-x                  &
     ** 2 / L ** 2) ** 2 * sinh(y) / L ** 4 / (cosh(y) + exp(-x ** 2                    &
     / L ** 2)) ** 3 - lambda * sinh(y) / (cosh(y) + exp(-x ** 2 / L                    &
     ** 2)) + 0.3D1 * lambda * cosh(y) * sinh(y) / (cosh(y) + exp(-x                    &
      ** 2 / L ** 2)) ** 2 - 0.2D1 * lambda * sinh(y) ** 3 / (cosh(y)                   &
      + exp(-x ** 2 / L ** 2)) ** 3) / cosh(x / Letax) ** 2 / cosh(Y /                  &
     Letay) ** 2 + 0.4D1 * vyinf * Lvx * tanh(x / Lvx) * lambda * x *                   &
     exp(-x ** 2 / L ** 2) * sinh(y / Lvy) / Lvy ** 2 / cosh(y / Lvy)                   &
     ** 3 / L ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)) + 0.2D1 * vyinf                  &
      * Lvx * tanh(x / Lvx) * lambda * x * exp(-x ** 2 / L ** 2) * sinh                 &
     (y) / Lvy / cosh(y / Lvy) ** 2 / L ** 2 / (cosh(y) + exp(-x **                     &
     2 / L ** 2)) ** 2 - vyinf * (0.1D1 - tanh(y / Lvy) ** 2) * lambda                  &
      * sinh(y) / Lvy / cosh(x / Lvx) ** 2 / (cosh(y) + exp(-x ** 2 /                   &
      L ** 2)) - vyinf * tanh(y / Lvy) * lambda * cosh(y) / cosh(x /                    &
     Lvx) ** 2 / (cosh(y) + exp(-x ** 2 / L ** 2)) + vyinf * tanh(y /                   &
      Lvy) * lambda * sinh(y) ** 2 / cosh(x / Lvx) ** 2 / (cosh(y) +                    &
     exp(-x ** 2 / L ** 2)) ** 2
 !DEDY=DEDY/Escl

 dEdz(1)=0.0_num
 dEdz(2)=0.0_num
 dEdz(3)=0.0_num
 !DEDZ=dEdz/Escl

! no time dependence
 dEdt(1)=0.0_num
 dEdt(2)=0.0_num
 dEdt(3)=0.0_num
 DEDT= dEdt*Tscl/Escl

END SUBROUTINE TESTFIELDS

END MODULE test_fields

