! *****************************************************************************
module module_friction
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Sources module (friction)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
function friction_coef(h, Rh, q, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the Fanning friction coefficient depending on the friction law 
! (<=> CF in TELEMAC2D)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: frict_law, kappa, nu_w, eps, eps_num, g, write_listing
implicit none
real(kind=dp), intent(in)    :: h
real(kind=dp), intent(in)    :: Rh
real(kind=dp), intent(in)    :: fc
real(kind=dp), intent(in)    :: q
real(kind=dp)                :: friction_coef
real(kind=dp)                :: aux, aux1, aux2, r, Re
! Constant for the Colebrook-White, Haaland and Nikuradse
real(kind=dp), parameter     :: C3 = 12.d0 ! 14.8d0 for Barr equivalent
! Colebrook-White iteration variables
integer                      :: iter
integer, parameter           :: fp_nitermax = 20
real(kind=dp), parameter     :: fp_prec = 1.d-8
real(kind=dp)                :: x, x0, fp_residual
real(kind=dp)                :: a1, a2, usr
! -----------------------------------------------------------------------------
if (q.ne.0 .and. h.ge.eps .and. Rh.ge.eps_num) then
  if(frict_law.eq.0) then
    friction_coef = 0.d0
  ! *********
  ! Fanning :
  ! *********
  elseif(frict_law.eq.1) then
    friction_coef = fc
  ! *********
    ! Chezy :
  ! *********  
  elseif (frict_law.eq.2) then
    friction_coef = 2d0*g/(fc**2.d0)
  ! ***********
  ! Strickler :
  ! ***********
  elseif (frict_law.eq.3) then
    friction_coef = 2d0*g/(fc**2.d0*Rh**(1.d0/3.d0))
  ! ***********
  ! Nikuradse :
  ! ***********  
  elseif (frict_law.eq.4) then
    aux = max(1.001d0, C3*Rh/fc)
    friction_coef = 2.d0/(log(aux)/kappa)**2.d0
  ! *****************
  ! Colebrook-White :
  ! *****************
  elseif (frict_law.eq.5) then
    Re = 4d0*abs(q/h)*Rh/nu_w ! hydraulic Reynolds number
    Re = max(Re, eps_num)
    ! initital guess with Haaland
    x0 = -3.6d0*log10( 6.9d0/Re + (fc/(C3*Rh))**1.11d0 )
    x0 = max(x0, 1d0)
    ! iterative solution of Colebrook-White equation
    iter = 1
    fp_residual = 1.d0
    do while (iter.le.fp_nitermax .and. fp_residual.ge.fp_prec)
      x = -4.d0*log10( (2.51d0/(2d0*Re))*x0 + fc/(C3*Rh) )
      x = max(x, 1d0)
      fp_residual = abs(x-x0)
      x0 = x
      iter = iter + 1
    enddo
    friction_coef = x**(-2d0)
  ! *********
  ! Haaland :
  ! *********
  elseif (frict_law.eq.6) then
    Re = 4d0*abs(q/h)*Rh/nu_w ! hydraulic Reynolds number
    Re = max(Re, eps_num)
    friction_coef = (-3.6d0*log10(6.9d0/Re+(fc/(C3*Rh))**1.11d0))**(-2.d0)
  ! ******
  ! Barr :
  ! ******
  elseif (frict_law.eq.7) then
    Re = abs(q/h)*Rh/nu_w ! Reynolds number
    Re = max(Re, eps_num)
    r = fc/Rh
    aux1 = 1.1295d0*log10(Re/1.75d0)
    aux2 = Re*(1d0 + (Re**0.52d0*r**0.7d0)/37.22d0)
    friction_coef = (-4.d0*log10(r/14.8d0 + aux1/aux2 ))**(-2.d0)
  ! **********
  ! Bathurst :
  ! **********
  elseif (frict_law.eq.8) then
    friction_coef = (-1.987d0*2.d0*log10(fc/(Rh*5.15d0)))**(-2.d0)
  ! **********
  ! Machiels :
  ! **********
  elseif (frict_law.eq.9) then
    Re = abs(q/h)*Rh/nu_w ! Reynolds number
    Re = max(Re, eps_num)
    r = fc/Rh
    if (r.le.0.05d0) then
      ! Barr formula
      aux1 = 1.1295d0*log10(Re/1.75d0)
      aux2 = Re*(1d0 + (Re**0.52d0*r**0.7d0)/37.22d0)
      friction_coef = (-4.d0*log10(r/14.8d0 + aux1/aux2 ))**(-2.d0)
    else if (0.05d0.lt.r .and. r.le.0.15d0) then
      ! Machiels formula (transition zone)
      x = 1469.76d0*(r**3d0) -382.83d0*(r**2d0) + 9.89d0*r + 5.22d0
      friction_coef = 1./(2.*x)**2
    else
      ! Bathurst formula
      friction_coef = (-1.987d0*2.d0*log10(fc/(Rh*5.15d0)))**(-2.d0)
    end if
  ! **********
  ! Ferguson :
  ! **********
  elseif (frict_law.eq.10) then
    a1 = 6.5d0
    a2 = 2.5d0
    usr = 3d0*Rh/fc ! D84 = ks/3
    aux = a1*a2*usr/sqrt(a1**2d0 + (a2**2d0)*(usr**(5.d0/3.d0)))
    friction_coef = (aux/sqrt(2d0))**(-2.d0)
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown BOTTOM FRICTION LAW '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
else
  friction_coef = 0d0
endif ! q!=0, h>0, Rh>0
! -----------------------------------------------------------------------------
! WARNING for high friction coefficient
if ((friction_coef.gt.1.d1).and.(h.gt.1d-2).and.(write_listing)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING: HIGH FRICTION COEFFICIENT: ', friction_coef
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  friction_coef = 1.d1
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
subroutine source_friction(h, Rh, q, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! time splitting scheme (i.e. friction solved after advection)
! time schemes for friction: explicit / semi-implicit / implicit
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: frict_law, frict_scheme, dt, eps, eps_num, kappa, nu_w
implicit none
real(kind=dp), intent(in)    :: h
real(kind=dp), intent(in)    :: Rh
real(kind=dp), intent(in)    :: fc
real(kind=dp), intent(inout) :: q
real(kind=dp)                :: Cf, af
! -----------------------------------------------------------------------------
if ((frict_law.ne.0).and.(fc.gt.0)) then
  ! no friction if dry cell or near zero hydraulic radius
  if (q.ne.0 .and. h.ge.eps .and. Rh.ge.eps_num) then
    ! compute friction coefficient and af
    Cf = friction_coef(h, Rh, q, fc)
    af = dt*Cf/(2d0*h*Rh)
    ! apply friction term
    ! *********
    ! Explicit:
    ! *********
    if(frict_scheme.eq.1) then
      q = q*(1 - af*abs(q))
    ! **************  
    ! Semi-implicit:
    ! **************
    elseif(frict_scheme.eq.2) then
      q = q/(1.d0 + af*abs(q))
    ! *********
    ! Implicit:
    ! *********  
    elseif(frict_scheme.eq.3) then
      if(q.gt.0) then
        q = (-1.d0 + dsqrt(1.d0 + 4.d0*af*q))/(2.d0*af)
      else
        q = ( 1.d0 - dsqrt(1.d0 - 4.d0*af*q))/(2.d0*af)
      endif
    else
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unknown BOTTOM FRICTION SCHEME '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
  endif 
endif ! q!=0, h>0, Rh>0
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
function friction_velocity(h, u, L, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! friction velocity
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho
implicit none
real(kind=dp), intent(in) :: h, u, L, fc
real(kind=dp)             :: friction_velocity
! -----------------------------------------------------------------------------
friction_velocity = sqrt(tau(h, u, L, 0d0, fc)/rho)
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function friction_slope(h, Rh, u, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! friction slope J
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, eps, g
implicit none
real(kind=dp), intent(in) :: h, Rh, u, fc
real(kind=dp)             :: Cf
real(kind=dp)             :: friction_slope
! -----------------------------------------------------------------------------
! Fanning friction coefficient
Cf = friction_coef(h, Rh, h*u, fc)
if (Cf.gt.0d0) then
  ! Darcy-Weisbach formula
  friction_slope = (Cf*u*abs(u))/(2d0*g*Rh)
else
  friction_slope = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function hydraulic_radius(h, L)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! hydraulic radius
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: variable_width, eps, Rh_approx, Lcst
implicit none
real(kind=dp), intent(in) :: h, L
real(kind=dp)             :: hydraulic_radius
! -----------------------------------------------------------------------------
if (variable_width.eqv..True.) then
  hydraulic_radius = max(h*L/(2d0*h + L), eps)
else
  if (Rh_approx.eq.0) then
    hydraulic_radius = max(h, eps)
  else
    hydraulic_radius = max((h*Lcst)/(Lcst + 2.d0*h), eps)
  endif
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function tau(h, u, L, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Bottom shear stress: Tau = rho*g*h*J
! where J is the friction slope.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, eps, g
implicit none
real(kind=dp), intent(in) :: h, u, L, Jeq, fc
real(kind=dp)             :: Rh, J, tau
! -----------------------------------------------------------------------------
! Check if not dry cell
if (h.ge.eps) then
  Rh = hydraulic_radius(h, L)
  ! friction slope
  if (Jeq.eq.0d0) then
    J = friction_slope(h, Rh, u, fc)
  else
    J = Jeq
  endif
  tau = rho*g*Rh*J ! rho*g*h*J in 1D
! If dry cell, no shear stress
else
  tau = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



end module
! *****************************************************************************