! *****************************************************************************
module module_bedload
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Bedload module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g, eps
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                               BEDLOAD FLUXES
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine compute_bedload_solidfluxes(zf, zb, h, u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes bedload rate qs [m^2/s]
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: qsed, Lx, porosite, fric, bedload_coef, &
                       bed_fric_coef, eps_b, eps_b_damping, pi
implicit none
real(kind=dp), intent(in) :: zb (0:nx+1)
real(kind=dp), intent(in) :: zf (0:nx+1)
real(kind=dp), intent(in) :: h  (0:nx+1)
real(kind=dp), intent(in) :: u  (0:nx+1)
integer                   :: i
real(kind=dp)             :: b, rb, fc
! -----------------------------------------------------------------------------
! loop on FV cells
!$OMP parallel do  &
!$OMP private(b, rb, fc) &
!$OMP shared(qsed, h, u, Lx, zb, zf, fric, eps_b, eps_b_damping, &
!$OMP        bedload_coef, bed_fric_coef, porosite)
do i=0,nx+1
  ! solid flux computation
  fc = fric(i)*bed_fric_coef
  qsed(i) = qs(h(i), u(i), Lx(i), 0d0, fc)

  ! non-erodible bed layer damping zone
  if (eps_b_damping.eq.1) then
    b = max(zb(i)-zf(i), 0d0)
    if (b.le.eps_b) then
      rb = b/eps_b
      qsed(i) = qsed(i)*damping_coef(rb)
    endif
  endif

  ! qsed = xi*qs with xi = 1/(1-p)
  qsed(i) = bedload_coef*qsed(i)/(1.d0 - porosite)
enddo
!$OMP end parallel do
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
elemental function damping_coef(r)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Coefficient to reduce solid flux if zb < zf + eps_b
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps_num, pi
implicit none
real(kind=dp), intent(in) :: r
real(kind=dp)             :: damping_coef
! -----------------------------------------------------------------------------
if (r.ge.1d0) then
  damping_coef = 1.d0
else if ((r.gt.eps_num).and.(r.lt.1d0)) then
  damping_coef = 0.5d0*(1d0 + sin(pi*(r-0.5d0)))
else 
  damping_coef = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes solid discharge: qs [m^2/s]
! 
! Friction slope equilibrium option:
! - if Jeq = 0 => Tau = rho*g*h*J
! - if Jeq > 0 => Tau = rho*g*h*Jeq (for dirichlet bnd conditions)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: bedload_formula
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: qs
! -----------------------------------------------------------------------------
! check if not dry cell, else zero flux
if (h.gt.eps) then
  ! ***********
  ! Grass model
  ! ***********
  if (bedload_formula.eq.1) then
    qs = qs_grass(u)
  ! ********************
  ! Grass modified model
  ! ********************
  elseif (bedload_formula.eq.2) then
    qs = qs_grass_mod(h, u)
  ! *********
  ! MPM model (Meyer-Peter, E. and Muller, 1948)
  ! *********
  elseif (bedload_formula.eq.3) then
    qs = qs_mpm(h, u, Lx, Jeq, fc)
  ! *****************
  ! Standardizd model (Zanke U. and Roland A., 2020)
  ! 1: MPM modified
  ! *****************
  elseif (bedload_formula.eq.4) then
    qs = qs_std(h, u, Lx, Jeq, fc)
  ! ************
  ! Lefort model (Courlis)
  ! ************
  elseif (bedload_formula.eq.5) then
    qs = qs_Lefort(h, u, Lx, Jeq, fc)
  ! ******************
  ! Recking 2013 model (Courlis)
  ! ******************
  elseif (bedload_formula.eq.6) then
    qs = qs_Reck2013(h, u, Lx, Jeq, fc)
  ! ******************
  ! Recking 2015 model (Courlis)
  ! ******************
  elseif (bedload_formula.eq.7) then
    qs = qs_Reck2015(h, u, Lx, 1, 1, Jeq, fc)
  ! **************************
  ! Engelund and Hansen (1967)
  ! **************************
  elseif (bedload_formula.eq.8) then
    qs = qs_eh(h, u, Lx, Jeq, fc)
  ! **************************
  ! Van Rijn (1984)
  ! **************************
  elseif (bedload_formula.eq.9) then
    qs = qs_vr(h, u, Lx, Jeq, fc)
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown BEDLOAD FORMULA '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
else
  qs = 0.d0
endif     
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_grass(u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Grass (1981)
! Ref: A.J. Grass. Sediment transport by waves and currents. 
!      SERC, London, report fl-29 edition, 1981.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Ag, mg
implicit none
real(kind=dp), intent(in) :: u
real(kind=dp)             :: qs_grass
! -----------------------------------------------------------------------------
qs_grass = Ag*u*abs(u)**(mg-1.d0)
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_grass_mod(h, u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Ag
implicit none
real(kind=dp), intent(in) :: h, u
real(kind=dp)             :: qs_grass_mod
! -----------------------------------------------------------------------------
qs_grass_mod = Ag*h*u*abs(u)**2.d0
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_mpm(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Meyer-Peter & Muller (1948)
! Ref: Meyer-Peter E. Muller R. Formulas for bed-load transport. 
!      Report on 2nd meeting on international associatio on hydraulic structures
!      research, pages 39–64, Stockholm 1948.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: g, rho, rhos, d50, tau_c, tau_adim
use module_friction, only: tau
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: qs_mpm, taus, coef, R
!skin friction correction
real(kind=dp)             :: Kr, corr, d90
logical, parameter        :: skin_friction = .False.
! -----------------------------------------------------------------------------
! skin friction correction
if (skin_friction.eqv..True.) then
  d90 = 2.1d0*d50
  Kr = 26d0/(d90**(1d0/6d0))
  corr = (fc/Kr)**(3d0/2d0)
else
  corr = 1.d0
endif
! non-dimensional shear stress
taus = tau(h, u, Lx, Jeq, fc)*tau_adim*corr
! qs*
qs_mpm = 8.d0*max((taus-tau_c), 0.d0)**(3.d0/2.d0)
! qs
R = (rhos-rho)/rho
coef = dsqrt(R*g*d50**3d0)
qs_mpm = dsign(1d0,taus)*coef*qs_mpm
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_std(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Standardized formula 
! Ref: ZANKE, Ulrich et ROLAND, Aron. Sediment bed-load transport: 
!      A standardized notation. Geosciences, 2020
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: tau_c, rhos, rho, d50, g, eps, smpm_a, smpm_b, tau_adim
use module_friction, only: tau
implicit none 
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: qs_std, taus, coef, tausp, R
! -----------------------------------------------------------------------------
! non-dimensional shear stress
taus = tau(h, u, Lx, Jeq, fc)*tau_adim
! MPM modified
R = (rhos-rho)/rho
tausp = abs(taus)
! check tau+*>tau_c*
if (tausp.gt.tau_c) then
  ! qs*
  qs_std = smpm_a*(dsqrt(tausp)-smpm_b*dsqrt(tau_c))*(tausp-tau_c)
  ! qs
  coef = dsqrt(R*g*d50**3d0)
  qs_std = dsign(1d0,taus)*coef*qs_std      
else
  qs_std = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_eh(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Engelund and Hansen (1967)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rhos, rho, d50, g, eps, tau_adim
use module_friction, only: tau, hydraulic_radius
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: qs_eh, taus, Rh, R, coef
! -----------------------------------------------------------------------------
R = (rhos-rho)/rho
coef = dsqrt(R*g*d50**3d0)
Rh = hydraulic_radius(h, Lx)
taus = tau(h, u, Lx, Jeq, fc)*tau_adim
! solid flux
qs_eh = (0.05d0*coef*(fc**2)*(Rh**(1d0/3d0))/g) &
&       *(taus**(5d0/2d0))
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_vr(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Van Rijn (1984)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rhos, rho, d50, g, nu_w, eps, tau_c, tau_adim
use module_friction, only: tau, hydraulic_radius
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: qs_vr, taus, Rh, R, coef, Ds
! -----------------------------------------------------------------------------
R = (rhos-rho)/rho
coef = dsqrt(R*g*d50**3d0)
Rh = hydraulic_radius(h, Lx)
taus = tau(h, u, Lx, Jeq, fc)*tau_adim
! solid flux
if (abs(taus).gt.tau_c) then
  Ds = d50*(R*g/nu_w**2d0)**(1d0/3d0)
  qs_vr = 0.053d0*coef*(((taus/tau_c) - 1d0)**2.1d0)/(Ds**0.3d0)
else
  qs_vr = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_Lefort(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Lefort formula
! Ref: Telemac/Documentation/Courlis
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, eps, d50, d16, dm, d84, g, rhos
use module_friction, only: friction_slope, hydraulic_radius
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: Rh, J, q, qstar, kkp, dms, R, Cd
real(kind=dp)             :: m0, C_dune, C_M, mz, Fq, m
real(kind=dp)             :: qs_Lefort
real(kind=dp)             :: QV, Q0
! -----------------------------------------------------------------------------
R = (rhos-rho)/rho
q = h*u
QV = q*Lx
Rh = hydraulic_radius(h, Lx)
! friction slope
if (Jeq.eq.0d0) then
  J = friction_slope(h, Rh, u, fc)
else
  J = Jeq
endif
! -----------------------------------------------------------------------------
if (J .gt. 0.d0) then
  ! q*
  qstar = q/dsqrt(g*J*dm**3d0)
  ! Ks/Kp
  if (qstar .lt. 200.d0) then 
    kkp = 0.75d0*(qstar/200)**0.23d0
  else
    kkp = 0.75d0
  endif
  ! dm*
  dms = dm*(g*R/(1.d-6)**2d0)**(1d0/3d0)
  ! C(dm*)
  Cd = 0.0444d0*(1.d0 + 15d0/(1d0 + dms) - 1.5d0*exp(-dms / 75d0))
  ! m0
  m0 = 1.6 + 0.06*log10(J)
  ! q0
  Q0 = dsqrt(g*((R*dm)**3d0))*Cd*(dm/Lx)**(1d0/3d0)* &
&       (kkp)**(-0.5d0)*J**(-m0)*Lx
  ! Cdune
  if ((kkp .lt. 0.63d0) .AND. (dms .lt. 14.d0)) then 
    C_dune = 1d0 - 1.4d0*exp(-0.9d0*((kkp)**2d0*dsqrt(QV/Q0)))
  else 
    C_dune = 1d0
  endif
  ! CM
  if (qstar .lt. 200d0) then 
    C_M = (qstar + 2.5d0)/200d0
  else  
    C_M = 1d0
  ! mz
  endif
  if (QV/Q0 .gt. 3.4d0) then
    mz = 1d0 + (0.38d0/dms**0.45d0)*(q/dsqrt(g*dm**3d0))**0.192d0
  else 
    mz = 1d0
  endif
  ! F(q)
  if (QV .lt. Q0) then
    Fq = 0.06d0*C_M*(QV/Q0)
  else 
    Fq = (6.1d0*(1d0 - 0.938d0*(Q0/QV)**0.284d0)**1.66d0)**mz
  endif
  m = 1.8d0 + 0.08d0*log10(J)
  ! qs
  qs_Lefort = QV*1.7d0*J**m*(R+1d0)/R**1.65*(0.5d0* & 
&             (d84/d50+ d50/d16))**0.2d0*C_dune*Fq/((R+1d0))/Lx
else
  qs_Lefort = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_Reck2013(h, u, Lx, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Recking formula 
! source :  Recking2013 & Telemac/Documentation/Courlis
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, eps, d50, dm, d84, g, rhos
use module_friction, only: friction_slope, hydraulic_radius
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
real(kind=dp)             :: Rh, J, q, Phi, R, p
real(kind=dp)             :: qs_Reck2013
real(kind=dp)             :: q_star, taum_star, tau84_star 
! -----------------------------------------------------------------------------
R = (rhos-rho)/rho
q = h*u
Rh = hydraulic_radius(h, Lx)
! friction slope
if (Jeq.eq.0d0) then
  J = friction_slope(h, Rh, u, fc)
else
  J = Jeq
endif
! -----------------------------------------------------------------------------
if (J.gt.0.d0) then
  ! taum*
  if (dm.gt.0.002) then
    taum_star = (5d0*J + 0.06d0)*(d84/d50)**(4.4d0*sqrt(J)-1.5d0)
  else 
    taum_star = 0.045d0
  endif
  ! q*
  q_star = q/sqrt(g*J*d84**3d0)
  ! gamma
  if (q_star .lt. 100.) then
    p = 0.5d0*(1.d0-0.545d0)
  else
    p = 0.5d0*(1.d0-0.395d0)
  endif
  ! tau84*
  tau84_star = J/( R*d84*(2d0/Lx + 74d0*(p**2.6d0)*((g*J)**p)* & 
                  (q**(-2d0*p))*(d84**(3d0*p-1d0))) )
  ! phi
  Phi = 14.d0*tau84_star**2.5d0*(1d0/(1d0+(taum_star/tau84_star)**4d0))
  ! qs
  qs_Reck2013 = Phi*dsqrt(g*R*d84**3d0)
else
  qs_Reck2013 = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function qs_Reck2015(h, u, Lx, morpho, qs_option, Jeq, fc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Recking formula 
! source :  Recking2013 & Telemac/Documentation/Courlis
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, eps, d50, d84, g, rhos
use module_friction, only: friction_slope, hydraulic_radius
implicit none
real(kind=dp), intent(in) :: h, u, Lx, Jeq, fc
integer,       intent(in) :: morpho, qs_option
real(kind=dp)             :: Rh, J, q, Phi, R, p, q_star, taum_star, tau84_star
real(kind=dp)             :: qs_Reck2015
! -----------------------------------------------------------------------------
R = (rhos-rho)/rho
q = h*u
Rh = hydraulic_radius(h, Lx)
! friction slope
if (Jeq.eq.0d0) then
  J = friction_slope(h, Rh, u, fc)
else
  J = Jeq
endif
! -----------------------------------------------------------------------------
if (J.gt.0.d0) then
  ! taum*
  if (morpho .eq. 1) then 
    taum_star = (5*J + 0.06)*((d84/d50)**(4.4*dsqrt(J) - 1.5))
  elseif (morpho .eq. 2) then 
    taum_star = 1.5*J**0.75
  endif
  ! tau84*
  if (qs_option .eq. 1) then 
    q_star = q/sqrt(g*J*d84**3d0)
    if (q_star .lt. 100.d0) then 
      p = 0.23d0
    else
      p = 0.31d0
    endif
    tau84_star = 0.015d0*J*q_star**(2*p)/(R*p**2.5)
  elseif (qs_option .eq. 2) then
    tau84_star = Rh*J/(R*d84)
  endif          
  ! phi
  Phi = 14d0*tau84_star**2.5d0/(1d0 + (taum_star/tau84_star)**4d0)
  ! qs
  qs_Reck2015 = Phi*dsqrt(g*R*d84**3d0)
else
  qs_Reck2015 = 0.d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function critical_shields(d50)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute critical shear stress from d50
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: g, rho, rhos, nu_w
implicit none
real(kind=dp), intent(in) :: d50
real(kind=dp)             :: dadim, critical_shields
! -----------------------------------------------------------------------------
dadim = d50*(((rhos-rho)/rho)*g/nu_w**2d0)**(1d0/3d0)
if (dadim.le.4d0) then
  critical_shields = 0.24d0*dadim**(-1d0)
elseif (dadim.le.10d0) then
  critical_shields = 0.14d0*dadim**(-0.64d0)
elseif (dadim.le.20d0) then
  critical_shields = 0.04d0*dadim**(-0.10d0)
elseif (dadim.le.150d0) then
  critical_shields = 0.013d0*dadim**(0.29d0)
else
  critical_shields = 0.045d0
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                   BEDLOAD FLUXES DERIVATIVES
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
function dtausdh(h, Lx, taus)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the derivative dtaus+/dh
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps_num, frict_law, variable_width, Rh_approx
use module_friction, only: hydraulic_radius
implicit none
real(kind=dp), intent(in) :: h, Lx, taus
real(kind=dp)             :: dtausdh, hc, Rh, dRhdh
! -----------------------------------------------------------------------------
! to avoid division by zero on dry cells
hc = max(h, eps) 
Rh = hydraulic_radius(h, Lx)

if(frict_law.eq.2) then
  ! Chezy :
  dtausdh = -2.d0*taus/hc
elseif (frict_law.eq.3) then
  ! Strickler :
  if (variable_width.eqv..True.) then
    dRhdh = (Lx**2d0)/((Lx+2d0*hc)**2d0)
    dtausdh = taus*((-2d0/hc) - (dRhdh/(3d0*Rh**(1d0/3d0))))
  else
    if (Rh_approx.ne.0) then
      dtausdh = (-7.d0/3.d0)*taus/hc
    else
      dRhdh = 1d0/((1d0+2d0*hc)**2d0)
      dtausdh = taus*((-2d0/hc) - (dRhdh/(3d0*Rh**(1d0/3d0))))
    endif
  endif
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
function dtausdq(h, u, taus)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the derivative dtaus+/dq
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, eps
implicit none
real(kind=dp), intent(in) :: h, u, taus
real(kind=dp)             :: dtausdq, hc
! -----------------------------------------------------------------------------
! to avoid division by zero on dry cells
hc = max(h, eps)
dtausdq = 2.d0*taus/(hc*u)
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
subroutine compute_bedload_derivatives(i, h, u, dqsdh, dqsdq, m)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes qs derivatives: returns dqs/dh and dqs/dq
! -> on cell i (if m=0) 
! -> on interface i between cells i and i+1 (if m=1) (For Roe)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: bedload_formula, qs_derivatives, bed_fric_coef, &
                       Ag, mg, eps, eps_num, Rh, Lx, nu_w, &
                       fric, rhos, rho, d50, g, tau_c, &
                       tau_adim, variable_width, smpm_a, smpm_b, &
                       qs_derivatives_dh, qs_derivatives_du
use module_friction, only: tau
implicit none
integer, intent(in)          :: i, m
real(kind=dp), intent(in)    :: h, u
real(kind=dp), intent(inout) :: dqsdh, dqsdq
real(kind=dp)                :: hc, dqsdu, dqsdtaus, R, coef, taus
real(kind=dp)                :: Rhc, Lc, fc, dh, du, qd0, qdp, qdm
real(kind=dp)                :: aux0, aux1, aux2, dRhdh, Ds
! -----------------------------------------------------------------------------
! to avoid division by zero on dry cells
hc = max(h, eps)
! channel width and hydraulic radius
if (variable_width.eqv..True.) then
  if (m.eq.0) then
    ! taking the value of cells i
    Lc  = Lx(i)
    Rhc = Rh(i)
    fc  = fric(i)*bed_fric_coef
  else 
    ! taking the mean value of cells i and i+1
    Lc  = 0.5d0*(Lx(i) + Lx(i+1))
    Rhc = 0.5d0*(Rh(i) + Rh(i+1))
    fc  = 0.5d0*(fric(i) + fric(i+1))*bed_fric_coef
  endif
else
  if (m.eq.0) then
    Lc  = 1d0
    Rhc = Rh(i)
    fc = fric(i)*bed_fric_coef
  else 
    Lc  = 1d0
    Rhc = 0.5d0*(Rh(i) + Rh(i+1))
    fc  = 0.5d0*(fric(i) + fric(i+1))*bed_fric_coef
  endif
endif
! -----------------------------------------------------------------------------
! Exact derivatives
! -----------------------------------------------------------------------------
if (qs_derivatives.eq.0) then
  ! ***********
  ! Grass model
  ! ***********
  if (bedload_formula.eq.1) then
    dqsdu = Ag*(dabs(u)**(mg-1.d0)+sign(1.d0,u)*(mg-1.d0)*u*dabs(u)**(mg-2.d0))
    dqsdq = dqsdu/hc
    dqsdh = -u*dqsdq
    ! simpler formula
    !dqsdh = -mg*Ag*(((h*u)**mg)/(h**(mg+1)))
    !dqsdq =  mg*Ag*(((h*u)**(mg-1))/(h**mg))

  ! ********************
  ! Grass modified model
  ! ********************
  elseif (bedload_formula.eq.2) then
    dqsdu = Ag*(dabs(u)**(mg-1.d0)+sign(1.d0,u)*(mg-1.d0)*u*dabs(u)**(mg-2.d0))
    dqsdq = dqsdu
    dqsdh = -u*dqsdq

  ! *********
  ! MPM model
  ! *********
  elseif (bedload_formula.eq.3) then
    taus = tau(h, u, Lc, 0d0, fc)*tau_adim
    ! check tau*>tau_c*
    if (abs(taus).gt.tau_c) then
      R = (rhos-rho)/rho
      coef = dsqrt(R*g*d50**3d0) ! qs = coef*(qs*)
      ! dqs*/dtau*
      dqsdtaus = 12.d0*max((taus-tau_c), 0.d0)**(1.d0/2.d0)
      ! dqs/dh etqs/dq
      dqsdh = coef*dqsdtaus*dtausdh(h, Lc, taus)
      dqsdq = coef*dqsdtaus*dtausdq(h, u, taus)
    else
      dqsdh = 0.d0
      dqsdq = 0.d0
    endif

  ! *****************
  ! Standardizd model (Zanke U. and Roland A., 2020)
  ! 1: MPM modified
  ! *****************
  elseif (bedload_formula.eq.4) then
    taus = tau(h, u, Lc, 0d0, fc)*tau_adim
    ! check tau*>tau_c*
    if (abs(taus).gt.tau_c) then
      R = (rhos-rho)/rho
      coef = dsqrt(R*g*d50**3d0) ! qs = coef*(qs*)
      ! dqs*/dtau*
      dqsdtaus = (smpm_a/(2.d0*dsqrt(max(taus, eps_num)))) &
                *(taus - tau_c) + smpm_a*(dsqrt(taus) - smpm_b*dsqrt(tau_c))
      ! dqs/dh et dqs/dq
      dqsdh = coef*dqsdtaus*dtausdh(h, Lc, taus)
      dqsdq = coef*dqsdtaus*dtausdq(h, u, taus)
    else
      dqsdh = 0.d0
      dqsdq = 0.d0
    endif

  ! **************************
  ! Engelund and Hansen (1967)
  ! **************************
  elseif (bedload_formula.eq.8) then
    R = (rhos-rho)/rho
    coef = dsqrt(R*g*d50**3d0) ! qs = coef*(qs*)
    taus = tau(h, u, Lc, 0d0, fc)*tau_adim
    dRhdh = (Lc**2d0)/((Lc+2d0*hc)**2d0)
    aux0 = (0.05d0*fc**2)/g
    aux1 = dRhdh*(taus**(5d0/2d0))/(3*Rhc**(3d0/2d0))
    aux2 = 5d0*(Rhc**(1d0/3d0))*(taus**(3d0/2d0))/2d0
    dqsdh = coef*aux0*( aux1 + aux2*dtausdh(h, Lc, taus) )
    dqsdq = coef*aux0*aux2*dtausdq(h, u, taus)

  ! ***************
  ! Van Rijn (1967)
  ! ***************
  elseif (bedload_formula.eq.9) then
    R = (rhos-rho)/rho
    coef = dsqrt(R*g*d50**3d0) ! qs = coef*(qs*)
    taus = tau(h, u, Lc, 0d0, fc)*tau_adim
    ! solid flux
    if (abs(taus).gt.tau_c) then
      Ds = d50*(R*g/nu_w**2d0)**(1d0/3d0)
      dqsdtaus = 0.1113d0*(((taus/tau_c)-1d0)**1.1d0)/(tau_c*(Ds**0.3d0))
      dqsdh = coef*dqsdtaus*dtausdh(h, Lc, taus)
      dqsdq = coef*dqsdtaus*dtausdq(h, u, taus)
    else
      dqsdh = 0.d0
      dqsdq = 0.d0
    endif

  endif
! -----------------------------------------------------------------------------
! Approximate derivatives
! -----------------------------------------------------------------------------
else
  ! increments for finite differences
  dh = qs_derivatives_dh
  if (u.gt.0d0) then
    du = qs_derivatives_du
  else
    du =-1d0*qs_derivatives_du
  endif

  ! finite differences
  ! ******************
  ! forward FD
  if(qs_derivatives.eq.1) then
    if ((bedload_formula.eq.1).or.(bedload_formula.eq.2)) then
      dqsdu = (qs(h,u+du,Lc,0d0,fc) - qs(h,u,Lc,0d0,fc))/du
      dqsdh = 0d0
    else
      qd0 = qs(h,u,Lc,0d0,fc)
      dqsdu = (qs(h,u+du,Lc,0d0,fc) - qd0)/du
      dqsdh = (qs(h+dh,u,Lc,0d0,fc) - qd0)/dh
    endif

  ! backward FD
  else if(qs_derivatives.eq.2) then
    if ((bedload_formula.eq.1).or.(bedload_formula.eq.2)) then
      dqsdu = (qs(h,u,Lc,0d0,fc) - qs(h,u-du,Lc,0d0,fc))/du
      dqsdh = 0d0
    else
      qd0 = qs(h,u,Lc,0d0,fc)
      dqsdu = (qd0 - qs(h,u-du,Lc,0d0,fc))/du
      dqsdh = (qd0 - qs(h-dh,u,Lc,0d0,fc))/dh
    endif

  ! centered FD
  else if(qs_derivatives.eq.3) then
    if ((bedload_formula.eq.1).or.(bedload_formula.eq.2)) then
      dqsdu = (qs(h,u+du,Lc,0d0,fc) - qs(h,u-du,Lc,0d0,fc))/(2d0*du)
      dqsdh = 0d0
    else
      dqsdu = (qs(h,u+du,Lc,0d0,fc) - qs(h,u-du,Lc,0d0,fc))/(2d0*du)
      dqsdh = (qs(h+dh,u,Lc,0d0,fc) - qs(h-dh,u,Lc,0d0,fc))/(2d0*dh)
    endif

  ! optim FD: ensure overestimation of the spectrum
  else if(qs_derivatives.eq.4) then
    if ((bedload_formula.eq.1).or.(bedload_formula.eq.2)) then
      qd0 = qs(h,u,Lc,0d0,fc)
      qdp = qs(h,u+du,Lc,0d0,fc)
      qdm = qs(h,u-du,Lc,0d0,fc)
      dqsdu = max((qdp - qd0)/du, (qd0 - qdm)/du, (qdp - qdm)/(2d0*du))
      dqsdh = 0d0
    else
      qd0 = qs(h,u,Lc,0d0,fc)
      qdp = qs(h,u+du,Lc,0d0,fc)
      qdm = qs(h,u-du,Lc,0d0,fc)
      dqsdu = max((qdp - qd0)/du, (qd0 - qdm)/du, (qdp - qdm)/(2d0*du))
      qdp = qs(h+dh,u,Lc,0d0,fc)
      qdm = qs(h-dh,u,Lc,0d0,fc)
      dqsdh = min((qdp - qd0)/dh, (qd0 - qdm)/dh, (qdp - qdm)/(2d0*dh))
    endif
  endif

  dqsdq = dqsdu/hc
  dqsdh = dqsdh-u*dqsdq

endif ! approx
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************