! *****************************************************************************
module module_suspension
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Bedload module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, eps
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                   NON-COHESIVE SUSPENSION SOURCE TERM
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine compute_suspension_ceq(h, u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute suspended sediment equilirbium concentration
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rho, g, Lx, fric, rhos, eps_num, Ceq_susp, Csratio, &
                       suspension_formula, tau_adim, difu, susp_critdiff, &
                       nu_u, susp_nu_u, susp_taus_scheme, bed_fric_coef, &
                       ws_susp, rouse_profile, rouse_profile_opt
use module_bedload, only: qs, damping_coef
use module_friction, only: tau
implicit none
real(kind=dp), intent(in) :: h  (0:nx+1)
real(kind=dp), intent(in) :: u  (0:nx+1)
real(kind=dp)             :: qsed, taus, taub, taup, taum, fr
real(kind=dp)             :: ustar, Ro, B, csr, aux
real(kind=dp)             :: zref = 0d0
real(kind=dp)             :: term0, term1, term2, term3, term4
real(kind=dp)             :: poly0, poly1, poly2, poly3, poly4
real(kind=dp)             :: hb, Rob, coefconv
integer                   :: i
logical, parameter        :: debug = .False.
! -----------------------------------------------------------------------------
!$OMP parallel do  &
!$OMP private(qsed, taus, fr, taub, taup, taum) &
!$OMP shared(suspension_formula, Ceq_susp, h, u, Lx, fric, tau_adim, &
!$OMP        rhos, difu, nu_u, susp_nu_u, susp_critdiff, &
!$OMP        susp_taus_scheme)
do i=1,nx
  ! ********************
  ! compute shear stress
  ! ********************
  if (suspension_formula.gt.1) then
    taub = tau(h(i), u(i), Lx(i), 0d0, fric(i)*bed_fric_coef)
    ! centered scheme
    if (susp_taus_scheme.eq.0) then
      taus = taub*tau_adim
    ! upwind scheme
    else if (susp_taus_scheme.eq.1) then
      taup = tau(h(i+1), u(i+1), Lx(i+1), 0d0, fric(i+1)*bed_fric_coef)
      taum = tau(h(i-1), u(i-1), Lx(i-1), 0d0, fric(i-1)*bed_fric_coef)
      fr = u(i)/max(sqrt(g*h(i)), eps_num)
      if (fr.ge.1) then
        taus = 0.5d0*(taub + taup)*tau_adim
      else
        taus = 0.5d0*(taum + taub)*tau_adim
      endif
    else
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unknown SCHEME FOR SHEAR STRESS                    '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
  endif

  ! ****************************************
  ! additional diffusivity for critical flow
  ! ****************************************
  if (susp_critdiff.eqv..True.) then
    fr = u(i)/max(sqrt(g*h(i)), eps_num)
    difu(i) = nu_u + susp_nu_u*fr**2d0
  endif

  ! ***********
  ! compute Ceq
  ! ***********
  ! compute Ceq from qs, with qs from bedload module
  ! ------------------------------------------------
  if (suspension_formula.eq.1) then
    qsed = qs(h(i), u(i), Lx(i), 0d0, fric(i)*bed_fric_coef)
    Ceq_susp(i) = rhos*qsed/max(h(i)*u(i), eps_num)
    rouse_profile = .False.
  ! Engelund and Hansen (1967)
  ! --------------------------
  else if (suspension_formula.eq.2) then
    call compute_ceq_eh(h(i), u(i), Lx(i), fric(i)*bed_fric_coef, &
                        taus, Ceq_susp(i))
  ! Van Rijn (1984)
  ! ---------------
  else if (suspension_formula.eq.3) then
    call compute_ceq_vr(fric(i), taus, Ceq_susp(i), zref)
  ! Zyserman and Fredsoe (1994)
  ! ---------------------------
  else if (suspension_formula.eq.4) then
    call compute_ceq_zf(taus, Ceq_susp(i), zref)
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown SUSPENSION FORMULA                         '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif

  ! ***************************************
  ! compute Csrtaio (for C(zref)=C*Csratio)
  ! ***************************************
  ! Rouse number
  ustar = sqrt(taub/rho)
  Ro = ws_susp/(0.4d0*max(ustar, eps_num))
  Ro = min(Ro, 4d0)

  ! Rouse profile -> csratio
  if (rouse_profile.eqv..True.) then
    ! Option 1: 
    ! ~~~~~~~~~
    if (rouse_profile_opt.eq.1) then
      hb = max(h(i), zref)
      term0 = hb-zref
      term1 = hb*log(hb/zref)
      term2 = (hb/zref)*(hb-zref)
      term3 = (0.5D0*hb)*((hb**2-zref**2)/zref**2)
      term4 = (hb/3.0)  *((hb**3-zref**3)/zref**3)
      poly0 = term0
      poly1 = term1 - term0
      poly2 = term2 - 2.D0*term1 + term0
      poly3 = term3 - 3.D0*term2 + 3.D0*term1 - term0
      poly4 = term4 - 4.D0*term3 + 6.D0*term2 - 4.D0*term1 + term0
      Rob = int(Ro)
      coefconv = Ro-Rob
      poly0 = (hb-zref)/max(poly0,eps_num)
      poly1 = ((hb-zref)**2/zref)/max(poly1,eps_num)
      poly2 = ((hb-zref)**3/zref**2)/max(poly2,eps_num)
      poly3 = ((hb-zref)**4/zref**3)/max(poly3,eps_num)
      poly4 = ((hb-zref)**5/zref**4)/max(poly4,eps_num)
      if(Rob.eq.0) then
        csr = (1-coefconv)*poly0+coefconv*poly1
      elseif(Rob.eq.1) then
        csr = (1-coefconv)*poly1+coefconv*poly2
      elseif(Rob.eq.2) then
        csr = (1-coefconv)*poly2+coefconv*poly3
      elseif(Rob.eq.3) then
        csr = (1-coefconv)*poly3+coefconv*poly4
      elseif(Rob.eq.4) then
        csr = poly4
      endif
      csr = max(csr, eps_num)
    ! Option 2:
    ! ~~~~~~~~~ 
    else if (rouse_profile_opt.eq.2) then
      B = zref/max(h(i), zref)
      aux = Ro-1.d0
      if(abs(aux).gt.1.d-4) then
        aux = min(aux, 3.d0)
        csr = B*(1.d0-B**aux)/aux
      else
        csr = -B*log(B)
      endif
      csr = max(1.d0/max(csr, eps_num),1.d0)
    endif
  else
    csr = 1d0
  endif
  
  Csratio(i) = csr

  ! *****
  ! debug
  ! *****
  if (debug.eqv..True.) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) "taub    = ", taub
    write(*,*) "taus    = ", taus
    write(*,*) "Ceq     = ", Ceq_susp(i)
    write(*,*) "Ro      = ", Ro
    write(*,*) "Csratio = ", Csratio(i)
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
enddo
!$OMP end parallel do
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine source_suspension(zb, zf, h, Ceq, Csr, t, ht, EmD)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute suspended sediment erosion-deposition source term;
! time schemes for Tracer: explicit or implicit
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps_num, dt, susp_trac_scheme, suspension_coef, ws_susp,&
                       eps_b_damping, eps_b
use module_bedload, only: qs, damping_coef
real(kind=dp), intent(in)    :: zf, zb, h, Ceq, Csr, t
real(kind=dp), intent(inout) :: ht, EmD
real(kind=dp)                :: aux, b, rb
! -----------------------------------------------------------------------------
! Erosion-Deposition source term
! Here: Czref = t*Csratio
EmD = suspension_coef*ws_susp*(Ceq-t*Csr)

! Non-erodible bed layer damping zone
if (eps_b_damping.eq.1) then
  b = max(zb-zf, 0d0)
  if (b.le.eps_b) then
    rb = b/eps_b
  else
    rb = 1d0
  endif
endif

! Explicit time scheme with limiter
! *********************************
if (susp_trac_scheme.eq.0) then
  ! erosion limiter  
  if (EmD.gt.0d0) then 
    EmD = EmD*damping_coef(rb)
  endif
  ! deposition limiter
  EmD = max(EmD, -ht/dt)
  ! Update tracer
  ht = ht + dt*EmD

! Implicit time scheme
! ********************
else
  aux = dt*suspension_coef*ws_susp
  ! erosion limiter
  if (EmD.gt.0d0) then 
    aux = aux*damping_coef(rb)
  endif
  ! Update tracer
  ht = (ht + aux*Ceq)/(1d0 + aux*Csr/max(h, eps))
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                          SUSPENSION FUNCTIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
function settling_velocity(d50)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute suspended sediments settling velocity
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: g, rho, rhos, nu_w
implicit none
real(kind=dp), intent(in) :: d50
real(kind=dp)             :: settling_velocity
real(kind=dp)             :: s
! -----------------------------------------------------------------------------
s = rhos/rho
! Low Reynolds : Stokes' law
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
if (d50 .lt. 1d-4) then
  settling_velocity = abs(s-1d0)*(g*d50**2d0)/(18d0*nu_w)
! High Reynolds :  Van Rijn (1985)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif (d50 .gt. 1d-3) then
  settling_velocity =  1.1d0*sqrt((s-1d0)*g*d50)
! Intermediate : Ruby and Zanke (1977)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else 
  settling_velocity = (1d1*nu_w/d50)*(sqrt(1d0+ abs(s-1d0) &
                                    *g*d50**3d0/(1d2*nu_w**2d0))-1d0)
endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
subroutine compute_ceq_eh(h, u, Lx, fc, taus, ceq_eh)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Engelund and Hansen (1967)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rhos, rho, d50, g, eps, eps_num
use module_friction, only: hydraulic_radius, tau
implicit none
real(kind=dp), intent(in)    :: h, u, Lx, fc, taus
real(kind=dp), intent(inout) :: ceq_eh
real(kind=dp)                :: qs_eh, Rh, R, coef
! -----------------------------------------------------------------------------
R = (rhos-rho)/rho
coef = dsqrt(R*g*d50**3d0)
Rh = hydraulic_radius(h, Lx)
! solid flux
qs_eh = (0.05d0*coef*(fc**2)*(Rh**(1d0/3d0))/g)*(taus**(5d0/2d0))
! equilibium concentration
ceq_eh = rhos*qs_eh/max(h*u, eps_num)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_ceq_vr(fc, taus, ceq_vr, zref)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! van Rijn L.C. Sediment transport - Part II : suspended load. 
! J. of Hydraulic Division, HY11:1631–1641, 1984.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rhos, rho, d50, g, nu_w, eps, tau_c, eps_num, frict_law
implicit none
real(kind=dp), intent(in)    :: fc, taus
real(kind=dp), intent(inout) :: ceq_vr, zref
real(kind=dp)                :: R, Ds, ks, aux
! -----------------------------------------------------------------------------
! check friction formula and get correct ks value
if ((frict_law.ge.0).and.(frict_law.le.2)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : Friction law not compatible with Van Rijn  '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
else if (frict_law.eq.3) then
  ! converting Strickler coefficient to roughness heigth
  ks = (8.3d0*dsqrt(g)/fc)**6d0
else if (frict_law.gt.3) then
  ks = fc
endif
! compute Ceq
if (abs(taus).gt.tau_c) then
  R = (rhos-rho)/rho
  Ds = d50*(R*g/(nu_w**2d0))**(1d0/3d0)
  zref = 0.5d0*ks
  aux = (taus-tau_c)/tau_c
  ceq_vr = 0.015*d50*sqrt(aux**3d0)/(zref*Ds**0.3d0)
else
  ceq_vr = 0.d0
endif
! convert to kg.m-3
ceq_vr = min(ceq_vr, 0.6d0)*rhos
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_ceq_zf(taus, ceq_zf, zref)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Zyserman J.A. and Fredsoe J. 
! Data analysis of bed concentration of suspended sediment.
! Journal of Hydraulic Engineering, 120:1021–1042, 1994.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: rhos, rho, g, nu_w, eps, tau_c, eps_num, d50
implicit none
real(kind=dp), intent(in)    :: taus
real(kind=dp), intent(inout) :: ceq_zf, zref
real(kind=dp)                :: aux
! -----------------------------------------------------------------------------
if (abs(taus).gt.tau_c) then
  aux = (taus-tau_c)**1.75d0
  ceq_zf = 0.331d0*aux/(1d0 + 0.72d0*aux)
else
  ceq_zf = 0.d0
endif
! convert to kg.m-3
ceq_zf = min(ceq_zf, 0.6d0)*rhos
zref = 3d0*d50
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************