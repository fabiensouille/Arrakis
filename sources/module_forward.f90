! *****************************************************************************
module module_forward
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Time update module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, dt, ci, eps, nx, g
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                           TIME UPDATE FUNCTIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine forward
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! solution increment over one time step
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: time_scheme, exner
implicit none
! -----------------------------------------------------------------------------
! Solver for Saint-Venant equations
! -----------------------------------------------------------------------------
if (exner.eqv..False.) then
  ! *************
  ! Steady kernel
  ! *************
  if (time_scheme.eq.0) then
    call update_sv_steady
  ! ************
  ! Euler scheme
  ! ************
  elseif (time_scheme.eq.1) then
    call update_sv_euler
  ! ***********
  ! Heun scheme
  ! ***********
  elseif (time_scheme.eq.2) then
    call update_sv_heun
  ! *********************************
  ! Runge-Kutta (second order) scheme
  ! *********************************
  elseif (time_scheme.eq.3) then
    call update_sv_rk2
  ! ********************************
  ! Runge-Kutta (third order) scheme
  ! ********************************
  elseif (time_scheme.eq.4) then
    call update_sv_rk3
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown TIME SCHEME for SV '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
! -----------------------------------------------------------------------------
! Solver for Saint-Venant-Exner equations
! -----------------------------------------------------------------------------
else
  ! *************
  ! Steady kernel
  ! *************
  if (time_scheme.eq.0) then
    call update_sve_steady
  ! ************
  ! Euler scheme
  ! ************
  elseif (time_scheme.eq.1) then
    call update_sve_euler
  ! ***********
  ! Heun scheme
  ! ***********
  elseif (time_scheme.eq.2) then
    call update_sve_heun
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown TIME SCHEME for SVE '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
endif 
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
elemental subroutine update_u(h, q, u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! update of non-conservative variables: u from q
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(in)    :: h
real(kind=dp), intent(inout) :: q
real(kind=dp), intent(out)   :: u
! -----------------------------------------------------------------------------
if (h>eps) then
  u = q/h
else
  u = 0d0
  q = 0d0
endif
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
elemental subroutine update_t(h, ht, t)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! update of non-conservative variables: t from ht
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(in)  :: h, ht
real(kind=dp), intent(out) :: t
! -----------------------------------------------------------------------------
if (h>eps) then
  t = ht/h
else
  t = ht/eps
endif
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
pure subroutine update_Rh(h, Lx, Rh)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! update hydraulic radius
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: variable_width, eps, Rh_approx, Lcst
implicit none
! -----------------------------------------------------------------------------
real(kind=dp), intent(in)    :: h
real(kind=dp), intent(in)    :: Lx
real(kind=dp), intent(inout) :: Rh
! -----------------------------------------------------------------------------
if (variable_width.eqv..True.) then
  Rh = max((h*Lx)/(Lx + 2.d0*h), eps)
else
  if (Rh_approx.eq.0) then
      Rh = max(h, eps)
  else
      Rh = max((h*Lcst)/(Lcst + 2.d0*h), eps)
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                               STEADY KERNEL
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine update_sv_steady
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Steady kernel
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,     only: ntrac, Rh, Lx, zb, h, u, q, t, ht, dx, &
                           ha0, hb0, qa0, qb0, write_listing, &
                           std_scheme_opt, std_critical_opt, std_critical_fr
use module_timestep, only: init_time_step
use module_steady,   only: compute_water_line, compute_water_line_linear, &
                           compute_tracers
use module_bc,       only: compute_cdl_steady
implicit none
integer            :: i
real(kind=dp)      :: hd, qd, zd, hcri, hcri_clip, dxgd
logical            :: upstream_control, skip_update, critical_flag, conv_flag
! -----------------------------------------------------------------------------
call init_time_step(h, u)
call compute_cdl_steady(upstream_control, skip_update, ha0, hb0, qa0, qb0)
! update of h, q
if (skip_update.eqv..False.) then
  ! compute flowrate
  do i=1,nx+1
    q(i) = q(0)*Lx(0)/Lx(i)
  enddo
  ! ---------------------------------------------------------------------------
  ! compute water line
  critical_flag = .False.
  hd = h(nx+1)
  qd = q(nx+1)
  zd = zb(nx+1)
  dxgd = dx(nx)
  do i=nx,1,-1
    ! water line solver (linear/non-linear)
    if (std_scheme_opt.eq.0) then
      call compute_water_line(i, h(i), q(i), zb(i), hd, qd, zd, dxgd, conv_flag)
    elseif (std_scheme_opt.eq.1) then
      call compute_water_line_linear(i, h(i), zb(i), hd, qd, zd, dxgd, conv_flag)
    endif
    ! check for critical depth
    hcri = (q(i)/sqrt(g))**(2d0/3d0)
    hcri_clip = (q(i)/(std_critical_fr*sqrt(g)))**(2d0/3d0)
    if ((h(i).le.hcri).or.(conv_flag.eqv..False.)) then
      h(i) = hcri_clip
      critical_flag = .True.
      if (std_critical_opt.eq.1) then
        dxgd = dxgd + dx(i-1)
      else
        hd = h(i)
        qd = q(i)
        zd = zb(i)
        dxgd = dx(i-1)
      endif
    else
      hd = h(i)
      qd = q(i)
      zd = zb(i)
      dxgd = dx(i-1)
    endif
  enddo
  if ((critical_flag).and.(write_listing)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : CRITICAL FLOW DETECTED IN THE DOMAIN '
    write(*,*) '           CONSIDER SWITCHING TO UNSEADY SOLVER '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! ---------------------------------------------------------------------------
  ! update of Rh and u
  do i=1,nx
    call update_Rh(h(i), Lx(i), Rh(i))
    call update_u(h(i), q(i), u(i))
  enddo
endif
! update of ht
if (ntrac>0) then
  call compute_tracers(h, u, q, t, ht)
endif
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
subroutine update_sve_steady
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Steady kernel
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ntrac, Rh, Lx, zf, zb, h, u, q, t, ht, dx, qsed, &
                       bedload, suspension, porosite, rhos, flut_susp, &
                       fluz_susp, Ceq_susp, Csratio, flut_susp, &
                       ha0, hb0, qa0, qb0, std_scheme_opt, write_listing, &
                       std_critical_opt, std_critical_fr
use module_timestep,   only: init_time_step
use module_steady,     only: compute_water_line, compute_tracers, &
                             compute_water_line_linear
use module_bc,         only: compute_cdl_steady
use module_bedload,    only: compute_bedload_solidfluxes
use module_suspension, only: compute_suspension_ceq, source_suspension
use module_fluxes_sve, only: compute_slope_effect, init_fluxes_exner, &
                             compute_flused_limiter
implicit none
integer            :: i, k
real(kind=dp)      :: fact
real(kind=dp)      :: fluz_lim(1:nx)
real(kind=dp)      :: hd, qd, zd, hcri, hcri_clip, dxgd
logical            :: upstream_control, skip_update, critical_flag, conv_flag
! -----------------------------------------------------------------------------
call init_time_step(h, u)
call init_fluxes_exner()
! bedload solid flux
if (bedload.eqv..True.) then
  call compute_bedload_solidfluxes(zf, zb, h, u)
endif
! If suspension, compute fluz_susp
if (suspension.eqv..True.) then
  call compute_suspension_ceq(h, u)
endif
! slope effect
call compute_slope_effect(zb)
! boundary condition
call compute_cdl_steady(upstream_control, skip_update, ha0, hb0, qa0, qb0)
call compute_flused_limiter(zf, zb, fluz_lim)
! -----------------------------------------------------------------------------
! update of h, q
if (skip_update.eqv..False.) then
  ! compute flowrate
  do i=1,nx+1
    q(i) = q(0)*Lx(0)/Lx(i)
  enddo
  ! ---------------------------------------------------------------------------
  ! compute water line
  hd = h(nx+1)
  qd = q(nx+1)
  zd = zb(nx+1)
  dxgd = dx(nx)
  do i=nx,1,-1
    ! water line solver (linear/non-linear)
    if (std_scheme_opt.eq.0) then
      call compute_water_line(i, h(i), q(i), zb(i), hd, qd, zd, dxgd, conv_flag)
    elseif (std_scheme_opt.eq.1) then
      call compute_water_line_linear(i, h(i), zb(i), hd, qd, zd, dxgd, conv_flag)
    endif
    ! check for critical depth
    hcri = (q(i)/sqrt(g))**(2d0/3d0)
    hcri_clip = (q(i)/(std_critical_fr*sqrt(g)))**(2d0/3d0)
    if ((h(i).le.hcri).or.(conv_flag.eqv..False.)) then
      h(i) = hcri_clip
      critical_flag = .True.
      if (std_critical_opt.eq.1) then
        dxgd = dxgd + dx(i-1)
      else
        hd = h(i)
        qd = q(i)
        zd = zb(i)
        dxgd = dx(i-1)
      endif
    else
      hd = h(i)
      qd = q(i)
      zd = zb(i)
      dxgd = dx(i-1)
    endif
  enddo
  if ((critical_flag).and.(write_listing)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : CRITICAL FLOW DETECTED IN THE DOMAIN '
    write(*,*) '           CONSIDER SWITCHING TO UNSTEADY SOLVER '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! ---------------------------------------------------------------------------
  ! update of Rh and u
  do i=1,nx
    call update_Rh(h(i), Lx(i), Rh(i))
    call update_u(h(i), q(i), u(i))
  enddo
endif
! update of ht
if (ntrac>0) then
  call compute_tracers(h, u, q, t, ht)
  do i=1,nx
    do k=1,ntrac
        ! erosion-deposition source term for tracer
      if (k.eq.1) then
        call source_suspension(zb(i), zf(i), h(i), Ceq_susp(i), &
                               Csratio(i), t(i,k), ht(i,k), flut_susp(i))
      endif
    enddo
  enddo
endif
! Exner update (upwind)
do i=1,nx
  fact = dt/(Lx(i)*ci(i))
  ! erosion-deposition bed flux
  fluz_susp(i) = ci(i)*flut_susp(i)/((1.d0 - porosite)*rhos)
  ! update of zb
  zb(i) = zb(i) - fact*(min(qsed(i)*Lx(i) - qsed(i-1)*Lx(i-1) &
                     + fluz_susp(i)*Lx(i), fluz_lim(i)*Lx(i)))
enddo
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                      SAINT-VENANT FINITE VOLUME UPDATE
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine update_sv_euler
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Euler
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Rh, Lx, zb, h, u, q, t, ht, fric, & 
                       ntrac, ordre, ordre_trac, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r
use module_timestep,   only: init_time_step, adjust_time_step
use module_high_order, only: compute_reconstructions
use module_fluxes_sv,  only: init_fluxes, compute_flux_sw
use module_bc,         only: compute_cdl, update_gc
use module_friction,   only: source_friction
implicit none
real(kind=dp) :: fact
integer       :: i, k
! -----------------------------------------------------------------------------
call init_time_step(h, u)
call init_fluxes
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h, u, q, t)
endif
call compute_flux_sw(zb, h, u, t)
call adjust_time_step
call compute_cdl
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do  &
!$OMP private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        t, u, Rh, Lx, fric, ci, dt, ntrac, nx)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  ! update of h, q
  h(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
  q(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                                -flut_r(i-1,k)+flut_src_r(i-1,k))
      call update_t(h(i), ht(i,k), t(i,k))
    enddo
  endif
  ! Sources
  call update_Rh(h(i), Lx(i), Rh(i))
  call source_friction(h(i), Rh(i), q(i), fric(i))
  ! update u
  call update_u(h(i), q(i), u(i))
enddo
!$OMP end parallel do
call update_gc(h, q, u, ht, t)
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
subroutine update_sv_heun
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Heun 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Rh, Lx, zb, h, u, q, t, ht, fric, &
                       ntrac, ordre, ordre_trac, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r
use module_timestep,   only: init_time_step, adjust_time_step
use module_high_order, only: compute_reconstructions
use module_fluxes_sv,  only: init_fluxes, compute_flux_sw
use module_bc,         only: compute_cdl, update_gc
use module_friction ,    only: source_friction
implicit none
real(kind=dp) :: h0   (0:nx+1)
real(kind=dp) :: q0   (0:nx+1)
real(kind=dp) :: ht0  (0:nx+1,ntrac)
real(kind=dp) :: dt1, dt2, aux1, aux2, gamma, fact
integer       :: i, k, iord, mord
! -----------------------------------------------------------------------------
! time order
!   if mord==1: Euler scheme (1st order)
!   if mord==2: Heun scheme (second order)
mord = 2
! -----------------------------------------------------------------------------
if (mord==2) then
  do i=0,nx+1
    h0(i) = h(i)
    q0(i) = q(i)
    if (ntrac>0) then
      do k=1,ntrac
        ht0(i,k) = ht(i,k)
      enddo
    endif
  enddo
endif
! -----------------------------------------------------------------------------
! Heun iterations (iord=1: first, iord=2: second)
! -----------------------------------------------------------------------------
do iord=1,mord
  call init_time_step(h, u)
  call init_fluxes
  if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
    call compute_reconstructions(zb, h, u, q, t)
  endif
  call compute_flux_sw(zb, h, u, t)
  call adjust_time_step
  call compute_cdl
  if (iord==1) then
    dt1 = dt
  endif
  ! -----------------------------------------------------------------------------
  ! OMP DIRECTIVES
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  !$OMP parallel do  &
  !$OMP private(fact, k) &
  !$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
  !$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
  !$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
  !$OMP        t, u, Rh, Lx, fric, ci, dt, ntrac, nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  do i=1,nx
    fact = dt/ci(i)
    ! update of h, q
    h(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
    q(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
    ! update of ht
    if (ntrac>0) then
      do k=1,ntrac
        ht(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                                 - flut_r(i-1,k)+flut_src_r(i-1,k))
        call update_t(h(i), ht(i,k), t(i,k))
      enddo
    endif
    ! sources
    call update_Rh(h(i), Lx(i), Rh(i))
    call source_friction(h(i), Rh(i), q(i), fric(i))
    ! update u
    call update_u(h(i), q(i), u(i))
  enddo
  !$OMP end parallel do
  call update_gc(h, q, u, ht, t)
enddo
! -----------------------------------------------------------------------------
! Heun update:
! -----------------------------------------------------------------------------
if (mord==2) then
  dt2=dt
  aux1 = 2d0*dt1*dt2
  aux2 = dt1 + dt2
  dt = aux1/aux2
  gamma = aux1/aux2**2
  ! OMP DIRECTIVES
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  !$OMP parallel do  &
  !$OMP private(k) &
  !$OMP shared(gamma, h, h0, q, q0, ht, ht0, t, u, Rh, ntrac)
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  do i=1,nx
    ! update of h, q
    h(i) = (1d0-gamma)*h0(i) + gamma*h(i)
    q(i) = (1d0-gamma)*q0(i) + gamma*q(i)
    ! update of ht
    if (ntrac>0) then
      do k=1,ntrac
        ht(i,k) = (1d0-gamma)*ht0(i,k) + gamma*ht(i,k)
        call update_t(h(i), ht(i,k), t(i,k))
      enddo
    endif
    ! update u
    call update_u(h(i), q(i), u(i))
  enddo
  !$OMP end parallel do
  call update_gc(h, q, u, ht, t)
endif
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
subroutine update_sv_rk2
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Second order TVD Runge-Kutta
! Gottlieb, S. & Shu, C.-W. Total variation diminishing Runge-Kutta
! schemes Mathematics of computation, 1998, 67, 73-85
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Rh, Lx, zb, h, u, q, t, ht, fric, &
                       ntrac, ordre, ordre_trac, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r
use module_timestep,   only: init_time_step, adjust_time_step
use module_high_order, only: compute_reconstructions
use module_fluxes_sv,  only: init_fluxes, compute_flux_sw
use module_bc,         only: compute_cdl, update_gc
use module_friction ,    only: source_friction
implicit none
real(kind=dp) :: h1  (0:nx+1)
real(kind=dp) :: u1  (0:nx+1)
real(kind=dp) :: q1  (0:nx+1)
real(kind=dp) :: t1  (0:nx+1,1:ntrac)
real(kind=dp) :: ht1 (0:nx+1,1:ntrac)
real(kind=dp) :: fact
integer       :: i, k
! -----------------------------------------------------------------------------
! First step RK2
! -----------------------------------------------------------------------------
! inititalization of h1, q1, ht1
do i=0,nx+1
  h1(i) = h(i)
  q1(i) = q(i)
  u1(i) = u(i)
  if (ntrac>0) then
    do k=1,ntrac
      ht1(i,k) = ht(i,k)
      t1(i,k) = t(i,k)
    enddo
  endif
enddo
call init_time_step(h, u)   ! Only at first step
call init_fluxes
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h, u, q, t)
endif
call compute_flux_sw(zb, h, u, t)
call adjust_time_step       ! Only at first step
call compute_cdl
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        t, u, Rh, Lx, fric)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  h1(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
  q1(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht1(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                                - flut_r(i-1,k)+flut_src_r(i-1,k))
      call update_t(h1(i), ht1(i,k), t1(i,k))
    enddo
  endif
  ! sources
  call update_Rh(h1(i), Lx(i), Rh(i))
  call source_friction(h1(i), Rh(i), q1(i), fric(i))
  ! update u
  call update_u(h1(i), q1(i), u1(i))
enddo
!$OMP end parallel do
call update_gc(h1, q1, u1, ht1, t1)
! -----------------------------------------------------------------------------
! Second step RK2
! -----------------------------------------------------------------------------
call init_fluxes
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h1, u1, q1, t1)
endif
call compute_flux_sw(zb, h1, u1, t1)
call compute_cdl
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        h1, q1, ht1, t, u, Rh, Lx, fric)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  h(i) = 0.5d0*(h(i)+h1(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i) &
                                                         + fluh_src_r(i-1)))
  q(i) = 0.5d0*(q(i)+q1(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i) &
                                                         + fluq_src_r(i-1)))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht(i,k) = 0.5d0*(ht(i,k) + ht1(i,k) &
              - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                     - flut_r(i-1,k)+flut_src_r(i-1,k)))
      call update_t(h(i), ht(i,k), t(i,k))
    enddo
  endif
  ! sources
  call update_Rh(h(i), Lx(i), Rh(i))
  call source_friction(h(i), Rh(i), q(i), fric(i))
  ! update u
  call update_u(h(i), q(i), u(i))
enddo
!$OMP end parallel do
call update_gc(h, q, u, ht, t)
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
subroutine update_sv_rk3
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Third order TVD Runge-Kutta
! Gottlieb, S. & Shu, C.-W. Total variation diminishing Runge-Kutta
! schemes Mathematics of computation, 1998, 67, 73-85
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Rh, Lx, zb, h, u, q, t, ht, fric, &
                       ntrac, ordre, ordre_trac, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r
use module_timestep,   only: init_time_step, adjust_time_step
use module_high_order, only: compute_reconstructions
use module_fluxes_sv,  only: init_fluxes, compute_flux_sw
use module_bc,         only: compute_cdl, update_gc
use module_friction ,    only: source_friction
implicit none
real(kind=dp) :: h0  (0:nx+1)
real(kind=dp) :: q0  (0:nx+1)
real(kind=dp) :: h2  (0:nx+1)
real(kind=dp) :: u2  (0:nx+1)
real(kind=dp) :: q2  (0:nx+1)
real(kind=dp) :: t2  (0:nx+1,1:ntrac)
real(kind=dp) :: ht0 (0:nx+1,1:ntrac)
real(kind=dp) :: ht2 (0:nx+1,1:ntrac)
real(kind=dp) :: fact
integer       :: i, k
! -----------------------------------------------------------------------------
! First step RK3 -> W1 (strored in W)
! -----------------------------------------------------------------------------
do i=0,nx+1
  h0(i) = h(i)
  q0(i) = q(i)
  if (ntrac>0) then
    do k=1,ntrac
      ht0(i,k) = ht(i,k)
    enddo
  endif
enddo
call init_time_step(h, u)  ! Only at first step
call init_fluxes
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h, u, q, t)
endif
call compute_flux_sw(zb, h, u, t)
call adjust_time_step     ! Only at first step
call compute_cdl
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        t, u, Rh, Lx, fric)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  h(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
  q(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                               - flut_r(i-1,k)+flut_src_r(i-1,k))
      call update_t(h(i), ht(i,k), t(i,k))
    enddo
  endif
  ! sources
  call update_Rh(h(i), Lx(i), Rh(i))
  call source_friction(h(i), Rh(i), q(i), fric(i))
  ! update u
  call update_u(h(i), q(i), u(i))
enddo
!$OMP end parallel do
call update_gc(h, q, u, ht, t)
! -----------------------------------------------------------------------------
! Second step RK3 -> W2
! -----------------------------------------------------------------------------
do i=0,nx+1
  h2(i) = h(i)
  q2(i) = q(i)
  u2(i) = u(i)
  if (ntrac>0) then
    do k=1,ntrac
      ht2(i,k) = ht(i,k)
      t2(i,k) = t(i,k)
    enddo
  endif
enddo
call init_fluxes
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h, u, q, t)
endif
call compute_flux_sw(zb, h, u, t)
call compute_cdl
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        h2, q2, ht2, h0, q0, ht0, t2, t, u2, u, &
!$OMP        Rh, Lx, fric)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  h(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
  q(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                               - flut_r(i-1,k)+flut_src_r(i-1,k))
      call update_t(h(i), ht(i,k), t(i,k))
    enddo
  endif
  ! sources
  call update_Rh(h(i), Lx(i), Rh(i))
  call source_friction(h(i), Rh(i), q(i), fric(i))
  ! update u
  call update_u(h(i), q(i), u(i))
  ! -----------------------------------------------------------------------------
  ! maj RK step
  h2(i) = (3d0/4d0)*h0(i) + (1d0/4d0)*h(i)
  q2(i) = (3d0/4d0)*q0(i) + (1d0/4d0)*q(i)
  call update_u(h2(i), q2(i), u2(i))
  if (ntrac>0) then
    do k=1,ntrac
      ht2(i,k) = (3d0/4d0)*ht0(i,k) + (1d0/4d0)*ht(i,k)
      call update_t(h2(i), ht2(i,k), t2(i,k))
    enddo
  endif
enddo
!$OMP end parallel do
call update_gc(h2, q2, u2, ht2, t2)
! -----------------------------------------------------------------------------
! Third step RK3 -> Wn+1
! -----------------------------------------------------------------------------
call init_fluxes
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h2, u2, q2, t2)
endif
call compute_flux_sw(zb, h2, u2, t2)
call compute_cdl
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        h2, q2, ht2, h0, q0, ht0, t2, t, u2, u, &
!$OMP        Rh, Lx, fric)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  h2(i) = h2(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
  q2(i) = q2(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht2(i,k) = ht2(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                                 - flut_r(i-1,k)+flut_src_r(i-1,k))
      call update_t(h2(i), ht2(i,k), t2(i,k))
    enddo
  endif
  ! sources
  call update_Rh(h2(i), Lx(i), Rh(i))
  call source_friction(h2(i), Rh(i), q2(i), fric(i))
  ! update u
  call update_u(h2(i), q2(i), u2(i))
  ! -----------------------------------------------------------------------------
  ! maj RK step
  h(i) = (1d0/3d0)*h0(i) + (2d0/3d0)*h2(i)
  q(i) = (1d0/3d0)*q0(i) + (2d0/3d0)*q2(i)
  call update_u(h(i), q(i), u(i))
  if (ntrac>0) then
    do k=1,ntrac
      ht(i,k) = (1d0/3d0)*ht0(i,k) + (2d0/3d0)*ht2(i,k)
      call update_t(h(i), ht(i,k), t(i,k))
    enddo
  endif
enddo
!$OMP end parallel do
call update_gc(h, q, u, ht, t)
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                         SAINT-VENANT-EXNER UPDATE
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine update_sve_euler
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Euler for the SVE system
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Rh, Lx, zf, zb, h, u, q, t, ht, fric, ntrac, ordre, &
                       ordre_trac, bedload, suspension, porosite, rhos, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r, &
                       fluqsed_l, fluqsed_r, fluqsed_src_l, fluqsed_src_r, &
                       fluz_susp, Ceq_susp, Csratio, flut_susp
use module_timestep,   only: init_time_step, adjust_time_step
use module_high_order, only: compute_reconstructions
use module_fluxes_sv,  only: init_fluxes, compute_flux_sw
use module_fluxes_sve, only: init_fluxes_exner, compute_flused_limiter, &
                             compute_flux_swe, compute_slope_effect
use module_bedload,    only: compute_bedload_solidfluxes
use module_suspension, only: compute_suspension_ceq, source_suspension
use module_bc,         only: compute_cdl, update_gc
use module_friction ,    only: source_friction
implicit none
real(kind=dp) :: fluz_lim(1:nx)
real(kind=dp) :: fact
integer       :: i, k
! -----------------------------------------------------------------------------
call init_time_step(h, u)
call init_fluxes
call init_fluxes_exner()
! bedload solid flux
if (bedload.eqv..True.) then
  call compute_bedload_solidfluxes(zf, zb, h, u)
endif
! second order reconstructions
if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
  call compute_reconstructions(zb, h, u, q, t)
endif
! If bedload, solve coupled SVE system
if (bedload.eqv..True.) then
  call compute_flux_swe(zb, h, u, t)
! If not, solve standard SV system
else
  call compute_flux_sw(zb, h, u, t)
endif
! If suspension, compute fluz_susp
if (suspension.eqv..True.) then
  call compute_suspension_ceq(h, u)
endif
! Slope effect
call compute_slope_effect(zb)
call compute_cdl
call adjust_time_step
call compute_flused_limiter(zf, zb, fluz_lim)
! -----------------------------------------------------------------------------
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do  &
!$OMP private(fact, k) &
!$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        zb, fluqsed_l, fluqsed_r, zf, fluqsed_src_l, fluqsed_src_r, &
!$OMP        Ceq_susp, flut_susp, fluz_susp, fluz_lim, &
!$OMP        t, u, Rh, Lx, fric, ci, dt, ntrac, nx, porosite, rhos)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx
  fact = dt/ci(i)
  ! update of h, q
  h(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
  q(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
  ! update of ht
  if (ntrac>0) then
    do k=1,ntrac
      ht(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                               - flut_r(i-1,k)+flut_src_r(i-1,k))
      ! erosion-deposition source term for tracer
      if (k.eq.1) then
        call source_suspension(zb(i), zf(i), h(i), Ceq_susp(i), &
                               Csratio(i), t(i,k), ht(i,k), flut_susp(i))
      endif
      call update_t(h(i), ht(i,k), t(i,k))
    enddo
  endif
  ! erosion-deposition bed flux
  fluz_susp(i) = ci(i)*flut_susp(i)/((1.d0 - porosite)*rhos)
  ! update of zb
  zb(i) = zb(i) - fact*(min(fluqsed_l(i)   + fluqsed_src_l(i) &
                          - fluqsed_r(i-1) + fluqsed_src_r(i-1) &
                          + fluz_susp(i), fluz_lim(i)))
  ! sources
  call update_Rh(h(i), Lx(i), Rh(i))
  call source_friction(h(i), Rh(i), q(i), fric(i))
  ! update u
  call update_u(h(i), q(i), u(i))
enddo
!$OMP end parallel do
call update_gc(h, q, u, ht, t)
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
subroutine update_sve_heun
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update step - Heun for the SVE system
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: Rh, Lx, zf, zb, h, u, q, t, ht, fric, ntrac, ordre, &
                       ordre_trac, bedload, suspension, porosite, rhos, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r, &
                       fluqsed_l, fluqsed_r, fluqsed_src_l, fluqsed_src_r, &
                       fluz_susp, Ceq_susp, Csratio, flut_susp
use module_timestep,   only: init_time_step, adjust_time_step
use module_high_order, only: compute_reconstructions
use module_fluxes_sv,  only: init_fluxes, compute_flux_sw
use module_fluxes_sve, only: init_fluxes_exner, compute_flused_limiter, &
                             compute_flux_swe, compute_slope_effect
use module_bedload,    only: compute_bedload_solidfluxes
use module_suspension, only: compute_suspension_ceq, source_suspension
use module_bc,         only: compute_cdl, update_gc
use module_friction ,    only: source_friction
implicit none
real(kind=dp) :: h0   (0:nx+1)
real(kind=dp) :: q0   (0:nx+1)
real(kind=dp) :: zb0  (0:nx+1)
real(kind=dp) :: ht0  (0:nx+1,ntrac)
real(kind=dp) :: fluz_lim(1:nx)
real(kind=dp) :: dt1, dt2, aux1, aux2, gamma, fact
integer       :: i, k, iord, mord
! -----------------------------------------------------------------------------
! time order
!   if mord==1: Euler scheme (1st order)
!   if mord==2: Heun scheme (second order)
mord = 2
! -----------------------------------------------------------------------------
if (mord==2) then
  do i=0,nx+1
    h0(i)  = h(i)
    q0(i)  = q(i)
    zb0(i) = zb(i)
    if (ntrac>0) then
      do k=1,ntrac
        ht0(i,k) = ht(i,k)
      enddo
    endif
  enddo
endif
! -----------------------------------------------------------------------------
! Heun iterations (iord=1: first, iord=2: second)
! -----------------------------------------------------------------------------
do iord=1,mord
  call init_time_step(h, u)
  call init_fluxes
  call init_fluxes_exner()
  ! bedload solid flux
  if (bedload.eqv..True.) then
    call compute_bedload_solidfluxes(zf, zb, h, u)
  endif
  ! second order reconstructions
  if ((ordre.ge.2).or.(ordre_trac.ge.2)) then
    call compute_reconstructions(zb, h, u, q, t)
  endif
  ! If bedload, solve coupled SVE system
  if (bedload.eqv..True.) then
    call compute_flux_swe(zb, h, u, t)
  else 
    call compute_flux_sw(zb, h, u, t)
  endif
  ! If suspension, compute fluz_susp
  if (suspension.eqv..True.) then
    call compute_suspension_ceq(h, u)
  endif
  ! Slope effect
  call compute_slope_effect(zb)
  call compute_cdl
  call adjust_time_step
  if (iord==1) then
    dt1 = dt
  endif
  call compute_flused_limiter(zf, zb, fluz_lim)
  ! -----------------------------------------------------------------------------
  ! OMP DIRECTIVES
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  !$OMP parallel do  &
  !$OMP private(fact, k) &
  !$OMP shared(h, fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
  !$OMP        q, fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
  !$OMP        ht,flut_l, flut_r, flut_src_l, flut_src_r, &
  !$OMP        zb, fluqsed_l, fluqsed_r, zf, fluqsed_src_l, fluqsed_src_r, &
  !$OMP        Ceq_susp, flut_susp, fluz_susp, fluz_lim, &
  !$OMP        t, u, Rh, Lx, fric, ci, dt, ntrac, nx, porosite, rhos)
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  do i=1,nx
    fact = dt/ci(i)
    ! update of h, q
    h(i) = h(i) - fact*(fluh_l(i)-fluh_r(i-1) + fluh_src_l(i)+fluh_src_r(i-1))
    q(i) = q(i) - fact*(fluq_l(i)-fluq_r(i-1) + fluq_src_l(i)+fluq_src_r(i-1))
    ! update of ht
    if (ntrac>0) then
      do k=1,ntrac
        ht(i,k) = ht(i,k) - fact*( flut_l(i,k)  +flut_src_l(i,k) &
                                 - flut_r(i-1,k)+flut_src_r(i-1,k))
        ! erosion-deposition source term for tracer
        if (k.eq.1) then
          call source_suspension(zb(i), zf(i), h(i), Ceq_susp(i), &
                                 Csratio(i), t(i,k), ht(i,k), flut_susp(i))
        endif
      enddo
    endif
    ! erosion-deposition bed flux
    fluz_susp(i) = ci(i)*flut_susp(i)/((1.d0 - porosite)*rhos)
    ! update of zb
    zb(i) = zb(i) - fact*(min(fluqsed_l(i)   + fluqsed_src_l(i) &
                            - fluqsed_r(i-1) + fluqsed_src_r(i-1) &
                            + fluz_susp(i), fluz_lim(i)))
    ! sources
    call update_Rh(h(i), Lx(i), Rh(i))
    call source_friction(h(i), Rh(i), q(i), fric(i))
    ! update u
    call update_u(h(i), q(i), u(i))
  enddo
  !$OMP end parallel do
  call update_gc(h, q, u, ht, t)
enddo
! -----------------------------------------------------------------------------
! Heun update:
! -----------------------------------------------------------------------------
if (mord==2) then
  dt2=dt
  aux1 = 2d0*dt1*dt2
  aux2 = dt1 + dt2
  dt = aux1/aux2
  gamma = aux1/aux2**2
  ! OMP DIRECTIVES
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  !$OMP parallel do  &
  !$OMP private(k) &
  !$OMP shared(gamma, h, h0, q, q0, zb, zb0, ht, ht0, &
  !$OMP        t, u, Rh, ntrac )
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  do i=1,nx
    ! update of h, q
    h(i)  = (1d0-gamma)*h0(i)  + gamma*h(i)
    q(i)  = (1d0-gamma)*q0(i)  + gamma*q(i)
    zb(i) = (1d0-gamma)*zb0(i) + gamma*zb(i)
    ! update of ht
    if (ntrac>0) then
      do k=1,ntrac
        ht(i,k) = (1d0-gamma)*ht0(i,k) + gamma*ht(i,k)
        call update_t(h(i), ht(i,k), t(i,k))
      enddo
    endif
    ! update u
    call update_u(h(i), q(i), u(i))
  enddo
  !$OMP end parallel do
  call update_gc(h, q, u, ht, t)
endif
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



end module
! *****************************************************************************