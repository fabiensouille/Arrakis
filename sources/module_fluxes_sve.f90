! *****************************************************************************
module module_fluxes_sve
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for Shallow Water Exner equations fluxes computation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g, eps, eps_num
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                             SVE NUMERICAL FLUX
!                            (FOR COUPLED SCHEMES)
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine init_fluxes_exner()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initialisation of fluxes
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: fluqsed_l, fluqsed_r, fluqsed_src_l, fluqsed_src_r, &
                       fluz_susp, flut_susp, Ceq_susp
implicit none
! -----------------------------------------------------------------------------
fluqsed_l     = 0d0
fluqsed_r     = 0d0
fluqsed_src_l = 0d0
fluqsed_src_r = 0d0
fluz_susp     = 0d0
flut_susp     = 0d0
Ceq_susp      = 0d0
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_flux_swe(zb, h, u, t)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: hydrostarec, ntrac, eps, mlamb, &
                       ordre, ordre_trac, rec_vel, variable_width, dx, &
                       diffusion_u, difu, diffusion_t, nu_t, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r, &
                       corr_bp, corr_bm, corr_hp, corr_hm, &
                       corr_up, corr_um, corr_tp, corr_tm, &
                       qsed, fluqsed_l, fluqsed_r, fluqsed_src_l, fluqsed_src_r
use module_fluxes_sv,  only: flux_src_zb
implicit none
real(kind=dp), intent(in) :: zb (0:nx+1)
real(kind=dp), intent(in) :: h  (0:nx+1)
real(kind=dp), intent(in) :: u  (0:nx+1)
real(kind=dp), intent(in) :: t  (0:nx+1,ntrac)
real(kind=dp)             :: hg, hd, zbg, zbd, ug, ud, qg, qd, tg, td
real(kind=dp)             :: qsg, qsd, dzp, dzm, hip, him, etag, etad, zip12
real(kind=dp)             :: fhg, fhd, fqg, fqd, ft, ft_diff, fu_diff
real(kind=dp)             :: fhg_src, fhd_src, fqg_src, fqd_src
real(kind=dp)             :: fqsg_src, fqsd_src, fqsg, fqsd, mc
real(kind=dp)             :: lambda_roe (3), A_roe (3,3), D_roe (3,3)
integer                   :: i, k
logical                   :: rec_flag
! -----------------------------------------------------------------------------
! Loop on cell interfaces
! ~~~~~~~~~~~~~~~~~~~~~~~
! OMP DIRECTIVES
! ~~~~~~~~~~~~~~~~~~~~~~~
!$OMP parallel do  &
!$OMP private(k, rec_flag, zbg, zbd, hg, hd, ug, ud, qg, qd, &
!$OMP         him, hip, dzm, dzp, etag, etad, zip12, fu_diff, &
!$OMP         fhg, fhd, fqg, fqd, fhg_src, fhd_src, fqg_src, &
!$OMP         fqd_src, qsg, qsd, fqsg, fqsd, fqsg_src, fqsd_src, &
!$OMP         lambda_roe, A_roe, D_roe, mc, tg, td, ft, ft_diff) &
!$OMP shared(fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        fluqsed_l, fluqsed_r, fluqsed_src_l, fluqsed_src_r, &
!$OMP        corr_bp, corr_bm, corr_hp, corr_hm, &
!$OMP        corr_up, corr_um, corr_tp, corr_tm, &
!$OMP        diffusion_u, difu, diffusion_t, nu_t, &
!$OMP        ntrac, rec_vel, variable_width, ordre, ordre_trac, &
!$OMP        hydrostarec, zb, h, u, t, dx, qsed) &
!$OMP reduction(max:mlamb)
! ~~~~~~~~~~~~~~~~~~~~~~~
do i=1,nx-1

  ! left and right states
  ! ---------------------
  zbg = zb(i)
  zbd = zb(i+1)
  hg = h(i)
  hd = h(i+1)
  ug = u(i)
  ud = u(i+1)
  qsg = qsed(i)
  qsd = qsed(i+1)

  ! second order reconstruction
  ! ---------------------------
  if (ordre.ge.2) then
    if ( zbg.lt.(hd+zbd) .and. zbd.lt.(hg+zbg) &
        .and. 2.d0*abs(corr_bp(i)).lt.hd   &
        .and. 2.d0*abs(corr_bp(i)).lt.hg   &
        .and. 2.d0*abs(corr_bm(i+1)).lt.hd &
        .and. 2.d0*abs(corr_bm(i+1)).lt.hg ) then
      ! not a wet dry interface: reconstructions
      rec_flag = .true.
      zbg = zbg + corr_bp(i)
      zbd = zbd + corr_bm(i+1)
      hg = hg + corr_hp(i)
      hd = hd + corr_hm(i+1)
      ! reconstruction of velocity or discharge
      if ((hd.gt.0.0d0).and.(hg.gt.0.0d0)) then
        if (rec_vel.eqv..true.) then
          ! velocity reconstruction
          ug = ug + corr_up(i)
          ud = ud + corr_um(i+1)
        else
          ! discharge reconstruction
          qg = h(i)*ug   + corr_up(i)
          qd = h(i+1)*ud + corr_um(i+1)
          ug = qg/hg
          ud = qd/hd
        endif
      else
        hg = 0.0d0
        hd = 0.0d0
        ug = 0.0d0
        ud = 0.0d0
      endif
      ! adapt time step
      ! recompute lambda max from rec variables
      mlamb = max(mlamb, max(abs(ug)+sqrt(g*hg), abs(ud)+sqrt(g*hd)))
    else
      ! wet/dry interface: stay first order
      rec_flag = .false.
    endif
  endif

  ! hydrostatic reconstruction
  ! --------------------------
  ! Audusse & al. A fast and stable well-balanced scheme
  ! with hydrostatic reconstruction for shallow water flows
  if (hydrostarec.eq.1) then
    dzm = max(0.d0, zbd- zbg)
    him = max(0.d0, hg - dzm)
    dzp = max(0.d0, zbg- zbd)
    hip = max(0.d0, hd - dzp)
  ! Chen and Noelle. A new hydrostatic reconstruction
  ! scheme based on subcell reconstructions
  elseif (hydrostarec.eq.2) then
    etag = hg + zbg
    etad = hd + zbd
    zip12 = min(max(zbg, zbd), min(etag, etad))
    him = min(etag - zip12, hg)
    hip = min(etad - zip12, hd)
  else
    him = hg
    hip = hd
  endif

  ! Flux computation SVE
  ! --------------------
  call flux_sve(mc, fhg, fhd, fqg, fqd, fqsg, fqsd, i, him, hip, &
                ug, ud, zbg, zbd, qsg, qsd, lambda_roe, A_roe, D_roe)
  mlamb = max(mlamb, mc)
  fluh_l(i) = fhg
  fluh_r(i) = fhd
  fluq_l(i) = fqg
  fluq_r(i) = fqd
  fluqsed_l(i) = fqsg
  fluqsed_r(i) = fqsd

  ! bathy source (scheme for Sb=-gh.grad(zb))
  ! ------------
  call flux_src_zb_sve(fhg_src, fhd_src, fqg_src, fqd_src, &
                       fqsg_src, fqsd_src, hg, hd, him, hip, zbg, zbd, &
                       lambda_roe, A_roe, D_roe)
  fluh_src_l(i) = fhg_src
  fluh_src_r(i) = fhd_src
  fluq_src_l(i) = fqg_src
  fluq_src_r(i) = fqd_src
  fluqsed_src_l(i) = fqsg_src
  fluqsed_src_r(i) = fqsd_src

  ! second order correction (for hydrostarec!=0)
  if (hydrostarec.ne.0 .and. ordre.ge.2 .and. rec_flag.eqv..true.) then
    fluq_src_l(i) = fluq_src_l(i) + 0.5d0*g*(hg+h(i)  )*corr_bp(i)
    fluq_src_r(i) = fluq_src_r(i) - 0.5d0*g*(hd+h(i+1))*corr_bm(i+1)
  endif
    
  ! width source (scheme for SL=-(hu, hu**2, qs).grad(L)/L)
  ! ------------
  if (variable_width.eqv..True.) then
    call flux_src_width_sve(fhg_src, fhd_src, fqg_src, fqd_src, &
                            fqsg_src, fqsd_src, i, hg, hd, ug, ud, &
                            qsg, qsd, lambda_roe, A_roe, D_roe)

    fluh_src_l(i) = fluh_src_l(i) + fhg_src
    fluh_src_r(i) = fluh_src_r(i) + fhd_src
    fluq_src_l(i) = fluq_src_l(i) + fqg_src
    fluq_src_r(i) = fluq_src_r(i) + fqd_src
    fluqsed_src_l(i) = fluqsed_src_l(i) + fqsg_src
    fluqsed_src_r(i) = fluqsed_src_r(i) + fqsd_src
  endif

  ! diffusion flux
  ! --------------
  if (diffusion_u.eqv..True.) then
    fu_diff = 0.25d0*(hg+hd)*(difu(i+1)+difu(i))*(u(i+1)-u(i))/dx(i)
    fluq_l(i) = fluq_l(i) - fu_diff
    fluq_r(i) = fluq_r(i) - fu_diff
  endif

  ! tracers flux
  ! ------------
  if (ntrac.gt.0) then
    do k=1,ntrac
      ! second order reconstructions
      if (ordre_trac.eq.1) then
        tg = t(i,k)
        td = t(i+1,k)
      elseif (ordre_trac.ge.2) then
        tg = t(i,k)   + corr_tp(i,k)
        td = t(i+1,k) + corr_tm(i+1,k)
      endif
      ! upwind flux
      if (fhg.ge.0d0) then
        ft = tg*fhg
      else
        ft = td*fhg
      endif
      ! advection flux
      flut_l(i,k) = flut_l(i,k) + ft
      flut_r(i,k) = flut_r(i,k) + ft

      ! width source term
      if (variable_width.eqv..True.) then
        flut_src_l(i,k) = flut_src_l(i,k) + fhg_src*tg
        flut_src_r(i,k) = flut_src_r(i,k) + fhd_src*td
      endif

      ! diffusion flux
      if (diffusion_t.eqv..True.) then
        ft_diff = 0.5d0*(hg+hd)*nu_t(k)*(t(i+1,k)-t(i,k))/dx(i)
        flut_l(i,k) = flut_l(i,k) - ft_diff
        flut_r(i,k) = flut_r(i,k) - ft_diff
      endif
    enddo
  endif

end do
!$OMP end parallel do 
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_slope_effect(zb)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dx, diffusion_z, nu_z, slope_cr, fluqsed_l, fluqsed_r
implicit none
real(kind=dp), intent(in) :: zb (0:nx+1)
real(kind=dp)             :: slope, fz_diff
integer                   :: i
! -----------------------------------------------------------------------------
! Loop on cell interfaces
do i=1,nx-1
  ! slope diffusion flux
  ! --------------------
  if (diffusion_z.eqv..True.) then
    slope = (zb(i+1)-zb(i))/dx(i)
    if (abs(slope).lt.slope_cr) then
      fz_diff = 0d0
    else
      fz_diff = nu_z*(abs(slope)-slope_cr)*sign(1d0, slope)
    endif
    fluqsed_l(i) = fluqsed_l(i) - fz_diff
    fluqsed_r(i) = fluqsed_r(i) - fz_diff
  endif
end do
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                     NUMERICAL SVE FLUXES COMPUTATION
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine flux_sve(mc, f1g, f1d, f2g, f2d, f3g, f3d, i, hg, hd, ug, ud, &
                    zg, zd, qsg, qsd, lambda_roe, A_roe, D_roe)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Flux selector : return the flux depending on the selected method
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: method, swe_scheme_opt
use module_fluxes_sv, only: flux_roe, flux_hll, flux_acu_sw, flux_cin
implicit none
real(kind=dp), intent(out)   :: mc       ! max wave speed
real(kind=dp), intent(inout) :: f1g, f1d, f2g, f2d, f3g, f3d
real(kind=dp), intent(inout) :: lambda_roe(3), A_roe(3,3), D_roe(3,3)
real(kind=dp), intent(in)    :: hg, hd, ug, ud, zg, zd, qsg, qsd
integer,       intent(in)    :: i
! -----------------------------------------------------------------------------
! Coupled schemes
! ---------------
if (swe_scheme_opt.eq.0) then
  ! *******
  ! Roe SVE
  ! *******
  if (method.eq.1) then
    call flux_roe_sve(mc, f1g, f1d, f2g, f2d, f3g, f3d, i, hg, hd, ug, ud, &
                      zg, zd, qsg, qsd, lambda_roe, A_roe, D_roe)
  ! *******
  ! HLL SVE
  ! *******
  elseif (method.eq.2) then
    call flux_hll_sve(mc, f1g, f2g, f3g, i, hg, hd, ug, ud, zg, zd, qsg, qsd)
    f1d = f1g
    f2d = f2g
    f3d = f3g
  ! *******
  ! ACU SVE
  ! *******
  elseif (method.eq.3) then
    call flux_acu_sve(mc, f1g, f1d, f2g, f2d, f3g, f3d, i, hg, hd, ug, ud, &
                      zg, zd, qsg, qsd)
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown BEDLOAD COUPLED SCHEME '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
! -----------------------------------------------------------------------------
! De-coupled schemes
! ------------------
elseif (swe_scheme_opt.ne.0) then
  ! 
  ! HYDRO FLUXES (SAINT-VENANT ONLY)
  ! 
  ! ******
  ! Roe SV 
  ! ******
  if (method.eq.1) then
    call flux_roe(mc, f1g, f1d, f2g, f2d, ug, ud, hg, hd)
  ! ******
  ! HLL SV
  ! ******
  elseif (method.eq.2) then
    call flux_hll(mc, f1g, f2g, ug, ud, hg, hd)
    f1d = f1g
    f2d = f2g
  ! ******
  ! ACU SV
  ! ******
  elseif (method.eq.3) then
    call flux_acu_sw(mc, f1g, f2g, f2d, ug, ud, hg, hd, zg, zd)
    f1d = f1g
  ! **********
  ! Kinetic SV
  ! **********
  elseif (method.eq.4) then
    call flux_cin(mc, f1g, f2g, ug, ud, hg, hd)
    f1d = f1g
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Unknown ADVECTION SCHEME '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  ! 
  ! BEDLOAD FLUXES (EXNER ONLY)
  ! 
  ! ********
  ! Centered
  ! ********
  if (swe_scheme_opt.eq.1) then
    call flux_bedload_centered(f3d, qsg, qsd)
    f3g = f3d
  ! *********************
  ! Upwind (with c=dq/dz)
  ! *********************
  elseif (swe_scheme_opt.eq.2) then
    call flux_bedload_upwind(f3d, qsg, qsd, zg, zd)
    f3g = f3d
  ! ********************
  ! Upwind (with Froude)
  ! ********************
  elseif (swe_scheme_opt.eq.3) then
    call flux_bedload_upwind_fr(f3d, qsg, qsd, ug, ud, hg, hd)
    f3g = f3d
  ! *********
  ! HLL Exner
  ! *********
  elseif (swe_scheme_opt.eq.4) then
    call flux_bedload_hll(f3g, i, qsg, qsd, ug, ud, hg, hd, zg, zd)
    f3d = f3g
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Unknown BEDLOAD SCHEME OPTION '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_hll_sve(mc, fh, fq, fqs, i, hg, hd, ug, ud, zg, zd, qsedg, qsedd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! HLL-SVE scheme 
! A 2D HLL-based weakly coupled model for transient flows on mobile beds.
! Meurice, Robin and Soares-Frazao, Sandra
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_fluxes_sv, only: Fluxh, Fluxq
implicit none
real(kind=dp), intent(in)    :: hg, hd, ug, ud, zg, zd, qsedg, qsedd
integer,       intent(in)    :: i
real(kind=dp), intent(inout) :: fh, fq, fqs
real(kind=dp), intent(out)   :: mc
real(kind=dp)                :: c1, c2, c3
real(kind=dp)                :: lambg(3), lambd(3)
! -----------------------------------------------------------------------------
! SVE polynomial roots
call sve_waves(i,   hg, ug, lambg, 0)
call sve_waves(i+1, hd, ud, lambd, 0)
c1 = min(lambg(1), lambd(1))
c2 = max(lambg(2), lambd(2))
c3 = max(lambg(3), lambd(3))
mc = max(c1, max(c2, c3))

! Hydro HLL flux
if (c1.ge.0d0) then
  fh = Fluxh(hg,ug)
  fq = Fluxq(hg,ug)
elseif ((c1.lt.0d0).and.(0d0.lt.c3)) then
  fh = (c3*Fluxh(hg,ug)-c1*Fluxh(hd,ud))/(c3-c1) + c1*c3*(hd-hg)/(c3-c1)
  fq = (c3*Fluxq(hg,ug)-c1*Fluxq(hd,ud))/(c3-c1) + c1*c3*(hd*ud-hg*ug)/(c3-c1)
else
  fh = Fluxh(hd,ud)
  fq = Fluxq(hd,ud)
endif

! Sediment HLL flux
if (c1.ge.0d0) then
  fqs = qsedg
elseif ((c1.lt.0d0).and.(0d0.lt.c2)) then
  fqs = (c2*qsedg -c1*qsedd)/(c2-c1) + c1*c2*(zd-zg)/(c2-c1)
else
  fqs = qsedd
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_acu_sve(mc, fhg, fhd, fqg, fqd, fqsg, fqsd, i, hg, hd, ug, ud,& 
                        zg, zd, qsedg, qsedd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ACU-SVE intermediate states
! A simple three-wave approximate Riemann solver for the SVE equations. 
! E. AUDUSSE, C. CHALONS AND P. UNG
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_fluxes_sv, only: Fluxh, Fluxq
implicit none
real(kind=dp), intent(in)    :: hg, hd, ug, ud, zg, zd
real(kind=dp), intent(in)    :: qsedg, qsedd
integer,       intent(in)    :: i
real(kind=dp), intent(inout) :: fhg, fhd, fqg, fqd, fqsg, fqsd
real(kind=dp), intent(out)   :: mc
real(kind=dp)                :: c1, c2, c3, hs1, hs2, qs1, qs2, zs1, zs2
real(kind=dp)                :: dqs, aux, hsa, hhll, qhll, src, dz
real(kind=dp)                :: lambg(3), lambd(3)
! -----------------------------------------------------------------------------
! SVE polynomial roots
call sve_waves(i,   hg, ug, lambg, 0)
call sve_waves(i+1, hd, ud, lambd, 0)
c1 = min(lambg(1), lambd(1))
c2 = max(lambg(2), lambd(2))
c3 = max(lambg(3), lambd(3))
c2 = max(c2, c3)
mc = max(c1, c2)

! bottom intermediate state
dqs = qsedd - qsedg
aux = max(c1**2.d0 + c2**2.d0, eps_num)
zs1 = zg + (c1/aux)*dqs
zs2 = zd - (c2/aux)*dqs

! water depth intermediate state
aux = max(c2 - c1, eps_num)
hhll = (c2*hd-c1*hg)/aux + (Fluxh(hg,ug)-Fluxh(hd,ud))/aux
hs1 = hhll + (c2/aux)*(zs2 - zs1)
hs2 = hhll + (c1/aux)*(zs2 - zs1)

! positivity of water depth
if ((zs2 - zs1).ge.0.d0) then
  hsa = hs2
  hs2 = max(hs2, 0.d0)
  hs1 = hs1 - (c2/max(c1, eps_num))*(hsa-hs2)
else
  hsa = hs1
  hs1 = max(hs1, 0.d0)
  hs2 = hs2 - (c1/max(c2, eps_num))*(hsa-hs1)
endif

! bottom gradient reconstruction
dz = zd - zg
if (dz.ge.0.d0) then
  src = -0.5d0*(hg+hd)*min(max(hg, hd), dz)  !Eq. 35a
else
  src = -0.5d0*(hg+hd)*max(-max(hg, hd), dz) !Eq. 35b
endif

! flowrate intermediate state
qhll = (c2*hd*ud-c1*hg*ug)/aux + (Fluxq(hg,ug)-Fluxq(hd,ud))/aux
qs1 = qhll + g*src/aux
qs2 = qhll + g*src/aux

! left and right fluxes
fhg  = c1*(hs1 - hg)
fhd  = c2*(hs2 - hd)
fqg  = c1*(qs1 - hg*ug)
fqd  = c2*(qs2 - hd*ud)
fqsg = c1*(zs1 - zg)
fqsd = c2*(zs2 - zd)
! -----------------------------------------------------------------------------
end subroutine  
! *****************************************************************************



! *****************************************************************************
subroutine flux_roe_sve(mc, f1g, f1d, f2g, f2d, f3g, f3d, i, hg, hd, ug, ud, &
                        zg, zd, qsg, qsd, lambda_roe, A_roe, D_roe)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Roe scheme for SVE
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_fluxes_sv, only: Fluxh, Fluxq
implicit none
real(kind=dp), intent(out)   :: mc
real(kind=dp), intent(inout) :: f1g, f1d, f2g, f2d, f3g, f3d
real(kind=dp), intent(inout) :: lambda_roe(3), A_roe(3,3), D_roe(3,3)
real(kind=dp), intent(in)    :: ug, ud, hg, hd, zg, zd, qsg, qsd
integer,       intent(in)    :: i
integer                      :: j
real(kind=dp)                :: jump_w(3)
integer, parameter           :: roe_method = 0
! -----------------------------------------------------------------------------
! variable difference vector
jump_w(1) = hd - hg
jump_w(2) = hd*ud - hg*ug
jump_w(3) = zd - zg

! **************************************************************
! Roe flux, with D_roe = |A|: f = 0.5*(fl+fr) - 0.5*|a|*(wr-wl))
! **************************************************************
if (roe_method.eq.0) then
  ! compute roe matrix D_roe=|A|=R*|lambda|*R1
  call aroe_sve(i, hg, hd, ug, ud, 0, lambda_roe, A_roe, D_roe, mc)
  ! evaluate roe's flux: f=(fl+fr-D_roe*(wr-wl))/2
  f1g = Fluxh(hg,ug) + Fluxh(hd,ud)
  f2g = Fluxq(hg,ug) + Fluxq(hd,ud)
  f3g = qsg + qsd
  do j=1,3
    f1g = f1g - D_roe(1,j)*jump_w(j)
    f2g = f2g - D_roe(2,j)*jump_w(j)
    f3g = f3g - D_roe(3,j)*jump_w(j)
  enddo
  f1g = 0.5d0*f1g
  f2g = 0.5d0*f2g
  f3g = 0.5d0*f3g
  f1d = f1g
  f2d = f2g
  f3d = f3g

! **************************************************
! Roe flux, with D_roe = A- : fg = (fl + a-*(wr-wl))
! **************************************************
elseif (roe_method.eq.1) then
  ! compute roe matrix D_roe=A-=R*(Lambda-)*R1
  call aroe_sve(i, hg, hd, ug, ud, -1, lambda_roe, A_roe, D_roe, mc)
  ! evaluate roe's flux
  f1g = Fluxh(hg,ug)
  f2g = Fluxq(hg,ug)
  f3g = qsg
  do j=1,3
    f1g = f1g + D_roe(1,j)*jump_w(j)
    f2g = f2g + D_roe(2,j)*jump_w(j)
    f3g = f3g + D_roe(3,j)*jump_w(j)
  enddo
  f1d = f1g
  f2d = f2g
  f3d = f3g

! **************************************************
! Roe flux, with D_roe = A+ : fd = (fr - a+*(wr-wl))
! **************************************************
elseif (roe_method.eq.2) then
  ! compute roe matrix D_roe=A+=R*(Lambda+)*R1
  call aroe_sve(i, hg, hd, ug, ud, +1, lambda_roe, A_roe, D_roe, mc)
  ! evaluate roe's flux
  f1d = Fluxh(hd,ud)
  f2d = Fluxq(hd,ud)
  f3d = qsd
  do j=1,3
    f1d = f1d - D_roe(1,j)*jump_w(j)
    f2d = f2d - D_roe(2,j)*jump_w(j)
    f3d = f3d - D_roe(3,j)*jump_w(j)
  enddo
  f1g = f1d
  f2g = f2d
  f3g = f3d
! -----------------------------------------------------------------------------
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Unknown Roe flux method, check flux_roe '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine aroe_sve(face, hg, hd, ug, ud, matrix_sign, &
                    lambda_roe, A_roe, D_roe, mc)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! computes Roe's matrix for the SVE system (A_roe and D_roe)
! matrix_sign =  0: D_roe = |A|
! matrix_sign = -1: D_roe =  A-
! matrix_sign = +1: D_roe =  A+
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: write_listing
use module_linalg, only: invert_3x3_matrix
implicit none
integer, intent(in)          :: face, matrix_sign
real(kind=dp), intent(in)    :: hg, hd, ug, ud
real(kind=dp), intent(inout) :: lambda_roe(3), A_roe(3,3), D_roe(3,3)
real(kind=dp), intent(out)   :: mc
real(kind=dp)                :: R(3,3), R1(3,3), Lamb(3)
real(kind=dp)                :: hc, uc, cc, cc2, hge, hde, alphat, betat
integer                      :: i, j, k
logical                      :: flag
! -----------------------------------------------------------------------------
! dry state correction
hge = max(eps_num, hg)
hde = max(eps_num, hd)
! Roe's average for shallow water
hc = (hge + hde)/2.d0
uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
cc2 = g*hc
cc = dsqrt(cc2)

! eigenvalues
! ***********
call sve_waves(face, hc, uc, lambda_roe, 1, alphat, betat)
! A_roe matrix (store for source terms)
A_roe(1,1) = 0d0
A_roe(1,2) = 1d0
A_roe(1,3) = 0d0
A_roe(2,1) = (cc**2)-(uc**2)
A_roe(2,2) = 2d0*uc
A_roe(2,3) = cc**2
A_roe(3,1) = alphat
A_roe(3,2) = betat
A_roe(3,3) = 0d0

! diagonal of Lambda : Lamb = |L|, L+ or L-
! *****************************************
do k=1,3
  if (matrix_sign.eq.-1) then
    Lamb(k) = min(lambda_roe(k), 0d0)
  else if (matrix_sign.eq.1) then
    Lamb(k) = max(lambda_roe(k), 0d0)
  else
    Lamb(k) = dabs(lambda_roe(k))
  endif
enddo

! diagonalization matrix R
! ************************
R(1,1) = 1.d0
R(1,2) = 1.d0
R(1,3) = 1.d0
R(2,1) = lambda_roe(1)
R(2,2) = lambda_roe(2)
R(2,3) = lambda_roe(3)
R(3,1) = (((uc - lambda_roe(1))**2)/cc2) - 1.d0
R(3,2) = (((uc - lambda_roe(2))**2)/cc2) - 1.d0
R(3,3) = (((uc - lambda_roe(3))**2)/cc2) - 1.d0
! inverse R-1 
call invert_3x3_matrix(R, R1, flag)
if ((flag.eqv..false.).and.(write_listing)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) "WARNING: could not compute inverse ", &
             "of eigen vectors matrix in Roe scheme "
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'             
endif

! matrix*diag*matrix multiply
! ***************************
do i=1,3
  do j=1,3
    D_roe(i,j) = 0.d0
    do k=1,3
      ! matrix D_roe=|A|=R|L|R-1, A+=R(L+)R-1 or A-=R(L-)R-1
      D_roe(i,j) = D_roe(i,j) + R(i,k)*Lamb(k)*R1(k,j)
    enddo
  enddo
enddo
! -----------------------------------------------------------------------------
mc = max(dabs(lambda_roe(1)), max(dabs(lambda_roe(2)), dabs(lambda_roe(3))))
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                     SOURCE FLUX AT INTERFACE  [Wg | Wd]
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine flux_src_zb_sve(f1g, f1d, f2g, f2d, f3g, f3d, hg, hd, him, hip, &
                           zbg, zbd, lambda_roe, A_roe, D_roe)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Scheme for the bathy source term Sb=-gh grad(zb)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,       only: g, hydrostarec
implicit none
real(kind=dp), intent(out) :: f1g, f1d, f2g, f2d, f3g, f3d
real(kind=dp), intent(in)  :: hg, hd, him, hip, zbg, zbd
real(kind=dp), intent(in)  :: lambda_roe(3), A_roe(3,3), D_roe(3,3)
real(kind=dp)              :: etag, etad, zip12, hc, src_centered
real(kind=dp)              :: Mp3(3,3), Mm3(3,3)
! -----------------------------------------------------------------------------
! initialization    
f1g = 0d0
f1d = 0d0
f2g = 0d0
f2d = 0d0
f3g = 0d0
f3d = 0d0
! -----------------------------------------------------------------------------
! Audusse & al. (Hydrostatic reconstruction)
! ******************************************
if (hydrostarec.eq.1) then
  f2g = + 0.5d0*g*(hg+him)*(hg-him)
  f2d = - 0.5d0*g*(hd+hip)*(hd-hip)

! Chen and Noelle
! ***************
elseif (hydrostarec.eq.2) then
  etag = hg + zbg
  etad = hd + zbd
  zip12 = min(max(zbg, zbd), min(etag, etad))
  f2g = + 0.5d0*g*(hg+him)*(zip12-zbg)
  f2d = - 0.5d0*g*(hd+hip)*(zip12-zbd)

! Centered (not well-ballanced)
! *****************************
elseif (hydrostarec.eq.3) then
  f2g = + 0.5d0*g*hg*(zbd-zbg)
  f2d = + 0.5d0*g*hd*(zbd-zbg)

! Centered Roe
! ************
elseif (hydrostarec.eq.4) then
  hc = 0.5d0*(hg+hd)
  src_centered = 0.5d0*g*hc*(zbd-zbg)
  f2g = src_centered
  f2d = src_centered

! Upwind Roe (Vazquez-Cendon)
! ***************************
elseif (hydrostarec.eq.5) then
  ! computes upwind matrices M+ and M-
  call src_upwind_roe_sve(lambda_roe, A_roe, D_roe, Mp3, Mm3)
  ! update source
  hc = 0.5d0*(hg+hd)
  src_centered = 0.5d0*g*hc*(zbd-zbg)
  f1g = src_centered*Mm3(1,2)
  f1d = src_centered*Mp3(1,2)
  f2g = src_centered*Mm3(2,2)
  f2d = src_centered*Mp3(2,2)
  f3g = src_centered*Mm3(3,2)
  f3d = src_centered*Mp3(3,2)

elseif (hydrostarec.gt.5) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : Unknown HYDROSTATIC RECONSTRUCTION '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_src_width_sve(f1g, f1d, f2g, f2d, f3g, f3d, i, hg, hd, ug, ud,& 
                              qsg, qsd, lambda_roe, A_roe, D_roe)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Scheme for the variable width source term (for SVE)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,      only: width_scheme, dx, Lx, dLx, nx
use module_fluxes_sv, only: src_upwind_roe_sv
implicit none
real(kind=dp), intent(out) :: f1g, f1d, f2g, f2d, f3g, f3d
real(kind=dp), intent(in)  :: hg, hd, ug, ud, qsg, qsd
real(kind=dp), intent(in)  :: lambda_roe(3), A_roe(3,3), D_roe(3,3)
integer,       intent(in)  :: i
real(kind=dp)              :: src_Lg, src_Ld, dLip12m, dLip12p
real(kind=dp)              :: Mp3(3,3), Mm3(3,3), Sg(3), Sd(3)
real(kind=dp)              :: hge, hde, hc, uc
integer, parameter         :: var_opt = 0
! -----------------------------------------------------------------------------
! initialization    
f1g = 0d0
f1d = 0d0
f2g = 0d0
f2d = 0d0
f3g = 0d0
f3d = 0d0
! -----------------------------------------------------------------------------
! computation of L'/L at i and i+1
! 
! 1st order FD
! ************
if(width_scheme.eq.1) then
  ! backward FD
  dLip12m = dLx(i-1)
  dLip12p = dLx(i)
  ! forward FD
  dLip12m = dLx(i)
  dLip12p = dLx(i+1)

! 2nd order centered FD
! *********************
elseif((width_scheme.eq.2).or.(width_scheme.eq.5)) then
  dLip12m = dLx(i)
  dLip12p = dLx(i)

! 3rd order FD
! ************
elseif(width_scheme.eq.3) then
  ! L'/L at i and i+1
  dLip12m = 0.5d0*(dLx(i-1) + dLx(i))
  dLip12p = 0.5d0*(dLx(i) + dLx(i+1))

! 4th order centered FD (work only for regular mesh i.e. dx=cst)
! *********************
elseif(width_scheme.eq.4) then
  dLip12m = (Lx(max(0,i-2))-8d0*Lx(max(0,i-1)) &
          + 8d0*Lx(min(i+1,nx+1))-Lx(min(i+2,nx+1)))/(12d0*dx(i))
  dLip12p = (Lx(max(0,i-2))-8d0*Lx(max(0,i-1)) &
          + 8d0*Lx(min(i+1,nx+1))-Lx(min(i+2,nx+1)))/(12d0*dx(i))

elseif(width_scheme.gt.5) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : Unknown VARIABLE WIDTH SCHEME '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
! Upwind Roe (Vazquez-Cendon)
! ***************************
if(width_scheme.eq.5) then
  ! compute 0.5(x_{i+1}-x_{i})*(L'/L)_i+1/2-
  ! and     0.5(x_{i+1}-x_{i})*(L'/L)_i+1/2+
  src_Lg = dx(i)*dLip12m/(Lx(i)+Lx(i+1))
  src_Ld = dx(i)*dLip12p/(Lx(i)+Lx(i+1))
  ! update source flux (theta * src)
  ! we use Roe's average states
  hge = max(eps_num, hg)
  hde = max(eps_num, hd)
  hc = (hge+hde)/2.d0
  uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
  f1g = src_Lg*hc*uc
  f1d = src_Ld*hc*uc
  f2g = src_Lg*hc*uc**2
  f2d = src_Ld*hc*uc**2
  f3g = src_Lg*qsg
  f3d = src_Ld*qsd
  ! computes upwind matrices M+ and M-
  call src_upwind_roe_sve(lambda_roe, A_roe, D_roe, Mp3, Mm3)
  Sg(1) = f1g
  Sg(2) = f2g
  Sg(3) = f3g
  Sd(1) = f1d
  Sd(2) = f2d
  Sd(3) = f3d
  ! upwinding the src
  f1g = Sg(1)*Mm3(1,1) + Sg(2)*Mm3(1,2) + Sg(3)*Mm3(1,3)
  f1d = Sd(1)*Mp3(1,1) + Sd(2)*Mp3(1,2) + Sd(3)*Mp3(1,3)
  f2g = Sg(1)*Mm3(2,1) + Sg(2)*Mm3(2,2) + Sg(3)*Mm3(2,3)
  f2d = Sd(1)*Mp3(2,1) + Sd(2)*Mp3(2,2) + Sd(3)*Mp3(2,3)
  f3g = Sg(1)*Mm3(3,1) + Sg(2)*Mm3(3,2) + Sg(3)*Mm3(3,3)
  f3d = Sd(1)*Mp3(3,1) + Sd(2)*Mp3(3,2) + Sd(3)*Mp3(3,3)
! ***************************
else
  ! compute 0.5(x_{i+1}-x_{i})*(L'/L)_i+1/2-
  ! and     0.5(x_{i+1}-x_{i})*(L'/L)_i+1/2+
  src_Lg = 0.5d0*dx(i)*dLip12m/Lx(i)
  src_Ld = 0.5d0*dx(i)*dLip12p/Lx(i+1)
  ! update source flux (theta * src)
  if(var_opt.eq.1) then
    ! we use Roe's average states
    hge = max(eps_num, hg)
    hde = max(eps_num, hd)
    hc = (hge+hde)/2.d0
    uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
    f1g = src_Lg*hc*uc
    f1d = src_Ld*hc*uc
    f2g = src_Lg*hc*uc**2
    f2d = src_Ld*hc*uc**2
  else
    f1g = src_Lg*hg*ug
    f1d = src_Ld*hd*ud
    f2g = src_Lg*hg*ug**2
    f2d = src_Ld*hd*ud**2
  endif
  f3g = src_Lg*qsg
  f3d = src_Ld*qsd
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine src_upwind_roe_sve(lambda_roe, A_roe, D_roe, Mp, Mm)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the upwind matrices for Roe scheme, returns
! M+- = I+-|A|A^-1 where A is the Roe matrix.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_bedload, only: compute_bedload_derivatives
use module_linalg,  only: matrix_product, invert_3x3_matrix
implicit none
real(kind=dp), intent(in)  :: lambda_roe(3)
real(kind=dp), intent(in)  :: A_roe(3,3), D_roe(3,3)
real(kind=dp), intent(out) :: Mm(3,3), Mp(3,3)
real(kind=dp)              :: A(3,3), A1(3,3)
integer                    :: i, j
logical                    :: flag
! ----------------------------------------------------------------------------- 
! centered scheme if null eigenvalue (M+-=I)
if ((lambda_roe(1).eq.0).or. &
&    (lambda_roe(2).eq.0).or. &
&    (lambda_roe(3).eq.0)) then
  do i=1,3
    do j=1,3
      if (i.eq.j) then
        Mm(i,j) = 1d0
        Mp(i,j) = 1d0
      else
        Mm(i,j) = 0d0
        Mp(i,j) = 0d0
      endif
    enddo
  enddo
! compute upwinding matrix: M+- = I+-|A|A^-1
else
  call invert_3x3_matrix(A_roe, A1, flag) ! A1 <- A_roe^-1
  call matrix_product(3, D_roe, A1, A)    ! A  <- |A|A1
  ! compute M+-
  do i=1,3
    do j=1,3
      Mm(i,j) = -A(i,j)
      Mp(i,j) = +A(i,j)
      if (i.eq.j) then
        Mm(i,j) = Mm(i,j) + 1d0
        Mp(i,j) = Mp(i,j) + 1d0
      endif
    enddo
  enddo
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                     SVE SYSTEM WAVE SPEEDS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine sve_waves(i, h, u, lambda, m, at, bt)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Returns SVE wave speeds 
! -> on cell i (if m=0) 
! -> on interface i between cells i and i+1 (if m=1) (For Roe)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: porosite, swe_roots_approx
use module_bedload, only: compute_bedload_derivatives
implicit none
integer, intent(in)                  :: i, m
real(kind=dp), intent(in)            :: u, h
real(kind=dp), intent(out)           :: lambda(3)
real(kind=dp), intent(out), optional :: at, bt
real(kind=dp)                        :: alphat, betat
! -----------------------------------------------------------------------------
! compute alpha and beta (sediment flux derivatives)
call compute_bedload_derivatives(i, h, u, alphat, betat, m)
alphat = alphat/(1.d0-porosite)
betat  = betat /(1.d0-porosite)
! optional outputs
if(present(at))then
  at = alphat
endif
if(present(bt))then
  bt = betat
endif
! ~~~~~~~~~
! Cardano's method for polynomial roots
! ~~~~~~~~~
if (swe_roots_approx.eq.1) then
  call cardano_waves(h, u, alphat, betat, lambda)
! ~~~~~~~~~~~~~~~~~~~~~~
! Soares-Frazao and Zech approximation of polynomial roots
! ~~~~~~~~~~~~~~~~~~~~~~
elseif (swe_roots_approx.eq.2) then
  call sfz_waves(h, u, betat, lambda)
! ~~~~~~~~
! Nickalls approximation of polynomial roots
! ~~~~~~~~
elseif (swe_roots_approx.eq.3) then
  call nickalls_waves(h, u, betat, lambda)
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : Unknown SVE WAVES APPROXIMATION                    '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine sfz_waves(h, u, betat, lambda)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Soares-Frazao and Zech approximation of SVE wave speeds.
! HLLC scheme with novel wave-speed estimators appropriate
! for two-dimensional shallow-water flow on erodible bed.
! S.Soares-Frazao and Y.Zech
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(in)    :: u, h, betat
real(kind=dp), intent(inout) :: lambda(3)
real(kind=dp)                :: c, aux
! -----------------------------------------------------------------------------
c = dsqrt(g*h)
aux = dsqrt((u-c)**2.d0  + 4.d0*betat*c**2.d0)
lambda(1) = 0.5d0*(u-c-aux)
lambda(2) = 0.5d0*(u-c+aux)
lambda(3) = u+c
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine nickalls_waves(h, u, betat, lambda)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Nickalls approximation of SVE wave speeds.
! A simple three-wave approximate Riemann solver for the
! Saint-Venant-Exner equations. E. AUDUSSE, C. CHALONS AND P. UNG
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(in)    :: u, h, betat
real(kind=dp), intent(inout) :: lambda(3)
real(kind=dp)                :: x0, Omega
real(kind=dp), parameter     :: tmp2s3 = 2.d0/3.d0
real(kind=dp), parameter     :: tmp1s3 = 1.d0/3.d0
! -----------------------------------------------------------------------------
x0 = tmp2s3*u
Omega = tmp1s3*dsqrt(u**2.d0 + 3.d0*g*h*(1.d0+betat))
lambda(1) = x0 - 2.d0*Omega  !Eq. 25
lambda(2) = 0.d0
lambda(3) = x0 + 2.d0*Omega  !Eq. 26
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine cardano_waves(h, u, alphat, betat, lambda)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! True SVE wave speeds determined with the Cardan method
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(in)    :: u, h, alphat, betat
real(kind=dp), intent(inout) :: lambda(3)
real(kind=dp)                :: cc2, a1, a2, a3
integer                      :: errflag
! -----------------------------------------------------------------------------
cc2 = g*h
a1 = -2.d0*u
a2 = u**2 - cc2*(1.d0 + betat)
a3 = -cc2*alphat
call cardan_roots(a1, a2, a3, lambda, errflag)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine cardan_roots(a1, a2, a3, lambda, flag)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! computes roots of third degree polynomial with cardan formula
! polynome: p(l) = l^3 + a1*l^2 + a2*l + a3
! https://fortran-lang.discourse.group/t/cardanos-solution-of-the-cubic-equation/111
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: pi, write_listing
implicit none
real(kind=dp), intent(in)  :: a1, a2, a3
real(kind=dp), intent(out) :: lambda(3)
integer,       intent(out) :: flag
real(kind=dp)              :: Q, R, D, theta, s, temp
complex(kind=dp)           :: omega, z2, z3
integer                    :: k
logical, parameter         :: verbose = .True.
! -----------------------------------------------------------------------------
Q = (a1**2 - 3.d0*a2)/9.d0
R = (2.d0*a1**3 - 9.d0*a1*a2 + 27.d0*a3)/54.d0
D = Q**3 + R**2
! -----------------------------------------------------------------------------
! real roots
if (D.ge.0.d0) then
  flag = 0
  theta = dacos(R/dsqrt(Q**3))
  do k = 1,3
      lambda(k) = (-2.d0*dsqrt(Q))*dcos( (theta+2.d0*pi*(k-1.d0))/3.d0 ) &
                - (a1/3.d0)
  enddo
! -----------------------------------------------------------------------------
! complex roots
else
  flag = 1
  s = -1.d0*sign(1.d0, R)*(dabs(R)+dsqrt(R**2-Q**3))**(1.d0/3.d0)
  omega = cmplx(-0.5d0, dsqrt(3.d0)/2.d0, kind=dp)
  lambda(1) = s + Q/s - a1/3.d0
  z2 = (s*omega) + ((Q/s)*omega**2) - a1/3.d0
  z3 = (s*omega**2) + ((Q/s)*omega) - a1/3.d0
  lambda(2) = real(z2)
  lambda(3) = real(z3)
  if ((verbose.eqv..True.).and.(write_listing)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING: complex eigenvalues in characteristic polynomial  '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
endif
! -----------------------------------------------------------------------------
! tri
if (lambda(1) > lambda(2)) then
  temp = lambda(2)
  lambda(2) = lambda(1)
  lambda(1) = temp
endif
if (lambda(2) > lambda(3)) then
  temp = lambda(3)
  lambda(3) = lambda(2)
  lambda(2) = temp
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!           SEDIMENT FLUX LIMITATION FOR NON ERODIBLE BED LAYER
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine compute_flused_limiter(zf, zb, fluqsed_lim)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Sediment flux limitation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ci, dt, qs_lim
implicit none
real(kind=dp), intent(in)    :: zf          (0:nx+1)
real(kind=dp), intent(in)    :: zb          (0:nx+1)
real(kind=dp), intent(inout) :: fluqsed_lim (1:nx)
real(kind=dp)                :: sedmass
integer                      :: i
! -----------------------------------------------------------------------------
do i=1,nx
  if (qs_lim.eq.1) then
    sedmass = max(abs(zb(i)-zf(i)), 0d0)
    fluqsed_lim(i) = sedmass*(ci(i)/dt)
  else 
    fluqsed_lim(i) = 0d0
  endif
enddo
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                     NUMERICAL BEDLOAD FLUXES COMPUTATION
!                             (DE-COUPLED SCHEMES)
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine flux_bedload_centered(fqsed, qsg, qsd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Exner centered scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: fqsed
real(kind=dp), intent(in)    :: qsg, qsd
! -----------------------------------------------------------------------------
fqsed = 0.5d0*(qsg + qsd)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bedload_upwind(fqsed, qsg, qsd, zg, zd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Exner Upwind scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: fqsed
real(kind=dp), intent(in)    :: qsg, qsd, zg, zd
real(kind=dp)                :: cc
! -----------------------------------------------------------------------------
! dry state correction
cc = (qsd-qsg)/max(zd-zg, eps_num)
! Upwind flux
if (cc >= 1.d0) then
  fqsed = qsg
else
  fqsed = qsd
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bedload_upwind_fr(fqsed, qsg, qsd, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Exner Upwind scheme (Froude)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: fqsed
real(kind=dp), intent(in)    :: qsg, qsd
real(kind=dp), intent(in)    :: ug, ud, hg, hd
real(kind=dp)                :: hge, hde, hc, uc, cc, Fr
! -----------------------------------------------------------------------------
! dry state correction
hge = max(eps, hg)
hde = max(eps, hd)
! Roe states and Froude number at interfaces: 
!hc = sqrt(hge*hde)
hc = (hge + hde)/2.d0
uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
cc = dsqrt(g*hc)
Fr = uc/cc
! Upwind flux
if (Fr <= 1.d0) then
  !fqsed = qsd  ! PU notes 2017
  fqsed = qsg  ! Courlis
else
  !fqsed = qsg  ! PU notes 2017
  fqsed = qsd  ! Courlis
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bedload_acu(fqsedg, fqsedd, i, qsg, qsd, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Exner ACU scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: fqsedg, fqsedd
real(kind=dp), intent(in)    :: ug, ud, hg, hd, qsg, qsd
integer,       intent(in)    :: i
real(kind=dp)                :: c1, c2, cc
real(kind=dp)                :: lambg(3), lambd(3)
! -----------------------------------------------------------------------------
! SVE wave speeds : retaining l1 and l2 for Exner
call sve_waves(i,   hg, ug, lambg, 0)
call sve_waves(i+1, hd, ud, lambd, 0)
c1 = min(lambg(1), lambd(1))
c2 = max(lambg(2), lambd(2))
! -----------------------------------------------------------------------------
! Centered flux version
! ~~~~~~~~~~~~~~~~~~~~~
!if ((c1.eq.0.d0).and.(c2.eq.0.d0)) then
!  cc = 0.d0
!else
!  cc = (c2**2.d0 - c1**2.d0)/(c2**2.d0 + c1**2.d0)
!endif
!fqsedg = 0.5d0*(qsg + qsd) - 0.5d0*cc*(qsd - qsg)
!fqsedg = fqsedd
! -----------------------------------------------------------------------------
! ACU version
! ~~~~~~~~~~~
cc = c2**2.d0 + c1**2.d0
fqsedg = qsg + (c1**2/cc)*(qsd - qsg)
fqsedd = qsd - (c2**2/cc)*(qsd - qsg)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bedload_hll(fqs, i, qsg, qsd, ug, ud, hg, hd, zg, zd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Exner HLL flux with lambda1 and lambda2 of the SVE system
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: fqs
real(kind=dp), intent(in)    :: qsg, qsd
real(kind=dp), intent(in)    :: ug, ud, hg, hd, zg, zd
integer,       intent(in)    :: i
real(kind=dp)                :: c1, c2
real(kind=dp)                :: lambg(3), lambd(3)
! -----------------------------------------------------------------------------
! SVE wave speeds : retaining l1 and l2 for Exner
call sve_waves(i,   hg, ug, lambg, 0)
call sve_waves(i+1, hd, ud, lambd, 0)
c1 = min(lambg(1), lambd(1))
c2 = max(lambg(2), lambd(2))
! -----------------------------------------------------------------------------
if (c1.ge.0d0) then
  fqs = qsg
elseif ((c1.lt.0d0).and.(0d0.lt.c2)) then
  fqs = (c2*qsg -c1*qsd)/(c2-c1) + c1*c2*(zd - zg)/(c2-c1)
else
  fqs = qsd
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bedload_ung_2017(fqsed1, fqsed2, qsg, qsd, ug , ud , hg , hd, &
                                 ug2, ud2, hg2, hd2)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! P. Ung (2017) scheme: 
! * upwind with Fr when Fluvial/Fluvial or torrential/torrential 
! * centered when Fluvial/Torrential
!
!       Ii-1  Ii  Ii+1
!       fr1  fr2  fr3
!   +----+----+----+----+  stencil
!     hg2  hg   hd  hd2 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: fqsed1, fqsed2
real(kind=dp), intent(in)    :: qsg, qsd, ug, ud, hg, hd, ug2, ud2, hg2, hd2
real(kind=dp)                :: hge, hde, hge2, hde2, hc, uc, fr1, fr2, fr3
! -----------------------------------------------------------------------------
hge = max(eps, hg)
hde = max(eps, hd)
hge2 = max(eps, hg2)
hde2 = max(eps, hd2)
! Froude at interface i-1/2
hc = (hge2 + hge)/2.d0
uc = (dsqrt(hg)*ug+dsqrt(hg2)*ug2)/(dsqrt(hge2)+dsqrt(hge))
fr1 = uc/dsqrt(g*hc)
! Froude at interface i+1/2
hc = (hge + hde)/2.d0
uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
fr2 = uc/dsqrt(g*hc)
! Froude at interface i+3/2
hc = (hde + hde2)/2.d0
uc = (dsqrt(hd2)*ud2+dsqrt(hd)*ud)/(dsqrt(hde)+dsqrt(hde2))
fr3 = uc/dsqrt(g*hc)
! flux i+1/2 left
if ((fr1.lt.1.d0).and.(fr2.lt.1.d0)) then
  !fqsed1 = qsd ! PU notes 2017
  fqsed1 = qsg ! Courlis
elseif ((fr1.ge.1.d0).and.(fr2.ge.1.d0)) then
  !fqsed1 = qsg ! PU notes 2017
  fqsed1 = qsd ! Courlis
elseif ((fr1.lt.1.d0).and.(fr2.ge.1.d0)) then
  fqsed1 = 0.5d0*qsd
elseif ((fr1.ge.1.d0).and.(fr2.lt.1.d0)) then
  fqsed1 = 0.d0
endif
! flux i+1/2 right
if ((fr2.lt.1.d0).and.(fr3.lt.1.d0)) then
  !fqsed2 = qsd ! PU notes 2017
  fqsed2 = qsg ! Courlis
elseif ((fr2.ge.1.d0).and.(fr3.ge.1.d0)) then
  !fqsed2 = qsg ! PU notes 2017
  fqsed2 = qsd ! Courlis
elseif ((fr2.lt.1.d0).and.(fr3.ge.1.d0)) then
  fqsed2 = 0.5d0*qsg
elseif ((fr2.ge.1.d0).and.(fr3.lt.1.d0)) then
  fqsed2 = 0.d0
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************