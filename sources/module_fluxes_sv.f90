! *****************************************************************************
module module_fluxes_sv
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for Shallow Water equations fluxes computation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g, eps, eps_num
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                             SHALLOW WATER FLUXES
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
elemental function Fluxh(h,u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none  
real(kind=dp),intent(in) :: h,u
real(kind=dp)            :: Fluxh
! -----------------------------------------------------------------------------
Fluxh = h*u
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
elemental function Fluxq(h,u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none  
real(kind=dp),intent(in) :: h,u
real(kind=dp)            :: Fluxq
! -----------------------------------------------------------------------------
Fluxq = h*u**2 + 0.5d0*g*h**2
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                        NUMERICAL FLUXES COMPUTATION
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine init_fluxes()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initialisation of fluxes
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ntrac, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r
implicit none
! -----------------------------------------------------------------------------
! hydro fluxes
fluh_l = 0d0
fluh_r = 0d0
fluq_l = 0d0
fluq_r = 0d0
fluh_src_l = 0d0
fluh_src_r = 0d0
fluq_src_l = 0d0
fluq_src_r = 0d0
! tracers
if (ntrac.gt.0) then
  flut_l = 0d0
  flut_r = 0d0
  flut_src_l = 0d0
  flut_src_r = 0d0
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_flux_sw(zb, h, u, t)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: hydrostarec, ntrac, eps, mlamb, &
                       ordre, ordre_trac, rec_vel, variable_width, dx, &
                       diffusion_u, difu, diffusion_t, nu_t, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r, &
                       corr_bp, corr_bm, corr_hp, corr_hm, &
                       corr_up, corr_um, corr_tp, corr_tm
implicit none
real(kind=dp), intent(in) :: zb (0:nx+1)
real(kind=dp), intent(in) :: h  (0:nx+1)
real(kind=dp), intent(in) :: u  (0:nx+1)
real(kind=dp), intent(in) :: t  (0:nx+1,ntrac)
real(kind=dp)             :: fhg, fhd, fqg, fqd, ft, ft_diff, fu_diff
real(kind=dp)             :: fhg_src, fhd_src, fqg_src, fqd_src
real(kind=dp)             :: hg, hd, zbg, zbd, ug, ud, qg, qd, tg, td
real(kind=dp)             :: dzp, dzm, hip, him, etag, etad, zip12, mc
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
!$OMP         fqd_src, tg, td, ft, ft_diff, mc) &
!$OMP shared(fluh_l, fluh_r, fluh_src_l, fluh_src_r, &
!$OMP        fluq_l, fluq_r, fluq_src_l, fluq_src_r, &
!$OMP        flut_l, flut_r, flut_src_l, flut_src_r, &
!$OMP        corr_bp, corr_bm, corr_hp, corr_hm, &
!$OMP        corr_up, corr_um, corr_tp, corr_tm, &
!$OMP        diffusion_u, difu, diffusion_t, nu_t, &
!$OMP        ntrac, rec_vel, variable_width, ordre, ordre_trac, &
!$OMP        hydrostarec, zb, h, u, t, dx) &
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

  ! second order reconstruction
  ! ---------------------------
  if (ordre.ge.2) then
    if ( zbg.lt.(hd+zbd) .and. zbd.lt.(hg+zbg) &
        .and.2.d0*abs(corr_bp(i)).lt.hd &  
        .and.2.d0*abs(corr_bp(i)).lt.hg &
        .and.2.d0*abs(corr_bm(i+1)).lt.hd &
        .and.2.d0*abs(corr_bm(i+1)).lt.hg ) then
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

  ! shallow water flux
  ! ------------------
  call flux_sv(mc, fhg, fhd, fqg, fqd, him, hip, ug, ud, zbg, zbd)
  mlamb = max(mlamb, mc)
  fluh_l(i) = fhg
  fluh_r(i) = fhd
  fluq_l(i) = fqg
  fluq_r(i) = fqd

  ! bathy source (scheme for Sb=-gh.grad(zb))
  ! ------------
  call flux_src_zb(fhg_src, fhd_src, fqg_src, fqd_src, hg, hd, him, hip, &
                   zbg, zbd, ug, ud)
  fluh_src_l(i) = fhg_src
  fluh_src_r(i) = fhd_src
  fluq_src_l(i) = fqg_src
  fluq_src_r(i) = fqd_src

  ! second order correction
  if (hydrostarec.ne.0 .and. ordre.ge.2 .and. rec_flag.eqv..true.) then
    fluq_src_l(i) = fluq_src_l(i) + 0.5d0*g*(hg+h(i)  )*corr_bp(i)
    fluq_src_r(i) = fluq_src_r(i) - 0.5d0*g*(hd+h(i+1))*corr_bm(i+1)
  endif

  ! width source (scheme for SL=-(hu, hu**2).grad(L)/L)
  ! ------------
  if (variable_width.eqv..True.) then
    call flux_src_width(i, fhg_src, fhd_src, fqg_src, fqd_src, &
                        hg, hd, ug, ud)
    fluh_src_l(i) = fluh_src_l(i) + fhg_src
    fluh_src_r(i) = fluh_src_r(i) + fhd_src
    fluq_src_l(i) = fluq_src_l(i) + fqg_src
    fluq_src_r(i) = fluq_src_r(i) + fqd_src
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
! *****************************************************************************
!                      NUMERICAL FLUX AT INTERFACE  [Wg | Wd]
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine flux_sv(mc, f1g, f1d, f2g, f2d, hg, hd, ug, ud, zg, zd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Flux selector : return the flux depending on the selected method
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: method
implicit none
real(kind=dp), intent(inout) :: mc
real(kind=dp), intent(out)   :: f1g, f1d, f2g, f2d
real(kind=dp), intent(in)    :: hg, hd, ug, ud, zg, zd
! -----------------------------------------------------------------------------
! initialization    
f1g = 0d0
f1d = 0d0
f2g = 0d0
f2d = 0d0
! 
! Roe
! ***
if (method.eq.1) then
  call flux_roe(mc, f1g, f1d, f2g, f2d, ug, ud, hg, hd)
! HLL (Harten, Lax, Van-Leer)
! ***************************
elseif (method.eq.2) then
  call flux_hll(mc, f1g, f2g, ug, ud, hg, hd)
  f1d = f1g
  f2d = f2g
! ACU (Audusse, Chalon, Ung)
! **************************
elseif (method.eq.3) then
  call flux_acu_sw(mc, f1g, f2g, f2d, ug, ud, hg, hd, zg, zd)
  f1d = f1g
! Kinetic
! *******
elseif (method.eq.4) then
  call flux_cin(mc, f1g, f2g, ug, ud, hg, hd)
  f1d = f1g
  f2d = f2g
! TEST
! ****
elseif (method.eq.5) then
  ! Centered
  !call flux_centered(f1g, f2g, ug, ud, hg, hd)
  ! Lax-Friedrich
  !call flux_laxfriedrich(f1g, f2g, ug, ud, hg, hd)
  ! Gallouet
  !call flux_gallouet(mc, f1g, f2g, ug, ud, hg, hd)
  ! Rusanov
  call flux_rusanov(mc, f1g, f2g, ug, ud, hg, hd)
  f1d = f1g
  f2d = f2g
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : Unknown ADVECTION SCHEME '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************
      


! *****************************************************************************
subroutine flux_centered(f1, f2, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Centered fluxes scheme / WARNING : THIS SCHEME IS UNSTABLE
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: f1, f2
real(kind=dp), intent(in)    :: ug, ud, hg, hd
! -----------------------------------------------------------------------------
f1 = 0.5d0*( Fluxh(hg,ug) + Fluxh(hd,ud) )
f2 = 0.5d0*( Fluxq(hg,ug) + Fluxq(hd,ud) )
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_laxfriedrich(f1, f2, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Lax-Friedrich scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dt, minci
implicit none
real(kind=dp), intent(inout) :: f1, f2
real(kind=dp), intent(in)    :: ug, ud, hg, hd
real(kind=dp)                :: c
! -----------------------------------------------------------------------------
c = minci/dt ! c=celerity <=> upwind scheme
f1 = 0.5d0*( Fluxh(hg,ug) + Fluxh(hd,ud) -c*(hd - hg) )
f2 = 0.5d0*( Fluxq(hg,ug) + Fluxq(hd,ud) -c*(hd*ud - hg*ug) )
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_rusanov(c, f1, f2, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Rusanov scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out)   :: c
real(kind=dp), intent(inout) :: f1, f2
real(kind=dp), intent(in)    :: ug, ud, hg, hd
! -----------------------------------------------------------------------------
c = max( abs(ug)+sqrt(g*hg), abs(ud)+sqrt(g*hd) )
f1 = 0.5d0*( Fluxh(hg,ug) + Fluxh(hd,ud) -c*(hd - hg) )
f2 = 0.5d0*( Fluxq(hg,ug) + Fluxq(hd,ud) -c*(hd*ud - hg*ug) )
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_roe(c, f1g, f1d, f2g, f2d, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Roe scheme (J.Comp.Physics, vol.43,1981)(lin. Riemann solver,flux-diff-split)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out)   :: c
real(kind=dp), intent(inout) :: f1g, f1d, f2g, f2d
real(kind=dp), intent(in)    :: ug, ud, hg, hd
real(kind=dp)                :: a(2,2), jump_w(2)
integer                      :: i
integer, parameter           :: roe_method = 3
! -----------------------------------------------------------------------------
! variable difference vector
jump_w(1) = hd - hg
jump_w(2) = hd*ud - hg*ug
! **************************************************************
! Roe flux, method with |A| : f = 0.5*(fl+fr) - 0.5*|a|*(wr-wl))
! **************************************************************
if (roe_method.eq.0) then
  ! compute roe matrix |A|=R*|lambda|*R1
  call aroe_pm(c, a, hd, hg, ud, ug, 0)
  ! evaluate roe's flux
  f1g = Fluxh(hg,ug) + Fluxh(hd,ud)
  f2g = Fluxq(hg,ug) + Fluxq(hd,ud)
  do i=1,2
    f1g = f1g - a(1,i)*jump_w(i)
    f2g = f2g - a(2,i)*jump_w(i)
  enddo
  f1g = f1g/2.d0
  f2g = f2g/2.d0
  f1d = f1g
  f2d = f2g
! *************************************************
! Roe flux, method with A- : fg = (fl + a-*(wr-wl))
! *************************************************
elseif (roe_method.eq.1) then
  ! compute roe matrix A-=R*(Lambda-)*R1
  call aroe_pm(c, a, hd, hg, ud, ug, -1)
  ! evaluate roe's flux
  f1g = Fluxh(hg,ug)
  f2g = Fluxq(hg,ug)
  do i=1,2
    f1g = f1g + a(1,i)*jump_w(i)
    f2g = f2g + a(2,i)*jump_w(i)
  enddo
  f1d = f1g
  f2d = f2g
! *************************************************
! Roe flux, method with A+ : fd = (fr - a+*(wr-wl))
! *************************************************
elseif (roe_method.eq.2) then
  ! compute roe matrix A+=R*(Lambda+)*R1
  call aroe_pm(c, a, hd, hg, ud, ug, +1)
  ! evaluate roe's flux
  f1d = Fluxh(hd,ud)
  f2d = Fluxq(hd,ud)
  do i=1,2
    f1d = f1d - a(1,i)*jump_w(i)
    f2d = f2d - a(2,i)*jump_w(i)
  enddo
  f1g = f1d
  f2g = f2d
! ***********************
! Roe flux: Mascaret like
! ***********************
elseif (roe_method.eq.3) then
  ! check left cell dry state         
  if (hg.lt.eps) then
    ! compute roe matrix A+=R*(Lambda+)*R1
    call aroe_pm(c, a, hd, hg, ud, ug, +1)
    ! evaluate roe's flux
    f1d = Fluxh(hd,ud)
    f2d = Fluxq(hd,ud)
    do i=1,2
      f1d = f1d - a(1,i)*jump_w(i)
      f2d = f2d - a(2,i)*jump_w(i)
    enddo
    f1g = f1d
    f2g = f2d
  else
    ! compute roe matrix A-=R*(Lambda-)*R1
    call aroe_pm(c, a, hd, hg, ud, ug, -1)
    ! evaluate roe's flux
    f1g = Fluxh(hg,ug)
    f2g = Fluxq(hg,ug)
    do i=1,2
      f1g = f1g + a(1,i)*jump_w(i)
      f2g = f2g + a(2,i)*jump_w(i)
    enddo
    f1d = f1g
    f2d = f2g
  endif
! -----------------------------------------------------------------------------
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : Unknown Roe flux method, check flux_roe '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine aroe(hd, hg, ud, ug, a)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute roe's matrix A
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out) :: a(2,2)
real(kind=dp), intent(in)  :: hd, hg, ud, ug
real(kind=dp)              :: hc, uc, cc, hge, hde
! -----------------------------------------------------------------------------
! dry state correction
hge = max(eps_num, hg)
hde = max(eps_num, hd)
! roe's average for shallow water
hc = (hge + hde)/2.d0
uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
cc = dsqrt(g*hc)
! Roe's matrix
a(1,1) = 0d0
a(1,2) = 1d0
a(2,1) = (cc**2)-(uc**2)
a(2,2) = 2d0*uc
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine aroe_pm(mc, a, hd, hg, ud, ug, matrix_sign)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute roe's matrix
! * matrix_sign =  0: returns |A|
! * matrix_sign = -1: returns  A-
! * matrix_sign = +1: returns  A+
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_linalg, only: invert_2x2_matrix
implicit none
integer,       intent(in)  :: matrix_sign
real(kind=dp), intent(out) :: mc
real(kind=dp), intent(out) :: a(2,2)
real(kind=dp), intent(in)  :: hd, hg, ud, ug
real(kind=dp)              :: r(2,2), r1(2,2), lambda(2), Lamb(2)
real(kind=dp)              :: hc, uc, cc, hge, hde
integer                    :: i, j, k
logical, parameter         :: entropy_fix=.True.
! -----------------------------------------------------------------------------
! dry state correction
hge = max(eps_num, hg)
hde = max(eps_num, hd)

! roe's average for shallow water
! *******************************
hc = (hge + hde)/2.d0
uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
cc = dsqrt(g*hc)

! eigenvalues
lambda(1) = uc-cc
lambda(2) = uc+cc

! entropy fix (Harten-Hyman, see Leveque 2002)
if (entropy_fix.eqv..True.) then
  call roe_entropy_fix(lambda, uc, ug, ud, hg, hd)
endif

! diagonal of Lambda (|L|, L+ or L-)
! **********************************
do k=1,2
  if (matrix_sign.eq.-1) then
    Lamb(k) = min(lambda(k), 0d0)
  else if (matrix_sign.eq.+1) then
    Lamb(k) = max(lambda(k), 0d0)
  else
    Lamb(k) = dabs(lambda(k))
  endif
enddo

! diagonalization matrix r
! ************************
r(1,1) = 1.d0
r(1,2) = 1.d0
r(2,1) = uc-cc
r(2,2) = uc+cc

! inverse r1
!call invert_2x2_matrix(r, r1, flag)
! **********************************
r1(1,1) = (uc+cc)/(2.d0*cc)
r1(1,2) = -1.d0/(2.d0*cc)
r1(2,1) =-(uc-cc)/(2.d0*cc)
r1(2,2) = 1.d0/(2.d0*cc)

! matrix*diag*matrix multiply
! ***************************
do i=1,2
  do j=1,2
    a(i,j) = 0.d0
    do k=1,2
      ! matrix |A|=R|L|R-1, A+=R(L+)R-1 or A-=R(L-)R-1
      a(i,j) = a(i,j)+r(i,k)*Lamb(k)*r1(k,j)
    enddo
  enddo
enddo
! estimate the maximum signal velocity
mc = max(dabs(lambda(1)), dabs(lambda(2)))
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine roe_entropy_fix(lambda, uc, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Harten-Hyman entropy fix for Roe scheme (See Leveque 2002)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(inout) :: lambda(2)
real(kind=dp), intent(in)    :: ug, ud, hg, hd, uc
real(kind=dp)                :: cg, cd, frg, frd, l1d, l2d, l1g, l2g
! -----------------------------------------------------------------------------
cg = dsqrt(g*hg)
cd = dsqrt(g*hd)
frd = ud/max(eps_num, cd)
frg = ug/max(eps_num, cg)
if ((frd>1.d0 ).and.(frg<1.d0).and.(uc>0.d0).and.(hg>0.d0)) then
  l1d = ud-cd
  l1g = ug-cg
  lambda(1) = l1g*(lambda(1) - l1d)/(l1g - l1d)
endif
if ((frd<1.d0 ).and.(frg>1.d0).and.(uc<0.d0).and.(hg>0.d0)) then
  l2d = ud+cd
  l2g = ug+cg
  lambda(2) = l2d*(lambda(2) - l2g)/(l2d - l2g)
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_gallouet(mc, f1, f2, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Gallouet scheme 
! (Buffard, Gallouet, Herard, CRAcadSci, t.326,Serie I,p.385-390, 1998)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out)   :: mc
real(kind=dp), intent(inout) :: f1, f2
real(kind=dp), intent(in)    :: ug, ud, hg, hd
real(kind=dp)                :: hm, um, cm, cg, cd, hc, uc, cc
! -----------------------------------------------------------------------------
cg=dsqrt(g*hg)
cd=dsqrt(g*hd)
hm = (hg+hd)/2.d0
um = (ug+ud)/2.d0
cm = (cg+cd)/2.d0
! dry state
if ((ug-ud).ge.2.d0*(cg+cd)) then
  uc = 0.d0
  cc = 0.d0
! intermediate state
else
  uc = um-(cd-cg)
  cc = (cd+cg)/4.d0*(2.d0-(ud-ug)/max(eps, cd+cg))
  mc = max(dabs(um+cm), dabs(um-cm))
endif
hc=cc*cc/g
! fluxes
f1 = Fluxh(hc,uc)
f2 = Fluxq(hc,uc)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_hll(mc, f1, f2, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! HLL scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out)   :: mc
real(kind=dp), intent(inout) :: f1, f2
real(kind=dp), intent(in)    :: ug, ud, hg, hd
real(kind=dp)                :: c1, c2
! -----------------------------------------------------------------------------
c1 = min( ug-sqrt(g*hg), ud-sqrt(g*hd) )
c2 = max( ug+sqrt(g*hg), ud+sqrt(g*hd) )
mc = max(c1, c2)
if (c1.ge.0d0) then
  f1 = Fluxh(hg,ug)
  f2 = Fluxq(hg,ug)
elseif ((c1.lt.0d0).and.(0d0.lt.c2)) then
  f1 = (c2*Fluxh(hg,ug)-c1*Fluxh(hd,ud))/(c2-c1) + c1*c2*(hd-hg)/(c2-c1)
  f2 = (c2*Fluxq(hg,ug)-c1*Fluxq(hd,ud))/(c2-c1) + c1*c2*(hd*ud-hg*ug)/(c2-c1)
else
  f1 = Fluxh(hd,ud)
  f2 = Fluxq(hd,ud)
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_acu_sw(mc, f1, f2g, f2d, ug, ud, hg, hd, zg, zd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ACU scheme (SW).
! A SIMPLE WELL-BALANCED AND POSITIVE NUMERICAL SCHEME FOR THE SHALLOW-WATER 
! SYSTEM. E. AUDUSSE, C. CHALONS AND P. UNG
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out)   :: mc
real(kind=dp), intent(inout) :: f1, f2g, f2d
real(kind=dp), intent(in)    :: ug, ud, hg, hd, zg, zd
real(kind=dp)                :: c1, c2, mc0, tmpg, tmpd, dz, src
! -----------------------------------------------------------------------------
! HLL flux
call flux_hll(mc0, f1, f2g, ug, ud, hg, hd)
f2d = f2g
! -----------------------------------------------------------------------------
c1 = min( ug-sqrt(g*hg), ud-sqrt(g*hd) )  !Eq. 2.20
c2 = max( ug+sqrt(g*hg), ud+sqrt(g*hd) )  !Eq. 2.21
mc = max(mc0, max(c1, c2))
if (c2-c1.ne.0) then
  tmpg = c1/(c2-c1)
  tmpd = c2/(c2-c1)
else
  tmpg = 0
  tmpd = 0
endif
dz = zd-zg
if (dz.ge.0) then
  src = -0.5d0*(hg+hd)*min(hg, dz)  !Eq. 2.27a
else
  src = -0.5d0*(hg+hd)*max(-hd, dz) !Eq. 2.27b
endif
! -----------------------------------------------------------------------------
f1 = f1 + tmpg*c2*dz    !Eq. 2.17a,b
f2g = f2g + tmpg*g*src  !Eq. 2.17a
f2d = f2d + tmpd*g*src  !Eq. 2.17b
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_cin(mc, f1, f2, ug, ud, hg, hd)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Kinetic scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(out)   :: mc
real(kind=dp), intent(inout) :: f1, f2
real(kind=dp), intent(in)    :: ug, ud, hg, hd
real(kind=dp)                :: cg, cd, sqrt3, aux1, aux2, extp, extd
real(kind=dp)                :: f1p, f1m, f2p, f2m
integer, parameter           :: kin_type = 1
! -----------------------------------------------------------------------------
cd = dsqrt(0.5d0*g*hd)
cg = dsqrt(0.5d0*g*hg)
! ************************************************************
! Kinetic scheme A: corresponds to the following chi function:
!   chi = (1/2sqrt(3)) * Ind_{abs{sqrt(3)}}
! ************************************************************
if (kin_type.eq.1) then
  sqrt3 = dsqrt(3.0d0)
  if (hg.gt.eps) then
    extp = max(0.d0, ug + cg*sqrt3)
    extd = max(0.d0, ug - cg*sqrt3)
    mc = max(extp, extd)
    aux1 = hg/(cg*2.d0*sqrt3)            !2D: = sqrt3*hg/cg
    f1p = aux1*(extp**2 - extd**2)/2.d0  !2D: = .../12.d0
    f2p = aux1*(extp**3 - extd**3)/3.d0  !2D: = .../18.d0
  else
    f1p = 0d0
    f2p = 0d0
  endif
  if (hd.gt.eps) then
    extp = min(0.d0, ud + cd*sqrt3)
    extd = min(0.d0, ud - cd*sqrt3)
    mc = max(extp, extd)
    aux2 = hd/(cd*2.d0*sqrt3)
    f1m = aux2*(extp**2 - extd**2)/2.d0
    f2m = aux2*(extp**3 - extd**3)/3.d0
  else
    f1m = 0d0
    f2m = 0d0
  endif
! ************************************************************
! Kinetic scheme B: corresponds to the following chi function:
!   chi =  a*sqrt(1-x^2/b^2)) with b=2 and a=1/pi
! ************************************************************
elseif (kin_type.eq.2) then
  f1p = fchi4( 1, 1, hg, ug, cg)
  f2p = fchi4( 1, 2, hg, ug, cg)
  f1m = fchi4(-1, 1, hd, ud, cd)
  f2m = fchi4(-1, 2, hd, ud, cd)
endif
! -----------------------------------------------------------------------------
f1 = f1p + f1m
f2 = f2p + f2m
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
function fchi4(idir, iorder, h, u, c)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Flux computeus for a*sqrt(1-x^2/b^2)) with b=2 and a=1/pi
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: pi, eps_num
implicit none
integer, intent(in)       :: idir
integer, intent(in)       :: iorder
real(kind=dp), intent(in) :: h, u, c
real(kind=dp)             :: fchi4, a, b, xi0, zi, zs, zs2, zi2, b2, c2, u2
real(kind=dp)             :: aux0, aux1, aux2, aux3, aux4, aux5, aux6
real(kind=dp)             :: asin_zs, asin_zi, sqrt_zs, sqrt_zi
! -----------------------------------------------------------------------------
b = 2.0d0
a = 1.0d0/pi
fchi4 = 0.0d0
! -----------------------------------------------------------------------------
if (c.ge.eps_num) then
  xi0 = -u/c
else
  xi0 = 0.0d0
endif
! -----------------------------------------------------------------------------
if (idir.eq.1.and.xi0.le.b) then
  ! idir=1 corresponds to the upwind xi>=0
  zs = max(xi0, b)
  zi = max(xi0, -1.0d0*b)
elseif (idir.eq.-1.and.xi0.ge.-b) then
  zs = min(xi0, b)
  zi = min(xi0, -1.0d0*b)
else
  zs = 0.d0
  zi = 0.d0
endif
! -----------------------------------------------------------------------------
zs2 = zs**2
zi2 = zi**2
b2 = b**2
asin_zs = dasin(zs/b)
asin_zi = dasin(zi/b)
sqrt_zs = dsqrt(b2-zs2)
sqrt_zi = dsqrt(b2-zi2)
aux0 = a*h/b
! -----------------------------------------------------------------------------
if (iorder.eq.0) then
  fchi4 = aux0*(0.5d0*zs*sqrt_zs + 0.5d0*b2*asin_zs) &
        - aux0*(0.5d0*zi*sqrt_zi + 0.5d0*b2*asin_zi)
elseif (iorder.eq.1) then
  fchi4 = aux0*(-c/3.0d0*(b2-zs2)**(1.5d0) + u/2.0d0*(zs*sqrt_zs+b2*asin_zs)) &
        - aux0*(-c/3.0d0*(b2-zi2)**(1.5d0) + u/2.0d0*(zi*sqrt_zi+b2*asin_zi))
elseif (iorder.eq.2) then
  c2 = c**2
  u2 = u**2
  aux1 = 0.125d0*b2
  aux2 = 2.0d0/3.0d0*c*u
  aux3 = (b2-zs2)**1.5d0
  aux4 = (b2-zi2)**1.5d0
  aux5 = zs*sqrt_zs+b2*asin_zs
  aux6 = zi*sqrt_zi+b2*asin_zi
  fchi4 = aux0*(c2*(-0.25d0*zs*aux3 + aux1*aux5) - aux2*aux3 + 0.5d0*u2*aux5) &
        - aux0*(c2*(-0.25d0*zi*aux4 + aux1*aux6) - aux2*aux4 + 0.5d0*u2*aux6)
else
    print*,'Wrong value for iorder in the definition of fchi4'
endif
! -----------------------------------------------------------------------------
end function fchi4
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                   SOURCE FLUX AT INTERFACE  [Wg | Wd]
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine flux_src_zb(f1g, f1d, f2g, f2d, hg, hd, him, hip, zbg, zbd, ug, ud)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Scheme for the bathy source term Sb=-gh grad(zb)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,   only: g, hydrostarec
implicit none
real(kind=dp), intent(out) :: f1g, f1d, f2g, f2d
real(kind=dp), intent(in)  :: hg, hd, him, hip, zbg, zbd, ug, ud
real(kind=dp)              :: etag, etad, zip12, hc, src_centered
real(kind=dp)              :: Mp2(2,2), Mm2(2,2)
! -----------------------------------------------------------------------------
! initialization
f1g = 0d0
f1d = 0d0
f2g = 0d0
f2d = 0d0

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
  call src_upwind_roe_sv(Mp2, Mm2, hg, hd, ug, ud)
  ! update source
  hc = 0.5d0*(hg+hd)
  src_centered = 0.5d0*g*hc*(zbd-zbg)
  f1g = src_centered*Mm2(1,2)
  f1d = src_centered*Mp2(1,2)
  f2g = src_centered*Mm2(2,2)
  f2d = src_centered*Mp2(2,2)

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
subroutine flux_src_width(i, f1g, f1d, f2g, f2d, hg, hd, ug, ud)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Scheme for the variable width source term (for Saint-Venant)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: width_scheme, dx, Lx, dLx
implicit none
integer, intent(in)        :: i
real(kind=dp), intent(out) :: f1g, f1d, f2g, f2d
real(kind=dp), intent(in)  :: hg, hd, ug, ud
real(kind=dp)              :: src_Lg, src_Ld, dLip12m, dLip12p
real(kind=dp)              :: Mp2(2,2), Mm2(2,2), Sg(2), Sd(2)
real(kind=dp)              :: hge, hde, hc, uc
integer, parameter         :: var_opt = 0
! -----------------------------------------------------------------------------
! initialization    
f1g = 0d0
f1d = 0d0
f2g = 0d0
f2d = 0d0
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
&          + 8d0*Lx(min(i+1,nx+1))-Lx(min(i+2,nx+1)))/(12d0*dx(i))
  dLip12p = (Lx(max(0,i-2))-8d0*Lx(max(0,i-1)) &
&          + 8d0*Lx(min(i+1,nx+1))-Lx(min(i+2,nx+1)))/(12d0*dx(i))

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
  ! computes upwind matrices M+ and M-
  call src_upwind_roe_sv(Mp2, Mm2, hg, hd, ug, ud)
  Sg(1) = f1g
  Sg(2) = f2g
  Sd(1) = f1d
  Sd(2) = f2d
  ! upwinding the src
  f1g = Sg(1)*Mm2(1,1) + Sg(2)*Mm2(1,2)
  f1d = Sd(1)*Mp2(1,1) + Sd(2)*Mp2(1,2)
  f2g = Sg(1)*Mm2(2,1) + Sg(2)*Mm2(2,2)
  f2d = Sd(1)*Mp2(2,1) + Sd(2)*Mp2(2,2)
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
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine src_upwind_roe_sv(Mp, Mm, hg, hd, ug, ud)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the upwind matrices for Roe scheme, returns
! M+- = I+-|A|A^-1 where A is the Roe matrix.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,   only: g
use module_linalg, only: matrix_product, invert_2x2_matrix
implicit none
real(kind=dp), intent(out) :: Mm(2,2), Mp(2,2)
real(kind=dp), intent(in)  :: hg, hd, ug, ud
real(kind=dp)              :: A(2,2), A1(2,2), absA(2,2)
real(kind=dp)              :: hge, hde, hc, uc, cc, mlamb
integer                    :: i, j
logical                    :: flag
! -----------------------------------------------------------------------------
! Roe's average
hge = max(eps_num, hg)
hde = max(eps_num, hd)
hc = (hge+hde)/2.d0
uc = (dsqrt(hd)*ud+dsqrt(hg)*ug)/(dsqrt(hge)+dsqrt(hde))
cc = dsqrt(g*hc)

! centered scheme if null eigenvalue (M+-=I)
! ******************************************
if ((uc-cc.eq.0).or.(uc+cc.eq.0)) then
  do i=1,2
    do j=1,2
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
! ******************************************
else
  A(1,1) = 0d0
  A(1,2) = 1d0
  A(2,1) = (cc**2)-(uc**2)
  A(2,2) = 2d0*uc
  call invert_2x2_matrix(A, A1, flag)
  call aroe_pm(mlamb, absA, hd, hg, ud, ug, 0)
  call matrix_product(2, absA, A1, A) ! A <- |A|A^-1
  ! compute M+-
  do i=1,2
    do j=1,2
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



end module
! *****************************************************************************