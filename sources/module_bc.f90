! *****************************************************************************
module module_bc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Boundary Conditions 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                              BOUNDARY CONDITIONS  
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine time_interp_bc(time, time_arr, val_arr, last_idx, val)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! boundary conditions
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: 
implicit none
real(kind=dp), intent(in)  :: time
real(kind=dp), intent(in)  :: time_arr(:)
real(kind=dp), intent(in)  :: val_arr(:)
integer,       intent(inout) :: last_idx
real(kind=dp), intent(out) :: val
real(kind=dp)              :: alpha
integer                    :: i, ni, low, high, mid
logical                    :: found_value
! -----------------------------------------------------------------------------
ni = size(val_arr)
val = 0.d0
found_value = .False.

if (ni.lt.2) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR READING LIQ BOUNDARY FILE, INVALID TIME ARRAY SIZE'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif

! Keep cached index in valid interval range [1, ni-1]
last_idx = max(1, min(last_idx, ni-1))

! Fast path: check cached interval first
if ((time_arr(last_idx) <= time).and.(time <= time_arr(last_idx+1))) then
  i = last_idx
  found_value = .True.
! Very common case for increasing time: next interval
elseif ((last_idx.lt.ni-1).and.(time_arr(last_idx+1) <= time).and. &
        (time <= time_arr(last_idx+2))) then
  i = last_idx + 1
  last_idx = i
  found_value = .True.
else
  ! Fallback: binary search on [last_idx, ni-1]
  low = last_idx
  high = ni-1
  do while (low <= high)
    mid = (low + high)/2
    if ((time_arr(mid) <= time).and.(time <= time_arr(mid+1))) then
      i = mid
      last_idx = mid
      found_value = .True.
      exit
    elseif (time < time_arr(mid)) then
      high = mid - 1
    else
      low = mid + 1
    endif
  enddo
endif

if (found_value) then
  alpha = (time_arr(i+1)-time)/(time_arr(i+1)-time_arr(i))
  val = alpha*val_arr(i) + (1.d0-alpha)*val_arr(i+1)
endif

if (found_value.eqv..False.) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR READING LIQ BOUNDARY FILE, CHECK TIME'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
end if
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine ghostcell_values(i, vars_bc, nvars_bc, time, time_bc, field_bc, &
                            bc_idx, &
                            zs, h, u, q, qs, t)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Prepare ghost cell values for boundary conditions when boundary file is used
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ntrac, clhorz, clqoru, bedload, Lx, trac_names
implicit none
integer,       intent(in)                   :: i
character(len=3), dimension(13), intent(in) :: vars_bc
integer,       intent(in)                   :: nvars_bc
real(kind=dp), intent(in)                   :: time
real(kind=dp), intent(in)                   :: time_bc(:)
real(kind=dp), intent(in)                   :: field_bc(:,:)
integer,       intent(inout)                :: bc_idx
real(kind=dp), intent(inout)                :: zs, h, u, q, qs
real(kind=dp), intent(inout)                :: t(1:ntrac)
real(kind=dp)                               :: val
integer                                     :: k, l
character(len=4)                            :: char_trac
! -----------------------------------------------------------------------------
! loop on fields in bc file
do k=2,nvars_bc ! skip time
  ! time interpolation
  call time_interp_bc(time, time_bc, field_bc(:,k-1), bc_idx, val)
  ! replace ghost cell imposed value with interpolated one
  if (vars_bc(k).eq.'H') then
    h = val
    clhorz = 1
  endif
  if (vars_bc(k).eq.'ZS') then
    zs = val
    clhorz = 0
  endif
  if (vars_bc(k).eq.'Q') then
    q = val/Lx(i)
    clqoru = 1
  endif
  if (vars_bc(k).eq.'U') then
    u = val
    clqoru = 0
  endif
  if ((bedload.eqv..True.).and.(vars_bc(k).eq.'QS')) then
    qs = val/Lx(i)
  endif
  ! tracers (tracer index up to 99 here)
  if (ntrac.gt.0) then
    do l=1,ntrac
      char_trac = trac_names(l)
      if ((vars_bc(k).eq."T"//char_trac(4:4)).or. &
          (vars_bc(k).eq."T"//char_trac(3:4))) then
        t(l) = val
      endif
    enddo
  endif ! end if trac
enddo ! end of field loop
! check values
if (h.lt.0d0) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR IN LIQ BOUNDARY FILE, IMPOSED WATER DEPTH H < 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_cdl
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! boundary conditions
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: time, use_bcfile_l, use_bcfile_r, &
                       bctype_l, bctypet_l, zsa, ha, ua, qa, ta, &
                       bctype_r, bctypet_r, zsb, hb, ub, qb, tb, &
                       time_bcl, time_bcr, field_bcl, field_bcr, &
                       bctypeqs_l, bctypeqs_r, vars_bcl, vars_bcr, &
                       nvars_bcl, nvars_bcr, qsa, qsb, Jeqa, Jeqb, &
                       bc_idx_l, bc_idx_r
use module_io, only: integer_to_string
implicit none
! -----------------------------------------------------------------------------
! LEFT BOUNDARY CONDITION
! -----------------------------------------------------------------------------
! ghost cell data from liq boundary file (time interpolation)
if (use_bcfile_l.eqv..True.) then
  call ghostcell_values(1, vars_bcl, nvars_bcl, time, time_bcl, &
                        field_bcl, bc_idx_l, &
                        zsa, ha, ua, qa, qsa, ta)
end if
! compute fluxes of left boundary
call flux_cdl_ghostcell(1, bctype_l, bctypet_l, bctypeqs_l, zsa, ha, ua, qa, &
                        ta, qsa, Jeqa)
! -----------------------------------------------------------------------------
! RIGHT BOUNDARY CONDITION
! -----------------------------------------------------------------------------
! ghost cell data from liq boundary file (time interpolation)
if (use_bcfile_r.eqv..True.) then
  call ghostcell_values(nx+1, vars_bcr, nvars_bcr, time, time_bcr, &
                        field_bcr, bc_idx_r, &
                        zsb, hb, ub, qb, qsb, tb)
end if
! compute fluxes of right boundary
call flux_cdl_ghostcell(nx, bctype_r, bctypet_r, bctypeqs_r, zsb, hb, ub, qb,&
                        tb, qsb, Jeqb)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_cdl_ghostcell(i, typcl, typclt, typclqs, zsu, hu, uu, qu, tu, &
                              qsu, Jeq)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! boundary conditions
! 
! Ghost cells:
!  u : user imposed values
!  i : internal
!  e : external
! 
!   Wu=ha,qa                             Wu=hb,qb <- user defined
!    o---|---o---     ...      ---o---|---o
!    0   A   1                   nx   B  nx+1
!   We       Wi                  Wi       We
! 
! i==1 : left bnd flux computation
! i==nx: right bnd flux computation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps, eps_num, h, u, t, ntrac, clhorz, clqoru, fric, &
                       flut_l, flut_r, exner, bedload, qsed, bedload_formula, &
                       zb, zf, dx, Lx, Rh, fric, bnd_extrapolation, bnd_flusrc, &
                       bnd_highord, ordre, corr_bp, corr_bm, bnd_warnings, &
                       mass_balance, mb_fh_l, mb_fq_l, mb_fqsed_l, &
                       mb_fh_r, mb_fq_r, mb_fqsed_r, &
                       bedload_coef, porosite, eps_b, eps_b_damping, bed_fric_coef
use module_bedload, only: qs, damping_coef
implicit none
integer, intent(in)       :: i, typcl, typclqs
integer, intent(in)       :: typclt(:)
real(kind=dp), intent(in) :: zsu, hu, uu, qu, qsu, Jeq
real(kind=dp), intent(in) :: tu(:)
real(kind=dp)             :: he, ue, qe, ce, te, ze, qsede, ueq, hec
real(kind=dp)             :: hi, ui, ci, zi, zip, zim, qsedi
real(kind=dp)             :: fh_bord, fq_bord, ft_bord, fqsed_bord
real(kind=dp)             :: b, rb
real(kind=dp)             :: lamb1, lamb2, regime, sgn, gradz, dz, dz_rec
real(kind=dp)             :: hcri
integer                   :: k, j, shift, bnd_highord_loc
character*5               :: bnd_side
logical                   :: force_subcritical = .True.
! -----------------------------------------------------------------------------
bnd_highord_loc = bnd_highord

! sign for left/right bnd
if (i.eq.1) then
  sgn = -1d0
  bnd_side = "LEFT "
else
  sgn = 1d0
  bnd_side = "RIGHT"
endif

! internal state
hi = h(i)
ui = u(i)
zi = zb(i)
zip = zb(i+1)
zim = zb(i-1)

! ------------------
! bottom ghost cell:
! ------------------
! SVE
! ~~~
if (exner.eqv..True.) then
  ! no extrapolation
  if (bnd_extrapolation.eq.0) then
    ze = zi
  ! recompute extrapolation because zb variable
  else
    ! dz of second order rec
    if (ordre.gt.1) then
      if (i.eq.1) then
        dz_rec = corr_bm(i+1)-corr_bp(i)
      else
        dz_rec = corr_bm(i)-corr_bp(i-1)
      endif
    else
      dz_rec = 0d0
    endif
    dz_rec = 0d0 ! TODO : test 
    ! extrapolation of zb        
    if (i.eq.1) then
      gradz = (zip-zi+dz_rec)/dx(i)
      ze = zi-gradz*dx(i-1)
    else
      gradz = (zi-zim+dz_rec)/dx(i-1)
      ze = zi+gradz*dx(i)
    endif
  endif
! SV
! ~~
else
  ! extrapolation computed in init_geometry
  if (i.eq.1) then
    ze = zim
  else
    ze = zip
  endif
endif
! boundary zb step
dz = zi-ze

! --------------------------------
! hydrodynamic boundary conditions
! --------------------------------
! ~~~~~~~~~~~~~
! wall boundary
! ~~~~~~~~~~~~~
if (typcl == 1) then
  bnd_highord_loc = 0
  he = max(0d0, hi + dz) ! +dz to get constant free surface
  ue = -ui
! ~~~~~~~~~
! h imposed
! ~~~~~~~~~
elseif (typcl == 2) then
  ! free surface imposed ~> compute water depth
  if (clhorz.eq.0) then
    he = zsu-ze
  else
    he = hu
  endif
  ! correction if extrapolation of zb
  if (bnd_extrapolation.eq.1) then
    hec = max(0d0, he-dz)
  else
    hec = he
  endif
  ! compute bc :
  ci = sqrt(g*hi)
  ce = sqrt(g*hec)
  lamb1 = ui - ci
  lamb2 = ui + ci
  regime = lamb1*lamb2
  ! fluvial
  if (regime.lt.0d0) then
    ! fluvial: conservation of Riemann invariant
    ue = ui + sgn*2.d0*(ci-ce)
  ! torrential
  else
    if (sgn*ui.le.0d0) then
      ! torrential inflow: imposition of h and u
      if (clqoru.eq.1) then
        qe = qu
        ue = qe/max(eps, he)
      else
        ue = uu
      endif
      if (bnd_warnings.eqv..True.) then
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*) ' WARNING : '//bnd_side//' BC'
        write(*,*) ' h imposed but torrential inflow : ' 
        write(*,*) ' he = ', he
        write(*,*) ' ue = ', ue
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      endif
    else
      ! torrential outflow: no need to impose h nor u
      he = hi
      ue = ui
      if (bnd_warnings.eqv..True.) then
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*) ' WARNING : '//bnd_side//' BC'
        write(*,*) ' h imposed but torrential outflow : ' 
        write(*,*) ' he = ', he
        write(*,*) ' ue = ', ue
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      endif
      ! force a switch back to h imposed if wanted by user (hu > hcri)
      if (force_subcritical.eqv..True.) then
        hcri = (hi*ui/sqrt(g))**(2d0/3d0)
        ! if imposed h > hcritical -> wall until sub-critical
        if (hu.gt.hcri) then
          he = hi
          ue =-ui
        endif
      endif
    endif
  endif
! ~~~~~~~~~
! u imposed
! ~~~~~~~~~
elseif (typcl == 3) then
  ! flowrate imposed ~> compute velocity
  if (clqoru.eq.1) then
    qe = qu
    ue = qe/max(eps, hi)
  ! velocity imposed
  else
    ue = uu
  endif
  ci = sqrt(g*hi)
  lamb1 = ui - ci
  lamb2 = ui + ci
  regime = lamb1*lamb2
  ! check that imposed velocity is not zero
  if (abs(ue).lt.eps_num) then
    bnd_highord_loc = 0
    ! Neuman
    he = hi + dz
    ue = ui
  else
    ! fluvial
    if (regime.lt.0d0) then
      ! fluvial: conservation of Riemann invariant
      he = ((ui - ue + sgn*2d0*sqrt(g*hi))**2)/(4d0*g) 
    ! torrential
    else
      if (sgn*ui.le.0d0) then
        ! torrential inflow: imposition of h and u
        if (clhorz.eq.0) then
          he = max(zsu - ze, eps)
        else
          he = max(hu, eps)
        endif
        if (bnd_warnings.eqv..True.) then
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(*,*) ' WARNING : '//bnd_side//' BC'
          write(*,*) ' u imposed but torrential inflow : ' 
          write(*,*) ' he = ', he
          write(*,*) ' ue = ', ue
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        endif
      else
        ! torrential outflow: no need to impose h nor u
        he = hi
        ue = ui
        if (bnd_warnings.eqv..True.) then
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(*,*) ' WARNING : '//bnd_side//' BC'
          write(*,*) ' u imposed but torrential outflow : ' 
          write(*,*) ' he = ', he
          write(*,*) ' ue = ', ue
          write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        endif
      endif
    endif
  endif
! ~~~~~~~~~~~~~~~~~~~~~~~~
! torrential inlet -> h, u 
! ~~~~~~~~~~~~~~~~~~~~~~~~
elseif (typcl == 4) then
  bnd_highord_loc = 0
  if (clhorz.eq.0) then
    he = zsu - ze
  else
    he = hu
  endif
  if (clqoru.eq.1) then
    ue = qu/max(eps, he)
  else
    ue = uu
  endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~
! torrential outlet -> none 
! ~~~~~~~~~~~~~~~~~~~~~~~~~
elseif (typcl == 5) then
  bnd_highord_loc = 0
  he = 0d0 ! hi
  ue = ui
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Neuman on h and u (strong imposition)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif (typcl == 6) then
  bnd_highord_loc = 0
  he = max(0d0, hi + dz)
  ue = ui
! ~~~~~~~~
! Periodic
! ~~~~~~~~
elseif (typcl == 7) then
  bnd_highord_loc = 0
  if (i.eq.1) then
    he = h(nx)
    ue = u(nx)
  elseif (i.eq.nx) then
    he = h(1)
    ue = u(1)
  endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Unknown LEFT or RIGHT BOUNDARY CONDITION   '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif

! ----------------------------
! bedload boundary conditions:
! ----------------------------
if (bedload.eqv..True.) then
  ! bedload flux on last cell
  qsedi = qsed(i)
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! Dirichlet, qs imposed
  ! ~~~~~~~~~~~~~~~~~~~~~
  if(typclqs.eq.1) then
    qsede = qsu
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Dirichlet, equilibrium slope
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if(typclqs.eq.4) then
    ! shift cell for computation of qsede
    shift = 0
    if (i.eq.1) then
      j = i+shift
    elseif (i.eq.nx) then
      j = i-shift
    endif
    ! Grass formula -> Jeq replaced by equilibrium velocity
    if (bedload_formula.le.2) then
      ueq = sqrt(Jeq*fric(j)**2*Rh(j)**(4d0/3d0))
      qsede = qs(h(j), ueq, Lx(j), 0d0, fric(j)*bed_fric_coef)
    ! Other formulas (which depend on J)
    else
      qsede = qs(h(j), u(j), Lx(j), Jeq, fric(j)*bed_fric_coef)
    endif
  ! ~~~~~~
  ! Dirichlet, qs computed from hydro ghost cell values
  ! ~~~~~~
  elseif(typclqs.eq.5) then
    qsede = qs(he, ue, Lx(i), 0d0, fric(i)*bed_fric_coef)
  ! ~~~~~~
  ! Neuman
  ! ~~~~~~
  elseif(typclqs.eq.2) then
    qsede = qsed(i)
  ! ~~~~~~~~
  ! Periodic
  ! ~~~~~~~~
  elseif(typclqs.eq.3) then
    if (i.eq.1) then
      qsede = qsed(nx)
    elseif (i.eq.nx) then
      qsede = qsed(1)
    endif
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Unknown LEFT or RIGHT BEDLOAD BOUNDARY CONDITION   '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif

  ! Apply damping zone, bedload_coef, porosity 
  ! (Only for Dirichlet qs, consistent with compute_bedload_solidfluxes)
  if ((typclqs.eq.1).or.(typclqs.eq.4).or.(typclqs.eq.5)) then
    ! damping zone
    if (eps_b_damping.eq.1) then
      b = max(zb(i)-zf(i), 0d0)
      if (b.le.eps_b) then
        rb = b/eps_b
        qsede = qsede*damping_coef(rb)
      endif
    endif
    ! porosity and bedload_coef
    qsede = bedload_coef*qsede/(1.d0 - porosite)
  endif
endif

! ----------------
! Flux computation
! ----------------
! Flux Shallow Water system
! *************************
if (bedload.eqv..False.) then
  ! left bnd
  if (i.eq.1) then
    call flux_bc_sv(i, fh_bord, fq_bord, he, hi, ue, ui, ze, zi, &
                    bnd_flusrc, bnd_highord_loc)
    if (mass_balance) then
      mb_fh_l    = fh_bord
      mb_fq_l    = fq_bord
      mb_fqsed_l = 0d0
    endif
  ! right bnd
  else
    call flux_bc_sv(i, fh_bord, fq_bord, hi, he, ui, ue, zi, ze, &
                    bnd_flusrc, bnd_highord_loc)
    if (mass_balance) then
      mb_fh_r    = fh_bord
      mb_fq_r    = fq_bord
      mb_fqsed_r = 0d0
    endif
  endif

! Flux Shallow Water Exner system
! *******************************
else 
  ! left bnd
  if (i.eq.1) then
    call flux_bc_sve(i, fh_bord, fq_bord, fqsed_bord, he, hi, ue, ui, &
                     ze, zi, qsede, qsedi, bnd_flusrc, bnd_highord_loc)
    if (mass_balance) then
      mb_fh_l    = fh_bord
      mb_fq_l    = fq_bord
      mb_fqsed_l = fqsed_bord
    endif
  ! right bnd
  else
    call flux_bc_sve(i, fh_bord, fq_bord, fqsed_bord, hi, he, ui, ue, &
                     zi, ze, qsedi, qsede, bnd_flusrc, bnd_highord_loc)
    if (mass_balance) then
      mb_fh_r    = fh_bord
      mb_fq_r    = fq_bord
      mb_fqsed_r = fqsed_bord
    endif
  endif
endif

! ----------------------------
! tracers boundary conditions:
! ----------------------------
if (ntrac.gt.0) then
  do k=1,ntrac
    ! ~~~~~~~~~
    ! Dirichlet
    ! ~~~~~~~~~
    if (typclt(k) == 1) then
      te = tu(k)
    ! ~~~~~~
    ! Neuman
    ! ~~~~~~
    elseif (typclt(k) == 2) then
      te = t(i,k)
    ! ~~~~~~~~
    ! Periodic
    ! ~~~~~~~~
    elseif (typclt(k) == 3) then
      if (i.eq.1) then
        te = t(nx,k)
      elseif (i.eq.nx) then
        te = t(1,k)
      endif
    else
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' Unknown LEFT or RIGHT TRACERBOUNDARY CONDITION   '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    ! ~~~~~~~
    ! outflow
    ! ~~~~~~~
    if (fh_bord*sgn.ge.0d0) then
      ft_bord = t(i,k)*fh_bord
    ! ~~~~~~
    ! inflow
    ! ~~~~~~
    else
      ft_bord = te*fh_bord
    endif
    if (i.eq.1) then
      flut_r(i-1,k) = flut_r(i-1,k) + ft_bord
    else
      flut_l(i,k) = flut_l(i,k) + ft_bord
    endif
  enddo
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bc_sv(i, fh_bord, fq_bord, hgi, hdi, ugi, udi, zbgi, zbdi, &
                      bnd_flusrc, bnd_highord)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute boundary fluxes for the Shallow Water system
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: hydrostarec, variable_width, ordre, rec_vel, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       corr_bp, corr_bm, corr_hp, corr_hm, corr_up, corr_um
use module_fluxes_sv, only: flux_sv, flux_src_zb, flux_src_width
implicit none
integer, intent(in)          :: i, bnd_flusrc, bnd_highord
real(kind=dp), intent(inout) :: fh_bord, fq_bord
real(kind=dp)                :: fhg, fqg, fhd, fqd
real(kind=dp)                :: fhg_src, fhd_src, fqg_src, fqd_src
real(kind=dp), intent(in)    :: hgi, hdi, ugi, udi, zbgi, zbdi
real(kind=dp)                :: hg, hd, ug, ud, zbg, zbd, dzp, dzm, hip, him
real(kind=dp)                :: qg, qd, etag, etad, zip12, mc
integer                      :: hydrec = 0
logical                      :: rec_flag
! -----------------------------------------------------------------------------
! compute flu src in bc
if (bnd_flusrc.eq.1) then
  hydrec = hydrostarec
else
  hydrec = 0
endif
! Left and right states
hg = hgi
hd = hdi
ug = ugi
ud = udi
zbg = zbgi
zbd = zbdi
! -----------------------------------------------------------------------------
! second order reconstruction
! ---------------------------
if ((ordre.ge.2).and.(bnd_highord.ne.0)) then
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
        qg = hgi*ug + corr_up(i)
        qd = hdi*ud + corr_um(i+1)
        ug = qg/hg
        ud = qd/hd
      endif
    else
      hg = 0.0d0
      hd = 0.0d0
      ug = 0.0d0
      ud = 0.0d0
    endif
  else
    ! wet/dry interface: stay first order
    rec_flag = .false.
  endif
endif

! hydrostatic reconstruction
! --------------------------
! Audusse & al. A fast and stable well-balanced scheme
! with hydrostatic reconstruction for shallow water flows
if (hydrec.eq.1) then
  dzm = max(0.d0, zbd- zbg)
  him = max(0.d0, hg - dzm)
  dzp = max(0.d0, zbg- zbd)
  hip = max(0.d0, hd - dzp)
! Chen and Noelle. A new hydrostatic reconstruction
! scheme based on subcell reconstructions
elseif (hydrec.eq.2) then
  etag = hg + zbg
  etad = hd + zbd
  zip12 = min(max(zbg, zbd), min(etag, etad))
  him = min(etag - zip12, hg)
  hip = min(etad - zip12, hd)
else
  him = hg
  hip = hd
endif
! 
! shallow water flux
! ------------------
call flux_sv(mc, fhg, fhd, fqg, fqd, him, hip, ug, ud, zbg, zbd)
! left bnd (i=1 <-> fluh_cell(i+1) - flu_rigth)
if (i.eq.1) then
  fluh_r(i-1) = fluh_r(i-1) + fhd
  fluq_r(i-1) = fluq_r(i-1) + fqd
  fh_bord = fhd
  fq_bord = fqd
! right bnd (i=nx <-> fluh_cell(i) + flu_left)
else
  fluh_l(i) = fluh_l(i) + fhg
  fluq_l(i) = fluq_l(i) + fqg
  fh_bord = fhg
  fq_bord = fqg
endif
! 
! bathy source
! ------------
if (bnd_flusrc.eq.1) then
  call flux_src_zb(fhg_src, fhd_src, fqg_src, fqd_src, &
                   hg, hd, him, hip, zbg, zbd, ug, ud)
  ! left bnd 
  if (i.eq.1) then
    fluh_src_r(i-1) = fluh_src_r(i-1) + fhd_src
    fluq_src_r(i-1) = fluq_src_r(i-1) + fqd_src
    ! second order correction
    if (hydrostarec.ne.0 .and. ordre.ge.2 .and. rec_flag.eqv..true.) then
      fluq_src_r(i-1) = fluq_src_r(i-1) - 0.5d0*g*(hd+hdi)*corr_bm(i)
    endif
  ! right bnd
  else
    fluh_src_l(i) = fluh_src_l(i) + fhg_src
    fluq_src_l(i) = fluq_src_l(i) + fqg_src
    ! second order correction
    if (hydrostarec.ne.0 .and. ordre.ge.2 .and. rec_flag.eqv..true.) then
      fluq_src_l(i) = fluq_src_l(i) + 0.5d0*g*(hg+hgi)*corr_bp(i)
    endif
  endif
endif ! endif flusrc
!
! width source
! ------------
if ((bnd_flusrc.eq.1).and.(variable_width.eqv..True.)) then
  call flux_src_width(i, fhg_src, fhd_src, fqg_src, fqd_src, hg, hd, ug, ud)
  ! left bnd
  if (i.eq.1) then
    fluh_src_r(i-1) = fluh_src_r(i-1) + fhd_src
    fluq_src_r(i-1) = fluq_src_r(i-1) + fqd_src
  ! right bnd
  else
    fluh_src_l(i) = fluh_src_l(i) + fhg_src
    fluq_src_l(i) = fluq_src_l(i) + fqg_src
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine flux_bc_sve(i, fh_bord, fq_bord, fqsb_bord, hgi, hdi, ugi, udi, &
                       zbgi, zbdi, qsgi, qsdi, bnd_flusrc, bnd_highord)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute boundary fluxes for the Shallow Water system
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: hydrostarec, variable_width, ordre, rec_vel, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       fluqsed_l, fluqsed_r, fluqsed_src_l, fluqsed_src_r, &
                       corr_bp, corr_bm, corr_hp, corr_hm, corr_up, corr_um
use module_fluxes_sv,  only: flux_src_zb
use module_fluxes_sve, only: flux_sve, flux_src_zb_sve, flux_src_width_sve
implicit none
integer, intent(in)          :: i, bnd_flusrc, bnd_highord
real(kind=dp), intent(inout) :: fh_bord, fq_bord, fqsb_bord
real(kind=dp), intent(in)    :: hgi, hdi, ugi, udi, zbgi, zbdi, qsgi, qsdi
real(kind=dp)                :: fhd, fqd, fqsd, fhg, fqg, fqsg
real(kind=dp)                :: fhg_src, fhd_src
real(kind=dp)                :: fqg_src, fqd_src, fqsg_src, fqsd_src
real(kind=dp)                :: hg, hd, ug, ud, zbg, zbd, qsg, qsd, dzp, dzm
real(kind=dp)                :: hip, him, qg, qd, etag, etad, zip12, mc
real(kind=dp)                :: lambda_roe (3), A_roe (3,3), D_roe (3,3)
integer                      :: hydrec = 0
logical                      :: rec_flag
! -----------------------------------------------------------------------------
! compute flu src in bc
if (bnd_flusrc.eq.1) then
  hydrec = hydrostarec
else
  hydrec = 0
endif
! Left and right states
hg = hgi
hd = hdi
ug = ugi
ud = udi
zbg = zbgi
zbd = zbdi
qsg = qsgi
qsd = qsdi
! -----------------------------------------------------------------------------
! second order reconstruction
! ---------------------------
if ((ordre.ge.2).and.(bnd_highord.ne.0)) then
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
        qg = hgi*ug + corr_up(i)
        qd = hdi*ud + corr_um(i+1)
        ug = qg/hg
        ud = qd/hd
      endif
    else
      hg = 0.0d0
      hd = 0.0d0
      ug = 0.0d0
      ud = 0.0d0
    endif
  else
    ! wet/dry interface: stay first order
    rec_flag = .false.
  endif
endif

! hydrostatic reconstruction
! --------------------------
! Audusse & al. A fast and stable well-balanced scheme
! with hydrostatic reconstruction for shallow water flows
if (hydrec.eq.1) then
  dzm = max(0.d0, zbd- zbg)
  him = max(0.d0, hg - dzm)
  dzp = max(0.d0, zbg- zbd)
  hip = max(0.d0, hd - dzp)
! Chen and Noelle. A new hydrostatic reconstruction
! scheme based on subcell reconstructions
elseif (hydrec.eq.2) then
  etag = hg + zbg
  etad = hd + zbd
  zip12 = min(max(zbg, zbd), min(etag, etad))
  him = min(etag - zip12, hg)
  hip = min(etad - zip12, hd)
else
  him = hg
  hip = hd
endif

! Flux SVE
! --------
call flux_sve(mc, fhg, fhd, fqg, fqd, fqsg, fqsd, i, him, hip, ug, ud, &
              zbg, zbd, qsg, qsd, lambda_roe, A_roe, D_roe)

! Flux update
! -----------
! left bnd (i=1 <-> fluh_cell(i+1) - flu_rigth)
if (i.eq.1) then
  fluh_r(i-1) = fluh_r(i-1) + fhd
  fluq_r(i-1) = fluq_r(i-1) + fqd
  fluqsed_r(i-1) = fluqsed_r(i-1) + fqsd
  fh_bord = fhd
  fq_bord = fqd
  fqsb_bord = fqsd
! right bnd (i=nx <-> fluh_cell(i) + flu_left)
else
  fluh_l(i) = fluh_l(i) + fhg
  fluq_l(i) = fluq_l(i) + fqg
  fluqsed_l(i) = fluqsed_l(i) + fqsg
  fh_bord = fhg
  fq_bord = fqg
  fqsb_bord = fqsg
endif

! bathy source
! ------------
if (bnd_flusrc.eq.1) then
  call flux_src_zb_sve(fhg_src, fhd_src, fqg_src, fqd_src, fqsg_src, fqsd_src,&
                       hg, hd, him, hip, zbg, zbd, lambda_roe, A_roe, D_roe)
  ! left bnd
  if (i.eq.1) then
    fluh_src_r(i-1) = fluh_src_r(i-1) + fhd_src
    fluq_src_r(i-1) = fluq_src_r(i-1) + fqd_src
    fluqsed_src_r(i-1) = fluqsed_src_r(i-1) + fqsd_src
    ! second order correction
    if (hydrostarec.ne.0 .and. ordre.ge.2 .and. rec_flag.eqv..true.) then
      fluq_src_r(i-1) = fluq_src_r(i-1) - 0.5d0*g*(hd+hdi)*corr_bm(i)
    endif
  ! right bnd
  else
    fluh_src_l(i) = fluh_src_l(i) + fhg_src
    fluq_src_l(i) = fluq_src_l(i) + fqg_src
    fluqsed_src_l(i) = fluqsed_src_l(i) + fqsg_src
    ! second order correction
    if (hydrostarec.ne.0 .and. ordre.ge.2 .and. rec_flag.eqv..true.) then
      fluq_src_l(i) = fluq_src_l(i) + 0.5d0*g*(hg+hgi)*corr_bp(i)
    endif
  endif
endif
!
! width source
! ------------
if ((bnd_flusrc.eq.1).and.(variable_width.eqv..True.)) then
  call flux_src_width_sve(fhg_src, fhd_src, fqg_src, fqd_src, &
                          fqsg_src, fqsd_src, i, hg, hd, ug, ud, qsg, qsd, &
                          lambda_roe, A_roe, D_roe)
  ! left bnd
  if (i.eq.1) then
    fluh_src_r(i-1) = fluh_src_r(i-1) + fhd_src
    fluq_src_r(i-1) = fluq_src_r(i-1) + fqd_src
    fluqsed_src_r(i-1) = fluqsed_src_r(i-1) + fqsd_src
  ! right bnd
  else
    fluh_src_l(i) = fluh_src_l(i) + fhg_src
    fluq_src_l(i) = fluq_src_l(i) + fqg_src
    fluqsed_src_l(i) = fluqsed_src_l(i) + fqsg_src
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine update_gc(h, q, u, ht, t)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Update ghost cells :
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ntrac
implicit none
real(kind=dp), intent(inout) :: h(0:nx+1), q(0:nx+1), u(0:nx+1)
real(kind=dp), intent(inout) :: ht(0:nx+1,ntrac), t(0:nx+1,ntrac)
integer                      :: k
! -----------------------------------------------------------------------------
! left ghost cell
! -----------------------------------------------------------------------------
h(0) = h(1) 
q(0) = q(1)
u(0) = u(1)
! tracers
if (ntrac.gt.0) then
  do k=1,ntrac
    ht(0,k) = ht(1,k)
    t (0,k) = t (1,k)
  enddo
endif
! -----------------------------------------------------------------------------
! right ghost cell
! -----------------------------------------------------------------------------
h(nx+1) = h(nx)
q(nx+1) = q(nx)
u(nx+1) = u(nx)
! tracers
if (ntrac.gt.0) then
  do k=1,ntrac
    ht(nx+1,k) = ht(nx,k)
    t (nx+1,k) = t (nx,k)
  enddo
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine compute_cdl_steady(upstream_control, skip_update, &
                              ha0, hb0, qa0, qb0)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! boundary conditions
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: time, niter, use_bcfile_l, use_bcfile_r, &
                       bctype_l, bctypet_l, zsa, ha, ua, qa, ta, &
                       bctype_r, bctypet_r, zsb, hb, ub, qb, tb, &
                       time_bcl, time_bcr, field_bcl, field_bcr, &
                       vars_bcl, vars_bcr, nvars_bcl, nvars_bcr, qsa, qsb, &
                       eps_evol_h, eps_evol_q, exner, &
                       bctypeqs_l, bctypeqs_r, qsa, qsb, Jeqa, Jeqb, &
                       bc_idx_l, bc_idx_r
implicit none
real(kind=dp), intent(inout)   :: ha0, hb0, qa0, qb0
logical,       intent(out)     :: upstream_control, skip_update
! -----------------------------------------------------------------------------
! check control side
if ((bctype_l==2).and.(bctype_r==3)) then
  upstream_control = .True.
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' UPSTREAM CONTROL NOT IMPLEMENTED YET  '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
else if ((bctype_l==3).and.(bctype_r==2)) then
  upstream_control = .False.
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' BOUNDARY CONDITION TYPE NOT COMPATIBLE WITH STEADY KERNEL  '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
! LEFT BOUNDARY CONDITION
! -----------------------------------------------------------------------------
! ghost cell data from liq boundary file (time interpolation)
if (use_bcfile_l.eqv..True.) then
  call ghostcell_values(1, vars_bcl, nvars_bcl, time, time_bcl, &
                        field_bcl, bc_idx_l, &
                        zsa, ha, ua, qa, qsa, ta)
end if
! compute fluxes of left boundary
call cdl_steady(1, bctype_l, bctypet_l, bctypeqs_l, zsa, ha, ua, qa, ta, &
                qsa, Jeqa)
! -----------------------------------------------------------------------------
! RIGHT BOUNDARY CONDITION
! -----------------------------------------------------------------------------
! ghost cell data from liq boundary file (time interpolation)
if (use_bcfile_r.eqv..True.) then
  call ghostcell_values(nx+1, vars_bcr, nvars_bcr, time, time_bcr, &
                        field_bcr, bc_idx_r, &
                        zsb, hb, ub, qb, qsb, tb)
end if
! compute fluxes of right boundary
call cdl_steady(nx, bctype_r, bctypet_r, bctypeqs_r, zsb, hb, ub, qb, tb, &
                qsb, Jeqb)
! -----------------------------------------------------------------------------
! check if imposed values change, else skip update
if (niter.eq.1) then
  ha0 = 0d0
  hb0 = 0d0
  qa0 = 0d0
  qb0 = 0d0
endif
if ((abs(ha-ha0).ge.eps_evol_h).or.&
    (abs(hb-hb0).ge.eps_evol_h).or.&
    (abs(qa-qa0).ge.eps_evol_q).or.&
    (abs(qb-qb0).ge.eps_evol_q)) then
  skip_update = .False.
  ha0 = ha
  hb0 = hb
  qa0 = qa
  qb0 = qb
else 
  skip_update = .True.
endif
if ((niter==1).or.(exner)) then
  skip_update = .False.
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine cdl_steady(i, typcl, typclt, typclqs, zsu, hu, uu, qu, tu, qsu, Jeq)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! boundary conditions for the steady kernel
! 
! Ghost cells:
!  u : user imposed values
!  i : internal
!  e : external
! 
!   Wu=ha,qa                             Wu=hb,qb <- user defined
!    o---|---o---     ...      ---o---|---o
!    0   A   1                   nx   B  nx+1
!   We       Wi                  Wi       We
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps, eps_num, Rh, Lx, zb, zf, h, u, q, fric, ntrac, ht, t, & 
                       clhorz, clqoru, &
                       bedload, bedload_formula, qsed,  mass_balance, &
                       mb_fh_l, mb_fq_l, mb_fqsed_l, &
                       mb_fh_r, mb_fq_r, mb_fqsed_r, &
                       bedload_coef, porosite, eps_b, eps_b_damping, bed_fric_coef
use module_bedload, only: qs, damping_coef
implicit none
integer, intent(in)       :: i, typcl, typclqs
integer, intent(in)       :: typclt(:)
real(kind=dp), intent(in) :: zsu, hu, uu, qu, qsu, Jeq
real(kind=dp), intent(in) :: tu(:)
real(kind=dp)             :: he, qe, ze, te, qsede, ueq, ue
real(kind=dp)             :: hcri
real(kind=dp)             :: b, rb
integer                   :: j, k, shift
! -----------------------------------------------------------------------------
he = h(i)
qe = q(i)

! bottom
if (i.eq.1) then
  ze = zb(i-1)
else
  ze = zb(i+1)
endif

! ~~~~~~~~~
! h imposed
! ~~~~~~~~~
if (typcl == 2) then
! free surface imposed ~> compute water depth
  if (clhorz.eq.0) then
    he = zsu-ze
  else
    he = hu
  endif
  ! prevent torrential on boundary
  hcri = (q(i)/sqrt(g))**(2d0/3d0)
  if (he.le.hcri) then
    he = hcri + eps
  endif
! ~~~~~~~~~
! u imposed
! ~~~~~~~~~
elseif (typcl == 3) then
  ! flowrate imposed
  if (clqoru.eq.1) then
    qe = qu
  ! velocity imposed ~> compute flowrate
  else
    qe = uu*max(h(i), eps)
  endif
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : BC TYPE NOT COMPATIBLE WITH STEADY SOLVER'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif

! Ghost cells
if (i.eq.1) then
  ! Left
  h(0) = he
  q(0) = qe
else
  ! right
  h(nx+1) = he
  q(nx+1) = qe
endif

! -----------------------------------------------------------------------------
! bedload boundary conditions:
! ----------------------------
if (bedload.eqv..True.) then
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! Dirichlet, qs imposed
  ! ~~~~~~~~~~~~~~~~~~~~~
  if(typclqs.eq.1) then
    qsede = qsu
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Dirichlet, equilibrium slope
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if(typclqs.eq.4) then
    ! shift cell for computation of qsede
    shift = 0
    if (i.eq.1) then
      j = i+shift
    elseif (i.eq.nx) then
      j = i-shift
    endif
    ! Grass formula -> Jeq replaced by equilibrium velocity
    if (bedload_formula.le.2) then
      ueq = sqrt(Jeq*fric(j)**2*Rh(j)**(4d0/3d0))
      qsede = qs(h(j), ueq, Lx(j), 0d0, fric(j)*bed_fric_coef)
    ! Other formulas (which depend on J)
    else
      qsede = qs(h(j), u(j), Lx(j), Jeq, fric(j)*bed_fric_coef)
    endif
  ! ~~~~~~
  ! Dirichlet, qs computed from hydro ghost cell values
  ! ~~~~~~
  elseif(typclqs.eq.5) then
    ue = qe/max(he, eps)
    qsede = qs(he, ue, Lx(i), 0d0, fric(i)*bed_fric_coef)
  ! ~~~~~~
  ! Neuman
  ! ~~~~~~
  elseif(typclqs.eq.2) then
    qsede = qsed(i)
  ! ~~~~~~~~
  ! Periodic
  ! ~~~~~~~~
  elseif(typclqs.eq.3) then
    if (i.eq.1) then
      qsede = qsed(nx)
    elseif (i.eq.nx) then
      qsede = qsed(1)
    endif
  else
    write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unknown LEFT or RIGHT BEDLOAD BOUNDARY CONDITION   '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif

  ! Apply damping zone, bedload_coef, porosity 
  ! (Only for Dirichlet qs, consistent with compute_bedload_solidfluxes)
  if ((typclqs.eq.1).or.(typclqs.eq.4).or.(typclqs.eq.5)) then
    ! damping zone
    if (eps_b_damping.eq.1) then
      b = max(zb(i)-zf(i), 0d0)
      if (b.le.eps_b) then
        rb = b/eps_b
        qsede = qsede*damping_coef(rb)
      endif
    endif
    ! bedload_coef and porosity
    qsede = bedload_coef*qsede/(1.d0 - porosite)
  endif

  ! Ghost cells
  if (i.eq.1) then
    ! Left
    qsed(0) = qsede
  else
    ! right
    qsed(nx+1) = qsede
  endif
endif

! -----------------------------------------------------------------------------
! mass balance:
! -------------
if (i.eq.1) then
  if (mass_balance) then
    mb_fh_l    = qe
    mb_fq_l    = 0d0
    mb_fqsed_l = qsede
  endif
elseif (i.eq.nx) then
  if (mass_balance) then
    mb_fh_r    = qe
    mb_fq_r    = 0d0
    mb_fqsed_r = qsede
  endif
endif

! -----------------------------------------------------------------------------
! tracers boundary conditions:
! ----------------------------
if (ntrac.gt.0) then
  do k=1,ntrac
    ! ~~~~~~~~~
    ! Dirichlet
    ! ~~~~~~~~~
    if (typclt(k) == 1) then
      te = tu(k)
    ! ~~~~~~
    ! Neuman
    ! ~~~~~~
    elseif (typclt(k) == 2) then
      te = t(i,k)
    ! ~~~~~~~~
    ! Periodic
    ! ~~~~~~~~
    elseif (typclt(k) == 3) then
      if (i.eq.1) then
        te = t(nx,k)
      elseif (i.eq.nx) then
        te = t(1,k)
      endif
    else
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unknown LEFT or RIGHT TRACER BOUNDARY CONDITION   '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif

    ! Ghost cells
    if (i.eq.1) then
      ! Left
      t(0,k) = te
      ht(0,k) = h(0)*te
    else
      ! right
      t(nx+1,k) = te
      ht(nx+1,k) = h(nx+1)*te
    endif
  enddo
endif

! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************