! *****************************************************************************
module module_alloc

! Allocation of arrays
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
subroutine allocation_geom
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Tables allocation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: nx, dx, ci, x, zb, zf, Lx, dLx, Rh, exner, &
                       variable_width
implicit none      
! -----------------------------------------------------------------------------
! allocation of geometrical parameters
! -----------------------------------------------------------------------------
allocate(dx (0:nx))
allocate(ci (1:nx))       ! ghost cells not included
allocate(x  (0:nx+1))     ! ghost cells included
allocate(zb (0:nx+1))
allocate(Rh (0:nx+1))
allocate(Lx (0:nx+1))
if (exner.eqv..True.) then
  allocate(zf (0:nx+1))
endif
if (variable_width.eqv..True.) then
  allocate(dLx(0:nx+1))
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine allocation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Tables allocation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: nx, h, u, q, t, ht, time_scheme, &
                       ntrac, ordre, aux1, aux2, fric, difu, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r, &
                       corr_bp, corr_bm, corr_hp, corr_hm, &
                       corr_up, corr_um, corr_tp, corr_tm, &
                       exner, qsed, fluqsed_l, fluqsed_r, &
                       fluqsed_src_l, fluqsed_src_r, taub, taus, &
                       fluz_susp, flut_susp, Ceq_susp, Csratio
implicit none
! -----------------------------------------------------------------------------
! main variables
! --------------
allocate(h          (0:nx+1))
allocate(u          (0:nx+1))
allocate(q          (0:nx+1))
allocate(fric       (0:nx+1))
allocate(difu       (0:nx+1))
if (time_scheme.ne.0) then
  allocate(fluh_l     (0:nx))
  allocate(fluh_r     (0:nx))
  allocate(fluq_l     (0:nx))
  allocate(fluq_r     (0:nx))
  allocate(fluh_src_l (0:nx))
  allocate(fluh_src_r (0:nx))
  allocate(fluq_src_l (0:nx))
  allocate(fluq_src_r (0:nx))
endif

! auxiliary arrays
! ----------------
allocate(aux1       (0:nx+1))
allocate(aux2       (0:nx+1))

! tracers
! -------
if (ntrac.gt.0) then
  allocate(t          (0:nx+1,ntrac))
  allocate(ht         (0:nx+1,ntrac))
  if (time_scheme.ne.0) then
    allocate(flut_l     (0:nx,ntrac))
    allocate(flut_r     (0:nx,ntrac))
    allocate(flut_src_l (0:nx,ntrac))
    allocate(flut_src_r (0:nx,ntrac))
  endif
endif

! hydro second order
! ------------------
if ((ordre.ge.2).and.((time_scheme.ne.0))) then
  allocate(corr_bp    (0:nx+1))
  allocate(corr_bm    (0:nx+1))
  allocate(corr_hp    (0:nx+1))
  allocate(corr_hm    (0:nx+1))
  allocate(corr_up    (0:nx+1))
  allocate(corr_um    (0:nx+1))
  if (ntrac.gt.0) then
    allocate(corr_tp  (0:nx+1,ntrac))
    allocate(corr_tm  (0:nx+1,ntrac))
  endif
endif

! Sediments
! ---------
if (exner.eqv..True.) then
  allocate(qsed         (0:nx+1))
  allocate(taub         (0:nx+1))
  allocate(taus         (0:nx+1))
  allocate(fluqsed_l    (0:nx))
  allocate(fluqsed_r    (0:nx))
  allocate(fluqsed_src_l(0:nx))
  allocate(fluqsed_src_r(0:nx))
  allocate(fluz_susp    (1:nx))
  allocate(flut_susp    (1:nx))
  allocate(Ceq_susp     (1:nx))
  allocate(Csratio      (1:nx))
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************




! *****************************************************************************
subroutine deallocation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Tables allocation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: x, ci, dx, zb, h, u, q, t, ht, Lx, dLx, Rh, ntrac,  &
                       ordre, aux1, aux2, fric, difu, variable_width, &
                       fluh_l, fluh_src_l, fluh_r, fluh_src_r, &
                       fluq_l, fluq_src_l, fluq_r, fluq_src_r, &
                       flut_l, flut_src_l, flut_r, flut_src_r, &
                       corr_bp, corr_bm, corr_hp, corr_hm, &
                       corr_up, corr_um, corr_tp, corr_tm, &
                       exner, zf, qsed, fluqsed_l, fluqsed_r, &
                       fluqsed_src_l, fluqsed_src_r, taub, taus, &
                       fluz_susp, flut_susp, Ceq_susp, Csratio, time_scheme
implicit none
! -----------------------------------------------------------------------------
! Geometry
! --------
deallocate(x)
deallocate(ci)
deallocate(dx)
deallocate(zb)
deallocate(Rh)
deallocate(Lx)
if (exner.eqv..True.) then
  deallocate(zf)
endif
if (variable_width.eqv..True.) then
  deallocate(dLx)
endif

! main variables
! --------------
deallocate(h)
deallocate(u)
deallocate(q)
deallocate(fric)
deallocate(difu)
if (time_scheme.ne.0) then
  deallocate(fluh_l)
  deallocate(fluh_r)
  deallocate(fluq_l)
  deallocate(fluq_r)
  deallocate(fluh_src_l)
  deallocate(fluh_src_r)
  deallocate(fluq_src_l)
  deallocate(fluq_src_r)
endif

! auxiliary arrays
! ----------------
deallocate(aux1)
deallocate(aux2)

! tracers
! -------
if (ntrac.gt.0) then
  deallocate(t)
  deallocate(ht)
  if (time_scheme.ne.0) then
    deallocate(flut_l)
    deallocate(flut_r)
    deallocate(flut_src_l)
    deallocate(flut_src_r)
  endif
endif

! hydro second order
! ------------------
if (ordre.ge.2) then
  deallocate(corr_bp)
  deallocate(corr_bm)
  deallocate(corr_hp)
  deallocate(corr_hm)
  deallocate(corr_up)
  deallocate(corr_um)
  if (ntrac.gt.0) then
    deallocate(corr_tp)
    deallocate(corr_tm)
  endif
endif

! Sediments
! ---------
if (exner.eqv..True.) then
  deallocate(qsed)
  deallocate(taub)
  deallocate(taus)
  deallocate(fluqsed_l)
  deallocate(fluqsed_r)
  deallocate(fluqsed_src_l)
  deallocate(fluqsed_src_r)
  deallocate(fluz_susp)
  deallocate(flut_susp)
  deallocate(Ceq_susp)
  deallocate(Csratio)
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************




end module
! *****************************************************************************