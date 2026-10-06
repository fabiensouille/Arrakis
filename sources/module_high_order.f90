! *****************************************************************************
module module_high_order
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Second order reconstruction module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, eps_num, nx, ntrac, g
implicit none
! ----------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                         HIGH ORDER RECONSTRUCTIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine compute_reconstructions(zb, h, u, q, t)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ordre, ordre_trac, rec_method, rec_vel, ntrac, &
                       muscl_ilim_h, muscl_ilim_u, muscl_ilim_t, &
                       corr_bp, corr_bm, corr_hp, corr_hm, &
                       corr_up, corr_um, corr_tp, corr_tm
implicit none
real(kind=dp), intent(in) :: zb (0:nx+1)
real(kind=dp), intent(in) :: h  (0:nx+1)
real(kind=dp), intent(in) :: u  (0:nx+1)
real(kind=dp), intent(in) :: q  (0:nx+1)
real(kind=dp), intent(in) :: t  (0:nx+1,ntrac)
integer                   :: k
! -----------------------------------------------------------------------------
! initialize corr arrays
! -----------------------------------------------------------------------------
corr_bp = 0d0
corr_bm = 0d0
corr_hp = 0d0
corr_hm = 0d0
corr_up = 0d0
corr_um = 0d0
if (ntrac.gt.0) then
  corr_tp = 0d0
  corr_tm = 0d0
endif

! Gradient reconstructions (see INRIA M-O.B. and E.A.)
! ----------------------------------------------------
if (rec_method.eq.1) then
  ! ***************
  ! ORDER 2 / HYDRO
  ! ***************
  if (ordre.eq.2) then
    call corr_grad_ord2(zb, muscl_ilim_h, corr_bp, corr_bm)
    call corr_grad_ord2(h , muscl_ilim_h, corr_hp, corr_hm)
    if (rec_vel.eqv..true.) then
      call corr_grad_ord2(u, muscl_ilim_u, corr_up, corr_um)
    else
      call corr_grad_ord2(q, muscl_ilim_u, corr_up, corr_um)
    endif
  elseif (ordre.gt.2) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' RECONSTRUCTION METHOD = 1, ONLY SECOND ORDER AVAILABLE '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif

  ! *****************
  ! ORDER 2 / TRACERS
  ! *****************
  if ((ntrac.gt.0).and.(ordre_trac.eq.2)) then
    do k=1,ntrac
      call corr_grad_ord2(t(:,k), muscl_ilim_t, corr_tp(:,k), corr_tm(:,k))
    enddo
  elseif ((ntrac.gt.0).and.(ordre_trac.gt.2)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' RECONSTRUCTION METHOD = 1, ONLY SECOND ORDER AVAILABLE '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif 

! MUSCL reconstructions
! ---------------------
elseif (rec_method.eq.2) then
  ! ***************
  ! ORDER 2 / HYDRO
  ! ***************
  if (ordre.eq.2) then
    call corr_var_ord2(zb, muscl_ilim_h, corr_bp, corr_bm)
    call corr_var_ord2(h , muscl_ilim_h, corr_hp, corr_hm)
    if (rec_vel.eqv..true.) then          
      call corr_var_ord2(u, muscl_ilim_u, corr_up, corr_um)
    else
      call corr_var_ord2(q, muscl_ilim_u, corr_up, corr_um)
    endif
  ! *******************
  ! ORDER 3 ENO / HYDRO
  ! *******************
  elseif (ordre.eq.3) then
    call corr_var_ord3(zb, muscl_ilim_h, corr_bp, corr_bm)
    call corr_var_ord3(h , muscl_ilim_h, corr_hp, corr_hm)
    if (rec_vel.eqv..true.) then
      call corr_var_ord3(u, muscl_ilim_u, corr_up, corr_um)
    else
      call corr_var_ord3(q, muscl_ilim_u, corr_up, corr_um)
    endif
  ! *********************
  ! ORDER 5: WENO / HYDRO
  ! *********************
  elseif (ordre.eq.5) then
    call corr_weno_ord5(zb, corr_bp, corr_bm)
    call corr_weno_ord5(h , corr_hp, corr_hm)
    if (rec_vel.eqv..true.) then
      call corr_weno_ord5(u, corr_up, corr_um)
    else
      call corr_weno_ord5(q, corr_up, corr_um)
    endif
  else
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ORDER NOT AVAILABLE WITH, RECONSTRUCTION METHOD = 2 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  ! TRACERS
  if (ntrac.gt.0) then
    ! *****************
    ! ORDER 2 / TRACERS
    ! *****************
    if (ordre_trac.eq.2) then
      do k=1,ntrac
        call corr_var_ord2(t(:,k), muscl_ilim_t, corr_tp(:,k), corr_tm(:,k))
      enddo
    ! *****************
    ! ORDER 3 / TRACERS
    ! *****************
    elseif (ordre_trac.eq.3) then
      do k=1,ntrac
        call corr_var_ord3(t(:,k), muscl_ilim_t, corr_tp(:,k), corr_tm(:,k))
      enddo
    ! *****************
    ! ORDER 5 / TRACERS
    ! *****************
    elseif (ordre_trac.eq.5) then            
      do k=1,ntrac
        call corr_weno_ord5(t(:,k), corr_tp(:,k), corr_tm(:,k))
      enddo
    else
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ORDER NOT AVAILABLE WITH, RECONSTRUCTION METHOD = 2 '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
  endif
! -----------------------------------------------------------------------------
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Unknown RECONSTRUCTION METHOD '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                  SECOND ORDER GRADIENTS RECONSTRUCTIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine corr_grad_ord2(f, ilim, corr_fp, corr_fm)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Returns correction arrays for field f:
! corr_fp         : f(i+1/2)- = f(i)   + corr_fp(i)
! corr_fm         : f(i+1/2)+ = f(i+1) + corr_fm(i+1)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: x, bnd_highord
implicit none
real(kind=dp), intent(inout) :: corr_fp (0:nx+1)
real(kind=dp), intent(inout) :: corr_fm (0:nx+1)
real(kind=dp), intent(in)    :: f (0:nx+1)
integer,       intent(in)    :: ilim
integer                      :: i
real(kind=dp)                :: gradf, dx
! -----------------------------------------------------------------------------
! loop on cell interfaces
!$OMP parallel do  &
!$OMP private(gradf) &
!$OMP shared(x, corr_fp, corr_fm, f, ilim)
do i=1,nx
  gradf = gradient(i, f, ilim)
  corr_fp(i) = +0.5d0*(x(i+1) - x(i))*gradf
  corr_fm(i) = -0.5d0*(x(i) - x(i-1))*gradf
enddo
!$OMP end parallel do
! -----------------------------------------------------------------------------
! addapted reconstructions on boundary edges
! (displaced gradients evaluations)
if (bnd_highord.eq.1) then
  ! left bnd:
  gradf = gradient(1, f, ilim)
  dx = x(1)-x(0)
  corr_fp(0) = +0.5d0*dx*gradf
  corr_fm(0) = -0.5d0*dx*gradf
  ! right bnd:
  gradf = gradient(nx, f, ilim)
  dx = x(nx+1)-x(nx)
  corr_fp(nx+1) = +0.5d0*dx*gradf
  corr_fm(nx+1) = -0.5d0*dx*gradf
else
  corr_fp(0) = 0.d0
  corr_fm(0) = 0.d0
  !corr_fm(1) = 0.d0
  !corr_fp(nx) = 0d0
  corr_fp(nx+1) = 0d0
  corr_fm(nx+1) = 0d0
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
function gradient(i, f, ilimitor)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: x, eps_num
implicit none
integer, intent(in)       :: i
integer, intent(in)       :: ilimitor
real(kind=dp), intent(in) :: f(0:nx+1)
real(kind=dp)             :: gradp, gradm, gradp2, gradm2, aux, gradient
real(kind=dp), parameter  :: beta = 0.33
! -----------------------------------------------------------------------------
gradp = (f(i+1)-f(i))/(x(i+1)-x(i))
gradm = (f(i)-f(i-1))/(x(i)-x(i-1))

! ponderation (see 2D)
!gradm = (1.d0+beta)*gradm - beta*gradp

! Minmod flux limitor
! -------------------
if (ilimitor.eq.1) then
  gradient = 0.5d0*(dsign(1.0d0,gradp)+dsign(1.0d0,gradm)) &
                  * min(dabs(gradp),dabs(gradm))

! Van Albala flux limitor
! -----------------------
elseif (ilimitor.eq.2) then
  if (dsign(1.0d0,gradp).ne.dsign(1.0d0,gradm)) then
    gradient = 0.0d0
  else
    aux = 0.5d0*(1.0d0+dsign(1.d0, gradp*gradm))
    gradp2 = max(gradp**2, eps_num)
    gradm2 = max(gradm**2, eps_num)
    gradient = aux*(gradp2*gradm + gradm2*gradp)/(gradm2+gradp2)
  endif

! Other flux limitors
! -------------------
elseif (ilimitor.gt.2) then
  gradient = gradp*Limitor(gradm/max(gradp, eps_num), ilimitor)

endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                         VARIABLE RECONSTRUCTIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine corr_var_ord2(f, ilim, corr_fp, corr_fm)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Returns correction arrays for field f:
! corr_fp         : f(i+1/2)- = f(i)   + corr_fp(i)
! corr_fm         : f(i+1/2)+ = f(i+1) + corr_fm(i+1)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: bnd_highord
implicit none
real(kind=dp), intent(inout) :: corr_fp (0:nx+1)
real(kind=dp), intent(inout) :: corr_fm (0:nx+1)
real(kind=dp), intent(in)    :: f (0:nx+1)
integer,       intent(in)    :: ilim
integer                      :: i
real(kind=dp)                :: rif, ff, dfp
! -----------------------------------------------------------------------------
! loop on cell interfaces
!$OMP parallel do  &
!$OMP private(dfp, rif, ff) &
!$OMP shared(corr_fp, corr_fm, f, ilim)
do i=1,nx
  dfp = f(i+1)-f(i)
  if (abs(dfp).lt.eps_num) then
    corr_fp(i) = 0.d0
    corr_fm(i) = 0.d0
  else
    rif = (f(i)-f(i-1))/dfp
    ff = Limitor(rif, ilim)
    corr_fp(i) = +0.5d0*dfp*ff
    corr_fm(i) = -0.5d0*dfp*ff
  endif
enddo
!$OMP end parallel do
! -----------------------------------------------------------------------------
! addapted reconstructions on boundary edges
! (displaced gradients evaluations)
if (bnd_highord.eq.1) then
  ! left bnd:
  dfp = f(2)-f(1)
  if (abs(dfp).lt.eps_num) then
    corr_fp(0) = 0.d0
    corr_fm(0) = 0.d0
  else
    rif = (f(1)-f(0))/dfp
    ff = Limitor(rif, ilim)
    corr_fp(0) = +0.5d0*dfp*ff
    corr_fm(0) = -0.5d0*dfp*ff
  endif
  ! right bnd:
  dfp = f(nx+1)-f(nx)
  if (abs(dfp).lt.eps_num) then
    corr_fp(nx+1) = 0.d0
    corr_fm(nx+1) = 0.d0
  else
    rif = (f(nx)-f(nx-1))/dfp
    ff = Limitor(rif, ilim)
    corr_fp(nx+1) = +0.5d0*dfp*ff
    corr_fm(nx+1) = -0.5d0*dfp*ff
  endif
else
  corr_fp(0) = 0.d0
  corr_fm(0) = 0.d0
  corr_fp(nx+1) = 0d0
  corr_fm(nx+1) = 0d0
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine corr_var_ord3(f, ilim, corr_fp, corr_fm)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Returns correction arrays for field f:
! corr_fp         : f(i+1/2)- = f(i)   + corr_fp(i)
! corr_fm         : f(i+1/2)+ = f(i+1) + corr_fm(i+1)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: bnd_highord
implicit none
real(kind=dp), intent(inout) :: corr_fp (0:nx+1)
real(kind=dp), intent(inout) :: corr_fm (0:nx+1)
real(kind=dp), intent(in)    :: f (0:nx+1)
integer,       intent(in)    :: ilim
integer                      :: i
real(kind=dp)                :: ff, dfp, dfm, ap, am
real(kind=dp), parameter     :: kappa = 1.d0/3.d0
! -----------------------------------------------------------------------------
am = 1.d0 - kappa
ap = 1.d0 + kappa
! loop on cell interfaces
!$OMP parallel do  &
!$OMP private(dfp, dfm, ff) &
!$OMP shared(corr_fp, corr_fm, f, ilim, am, ap)
do i=1,nx
  dfp = (f(i+1)-f(i))
  dfm = (f(i)-f(i-1))
  if (abs(dfp).lt.eps_num) then
    corr_fp(i) = 0.d0
    corr_fm(i) = 0.d0
  else        
    ff = Limitor(dfm/dfp, ilim)
    corr_fp(i) = +0.25d0*ff*(am*dfm + ap*dfp)
    corr_fm(i) = -0.25d0*ff*(am*dfp + ap*dfm)
  endif
enddo
!$OMP end parallel do
! -----------------------------------------------------------------------------
! addapted reconstructions on boundary edges
! (displaced gradients evaluations)
if (bnd_highord.eq.1) then
  ! left bnd:
  dfp = (f(2)-f(1))
  dfm = (f(1)-f(0))
  if (abs(dfp).lt.eps_num) then
    corr_fp(0) = 0.d0
    corr_fm(0) = 0.d0
  else        
    ff = Limitor(dfm/dfp, ilim)
    corr_fp(0) = +0.25d0*ff*(am*dfm + ap*dfp)
    corr_fm(0) = -0.25d0*ff*(am*dfp + ap*dfm)
  endif
  ! right bnd:
  dfp = (f(nx+1)-f(nx))
  dfm = (f(nx)-f(nx-1))
  if (abs(dfp).lt.eps_num) then
    corr_fp(nx+1) = 0.d0
    corr_fm(nx+1) = 0.d0
  else        
    ff = Limitor(dfm/dfp, ilim)
    corr_fp(nx+1) = +0.25d0*ff*(am*dfm + ap*dfp)
    corr_fm(nx+1) = -0.25d0*ff*(am*dfp + ap*dfm)
  endif
else
  corr_fp(0) = 0.d0
  corr_fm(0) = 0.d0
  corr_fp(nx+1) = 0d0
  corr_fm(nx+1) = 0d0
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
function Limitor(r, ilimitor)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! see https://en.wikipedia.org/wiki/Flux_limiter
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
integer, intent(in)       :: ilimitor
real(kind=dp), intent(in) :: r
real(kind=dp)             :: Limitor
real(kind=dp)             :: theta
! -----------------------------------------------------------------------------

! Minmod
! ------
if (ilimitor.eq.1) then
  Limitor = max(0.d0, min(1d0, r))

! Van Albala
! ----------
elseif (ilimitor.eq.2) then
  Limitor = (2.d0*r)/(1.d0+r**2)

! Van Leer
! --------
elseif (ilimitor.eq.3) then
  Limitor = (r+dabs(r)) / (1.d0+dabs(r))

! Monotonized central MC
! ----------------------
elseif (ilimitor.eq.4) then
  Limitor = max(0.d0, min(2d0*r, min(0.5d0*(1+r),2.d0)))

! Superbee
! --------
elseif (ilimitor.eq.5) then
  Limitor = max(0.d0, max(min(1d0, 2*r), min(2d0, r)))

! Generalized minmod
! ------------------
elseif (ilimitor.eq.6) then
  ! theta in [1, 2]
  ! theta = 2 less dissipative, theta = 1 more dissipative
  theta = 2.d0
  Limitor = max(0.d0, min(min(theta*r, (1.d0+r)/2.d0), theta))

! Generalized minmod
! ------------------
elseif (ilimitor.eq.7) then
  ! theta in [1, 2]
  ! theta = 2 less dissipative, theta = 1 more dissipative
  theta = 1.d0
  Limitor = max(0.d0, min(min(theta*r, (1.d0+r)/2.d0), theta))

endif
! -----------------------------------------------------------------------------
end function
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                           WENO RECONSTRUCTIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine corr_weno_ord5(f, corr_fp, corr_fm)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Fifth order weno reconstructions
! corr_fp         : f(i+1/2)- = f(i)   + corr_fp(i)
! corr_fm         : f(i+1/2)+ = f(i+1) + corr_fm(i+1)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none            
real(kind=dp), intent(inout) :: corr_fp (0:nx+1)
real(kind=dp), intent(inout) :: corr_fm (0:nx+1)
real(kind=dp), intent(in)    :: f (0:nx+1)
real(kind=dp)                :: fg, fd
real(kind=dp)                :: s1, s2, s3, f1, f2, f3, f4, f5
real(kind=dp)                :: a1, a2, a3, w1, w2, w3, usum
real(kind=dp)                :: td = 13d0/12d0
real(kind=dp)                :: uq = 1d0/4d0
real(kind=dp)                :: us = 1d0/6d0
real(kind=dp)                :: ss = 7d0/6d0
real(kind=dp)                :: cs = 5d0/6d0
real(kind=dp)                :: ut = 1d0/3d0
real(kind=dp)                :: os = 11d0/6d0
real(kind=dp), parameter     :: eps_weight = 1.d-16
integer                      :: i
! -----------------------------------------------------------------------------
! corr_fp
! ~~~~~~~
!$OMP parallel do  &
!$OMP private(f1, f2, f3, f4, f5, usum, fg, &
!$OMP         s1, s2, s3, a1, a2, a3, w1, w2, w3) &
!$OMP shared(corr_fp, f, td, uq, us, ss, cs, ut, os, nx)
do i=0,nx+1
  ! constant extrapolation on ghost cells
  f1 = f(max(i-2,0))
  f2 = f(max(i-1,0))
  f3 = f(i)
  f4 = f(min(i+1,nx+1))
  f5 = f(min(i+2,nx+1))
  ! rec
  s1 = td*(f1-2d0*f2+f3)**2d0+uq*(f1-4d0*f2+3d0*f3)**2d0
  s2 = td*(f2-2d0*f3+f4)**2d0+uq*(f2-f4)**2d0
  s3 = td*(f3-2d0*f4+f5)**2d0+uq*(3d0*f3-4d0*f4+f5)**2d0
  a1 = 0.1d0/max(eps_weight, s1**2)
  a2 = 0.6d0/max(eps_weight, s2**2)
  a3 = 0.3d0/max(eps_weight, s3**2)
  usum = 1.d0/(a1+a2+a3)
  w1 = a1*usum
  w2 = a2*usum
  w3 = a3*usum
  fg = w1*(ut*f1-ss*f2+os*f3)+w2*(-us*f2+cs*f3+ut*f4)+w3*(ut*f3+cs*f4-us*f5)
  corr_fp(i) = fg-f(i)
enddo
!$OMP end parallel do
! -----------------------------------------------------------------------------
! corr_fm
! ~~~~~~~
!$OMP parallel do  &
!$OMP private(f1, f2, f3, f4, f5, usum, fd, &
!$OMP         s1, s2, s3, a1, a2, a3, w1, w2, w3) &
!$OMP shared(corr_fm, f, td, uq, us, ss, cs, ut, os, nx)
do i=-1,nx
  ! constant extrapolation on ghost cells
  f1 = f(min(i+3,nx+1))
  f2 = f(min(i+2,nx+1))
  f3 = f(min(i+1,nx+1))
  f4 = f(max(i  ,0))
  f5 = f(max(i-1,0))
  ! rec
  s1 = td*(f1-2d0*f2+f3)**2d0+uq*(f1-4d0*f2+3d0*f3)**2d0
  s2 = td*(f2-2d0*f3+f4)**2d0+uq*(f2-f4)**2d0
  s3 = td*(f3-2d0*f4+f5)**2d0+uq*(3d0*f3-4d0*f4+f5)**2d0
  a1 = 0.1d0/max(eps_weight, s1**2)
  a2 = 0.6d0/max(eps_weight, s2**2)
  a3 = 0.3d0/max(eps_weight, s3**2)
  usum = 1.d0/(a1+a2+a3)
  w1 = a1*usum
  w2 = a2*usum
  w3 = a3*usum
  fd = w1*(ut*f1-ss*f2+os*f3)+w2*(-us*f2+cs*f3+ut*f4)+w3*(ut*f3+cs*f4-us*f5)
  corr_fm(i+1) = fd-f(i+1)
enddo
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                 FINITE DIFFERENCES GRADIENT
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine finite_difference(f, h, df, method)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute f'(x) with the finite difference method
! ~~~~~~~~~~~~~
! if method = 1: first order 
! if method = 2: second order centered
! if method = 4: fourth order centered
! if method = 6: sixth order centered
! if method = 9: ninth order centered
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), dimension(:), intent(in) :: f
real(kind=dp), intent(in)               :: h
real(kind=dp), intent(out)              :: df
integer, intent(in)                     :: method
integer                                 :: i, n
! -----------------------------------------------------------------------------
! first order forward
! ~~~~~~~~~~~~~~~~~~~~
if (method.eq.1) then
  n = size(f) ! should be = 2
  i = n-1
  df = (f(i+1)-f(i))/h
! second order centered
! ~~~~~~~~~~~~~~~~~~~~~
else if (method.eq.2) then
  n = size(f) ! should be = 3
  i = n-1
  df = (f(i+1)-f(i-1))/(2*h)
! fourth order centered
! ~~~~~~~~~~~~~~~~~~~~~
else if (method.eq.4) then
  n = size(f) ! should be = 5
  i = n-2     !           = 3
  df = (f(i-2)-8d0*f(i-1)+8d0*f(i+1)-f(i+2))/(12d0*h)
! sixth order centered
! ~~~~~~~~~~~~~~~~~~~~
else if (method.eq.6) then
  n = size(f) ! should be = 7
  i = n-3     !           = 4
  df = (-f(i-3)+9d0*f(i-2)-45d0*f(i-1)+f(i+3)-9d0*f(i+2)+45d0*f(i+1))/(60d0*h)
! ninth order centered
! ~~~~~~~~~~~~~~~~~~~~
else if (method.eq.9) then
  n = size(f) ! should be = 9
  i = n-4     !           = 5
  df = (-2d0*f(i-4)+25d0*f(i-3)-150d0*f(i-2)+600d0*f(i-1) &
        +2d0*f(i+4)-25d0*f(i+3)+150d0*f(i+2)-600d0*f(i+1)) /(840d0*h)
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Unknown FINITE DIFFERENCE GRADIENT METHOD '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************