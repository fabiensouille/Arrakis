! *****************************************************************************
module module_steady
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for Shallow Water equations fluxes computation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g, eps, eps_num
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                                STEADY SOLVER
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine compute_water_line(i, hg, qg, zbg, hd, qd, zbd, dx, conv_flag)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute water surface elevation at xg knowing hd and qg, qd.
! Solve non-linear equation f(zg) = 0
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         dx
!      o-------o
!   hg,qg  <- hd,qd
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,     only: Rh, Lx, fric, std_precision, std_max_iter, &
                           write_listing, std_newton_debug
use module_friction, only: friction_slope
implicit none
integer, intent(in)        :: i
real(kind=dp), intent(in)  :: dx, qg, zbg, hd, qd, zbd
real(kind=dp), intent(out) :: hg
logical,       intent(out) :: conv_flag
real(kind=dp)              :: Lg, Ld, fcg, fcd, zd
real(kind=dp)              :: z0, z, dz, J
real(kind=dp)              :: fz1, fz2, dfz, residual
integer                    :: iter
! -----------------------------------------------------------------------------
Lg = Lx(i)    ! width
Ld = Lx(i+1) 
fcg = fric(i) ! friction coeff
fcd = fric(i+1)
zd = hd + zbd ! free surface elevation

! Initial guess
! *************
J = friction_slope(hd, Rh(i+1), qd/max(hd, eps), fric(i+1))
z0 = hd + zbd + dx*J
z0 = z0 - f(z0, qg, zbg, Lg, fcg, zd, qd, zbd, Ld, fcd, dx)

! Non-linear solver for z = f(z, qg, zbg, hd, qd, zbd)  (Newton-Raphson)
! ****************************************************
iter = 0
residual = 1d0
dz = 1d-3 ! for FD derivative
conv_flag = .True.
do while (iter.le.std_max_iter .and. residual.ge.std_precision)
  fz1 = f(z0,    qg, zbg, Lg, fcg, zd, qd, zbd, Ld, fcd, dx)
  fz2 = f(z0+dz, qg, zbg, Lg, fcg, zd, qd, zbd, Ld, fcd, dx)
  dfz = (fz2 - fz1)/dz
  z = z0 - fz1/dfz
  residual = abs(z - z0)
  z0 = z
  iter = iter + 1
  if ((iter.ge.std_max_iter).and.&
      (residual.ge.std_precision).and.&
      (write_listing.or.std_newton_debug)) then
    conv_flag = .False.
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : STEADY SOLVER MAX ITER REACHED '
    write(*,*) ' -> CELL ', i
    write(*,*) ' -> RESIDUAL = ', residual
    write(*,*) ' -> FREE SURFACE ELEVATION [i+1] = ', zd
    write(*,*) ' -> FREE SURFACE ELEVATION [i]   = ', z
    write(*,*) ' -> WATER DEPTH [i] = ', max(z - zbd, 0d0)
    write(*,*) ' -> CONSIDER SWITCHING TO UNSTEADY SOLVER '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
end do

! Update h
! ********
hg = max(z - zbg, 0d0)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************




! *****************************************************************************
function f(zg, qg, zbg, Lg, fcg, &
           zd, qd, zbd, Ld, fcd, dx)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Implicit function for free surface update i.e. f(zg) = 0
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_friction, only: friction_slope, hydraulic_radius
implicit none
real(kind=dp), intent(in) :: zg, qg, zbg, zd, qd, zbd
real(kind=dp), intent(in) :: Lg, fcg, Ld, fcd, dx
real(kind=dp) :: f
real(kind=dp) :: hgc, hdc, ug, ud, Rhg, Rhd, Jg, Jd
real(kind=dp) :: dxu, Jdx
integer, parameter :: frict_scheme = 1
! -----------------------------------------------------------------------------
! compute water depths and velocities
hgc = max(zg - zbg, eps)
hdc = max(zd - zbd, eps)
ug = qg/hgc
ud = qd/hdc
! compute hydraulic radius
Rhg = hydraulic_radius(hgc, Lg)
Rhd = hydraulic_radius(hdc, Ld)
! compute friction slopes
Jg = friction_slope(hgc, Rhg, ug, fcg)
Jd = friction_slope(hdc, Rhd, ud, fcd)
! momentum term
dxu = (ud**2d0 - ug**2d0)/(2d0*g)
! friction term
if (Jg*Jd.gt.eps_num) then
  if (frict_scheme.eq.1) then
    Jdx = 2d0*dx*Jg*Jd/(Jg + Jd)
  else if (frict_scheme.eq.2) then
    Jdx = dx*Jd
  else if (frict_scheme.eq.3) then
    Jdx = dx*(Jd+Jg)/2d0
  endif
else 
  Jdx = 0d0
endif
! compute function
f = zg-zd - dxu - Jdx
! -----------------------------------------------------------------------------
end function f
! *****************************************************************************




! *****************************************************************************
subroutine compute_water_line_linear(i, hg, zbg, hd, qd, zbd, dx, conv_flag)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute water surface elevation at xg knowing hd and qg, qd.
! Solve linearized steady
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         dx
!      o-------o
!   hg,qg  <- hd,qd
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data,     only: Lx, fric
use module_friction, only: friction_slope, hydraulic_radius
implicit none
integer, intent(in)        :: i
real(kind=dp), intent(in)  :: dx, zbg, hd, qd, zbd
real(kind=dp), intent(out) :: hg
logical,       intent(out) :: conv_flag
real(kind=dp)              :: Lg, Ld, fcd
real(kind=dp)              :: Jg, Jd, Frd, ud, Rhd, Fr2
! -----------------------------------------------------------------------------
Lg = Lx(i)    ! width
Ld = Lx(i+1) 
fcd = fric(i+1) ! friction coeff
! Linearized algorithm
! ********************
Jg = (zbg - zbd)/dx
ud = qd/max(hd, eps)
Frd = ud/sqrt(g*max(hd, eps))
Fr2 = Frd*Frd
Rhd = hydraulic_radius(hd, Ld)
Jd = friction_slope(hd, Rhd, ud, fcd)
! Linearized water depth equation:
if (abs(1d0 - Fr2) > eps) then
  hg = hd + (Jd - Jg)*dx/(1d0 - Fr2)
else
  hg = hd
endif
! Ensure positive water depth
hg = max(hg, eps)
! Convergence flag (always true for this direct method)
conv_flag = .True.
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************




! *****************************************************************************
subroutine compute_tracers(h, u, q, t, ht)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: nx, ntrac, dt, dx, ci, Lx, nu_t, diffusion_t, &
                       write_listing
use module_linalg, only: invert_matrix, matrix_product
implicit none
real(kind=dp), intent(in)       :: h  (0:nx+1)
real(kind=dp), intent(in)       :: u  (0:nx+1)
real(kind=dp), intent(in)       :: q  (0:nx+1)
real(kind=dp), intent(inout)    :: t  (0:nx+1,ntrac)
real(kind=dp), intent(inout)    :: ht (0:nx+1,ntrac)
integer                         :: i, j, k
real(kind=dp)                   :: fact, fact0, max_dt
real(kind=dp)                   :: qip12, qim12, hip12, him12
real(kind=dp)                   :: ap, am, bp, bm, cc, cd, ce, cf, sum, bj
real(kind=dp), dimension(nx,nx) :: A, AINV
integer, parameter              :: tracer_time_scheme = 2
logical                         :: ok_flag
! -----------------------------------------------------------------------------
! explicit scheme
! ***************
if (tracer_time_scheme.eq.1) then
  do i=1,nx
    ! check CFL
    max_dt = q(i)/max(h(i), eps)*ci(i)
    if ((dt.gt.max_dt).and.(write_listing)) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' WARNING : CFL VIOLATION, CHECK YOUR TIME STEP < ', max_dt
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    endif
    fact0 = dt/ci(i)
    fact = fact0/Lx(i)
    do k=1,ntrac
      if (q(i).ge.0d0) then
        ht(i,k) = ht(i,k) - fact*(q(i)*Lx(i)*t(i,k) - q(i-1)*Lx(i-1)*t(i-1,k))
      else
        ht(i,k) = ht(i,k) - fact*(q(i+1)*Lx(i+1)*t(i+1,k) - q(i)*Lx(i)*t(i,k))
      endif
      ! diffusion term
      if (diffusion_t.eqv..True.) then
        hip12 = 0.5d0*(h(i)+h(i+1))
        him12 = 0.5d0*(h(i)+h(i-1))
        ht(i,k) = ht(i,k) + fact0*nu_t(k)*( hip12*(t(i+1,k)-t(i,k))/dx(i) &
                                          - him12*(t(i,k)-t(i-1,k))/dx(i-1))
      endif
      ! update t
      if (h(i)>eps) then
        t(i,k) = ht(i,k)/h(i)
      else
        t(i,k) = ht(i,k)/eps
      endif
    enddo
  enddo
! implicit scheme
! ***************
else if (tracer_time_scheme.eq.2) then
  do k=1,ntrac
    ! initialize matrices
    do i=1,nx
      do j=1,nx
        A(i,j) = 0d0
        AINV(i,j) = 0d0
      enddo
    enddo 
    ! assemble matrices
    do i=1,nx
      fact0 = dt/ci(i)
      fact = fact0/Lx(i)
      ! advection terms
      qip12 = 0.5d0*(q(i)+q(i+1))
      qim12 = 0.5d0*(q(i)+q(i-1))
      if (qip12.ge.0d0) then
        ap = -fact*Lx(i)*u(i)
      else
        ap = 0d0
      endif
      if (qip12.le.0d0) then
        am = -fact*Lx(i+1)*u(i+1)
      else
        am = 0d0
      endif
      if (qim12.ge.0d0) then
        bp = fact*Lx(i-1)*u(i-1)
      else
        bp = 0d0
      endif
      if (qim12.le.0d0) then
        bm = fact*Lx(i)*u(i)
      else
        bm = 0d0
      endif
      ! diffusion terms
      if (diffusion_t.eqv..True.) then
        hip12 = 0.5d0*(h(i)+h(i+1))
        him12 = 0.5d0*(h(i)+h(i-1))
        cc = fact0*nu_t(k)*hip12/(dx(i)*h(i+1))
        cd =-fact0*nu_t(k)*hip12/(dx(i)*h(i))
        ce =-fact0*nu_t(k)*him12/(dx(i-1)*h(i))
        cf = fact0*nu_t(k)*him12/(dx(i-1)*h(i-1))
      else
        cc = 0d0
        cd = 0d0
        ce = 0d0
        cf = 0d0
      endif
      ! assembly
      A(i,i) = 1 - (ap + bm + cd + ce)
      if(i.eq.1) then
        A(i,i+1) = -(am + cc)
      else if(i.eq.nx) then
        A(i,i-1) = -(bp + cf)
      else
        A(i,i-1) = -(bp + cf)
        A(i,i+1) = -(am + cc)
      endif 
    enddo
    ! invert A, to solve A*Yn+1 = Yn
    call invert_matrix(nx, A, AINV, ok_flag)
    ! update ht, t
    do i=1,nx
      ! update ht
      sum = 0d0
      do j=1,nx
        if(j.eq.1) then
          bj = ht(j,k) + (bp + cf)*ht(0,k) ! left bc
        elseif(j.eq.nx) then
          bj = ht(j,k) + (am + cc)*ht(nx+1,k) ! right bc
        else
          bj = ht(j,k)
        endif
        sum = sum + AINV(i,j)*bj
      enddo
      ht(i,k) = sum
      ! update t
      if (h(i)>eps) then
        t(i,k) = ht(i,k)/h(i)
      else
        t(i,k) = ht(i,k)/eps
      endif
    enddo
  enddo
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************