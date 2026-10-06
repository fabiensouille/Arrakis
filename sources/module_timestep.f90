! *****************************************************************************
module module_timestep
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Allocation of arrays
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, nx, g, eps, eps_num
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
subroutine init_time_step(h, u)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! time step initialization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: var_dt, dt, dt_const, max_dt, cfl, ci, minci, ordre, &
                       mlamb, ordre_trac, ntrac, difu, &
                       diffusion_u, nu_u, diffusion_z, nu_z,  diffusion_t, nu_t
implicit none
real(kind=dp), intent(in) :: h (0:nx+1)
real(kind=dp), intent(in) :: u (0:nx+1)
real(kind=dp)             :: advc, maxdifu
integer                   :: i, k
! -----------------------------------------------------------------------------
! variable time step
! ~~~~~~~~~~~~~~~~~~
if (var_dt.eqv..True.) then
  ! get min cell size
  minci = ci(1)
  do i=1,nx
    minci = min(ci(i), minci)
  enddo
  ! coef (half cfl condition)
  if((ordre.ge.2).or.(ordre_trac.ge.2)) then
    advc = 4d0
  else
    advc = 2d0
  endif
  ! advection
  mlamb = abs(maxval(u)) + sqrt(g*maxval(h)) ! initial guess
  mlamb = max(eps_num, mlamb)
  dt = cfl*minci/(advc*mlamb)
  ! diffusion of velocity
  if ((diffusion_u.eqv..True.).and.(nu_u.gt.0)) then
    maxdifu = maxval(difu)
    dt = min(dt, cfl*(minci**2./(4.*maxdifu)))
  endif
  ! diffusion of tracers
  if ((ntrac.gt.0).and.(diffusion_t.eqv..True.)) then
    do k=1,ntrac
      if (nu_t(k).gt.0) then
        dt = min(dt, cfl*(minci**2./(4.*nu_t(k))))
      endif
    enddo
  endif
  ! diffusion of zb
  if ((diffusion_z.eqv..True.).and.(nu_z.gt.0)) then
    dt = min(dt, cfl*(minci**2./(4.*nu_z)))
  endif
! constant time step
! ~~~~~~~~~~~~~~~~~~
else
  dt = dt_const
endif
! maximum time step allowed
! ~~~~~~~~~~~~~~~~~~~~~~~~~
dt = min(dt, max_dt)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine adjust_time_step
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Adjusting time step to mlamb (if variable time step)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: var_dt, dt, cfl, ordre, mlamb, minci, max_dt, ordre_trac
implicit none
real(kind=dp) :: advc
! -----------------------------------------------------------------------------
if (var_dt.eqv..True.) then
  ! coef (half cfl condition)
  if((ordre.ge.2).or.(ordre_trac.ge.2)) then
    advc = 4d0
  else
    advc = 2d0
  endif
  ! advection
  ! mlamb := maximum wave speed of the Approximate Riemann solver
  ! computed in compute_flux_sw or compute_flux_swe
  mlamb = max(eps_num, mlamb)
  dt = min(dt, cfl*minci/(advc*mlamb))
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************