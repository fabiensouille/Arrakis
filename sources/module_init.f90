! *****************************************************************************
module module_init
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Discretisation and Initialisation module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp, eps
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                              GEOMETRY AND MESH
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine init_geometry
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Discretisation parameters initialisation
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Convention:
! 
!  0   xA    1        i-1        i        i+1        n   xB   n+1
!  o----|----o-- ... --o----|----o----|----o-- ... --o----|----o
!  x0  I0   x1       xi-1  Ii-1  xi  Ii   xi+1      xn   In   xn+1
!  <-------> <------->           <--------> dx[i]    <--------> 
!     dx[0]    dx[1]         <-------->                 dx[n]
!                              ci[i] : cell size
! 
! Ghost cells:
! 
!  o--|--o--    ...      --o--|--o
!  0 xA  1                nx  xB nx+1
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: nx, dx, ci, xa, xb, x, zb, zf, zf_default, minci, &
                       use_mesh_file, mesh_file, exner, variable_width, &
                       Lx, Lcst, dLx, bnd_extrapolation
use module_io, only: split_string
use module_alloc, only: allocation_geom
implicit none
integer                         :: i
integer                         :: io, nlines, header_lines
integer                         :: ncols, col_x, col_zb, col_zf, col_l
integer                         :: j
real(kind=dp)                   :: dxcst, gradl, gradd
real(kind=dp), allocatable      :: row_data(:)
character(len=255)              :: var_line
character(len=3), dimension(10) :: vars
! -----------------------------------------------------------------------------
! No mesh file (constant dx, Lx)
! ------------------------------
if (use_mesh_file.eqv..False.) then
  call allocation_geom
  dxcst = (xb-xa)/real(nx)
  dx(0) = dxcst
  do i=1,nx
    dx(i) = dxcst
    ci(i) = dxcst
  enddo
  minci = dxcst
  do i=0,nx+1
    ! initialize node abscissa
    x(i) = xa + (dx(0)/2d0) + (i-1)*dx(0)
    ! initialize bottom elevation
    zb(i) = 0.d0
    ! initialize hard bottom elevation if sediment transport
    if (exner.eqv..True.) then
      zf(i) = min(zf_default, zb(i))
    endif
    ! initialize channel width
    Lx(i) = Lcst
    if (variable_width.eqv..True.) then
      dLx(i)= 0.d0
    endif
  enddo
! -----------------------------------------------------------------------------
! Mesh file
! ---------
else
  write(*,*) ''
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ' Reading mesh file   '
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ''
  ! ~~~~~~~~~~~~~~~~~~~
  ! read mesh file size
  ! ~~~~~~~~~~~~~~~~~~~
  header_lines = 2
  nlines = 0
  open(unit=1, file=trim(mesh_file)//'.geo')
  do
    read(1, *, iostat=io)
    if (io/=0) exit
    nlines = nlines + 1
  enddo
  nx = nlines - header_lines
  write(*,*) 'number of cells = ', nx
  call allocation_geom
  ! ~~~~~~~~~~~
  ! read header
  ! ~~~~~~~~~~~
  close(unit=1)
  open(unit=1, file=trim(mesh_file)//'.geo')
  read(1, *)
  read(1, '(A)') var_line
  ! parse header: identify column index for each known variable
  call split_string(var_line, vars)
  ncols = 0
  col_x  = 0
  col_zb = 0
  col_zf = 0
  col_l  = 0
  do j=1,10
    if (trim(vars(j)).eq.'') exit
    ncols = ncols + 1
    if (trim(vars(j)).eq.'X')  col_x  = j
    if (trim(vars(j)).eq.'ZB') col_zb = j
    if (trim(vars(j)).eq.'ZF') col_zf = j
    if (trim(vars(j)).eq.'L')  col_l  = j
  enddo
  ! check that mandatory columns are present
  if (col_x.eq.0) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Error reading the mesh file: missing column X'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  if (col_zb.eq.0) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Error reading the mesh file: missing column ZB'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  if ((col_zf.eq.0).and.(exner.eqv..True.)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Error reading the mesh file: missing column ZF'
    write(*,*) ' ZF is required when ACTIVATE BEDLOAD = yes'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  if ((col_l.eq.0).and.(variable_width.eqv..True.)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' Error reading the mesh file: missing column L'
    write(*,*) ' L is required when VARIABLE CHANNEL WIDTH = yes'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! read mesh file content
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! All columns are read generically into row_data.
  ! Only the required ones are assigned; extra columns are silently ignored.
  allocate(row_data(ncols))
  do i=1,nx
    read(1, *) row_data(1:ncols)
    x(i)  = row_data(col_x)
    zb(i) = row_data(col_zb)
    if (exner.eqv..True.) then
      zf(i) = row_data(col_zf)
      ! sanity check
      if(zf(i).gt.zb(i))then
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*) ' ERROR : ZF SHOULD BE <= ZB '
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        stop
      endif
    endif
    if (variable_width.eqv..True.) then
      Lx(i) = row_data(col_l)
    else
      Lx(i) = Lcst
    endif
  enddo
  deallocate(row_data)
  close(unit = 1)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! define discretization parameters
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do i=1,nx-1
    dx(i) = x(i+1)-x(i)
  enddo
  dx(0)  = dx(1)
  dx(nx) = dx(nx-1)
  x(0)    = x(1) -dx(0)
  x(nx+1) = x(nx)+dx(nx)
  do i=1,nx
    ci(i) = 0.5d0*(dx(i-1)+dx(i))
  enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! define values on ghost cells
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (bnd_extrapolation.eq.0) then
    ! set equal to value on neighboring cell
    zb(0)    = zb(1)
    zb(nx+1) = zb(nx)
    Lx(0)    = Lx(1)
    Lx(nx+1) = Lx(nx)
    if (exner.eqv..True.) then
      zf(0)    = zf(1)
      zf(nx+1) = zf(nx)
    endif
  else
    ! preserve longitudinal gradient
    ! -> topo
    gradl = (zb(2)-zb(1))/dx(1)
    gradd = (zb(nx)-zb(nx-1))/dx(nx-1)
    zb(0)    = zb(1) - gradl*dx(0)
    zb(nx+1) = zb(nx)+ gradd*dx(nx)
    ! -> width
    gradl = (Lx(2)-Lx(1))/dx(1)
    gradd = (Lx(nx)-Lx(nx-1))/dx(nx-1)
    Lx(0)    = Lx(1) - gradl*dx(0)
    Lx(nx+1) = Lx(nx)+ gradd*dx(nx)
    if (exner.eqv..True.) then
      gradl = (zf(2)-zf(1))/dx(1)
      gradd = (zf(nx)-zf(nx-1))/dx(nx-1)
      zf(0)    = zf(1) - gradl*dx(0)
      zf(nx+1) = zf(nx)+ gradd*dx(nx)
    endif
  endif
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  ! define xA and xB from x
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  xa = minval(x)
  xb = maxval(x)
  ! ~~~~~~~~~~~~~
  ! Initialize L'
  ! ~~~~~~~~~~~~~
  ! L' is defined at interfaces between cells dL(i) <- L'(i+1/2)
  if (variable_width.eqv..True.) then
    do i=0,nx
      dLx(i) = (Lx(i+1)-Lx(i))/dx(i)
    enddo
  end if
end if
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                            INITIAL CONDITIONS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine init_variables()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initialisation 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: use_ini_file, ini_file, exner, ordre, pi, nx, eps_num, &
                       x, dx, zf, zf_default, zb, Lx, Rh, h, u, q, t, ht, &
                       fric, frict_coef, qsed, taub, taus, aux1, aux2, ntrac, &
                       corr_bp, corr_bm, corr_hp, corr_hm, &
                       corr_up, corr_um, corr_tp, corr_tm, &
                       const_h, const_e, const_q,  const_t, &
                       bedload, g, rhos, rho, d50, tau_adim, &
                       type_init, bnd_extrapolation, ub, qa, qb, qsa, qsb, &
                       d16, dm, d84, eps_b_damping, eps_b, suspension, &
                       ws_susp, difu, nu_u, tau_c, trac_names
use module_io, only: split_string, count_vars
use module_forward, only: update_Rh
use module_bedload, only: critical_shields
use module_suspension, only: settling_velocity
implicit none
real(kind=dp), pointer           :: x_tmp (:)
real(kind=dp), pointer           :: f_tmp (:, :)
real(kind=dp), pointer           :: f_int (:, :)
real(kind=dp)                    :: initime, a0, a1
real(kind=dp)                    :: minb, minb_tmp
integer                          :: i, j, k, l, npoin
integer                          :: ndata, io
character(len=255)               :: line
character(len=255)               :: var_line
character(len=3), dimension(100) :: vars
character(len=4)                 :: char_trac
! -----------------------------------------------------------------------------
! initialization of arrays
! ------------------------
do i=0,nx+1
  ! CONSTANT ELEVATION and CONSTANT FLOWRATE
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (type_init==1) then
    ! main arrays
    h(i)  = max(0.d0, const_e - zb(i))
    q(i)  = const_q/Lx(i) ! volumic to lineic discharge
    u(i)  = q(i)/max(h(i), eps)
    ! tracers
    if (ntrac.gt.0) then
      do k=1,ntrac
        t(i,k) = const_t(k)
      enddo
    endif
  ! CONSTANT WATER DEPTH and CONSTANT FLOWRATE
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (type_init==2) then
    ! main arrays
    h(i)  = max(0.d0, const_h)
    q(i)  = const_q/Lx(i) ! volumic to lineic discharge
    u(i)  = q(i)/max(h(i), eps)
    ! tracers
    if (ntrac.gt.0) then
      do k=1,ntrac
        t(i,k) = const_t(k)
      enddo
    endif
  ! ALL NULL
  ! ~~~~~~~~
  else 
    ! main arrays
    h(i)  = 0.d0
    q(i)  = 0.d0
    u(i)  = 0.d0
    ! tracers
    if (ntrac.gt.0) then
      do k=1,ntrac
        t(i,k) = 0.d0
      enddo
    endif
  endif 
  ! initialise other arrays
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  ! initialize friction
  fric(i) = frict_coef
  ! initialize veolicty diffusion
  difu(i) = nu_u
  ! initialize bedload
  if (bedload.eqv..True.) then
    qsed(i) = 0d0
    taub(i) = 0d0
    taus(i) = 0d0
  endif
  ! hydraulic radius
  ! if variable_width init hydraulic radius
  call update_Rh(h(i), Lx(i), Rh(i))
  ! initialise auxiliary arrays
  aux1(i) = 0d0
  aux2(i) = 0d0
  ! initialise second order reconstructions
  if (ordre.ge.2) then
    corr_bp(i) = 0.d0
    corr_bm(i) = 0.d0
    corr_hp(i) = 0.d0
    corr_hm(i) = 0.d0
    corr_up(i) = 0.d0
    corr_um(i) = 0.d0
    if (ntrac.gt.0) then
      do k=1,ntrac
        corr_tp(i, k) = 0.0d0
        corr_tm(i, k) = 0.0d0
      enddo
    endif
  endif
enddo
! -----------------------------------------------------------------------------
! Continuation files (overwrite constant init)
! --------------------------------------------
if(use_ini_file.eqv..True.) then
  ! initialize variables of ini files to blank
  do i=1,size(vars)
    vars(i) = ''
  enddo
  ! Read continuation file
  ! ~~~~~~~~~~~~~~~~~~~~~~
  write(*,*) ''
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ' Reading ini file  '
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ''
  open(unit=2, file=trim(ini_file)//'.ini')
  read(2, *) initime, npoin
  read(2, '(A)') var_line
  call split_string(var_line, vars) ! read var names 
  call count_vars(vars, ndata) ! detect the number of var
  call check_ini_inputs(vars) ! checks var in file
  ! allocate tmp tables
  allocate(x_tmp(1:npoin))
  allocate(f_tmp(1:npoin, 1:ndata))
  allocate(f_int(0:nx+1,  1:ndata))
  ! read file content
  do i=1,npoin
    read(2, '(A)', iostat=io) line
    read(line, *, iostat=io) f_tmp(i,1:ndata)
    x_tmp(i) = f_tmp(i,1)
  enddo
  close(unit = 2)
  ! interpolate data on mesh (f_tmp->f_int)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~
  do l=1,ndata
    ! define values on internal cells
    do j=1,nx ! loop on internal cells
      if (x(j).le.x_tmp(1)) then
        f_int(j,l) = f_tmp(1,l)
      elseif (x(j).ge.x_tmp(npoin)) then
        f_int(j,l) = f_tmp(npoin,l)
      else
        do i=1,npoin-1 ! loop on temporary ini array
          if ((x(j).ge.x_tmp(i)).and.(x(j).le.x_tmp(i+1))) then
            a1 = (x(j)-x_tmp(i))/(x_tmp(i+1)-x_tmp(i))
            a0 = 1.d0 - a1
            f_int(j,l) = a0*f_tmp(i,l) + a1*f_tmp(i+1,l)
          endif
        enddo
      endif
    enddo
    ! define values on ghost cells
    if (bnd_extrapolation.eq.0) then
      f_int(0,l)    = f_int(1,l)
      f_int(nx+1,l) = f_int(nx,l)
    else
      f_int(0,l)    = f_int(1,l) -(f_int(2,l) -f_int(1,l))*(dx(0)/dx(1))
      f_int(nx+1,l) = f_int(nx,l)+(f_int(nx,l)-f_int(nx-1,l))*(dx(nx)/dx(nx-1))
    endif
  enddo
  ! assignment of variables depending on labels
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do l=1,ndata
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (trim(vars(l)).eq.'ZB') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zb(0:nx+1) = f_int(0:nx+1,l)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(vars(l)).eq.'ZS') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do j=0,nx+1
        h(j) = max(0.d0, f_int(j,l) - zb(j))
      enddo
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(vars(l)).eq.'H') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      h(0:nx+1) = f_int(0:nx+1,l)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(vars(l)).eq.'U') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      u(0:nx+1) = f_int(0:nx+1,l)
      do j=0,nx+1
        q(j) = h(j)*u(j)
      enddo
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(vars(l)).eq.'Q') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      q(0:nx+1) = f_int(0:nx+1,l)
      do j=0,nx+1
        if (h(j)>eps) then
          u(j) = q(j)/max(0.d0, h(j))
        else
          u(j) = 0d0
        endif
      enddo
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(vars(l)).eq.'FC') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      fric(0:nx+1) = f_int(0:nx+1,l)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(vars(l)).eq.'ZF') then
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zf(0:nx+1) = f_int(0:nx+1,l)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    endif
    ! ~~~~~~~
    ! tracers
    ! ~~~~~~~
    if (ntrac.gt.0) then
      do k=1,ntrac
        char_trac = trac_names(k)
        if ((trim(vars(l)).eq."T"//char_trac(4:4)).or. &
&            (trim(vars(l)).eq."T"//char_trac(3:4))) then
          t(0:nx+1,k) = f_int(0:nx+1, l)
        endif
      enddo
    endif ! end if trac
  enddo ! end assignment
  ! 
  ! deallocate tmp tables
  deallocate(x_tmp)
  deallocate(f_tmp)
  deallocate(f_int)
endif ! end if use_ini_file
! -----------------------------------------------------------------------------
do i=0,nx+1
  if (ntrac.gt.0) then
    do k=1,ntrac
      ht(i,k) = t(i,k)*h(i)
    enddo
  endif
enddo
! -----------------------------------------------------------------------------
! initialise other parameters (which depends on arrays)
! -----------------------------------------------------
! volumic to lineic discharge
qa = qa/Lx(0)
qb = qb/Lx(nx+1)
! flux sign convention for u,q in right BC
ub = -ub
qb = -qb
! Initialize bedload parameters
if (bedload.eqv..True.) then
  ! volumic to lineic discharge
  qsa = qsa/Lx(0)
  qsb = qsb/Lx(nx+1) 
  ! flux sign convention for qs in right BC
  qsb = -qsb
  ! size distribution (0. means auto-compute from d50)
  if (d16.eq.0d0) d16 = 0.5d0*d50 ! Only for Lefort bedload law
  if (dm .eq.0d0) dm  = 1.1d0*d50 ! Only for Lefort and Recking bedload laws
  if (d84.eq.0d0) d84 = 2.1d0*d50 ! Only for Lefort and Recking bedload laws
endif
if (suspension.eqv..True.) then
  ! suspension settling velocity
  if (ws_susp.eq.0.d0) then
    ws_susp = settling_velocity(d50)
  endif
endif
if (exner.eqv..True.) then
  ! init adim factor : tau_adim = taus/taub
  tau_adim = 1d0/((rhos-rho)*g*d50)
  ! init critical shear stress
  if (tau_c.eq.0d0) then
    tau_c = critical_shields(d50)
  endif
  ! erodible bed damping region
  if(eps_b_damping.eq.1)then
    ! compute min(zb-zf)
    minb = abs(zb(0)-zf(0))
    do i=1,nx+1
      minb_tmp = abs(zb(i)-zf(i))
      if (minb_tmp.lt.minb) then
        minb = minb_tmp
      endif
    enddo
    ! if damping detph set to 0.: init to 10% of min(zb-zf) 
    if (eps_b .eq. 0d0) then
      eps_b = max(eps_num, 0.1*minb)
      if (eps_b.eq.0.d0) then
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*) ' WARNING : "NON ERODIBLE BED DAMPING REGION DEPTH" = 0. '
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      endif
    endif
  endif 
endif ! end of exner init
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine check_ini_inputs(vars)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! checks variables in continuation file
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: use_mesh_file, exner
implicit none
character(len=3), dimension(:), intent(in) :: vars
logical    :: miss_h, miss_u, miss_q
logical    :: miss_zb, miss_zf, miss_zs
integer    :: k
! -----------------------------------------------------------------------------
miss_zs = .True.
miss_h = .True.
miss_u = .True.
miss_q = .True.
miss_zb = .True.
miss_zf = .True.
! -----------------------------------------------------------------------------
do k=1,size(vars)
  if (trim(vars(k)).eq.'ZB') then
    miss_zb = .False.
  else if (trim(vars(k)).eq.'ZS') then
    miss_zs = .False.
  else if (trim(vars(k)).eq.'H') then
    miss_h = .False.
  else if (trim(vars(k)).eq.'U') then
    miss_u = .False.
  else if (trim(vars(k)).eq.'Q') then
    miss_q = .False.
  else if (trim(vars(k)).eq.'ZF') then
    miss_zf = .False.
  endif
enddo
if ((miss_zs.eqv..True.).and.(miss_h.eqv..True.)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Error reading hydro ini file, must contain '
  write(*,*) ' either water depth H or free surface elevation ZS '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
if ((miss_u.eqv..True.).and.(miss_q.eqv..True.)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Error reading hydro ini file, must contain U or Q'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
if ((miss_zb.eqv..True.).and.(use_mesh_file.eqv..False.)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Hydro ini file do not contain ZB and no mesh file given, '
  write(*,*) ' ZB will be initiallized at 0. m'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
if ((miss_zf.eqv..True.).and.(use_mesh_file.eqv..False.) &
                        .and.(exner.eqv..True.)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' Hydro ini file do not contain ZF and no mesh file given, '
  write(*,*) ' ZF will be initiallized at -0.1 m'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************