! *****************************************************************************
module module_data
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Global parameters and variables module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use ISO_FORTRAN_ENV
implicit none
integer, parameter                      :: dp = REAL64
logical                                 :: use_omp

! -----------------------------------------------------------------------------
! main arrays
! ~~~~~~~~~~~
! x               : abscissa of mesh cells
! zb              : bottom elevation
! h               : water depth
! u               : water depth averaged velocity
! q               : water depth averaged discharge (q=h*u)
! fluh_l/r        : h fluxes on left and right side of interfaces
! fluq_l/r        : q fluxes on left and right side of interfaces
! fluh_src_l/r    : source fluxes for h
! fluq_src_l/r    : source fluxes for q
! fric            : friction coefficient
! aux1            : auxiliary array for debug
! aux2            : auxiliary array for debug
! -----------------------------------------------------------------------------
real(kind=dp), pointer                  :: x         (:)
real(kind=dp), pointer                  :: zb        (:)
real(kind=dp), pointer                  :: h         (:)
real(kind=dp), pointer                  :: u         (:)
real(kind=dp), pointer                  :: q         (:)
real(kind=dp), pointer                  :: fluh_l    (:)
real(kind=dp), pointer                  :: fluh_r    (:)
real(kind=dp), pointer                  :: fluq_l    (:)
real(kind=dp), pointer                  :: fluq_r    (:)
real(kind=dp), pointer                  :: fluh_src_l(:)
real(kind=dp), pointer                  :: fluh_src_r(:)
real(kind=dp), pointer                  :: fluq_src_l(:)
real(kind=dp), pointer                  :: fluq_src_r(:)
real(kind=dp), pointer                  :: fric      (:)
real(kind=dp), pointer                  :: difu      (:)
real(kind=dp), pointer                  :: aux1      (:)
real(kind=dp), pointer                  :: aux2      (:)

! -----------------------------------------------------------------------------
! spatial discretisation 
! ~~~~~~~~~~~~~~~~~~~~~~
! nx              : number of finite volume cells along x
! xa              : abscissa of left boundary interface
! xb              : abscissa of right boundary interface
! dx              : space step along x (dx(i) = x(i+1)-x(i))
! ci              : finite volume cell size
! minci           : size of the smallest cell (for dt computation)
! variable_width  : if the channel has variable width
! Lx              : channel width
! Rh              : channel hydraulic radius
! Rh_approx       : hydraulic radius approximation for constant width channels:
!                   if == 0: Rh = h*B/(B + 2*h)
!                   if != 0: Rh = h
! -----------------------------------------------------------------------------
integer                                 :: nx
real(kind=dp)                           :: xa
real(kind=dp)                           :: xb
real(kind=dp)                           :: Lcst
real(kind=dp), pointer                  :: dx (:)
real(kind=dp), pointer                  :: ci (:)
real(kind=dp)                           :: minci
logical                                 :: variable_width
real(kind=dp), pointer                  :: Lx (:)
real(kind=dp), pointer                  :: dLx(:)
real(kind=dp), pointer                  :: Rh (:)
integer                                 :: Rh_approx

! -----------------------------------------------------------------------------
! temporal discretisation 
! ~~~~~~~~~~~~~~~~~~~~~~~
! var_dt          : use variable time step (cfl condition)
! cfl             : cfl number
! dt_const        : time step specified by user in input file
! dt              : time step of the computation
! time            : current time
! time_max        : final time of the simulation
! mlamb           : maximum wave speed for time step limitation
! tstart          : initial time of the simulation
! niter           : current iteration
! nitermax        : maximum number of iteration
! when            : indicates initialization (0), temporal iterations (1) or finalization (2)
! eps_evol_h      : criteria on h evolution to trigger update [steady solver]
! eps_evol_h      : criteria on q evolution to trigger update [steady solver]
!                   Rmk: solver will not update if bnd condition do not change
! -----------------------------------------------------------------------------
logical                                 :: var_dt
real(kind=dp)                           :: cfl
real(kind=dp)                           :: dt_const
real(kind=dp)                           :: dt
real(kind=dp), parameter                :: max_dt = 3600.
real(kind=dp)                           :: time
real(kind=dp)                           :: time_max
real(kind=dp)                           :: mlamb
real(kind=dp), parameter                :: tstart = 0d0
integer                                 :: niter
integer                                 :: nitermax
integer                                 :: when
real(kind=dp), parameter                :: eps_evol_h = 1d-6
real(kind=dp), parameter                :: eps_evol_q = 1d-6

! -----------------------------------------------------------------------------
! physical parameters
! ~~~~~~~~~~~~~~~~~~~
! pi              : Pi = 3.14...
! g               : gravity acceleration
! kappa           : Von Karman constant
! rho             : density of water
! nu_w            : water molecular diffusivity
! frict_coef      : bottom friction coefficient
! frict_law       : bottom friction law
! diffusion_u     : activate diffusion of velocity
! nu_u            : water diffusivity
! -----------------------------------------------------------------------------
real(kind=dp), parameter                :: pi = 4.d0*datan(1.d0)
real(kind=dp), parameter                :: g = 9.80665
real(kind=dp), parameter                :: kappa = 0.41d0
real(kind=dp), parameter                :: rho = 1000.d0
real(kind=dp), parameter                :: nu_w = 1.d-6
real(kind=dp)                           :: frict_coef
integer                                 :: frict_law
logical                                 :: diffusion_u
real(kind=dp)                           :: nu_u

! -----------------------------------------------------------------------------
! io parameters
! ~~~~~~~~~~~~~
! ndisplay        : frequency of graphic printouts
! nrecord         : frequency of result recording
! irec            : index of written output file
! ficdon          : input file name
! res_file        : result file name
! pbr_file        : probe file name
! mesh_file       : mesh file name
! use_mesh_file   : activate mesh file
! nvars_outputs   : number of variables to write in result file
! vars_outputs    : variable names to write in result file
! probes          : number of probes
! probes_x        : coordinates of probes
! probe_buf_n     : number of time entries currently stored in buffer
! probe_buf_size  : currently allocated size of the buffer
! probe_buf_t     : times stored in buffer  (probe_buf_size)
! probe_buf_f     : probe values in buffer  (probe_buf_size, probes, nvars_probes)
! nvars_probes    : number of variables to write in probe file
! vars_probes     : variable names to write in probe file
! snaps           : number of snapshots
! snaps_t         : snapshots times
! dformat         : format specifier for writing doubles in files
!                   default : ES16.8E3 scientific notation with 8
!                   digits after coma and 3 in exponent ;
!                   use : ES25.17E3 for max precision
! bnd_warnings    : activate warnings for boundary conditions
! write_results   : activate writing of the results
! write_listing   : activate writing of the listing
! nrecord_probes  : frequency of probe and mass balance recording (0 = every time step)
! write_probes    : flag set by output_trigger to control probe/mb recording
! mass_balance    : activate mass balance boundary flux output file
! mb_option       : 0 = write FH, FQ, FQSED per unit width [m^2/s]
!                   1 = write volumic FH and FQSED only (multiplied by L) [m^3/s]
! mb_file         : mass balance file name
! mb_buf_n        : number of time entries stored in mass balance buffer
! mb_buf_size     : allocated size of mass balance buffer
! mb_buf_t        : times stored in buffer
! mb_buf_fh/fq/fqsed : boundary fluxes in buffer (size, 2): col1=left, col2=right
! mb_fh/fq/fqsed_l/r : current time step boundary fluxes (set in module_bc)
! -----------------------------------------------------------------------------
integer                                 :: ndisplay = 10
integer                                 :: nrecord = 10
integer                                 :: irec
character(len=255)                      :: ficdon
character(len=255)                      :: res_file
character(len=255)                      :: pbr_file
character(len=255)                      :: mesh_file
logical                                 :: use_mesh_file
integer                                 :: nvars_outputs
character(len=3), dimension(100)        :: vars_outputs
integer                                 :: probes
real(kind=dp), pointer                  :: probes_x (:)
integer                                 :: probe_buf_n
integer                                 :: probe_buf_size
real(kind=dp), pointer                  :: probe_buf_t  (:)
real(kind=dp), pointer                  :: probe_buf_f  (:,:,:)
integer                                 :: nvars_probes
character(len=3), dimension(100)        :: vars_probes
logical                                 :: mass_balance
integer                                 :: mb_option
integer                                 :: snaps
integer                                 :: snaps_i
real(kind=dp), pointer                  :: snaps_t (:)
character(len=8), parameter             :: dformat = 'ES16.8E3'
logical, parameter                      :: bnd_warnings = .False.
logical                                 :: write_results
logical                                 :: write_listing
integer                                 :: nrecord_probes
logical                                 :: write_probes
character(len=255)                      :: mb_file
integer                                 :: mb_buf_n
integer                                 :: mb_buf_size
real(kind=dp), pointer                  :: mb_buf_t     (:)
real(kind=dp), pointer                  :: mb_buf_fh    (:,:)
real(kind=dp), pointer                  :: mb_buf_fq    (:,:)
real(kind=dp), pointer                  :: mb_buf_fqsed (:,:)
real(kind=dp)                           :: mb_fh_l    = 0d0
real(kind=dp)                           :: mb_fq_l    = 0d0
real(kind=dp)                           :: mb_fqsed_l = 0d0
real(kind=dp)                           :: mb_fh_r    = 0d0
real(kind=dp)                           :: mb_fq_r    = 0d0
real(kind=dp)                           :: mb_fqsed_r = 0d0

! -----------------------------------------------------------------------------
! numerical parameters
! ~~~~~~~~~~~~~~~~~~~~
! method          : scheme for advection (Riemann solver)
! time_scheme     : time scheme
! ordre           : order for spatial reconstructions (MUSCL scheme)
! ordre_trac      : order for spatial reconstructions of tracers
! -----------------------------------------------------------------------------
integer                                 :: method
integer                                 :: time_scheme
integer                                 :: ordre
integer                                 :: ordre_trac

! -----------------------------------------------------------------------------
! second order reconstructions
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! corr_bp         : zb rec.: zb(i+1/2)- = zb(i)  + corr_bp(i)
! corr_bm         : zb rec.: zb(i+1/2)+ = zb(i+1)+ corr_bm(i+1)
! corr_hp         : h rec. : h(i+1/2)-  = h(i)   + corr_hp(i)
! corr_hm         : h rec. : h(i+1/2)+  = h(i+1) + corr_hm(i+1)
! corr_up         : u rec. : h(i+1/2)-  = u(i)   + corr_up(i)
! corr_um         : u rec. : h(i+1/2)+  = u(i+1) + corr_um(i+1)
! corr_tp         : t rec. : h(i+1/2)-  = t(i)   + corr_tp(i)
! corr_tm         : t rec. : h(i+1/2)+  = t(i+1) + corr_tm(i+1)
! -----------------------------------------------------------------------------
real(kind=dp), pointer                  :: corr_bp   (:)
real(kind=dp), pointer                  :: corr_bm   (:)
real(kind=dp), pointer                  :: corr_hp   (:)
real(kind=dp), pointer                  :: corr_hm   (:)
real(kind=dp), pointer                  :: corr_up   (:)
real(kind=dp), pointer                  :: corr_um   (:)
real(kind=dp), pointer                  :: corr_tp   (:,:)
real(kind=dp), pointer                  :: corr_tm   (:,:)

! -----------------------------------------------------------------------------
! initial conditions 
! ~~~~~~~~~~~~~~~~~~
! ini_file        : initial condition file name
! use_ini_file    : activate initial condition file
! type_init       : type of initial condition 
! const_h         : constant water depth initial condition
! const_e         : constant elevation initial condition
! const_q         : constant discharge initial condition
! -----------------------------------------------------------------------------
character(len=255)                      :: ini_file
logical                                 :: use_ini_file
integer                                 :: type_init
real(kind=dp)                           :: const_h
real(kind=dp)                           :: const_e
real(kind=dp)                           :: const_q

! -----------------------------------------------------------------------------
! boundary conditions
! ~~~~~~~~~~~~~~~~~~~
! clqoru          : u imposed(0) or q imposed(1) 
! bctype_l        : type of boundary condition at left boundary
! bctype_r        : type of boundary condition at right boundary
! ha              : left imposed water depth
! ua              : left imposed velocity
! qa              : left imposed discharge
! hb              : right imposed water depth
! ub              : right imposed velocity
! qb              : right imposed discharge
! ha0, hb0, qa0, qb0 : previous time step boundary values (for steady solver)
! use_bcfile_l    : use boundary condition file for left boundary
! use_bcfile_r    : use boundary condition file for right boundary
! bcfile_l        : boundary condition file name for left boundary
! bcfile_r        : boundary condition file name for left boundary
! time_bcl        : array of times read from left boundary file
! time_bcr        : array of times read from right boundary file
! field_bcl       : array of values read from left boundary file
! field_bcr       : array of values read from right boundary file
! nvars_bcl       : number of variables in left boundary condition file
! nvars_bcr       : number of variables in right boundary condition file
! vars_bcl        : name of variables in left boundary condition file
! vars_bcr        : name of variables in right boundary condition file
! bc_idx_l        : cached last interpolation interval index for left bc file
! bc_idx_r        : cached last interpolation interval index for right bc file
! bnd_extrapolation : extrapolation of geometry on ghost cells:
!                   if = 0: zb, zf, Lx equal to nearest
!                   if = 1: zb, zf, Lx preserve gradient
! -----------------------------------------------------------------------------
integer                                 :: clhorz
integer                                 :: clqoru
integer                                 :: bctype_l
integer                                 :: bctype_r
real(kind=dp)                           :: zsa
real(kind=dp)                           :: ha
real(kind=dp)                           :: ua
real(kind=dp)                           :: qa
real(kind=dp)                           :: zsb
real(kind=dp)                           :: hb
real(kind=dp)                           :: ub
real(kind=dp)                           :: qb
real(kind=dp)                           :: ha0, hb0, qa0, qb0
logical                                 :: use_bcfile_l
logical                                 :: use_bcfile_r
character(len=255)                      :: bcfile_l
character(len=255)                      :: bcfile_r
real(kind=dp), pointer                  :: time_bcl (:)
real(kind=dp), pointer                  :: time_bcr (:)
real(kind=dp), pointer                  :: field_bcl (:,:)
real(kind=dp), pointer                  :: field_bcr (:,:)
integer                                 :: nvars_bcl
integer                                 :: nvars_bcr
character(len=3), dimension(13)         :: vars_bcl
character(len=3), dimension(13)         :: vars_bcr
integer                                 :: bc_idx_l = 1
integer                                 :: bc_idx_r = 1
integer                                 :: bc_opt
integer                                 :: bnd_extrapolation
integer                                 :: bnd_flusrc
integer                                 :: bnd_highord

! -----------------------------------------------------------------------------
! tracers
! ~~~~~~~
! ntrac           : number of tracers
! t               : tracer values
! ht              : conservative variable for tracers
! flut_l/r        : ht fluxes on left and right side of interfaces
! flut_src_l/r    : ht source fluxes
! const_t         : initial tracer values
! bctypet_l       : type of boundary condition at left boundary
! bctypet_r       : type of boundary condition at right boundary
! ta              : left imposed tracer value
! tb              : right imposed tracer value
! trac_names      : precomputed tracer labels as 4-char strings (0001, ...)
! diffusion_t     : activate diffusion of tracers
! nu_t            : tracer diffusivity
! -----------------------------------------------------------------------------
integer                                 :: ntrac
real(kind=dp), pointer                  :: t         (:,:)
real(kind=dp), pointer                  :: ht        (:,:)
real(kind=dp), pointer                  :: flut_l    (:,:)
real(kind=dp), pointer                  :: flut_r    (:,:)
real(kind=dp), pointer                  :: flut_src_l(:,:)
real(kind=dp), pointer                  :: flut_src_r(:,:)
real(kind=dp), pointer                  :: const_t   (:)
integer,       pointer                  :: bctypet_l (:)
integer,       pointer                  :: bctypet_r (:)
real(kind=dp), pointer                  :: ta        (:)
real(kind=dp), pointer                  :: tb        (:)
character(len=4), pointer               :: trac_names(:)
logical                                 :: diffusion_t
real(kind=dp), pointer                  :: nu_t      (:)

! -----------------------------------------------------------------------------
! Sediments
! ~~~~~~~~~
! suspension        : activate suspension
! bedload           : activate bedload
! exner             : compute evolution of zb (bedload or suspension)
! zf_default        : default value for non erodible bed elevation
! zf                : non erodible bed elevation
!
! bedload
! ~~~~~~~
! qsed              : sediment solid flux
! taub              : bottom shear stress
! taus              : shields number
! taus_face         : shields number on cell interfaces
! fluqsed_l/r       : zb numerical fluxes
! fluqsed_src_l/r   : source fluxes for zb
! bctypeqs_l        : left bedload boundary condition type
! bctypeqs_r        : left bedload boundary condition type
! qsa               : imposed solid flux on left boundary
! qsb               : imposed solid flux on right boundary
! Jeqa              : imposed equilibrium slope on left boundary
! Jeqb              : imposed equilibrium slope on right boundary
! bedload_formula   : bedload formula
! bedload_coef      : bedload calibration coefficient
! porosite          : bed layer porosity
! rhos              : sediment density
! d50               : sediment median diameter
! Ag                : Ag coef for Grass law
! mg                : mg coef for Grass law
! tau_c             : critical shields for MPM law
! smpm_a            : standard MPM coefficient a
! smpm_b            : standard MPM coefficient b
!
! suspension
! ~~~~~~~~~~
! fluz_susp         : erosion-deposition flux for zb i.e. E-D
! Ceq_susp          : equilibrium concentration
! suspension_formula: formula for Ceq in E-D source term
! suspension_coef   : calibration coefficient for suspension
! ws_susp           : settling velocity
! susp_trac_scheme  : time scheme for tracer erosion-deposition

! slope effect
! ~~~~~~~~~~~~
! diffusion_z       : activate slope induced diffusion of zb
! nu_z              : slope diffusion coefficient
! slope_cr          : critical slope threshold 
! -----------------------------------------------------------------------------
logical                                 :: suspension
logical                                 :: bedload
logical                                 :: exner
real(kind=dp), parameter                :: zf_default = -0.1d0
real(kind=dp), pointer                  :: zf           (:)
! ~~> Bedload parameters
real(kind=dp), pointer                  :: qsed         (:)
real(kind=dp), pointer                  :: taub         (:)
real(kind=dp), pointer                  :: taus         (:)
real(kind=dp), pointer                  :: fluqsed_l    (:)
real(kind=dp), pointer                  :: fluqsed_r    (:)
real(kind=dp), pointer                  :: fluqsed_src_l(:)
real(kind=dp), pointer                  :: fluqsed_src_r(:)
integer                                 :: bctypeqs_l
integer                                 :: bctypeqs_r
real(kind=dp)                           :: qsa
real(kind=dp)                           :: qsb
real(kind=dp)                           :: Jeqa
real(kind=dp)                           :: Jeqb
integer                                 :: bedload_formula
real(kind=dp)                           :: bedload_coef
real(kind=dp)                           :: bed_fric_coef
real(kind=dp)                           :: porosite
real(kind=dp)                           :: rhos
real(kind=dp)                           :: tau_adim
real(kind=dp)                           :: d50
real(kind=dp)                           :: dm
real(kind=dp)                           :: d16
real(kind=dp)                           :: d84
real(kind=dp)                           :: Ag
real(kind=dp)                           :: mg
real(kind=dp)                           :: tau_c
real(kind=dp)                           :: smpm_a
real(kind=dp)                           :: smpm_b
! ~~> Suspension parameters
real(kind=dp), pointer                  :: flut_susp    (:)
real(kind=dp), pointer                  :: fluz_susp    (:)
real(kind=dp), pointer                  :: Ceq_susp     (:)
real(kind=dp), pointer                  :: Csratio      (:)
real(kind=dp)                           :: suspension_formula
real(kind=dp)                           :: suspension_coef
real(kind=dp)                           :: ws_susp
logical                                 :: susp_critdiff
logical                                 :: rouse_profile
integer, parameter                      :: rouse_profile_opt = 2
real(kind=dp)                           :: susp_nu_u
integer                                 :: susp_taus_scheme
integer                                 :: susp_trac_scheme
! ~~> Slope effect
logical                                 :: diffusion_z
real(kind=dp)                           :: nu_z
real(kind=dp)                           :: slope_cr

! -----------------------------------------------------------------------------
! Advanced numerical parameters
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! rec_method      : reconstruction of gradient(1) or variables(2)
! rec_vel         : reconstruction of velocity instead of discharge
! muscl_ilim_h    : MUSCL slope limitor on h
! muscl_ilim_u    : MUSCL slope limitor on u
! muscl_ilim_t    : MUSCL slope limitor on t
! hydrostarec     : scheme for bathy source term: Sb=-gh.grad(zb)
! frict_scheme    : friction scheme, for: SJ=-ghJ
! width_scheme    : width source term scheme
! eps             : minimum value of water depth, else dry cell
! eps_num         : to avoid division by zero
! std_precision   : precision for steady solver
! std_max_iter    : max iterations for steady solver
! std_critical_fr : critical h clipping for steady solver critical zones
! std_critical_opt: option for critical flow in steady solver
! 
! Bedload numerical parameters
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! swe_roots_approx  : approximation for the SVE wave speeds
! qs_derivatives    : approximation for dqs/dh and dqs/dq
! qs_derivatives_dh : increment dh for approx derivatives
! qs_derivatives_du : increment du for approx derivatives
! qs_lim            : use limitors for sediment flux
! eps_b_damping     : use damping region for sediment flux
! eps_b             : damping region thickness
! -----------------------------------------------------------------------------
logical                                 :: advanced_params
integer                                 :: rec_method
logical                                 :: rec_vel = .true.
integer                                 :: muscl_ilim_h
integer                                 :: muscl_ilim_u
integer                                 :: muscl_ilim_t
integer                                 :: hydrostarec
integer                                 :: frict_scheme
integer                                 :: width_scheme
integer                                 :: std_scheme_opt
real(kind=dp), parameter                :: eps = 1.d-6
real(kind=dp), parameter                :: eps_num = 1.d-16
real(kind=dp)                           :: std_precision
integer                                 :: std_max_iter
real(kind=dp)                           :: std_critical_fr
integer                                 :: std_critical_opt
logical, parameter                      :: std_newton_debug = .False.

! ~~> Bedload numerical parameters
integer                                 :: swe_scheme_opt
integer                                 :: swe_roots_approx
integer                                 :: qs_derivatives
real(kind=dp)                           :: qs_derivatives_dh
real(kind=dp)                           :: qs_derivatives_du
integer                                 :: qs_lim
integer                                 :: eps_b_damping
real(kind=dp)                           :: eps_b
! -----------------------------------------------------------------------------
end module
! *****************************************************************************