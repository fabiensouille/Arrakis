! *****************************************************************************
module module_checks
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Input / Output module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                                 INPUT CHECKS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine check_inputs()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - Checks values of inputs parameters if specified by user
! - Set default values of inputs parameters if not specified
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use OMP_LIB
use module_data
implicit none
logical                    :: stop_trigger=.FALSE.
integer                    :: k
! -----------------------------------------------------------------------------
! Check time parameters
! -----------------------------------------------------------------------------
! DURATION
if(time_max.le.0d0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "DURATION" SHOULD BE > 0. '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! MAXIMAL NUMBER OF ITERATIONS
if(nitermax.le.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "MAXIMAL NUMBER OF ITERATIONS" SHOULD BE >= 1 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! TIME STEP
if(dt_const.le.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "TIME STEP" SHOULD BE > 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! VARIABLE TIME STEP
if((var_dt.eqv..False.).and.(time_scheme.gt.0))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : "VARIABLE TIME STEP" = False ' 
  write(*,*) ' CFL CONDITION IS DISCARDED, YOU MAY ENCOUNTER INSTABILITIES' 
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
! CFL NUMBER
if(cfl.gt.1)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : "CFL NUMBER" SHOULD BE < 1 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
if(cfl.le.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "CFL NUMBER" SHOULD BE > 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! PARALLEL EXECUTION
if (use_omp.eqv..True.) then
  !$OMP parallel
  if (omp_get_thread_num()==0) then
  write(*,*) 'Running Arrakis in parallel mode (OpenMP) '
  write(*,*) 'Number of threads = ', omp_get_num_threads()
  end if
  !$OMP end parallel
else
  call omp_set_num_threads(1)
  write(*,*) 'Running Arrakis in sequential mode '
endif
! -----------------------------------------------------------------------------
! Check output parameters
! -----------------------------------------------------------------------------
! RESULT FILE
if(trim(res_file).eq.'')then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "RESULT FILE" NAME NOT VALID '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! LISTING FREQUENCY
if(ndisplay.lt.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "LISTING FREQUENCY" SHOULD BE >= 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! OUTPUT FREQUENCY
if(nrecord.lt.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "OUTPUT FREQUENCY" SHOULD BE >= 0  '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! VARIABLES
! TODO add list of possible names and error if unrecognized
! PROBES
if(probes.lt.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "PROBES" SHOULD BE >= 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! PROBES COORDINATES
do k=1,probes
  if((probes_x(k).lt.xa).or.(probes_x(k).gt.xb))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "PROBES COORDINATES SHOULD BE INSIDE DOMAIN        '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
enddo
! PROBES VARIABLES vars_probes
! TODO add list of possible names and error if unrecognized
! PROBES FILE NAME
if(trim(pbr_file).eq.'')then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "PROBES FILE NAME" NOT VALID '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! SNAPSHOTS
if(snaps.lt.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "SNAPSHOTS" SHOULD BE >= 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! SNAPSHOTS TIMES
do k=1,snaps
  if((snaps_t(k).lt.0.).or.(snaps_t(k).gt.time_max))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SNAPSHOTS TIMES" SHOULD BE 0. < . < DURATION '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
enddo
! -----------------------------------------------------------------------------
! Check mesh parameters
! -----------------------------------------------------------------------------
if(nx<=0 .or. xb<=xa)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : XB SHOULD BE > XA AND NX > 0 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! MESH FILE NAME
if((use_mesh_file.eqv..True.).and.(trim(mesh_file).eq.''))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "MESH FILE NAME" NOT VALID '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! -----------------------------------------------------------------------------
! Check physical parameters
! -----------------------------------------------------------------------------
! BOTTOM FRICTION LAW
if((frict_law.lt.0).or.(frict_law.gt.10))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "BOTTOM FRICTION LAW" SHOULD BE 0<=.<=10 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! BOTTOM FRICTION COEFFICIENT
if((frict_coef.lt.0.))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "BOTTOM FRICTION COEFFICIENT" SHOULD BE > 0. '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! VELOCITY DIFFUSIVITY
if(nu_u.lt.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "VELOCITY DIFFUSIVITY" SHOULD BE >= 0. '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! -----------------------------------------------------------------------------
! Check numerical parameters
! -----------------------------------------------------------------------------
! ADVECTION SCHEME
if((method.lt.1).or.(method.gt.4))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "ADVECTION SCHEME" SHOULD BE 1<=.<=4 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! TIME SCHEME
if((time_scheme.lt.0).or.(time_scheme.gt.4))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "TIME SCHEME" SHOULD BE 0<=.<=4 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
if((time_scheme.eq.3).or.(time_scheme.eq.4)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : TIME SCHEME = 3 or 4 IS EXPERIMENTAL'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
if((time_scheme.eq.0).and.(var_dt.eqv..True.)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : STEADY KERNEL REQUIRES CONSTANT TIME STEP          '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
if((time_scheme.eq.0).and.(ordre.ne.1)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : SECOND ORDER NOT AVAILABLE FOR STEADY KERNEL       '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
if((time_scheme.eq.0).and.(ordre_trac.ne.1)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : SECOND ORDER NOT AVAILABLE FOR STEADY KERNEL       '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! SPACE ORDER
if((ordre.lt.1).or.(ordre.eq.4).or.(ordre.gt.5))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "SPACE ORDER" SHOULD BE = 1,2,3 or 5  '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
if((ordre.gt.2)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : SPACE ORDER > 2 IS EXPERIMENTAL '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
if (ordre.gt.1) then
  bnd_extrapolation = 0
endif
! SPACE ORDER FOR TRACERS
if((ordre_trac.lt.1).or.(ordre_trac.eq.4).or.(ordre_trac.gt.5))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "SPACE ORDER FOR TRACERS" SHOULD BE = 1,2,3 or 5 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
if((ordre_trac.gt.2)) then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : SPACE ORDER FOR TRACERS > 2 IS EXPERIMENTAL '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
endif
! SVE SCHEME OPTION
if((swe_scheme_opt.lt.0).or.(swe_scheme_opt.gt.4))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "SVE SCHEME OPTION" SHOULD BE 0<=.<=4 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! -----------------------------------------------------------------------------
! Check initial condition
! -----------------------------------------------------------------------------
! INITIAL FILE NAME
if((use_ini_file.eqv..True.).and.(trim(ini_file).eq.''))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "INITIAL FILE NAME" NOT VALID '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! INITIAL CONDITION
if((type_init.lt.0).or.(type_init.gt.2))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "INITIAL CONDITION" SHOULD BE 0<=.<=2 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! CONSTANT WATER DEPTH
if(const_h.le.0d0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "CONSTANT WATER DEPTH" SHOULD BE > 0.'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
!CONSTANT ELEVATION no checks
!CONSTANT FLOWRATE no checks
! -----------------------------------------------------------------------------
! Check boundary conditions
! -----------------------------------------------------------------------------
! VELOCITY OR FLOWRATE FORMULATION
if((clqoru.ne.0).and.(clqoru.ne.1))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "VELOCITY OR FLOWRATE FORMULATION" SHOULD BE 0 OR 1'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! LEFT BOUNDARY CONDITION
if((bctype_l.lt.1).or.(bctype_l.gt.7))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "LEFT BOUNDARY CONDITION" SHOULD BE 1<=.<=7 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! LEFT DEPTH IMPOSED ha
if(ha.lt.0d0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "LEFT DEPTH IMPOSED" SHOULD BE >= 0.'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! LEFT VELOCITY IMPOSED no checks
! LEFT FLOWRATE IMPOSED no checks
! LEFT BOUNDARY FILE NAME bcfile_l
if((use_bcfile_l.eqv..True.).and.(trim(bcfile_l).eq.''))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "LEFT BOUNDARY FILE NAME" NOT VALID '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! RIGHT BOUNDARY CONDITION
if((bctype_r.lt.1).or.(bctype_r.gt.7))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "RIGHT BOUNDARY CONDITION" SHOULD BE 1<=.<=7 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! RIGHT DEPTH IMPOSED
if(hb.lt.0d0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "RIGHT DEPTH IMPOSED" SHOULD BE >= 0.'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! RIGHT VELOCITY IMPOSED no checks
! RIGHT FLOWRATE IMPOSED no checks
! RIGHT BOUNDARY FILE NAME bcfile_r
if((use_bcfile_r.eqv..True.).and.(trim(bcfile_r).eq.''))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "LEFT BOUNDARY FILE NAME" NOT VALID '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! BOUNDARY CONDITION OPTION
if((bc_opt.lt.0).or.(bc_opt.gt.3))then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "BOUNDARY CONDITION OPTION" SHOULD BE 0<=.<=3 '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
! -----------------------------------------------------------------------------
! Check tracer parameters
! -----------------------------------------------------------------------------
! NUMBER OF TRACERS
if(ntrac.lt.0)then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR : "NUMBER OF TRACERS" SHOULD BE >= 0'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop_trigger = .TRUE.
endif
do k=1,ntrac
  ! INITIAL TRACER VALUES no checks
  ! LEFT TRACER BOUNDARY CONDITION bctypet_l(1:ntrac)
  if((bctypet_l(k).lt.1).or.(bctypet_l(k).gt.3))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "LEFT TRACER BOUNDARY CONDITION" SHOULD BE 1<=.<=3 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! LEFT TRACER VALUE IMPOSED no checks
  ! RIGHT TRACER BOUNDARY CONDITION
  if((bctypet_r(k).lt.1).or.(bctypet_r(k).gt.3))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "RIGHT TRACER BOUNDARY CONDITION" SHOULD BE 1<=.<=3'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! RIGHT TRACER VALUE IMPOSED no checks
  ! TRACER DIFFUSIVITY
  if(nu_t(k).lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "TRACER DIFFUSIVITY" SHOULD BE >= 0. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
enddo
! -----------------------------------------------------------------------------
! Check suspension parameters
! -----------------------------------------------------------------------------
if (suspension.eqv..True.) then
  ! CHECK NUMBER OF TRACERS
  if (ntrac.ne.1) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "NUMBER OF TRACERS" SHOULD BE =1   '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif 
  ! SUSPENSION FORMULA
  if((suspension_formula.lt.1).or.(suspension_formula.gt.4))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SUSPENSION FORMULA" SHOULD BE 1<=.<=4 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SUSPENSION SETTLING VELOCITY
  if(ws_susp.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SUSPENSION SETTLING VELOCITY" SHOULD BE >= 0'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SUSPENSION CALIBRATION COEFFICIENT
  if(suspension_coef.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SUSPENSION CALIBRATION COEFFICIENT" SHOULD BE >= 0'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
endif
! -----------------------------------------------------------------------------
! Check bedload parameters
! -----------------------------------------------------------------------------
if (bedload.eqv..True.) then
  ! CHECK SVE SCHEME
  if (swe_scheme_opt.eq.0) then
    ! TODO
    if((method.lt.1).or.(method.gt.3))then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : "ADVECTION SCHEME" SHOULD BE 1<=.<=3'
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop_trigger = .TRUE.
    endif
    if((method.ne.1).and.(method.ne.3))then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' WARNING : ADVISED ADVECTION SCHEME FOR SVE ARE: ROE OR ACU'            
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    endif
  endif
  ! LEFT BEDLOAD BOUNDARY CONDITION
  if((bctypeqs_l.lt.1).or.(bctypeqs_l.gt.5))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "LEFT BEDLOAD BOUNDARY CONDITION" SHOULD BE 1<=.<=5'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! LEFT SOLID FLUX no checks
  ! LEFT EQUILIBRIUM SLOPE no checks
  ! RIGHT BEDLOAD BOUNDARY CONDITION
  if((bctypeqs_r.lt.1).or.(bctypeqs_r.gt.5))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "RIGHT BEDLOAD BOUNDARY CONDITION" SHOULD BE 1<=.<=5'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! RIGHT SOLID FLUX no checks
  ! RIGHT EQUILIBRIUM SLOPE no checks
  ! BEDLOAD FORMULA
  if((bedload_formula.lt.1).or.(bedload_formula.gt.9))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "BEDLOAD FORMULA" SHOULD BE 1<=.<=9 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if (bedload_formula.gt.2) then
    ! if not Grass, require Tau and friction law 2 or 3
    if ((frict_law.eq.0).or.(frict_law.eq.1)) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : FRICTION LAW INCOMPATIBLE WITH SVE '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop_trigger = .TRUE.
    endif
  endif
  ! BEDLOAD CALIBRATION COEFFICIENT
  if(bedload_coef.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "BEDLOAD CALIBRATION COEFFICIENT" SHOULD BE >= 0. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SEDIMENT POROSITY
  if(porosite.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SEDIMENT POROSITY" SHOULD BE >= 0.'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SEDIMENT DENSITY
  if(rhos.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SEDIMENT DENSITY" SHOULD BE >= 0. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SEDIMENT D50
  if(d50.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SEDIMENT D50" SHOULD BE >= 0. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif  ! SEDIMENT D16
  if(d16.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SEDIMENT D16" SHOULD BE >= 0. (0. = auto-compute from d50)'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SEDIMENT DM
  if(dm.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SEDIMENT DM" SHOULD BE >= 0. (0. = auto-compute from d50)'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SEDIMENT D84
  if(d84.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SEDIMENT D84" SHOULD BE >= 0. (0. = auto-compute from d50)'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif  
  ! GRASS AG
  if(Ag.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "GRASS AG" SHOULD BE >= 0. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! GRASS MG 
  if((mg.lt.1).or.(mg.gt.4)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "GRASS MG" SHOULD BE 1<=.<=4 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! MPM CRITICAL SHIELDS
  if(tau_c.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "MPM CRITICAL SHIELDS" SHOULD BE >= 0.'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! STD MPM COEF A
  if(smpm_a.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "STD MPM COEF A" SHOULD BE >= 0.'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! STD MPM COEF B
  if(smpm_b.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "STD MPM COEF B" SHOULD BE >= 0.'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
endif ! end if bedload
! -----------------------------------------------------------------------------
! Default advanced numercial parameters
! -----------------------------------------------------------------------------
! set default values if locked
if (advanced_params.eqv..False.) then
  ! RECONSTRUCTION METHOD
  rec_method = 2 ! force reconstruction of variables
  ! H FLUX LIMITER, U FLUX LIMITER, V FLUX LIMITER
  muscl_ilim_h = 1 ! force minmod
  muscl_ilim_u = 1
  muscl_ilim_t = 1
  ! HYDROSTATIC RECONSTRUCTION
  if (method.eq.3) then
    hydrostarec = 0 ! ACU: no rec
  else
    hydrostarec = 1 ! Other schemes
  endif
  if ((method.eq.1).and.(bedload.eqv..True.)) then
    hydrostarec = 4 ! ROE-SVE force hydrostarec=4
  endif
  ! BOTTOM FRICTION SCHEME
  frict_scheme = 3 ! implicit
  ! VARIABLE WIDTH SCHEME
  if ((method.eq.1).and.(bedload.eqv..True.)) then
    width_scheme = 5
  else
    width_scheme = 2
  endif
  ! SVE WAVES APPROXIMATION
  if (method.eq.1) then
    swe_roots_approx = 1 ! Cardan for ROE-SVE
  else if (method.eq.2) then
    swe_roots_approx = 2 ! SFZ for HLL-SVE
  else if (method.eq.3) then
    swe_roots_approx = 3 ! Nickalls for ACU-SVE
  endif
  ! BEDLOAD FLUX DERIVATIVES
  if ((qs_derivatives.eq.0).and.    & 
     ((bedload_formula.eq.5).or.    &
      (bedload_formula.eq.6).or.    &
      (bedload_formula.eq.7))) then
    qs_derivatives = 1
  endif
  ! if not Grass, require Tau and friction law 2 or 3
  if ((frict_law.ne.0).and. &
      (frict_law.ne.2).and. &
      (frict_law.ne.3).and. &
      (bedload_formula.ne.1)) then
    qs_derivatives = 1
  endif
  ! NON ERODIBLE BED FLUX LIMITER
  qs_lim = 1
  ! NON ERODIBLE BED DAMPING REGION
  eps_b_damping = 1
  ! NON ERODIBLE BED DAMPING REGION DEPTH
  !eps_b = 0.1 -> defined in init_variables

! -----------------------------------------------------------------------------
! Checks advanced numercial parameters
! -----------------------------------------------------------------------------
! unlocked advanced params
else
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : you have unlocked advanced numerical parameters  '
  write(*,*) '           make sure you know what you are doing !          '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  ! RECONSTRUCTION METHOD
  if((rec_method.lt.1).or.(rec_method.gt.2))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "RECONSTRUCTION METHOD" SHOULD BE 1 or 2 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! H FLUX LIMITER
  if((muscl_ilim_h.lt.1).or.(muscl_ilim_h.gt.6))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "H FLUX LIMITER" SHOULD BE 1<=.<=6 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! U FLUX LIMITER
  if((muscl_ilim_u.lt.1).or.(muscl_ilim_u.gt.6))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "U FLUX LIMITER" SHOULD BE 1<=.<=6 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! T FLUX LIMITER
  if((muscl_ilim_t.lt.1).or.(muscl_ilim_t.gt.6))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "T FLUX LIMITER" SHOULD BE 1<=.<=6 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! HYDROSTATIC RECONSTRUCTION
  if((hydrostarec.lt.0).or.(hydrostarec.gt.5))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "HYDROSTATIC RECONSTRUCTION" SHOULD BE 0<=.<=5 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if (hydrostarec.eq.3) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : Centered scheme for bottom is not well-balanced '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  if ((method.eq.3).and.(hydrostarec.ne.0)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : ACU scheme use its own reconstruction, '
    write(*,*) ' hydrostatic reconstruction forced to 0.    '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    hydrostarec = 0
  endif
  if ((method.eq.1).and.(hydrostarec.lt.4).and.(bedload.eqv..True.)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : ROE-SVE scheme is supposed to work '
    write(*,*) ' with either  HYDROSTATIC RECONSTRUCTION = 4 or 5 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  if ((method.eq.2).and.(hydrostarec.ne.1).and.(bedload.eqv..True.)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : HLL-SVE scheme is supposed to work '
    write(*,*) ' with HYDROSTATIC RECONSTRUCTION = 1 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! BOTTOM FRICTION SCHEME
  if((frict_scheme.lt.1).or.(frict_scheme.gt.3))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "BOTTOM FRICTION SCHEME" SHOULD BE 1<=.<=3 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! VARIABLE WIDTH SCHEME
  if((width_scheme.lt.1).or.(width_scheme.gt.5))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "VARIABLE WIDTH SCHEME" SHOULD BE 1<=.<=5 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if((width_scheme.eq.1).or.(width_scheme.eq.3) &
                        .or.(width_scheme.eq.4)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : VARIABLE WIDTH SCHEME = 1, 3 or 4 IS EXPERIMENTAL'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  if (bedload.eqv..True.) then
  ! SVE WAVES APPROXIMATION
  if((swe_roots_approx.lt.1).or.(swe_roots_approx.gt.3))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SVE WAVES APPROXIMATION" SHOULD BE 1<=.<=3 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if((method.eq.1).and.(swe_roots_approx.ne.1)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : SVE WAVES APPROXIMATION = 1 IS ADVISED WITH ROE '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  elseif((method.eq.2).and.(swe_roots_approx.ne.2)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : SVE WAVES APPROXIMATION = 2 IS ADVISED WITH HLL '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  elseif((method.eq.3).and.(swe_roots_approx.ne.3)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : SVE WAVES APPROXIMATION = 3 IS ADVISED WITH ACU '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! BEDLOAD FLUX DERIVATIVES
  if((qs_derivatives.lt.0).or.(qs_derivatives.gt.4))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "BEDLOAD FLUX DERIVATIVES" SHOULD BE 0<=.<=4 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if ((qs_derivatives.eq.0).and.    & 
     ((bedload_formula.eq.5).or.    &
      (bedload_formula.eq.6).or.    &
      (bedload_formula.eq.7))) then
    qs_derivatives = 1
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : Exact bedload flux derivatives  '
    write(*,*) '           not implemented for this bedload formula. '
    write(*,*) '           Switching to finite differences. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! if not Grass, require Tau and friction law 2 or 3
  if ((frict_law.ne.0).and. &
      (frict_law.ne.2).and. &
      (frict_law.ne.3).and. &
      (bedload_formula.ne.1)) then
    qs_derivatives = 1
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : Exact bedload flux derivatives  '
    write(*,*) '           not implemented for this friction formula. '
    write(*,*) '           Switching to finite differences. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! NON ERODIBLE BED FLUX LIMITER
  if((qs_lim.ne.0).and.(qs_lim.ne.1))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "NON ERODIBLE BED FLUX LIMITER" SHOULD BE 0 or 1 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if(qs_lim.eq.0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : "NON ERODIBLE BED FLUX LIMITER" IS 0,  '
    write(*,*) '           POSSIBLE VIOLATION OF NON ERODIBLE BED '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! NON ERODIBLE BED DAMPING REGION
  if((eps_b_damping.ne.0).and.(eps_b_damping.ne.1))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "NON ERODIBLE BED DAMPING REGION" SHOULD BE 0 or 1 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  if(eps_b_damping.eq.0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : "NON ERODIBLE BED DAMPING REGION" IS 0,'
    write(*,*) ' QS CAN BE POSITIVE EVEN IF NO MOBILE SEDIMENTS '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  ! NON ERODIBLE BED DAMPING REGION DEPTH
  if(eps_b.lt.0d0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "NON ERODIBLE BED DAMPING REGION DEPTH" '
    write(*,*) '         SHOULD BE >= 0.                    '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  endif ! bedload
  if (suspension.eqv..True.) then
  ! SUSPENSION TAU SCHEME
  if((susp_taus_scheme.ne.0).and.(susp_taus_scheme.ne.1))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SUSPENSION TAU SCHEME" SHOULD BE 0 or 1 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! SUSPENSION TRACER SCHEME
  if((susp_trac_scheme.ne.0).and.(susp_trac_scheme.ne.1))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "SUSPENSION TRACER SCHEME" SHOULD BE 0 or 1 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE. 
  endif
  endif ! suspension
  if (time_scheme.eq.0) then
  ! STEADY SOLVER PRECISION
  if((std_precision.lt.0d0))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "STEADY SOLVER PRECISION" SHOULD BE >= 0.'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! STEADY SOLVER MAX ITERATIONS
  if((std_max_iter.lt.0))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "STEADY SOLVER MAX ITERATIONS" SHOULD BE >= 0.'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif 
  ! STEADY SOLVER FROUDE IN CRITICAL ZONES
  if((std_critical_fr.lt.0d0).or.(std_critical_fr.gt.1d0))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "STEADY SOLVER FROUDE IN CRITICAL ZONES" '
    write(*,*) '         SHOULD BE >= 0. AND <= 1. '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  ! STEADY SOLVER CRITICAL OPTION
  if ((std_critical_opt.lt.0).or.(std_critical_opt.gt.1))then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : "STEADY SOLVER CRITICAL OPTION" SHOULD BE 0 OR 1 '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop_trigger = .TRUE.
  endif
  endif ! steady kernel
endif
! -----------------------------------------------------------------------------
if(stop_trigger.eqv..TRUE.) then
  write(*,*) '********************************************'
  write(*,*) ' Verify keywords, see arrakis.dico          '
  write(*,*) '********************************************'
  stop 
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine check_bc_inputs(vars_bc, bctype)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Check that given data at cdl is consistent with bc type
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: clqoru
implicit none
character(len=3), dimension(:), intent(in) :: vars_bc
integer,                        intent(in) :: bctype
! -----------------------------------------------------------------------------
! Check that first column is time
if (trim(vars_bc(1)).ne.'T') then
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' ERROR READING LIQ BOUNDARY FILE, MISSING TIME T '
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! Check imposed H or ZS
if ((bctype == 2).or.(bctype == 4)) then
  if ((trim(vars_bc(2)).ne.'H').and. &
      (trim(vars_bc(2)).ne.'ZS')) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR READING LIQ BOUNDARY FILE, PROVIDE  '
    write(*,*) ' EITHER WATER DEPTH H OR ELEVATION ZS '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
endif
! Check imposed Q or U
if ((bctype == 3).or.(bctype == 4)) then
  ! Check imposed Q
  if ((clqoru.eq.1).and. &
      (trim(vars_bc(2)).ne.'Q').and. &
      (trim(vars_bc(3)).ne.'Q')) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR READING LIQ BOUNDARY FILE, PROVIDE Q '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  ! Check imposed U
  if ((clqoru.eq.0).and. &
      (trim(vars_bc(2)).ne.'U').and. &
      (trim(vars_bc(3)).ne.'U')) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR READING LIQ BOUNDARY FILE, PROVIDE U '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
!******************************************************************************
