! *****************************************************************************
module module_io
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Input / Output module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
! *****************************************************************************
!                            INPUT PARAMETERS READER
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine read_inputs()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Read input file
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data
implicit none
character*201      :: line
logical            :: existence
character*20       :: version
character(len=255) :: cwd
character(len=255) :: tmp
! -----------------------------------------------------------------------------
! READ MAIN INPUT FILE
! -----------------------------------------------------------------------------
call getcwd(cwd)
ficdon = trim(cwd)//'/run.yml'
inquire(file=ficdon,exist=existence)
if(.not. existence) then
  write(*,*) ' ' 
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) ' WARNING : reading error' 
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '+ input data not available'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  stop
endif
! -----------------------------------------------------------------------------
open(unit=65,file=ficdon,status='unknown')
  read(65,'(a)') line
  read(65,'(a)') line
  read(65,'(a)') line
  read(65,'(a)') line
  if(index(line,'<ARRAKIS>')==0)then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING :incompatible data file'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
  if(index(line,'$2.0')/=0) then
    version = '2.0'
  else
    version = 'non'
  endif
close(65)
select case(trim(adjustl(version)))
  case('2.0')
    write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
    write(*,*) ' Input data version : ', trim(adjustl(version))
    write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
    write(*,*) ''
    call read_data()
  case default
    write(*,*) ' ' 
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' WARNING : reading error' 
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) '+ Unknown or unsupported version...'
    write(*,*) '+ Check your input data file (.yml)'
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop             
end select
! -----------------------------------------------------------------------------
! global path to files
tmp = trim(cwd)//'/RESU/'//trim(res_file)
res_file = tmp
tmp = trim(cwd)//'/'//trim(ini_file)
ini_file = tmp
tmp = trim(cwd)//'/RESU/'//trim(pbr_file)
pbr_file = tmp
tmp = trim(cwd)//'/RESU/'//trim(mb_file)
mb_file = tmp
! -----------------------------------------------------------------------------
! READ BOUNDARY CONDITION FILES      
call read_bc
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine read_data()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Read input file
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data 
implicit none
character*201      :: line
integer            :: i, err
character*200      :: s
character(len=255) :: var_line
! -----------------------------------------------------------------------------
open(unit=65,file=ficdon,status='unknown')
  s='debut'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)         
  ! ---------------------------------------------------------------------------
  s='Managment of computation'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  ! Managment of computation
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) nitermax ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) time_max ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) dt_const ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) var_dt   ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) cfl      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) use_omp  ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Managment of outputs'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  ! Managment of outputs
  ! ---------------------------------------------------------------------------
  do i=1,size(vars_outputs)
    vars_outputs(i) = ''
  enddo
  call tab55(line); read(65,*,iostat=err) res_file    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ndisplay    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) nrecord     ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,'(A)',iostat=err) var_line; if(err/=0) call lstop(s,line)
  call split_string(var_line, vars_outputs)
  call count_vars(vars_outputs, nvars_outputs)
  ! Probes
  do i=1,size(vars_outputs)
    vars_probes(i) = ''
  enddo
  call tab55(line); read(65,*,iostat=err) probes      ; if(err/=0) call lstop(s,line)
  allocate(probes_x (probes))
  call tab55(line); read(65,*,iostat=err) probes_x(1:probes) &
                                        & ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,'(A)',iostat=err) var_line; if(err/=0) call lstop(s,line)
  call split_string(var_line, vars_probes)
  call count_vars(vars_probes, nvars_probes)
  call tab55(line); read(65,*,iostat=err) pbr_file       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) nrecord_probes ; if(err/=0) call lstop(s,line)
  ! Snapshots
  call tab55(line); read(65,*,iostat=err) snaps       ; if(err/=0) call lstop(s,line)
  allocate(snaps_t (snaps))
  call tab55(line); read(65,*,iostat=err) snaps_t(1:snaps) &
                                        & ; if(err/=0) call lstop(s,line)
  snaps_i = 1
  ! Mass balance
  call tab55(line); read(65,*,iostat=err) mass_balance ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) mb_option    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) mb_file      ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Mesh'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  ! Mesh
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) nx ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) xa ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) xb ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) Lcst          ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) Rh_approx     ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) use_mesh_file ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) mesh_file     ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) variable_width; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Physical parameters'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  !  Physical parameters
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) frict_law    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) frict_coef   ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) diffusion_u  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) nu_u         ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Numerical parameters'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  !  Numerical parameters
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) method         ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) time_scheme    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ordre          ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ordre_trac     ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) swe_scheme_opt ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) std_scheme_opt ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Initial condition'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  ! Initial condition
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) use_ini_file ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ini_file     ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) type_init    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) const_h      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) const_e      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) const_q      ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Boundary conditions'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  !  Boundary conditions
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) clhorz    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) clqoru    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bctype_l  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) zsa       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ha        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ua        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qa        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) use_bcfile_l ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bcfile_l  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bctype_r  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) zsb       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) hb        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ub        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qb        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) use_bcfile_r ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bcfile_r  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bc_opt    ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  ! bounadry option
  if (bc_opt.eq.0) then
    bnd_extrapolation = 0
    bnd_flusrc = 0
    bnd_highord = 0
  else if (bc_opt.eq.1) then
    bnd_extrapolation = 0
    bnd_flusrc = 1
    bnd_highord = 0
  else if (bc_opt.eq.2) then
    bnd_extrapolation = 1
    bnd_flusrc = 1
    bnd_highord = 0
  else if (bc_opt.eq.3) then
    bnd_extrapolation = 1
    bnd_flusrc = 1
    bnd_highord = 1
  endif
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  s='Tracers'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  !  Tracers
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) ntrac     ; if(err/=0) call lstop(s,line)
  allocate(const_t    (ntrac))
  allocate(bctypet_l  (ntrac))
  allocate(bctypet_r  (ntrac))
  allocate(ta         (ntrac))
  allocate(tb         (ntrac))
  allocate(trac_names (ntrac))
  allocate(nu_t       (ntrac))
  call tab55(line); read(65,*,iostat=err) const_t(1:ntrac)  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bctypet_l(1:ntrac); if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ta(1:ntrac)       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bctypet_r(1:ntrac); if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) tb(1:ntrac)       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) diffusion_t       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) nu_t(1:ntrac)     ; if(err/=0) call lstop(s,line)
  do i=1,ntrac
    trac_names(i) = integer_to_string(i)
  enddo


  ! ---------------------------------------------------------------------------
  s='Bedload parameters'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  ! Sediments
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) porosite           ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) rhos               ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) d50                ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) d16                ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) dm                 ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) d84                ; if(err/=0) call lstop(s,line)
  ! suspension options
  call tab55(line); read(65,*,iostat=err) suspension         ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) suspension_formula ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) suspension_coef    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) ws_susp            ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) rouse_profile      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) susp_critdiff      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) susp_nu_u          ; if(err/=0) call lstop(s,line)
  ! bedload options
  call tab55(line); read(65,*,iostat=err) bedload            ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bctypeqs_l         ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qsa                ; if(err/=0) call lstop(s,line)        
  call tab55(line); read(65,*,iostat=err) Jeqa               ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bctypeqs_r         ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qsb                ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) Jeqb               ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bedload_formula    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bedload_coef       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) Ag                 ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) mg                 ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) tau_c              ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) smpm_a             ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) smpm_b             ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) diffusion_z        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) nu_z               ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) slope_cr           ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) eps_b              ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) bed_fric_coef      ; if(err/=0) call lstop(s,line)
  exner = suspension.or.bedload

  ! ---------------------------------------------------------------------------
  s='Advanced numerical parameters'
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  read(65,'(a)',iostat=err) line ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
  !  Advanced numerical parameters
  ! ---------------------------------------------------------------------------
  call tab55(line); read(65,*,iostat=err) advanced_params    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) rec_method         ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) muscl_ilim_h       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) muscl_ilim_u       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) muscl_ilim_t       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) hydrostarec        ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) frict_scheme       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) width_scheme       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) swe_roots_approx   ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qs_derivatives     ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qs_derivatives_dh  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qs_derivatives_du  ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) qs_lim             ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) eps_b_damping      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) susp_taus_scheme   ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) susp_trac_scheme   ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) std_precision      ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) std_max_iter       ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) std_critical_fr    ; if(err/=0) call lstop(s,line)
  call tab55(line); read(65,*,iostat=err) std_critical_opt   ; if(err/=0) call lstop(s,line)
  ! ---------------------------------------------------------------------------
close(65)
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************


! *****************************************************************************
subroutine lstop(c,l)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
character(*) :: c
character(*) :: l
integer n,m
n=len(trim(adjustl(c)))
m=len(trim(adjustl(l)))
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) ' ERROR : during reading of input data file. '
write(*,*) ' Check format.'
write(*,*) ' Marker : ',c(1:n)
write(*,*) ' Line   : ',l(1:m)
write(*,*) ' Stop.'
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
stop
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                           BOUNDARY CONDITION FILES
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine read_bc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! boundary conditions files
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: use_bcfile_l, bcfile_l, use_bcfile_r, bcfile_r, &
                       vars_bcl, nvars_bcl, vars_bcr, nvars_bcr, &
                       time_bcl, time_bcr, field_bcl, field_bcr, &
                       bctype_l, bctype_r
use module_checks, only: check_bc_inputs
implicit none
integer                                :: i, nlines, io
character(len=255)                     :: var_line
! -----------------------------------------------------------------------------
! initialize variables of BC files to blank
do i=1,size(vars_bcl)
  vars_bcl(i) = ''
enddo
do i=1,size(vars_bcr)
  vars_bcr(i) = ''
enddo
! read left liquid boundary file (default: BCL.liq)
if (use_bcfile_l.eqv..True.) then
  write(*,*) ''
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ' Reading left boundary condition file       '
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ''
  ! get the number of lines
  nlines = 0
  open(unit = 101, file = trim(bcfile_l) // '.liq')
  do
    read(101, *, iostat=io)
    if (io/=0) exit
    nlines = nlines + 1
  end do
  close(unit = 101)
  nlines = nlines - 3 ! remove header count
  ! read content
  open(unit = 101, file = trim(bcfile_l) // '.liq')
  read(101, *)
  read(101, '(A)') var_line
  read(101, *)
  call split_string(var_line, vars_bcl)
  call count_vars(vars_bcl, nvars_bcl) ! number of fiels to read
  ! allocate tmp tables
  allocate(time_bcl  (1:nlines))
  allocate(field_bcl (1:nlines, nvars_bcl-1))
  ! read file content
  do i = 1, nlines
    read(101, *) time_bcl(i), field_bcl(i,:)
  end do
  close(unit = 101)
  ! check that given data at cdl is consistent with bc type
  call check_bc_inputs(vars_bcl, bctype_l)
endif
! -----------------------------------------------------------------------------
! read right liquid boundary file (default: BCR.liq)
if (use_bcfile_r.eqv..True.) then
  write(*,*) ''
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ' Reading right boundary condition file      '
  write(*,*) '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'
  write(*,*) ''
  ! get the number of lines
  nlines = 0
  open(unit = 201, file = trim(bcfile_r) // '.liq')
  do
    read(201, *, iostat=io)
    if (io/=0) exit
    nlines = nlines + 1
  end do
  close(unit = 201)
  nlines = nlines - 3 ! remove header count
  ! read content
  open(unit = 201, file = trim(bcfile_r) // '.liq')
  read(201, *)
  read(201, '(A)') var_line
  read(201, *)
  call split_string(var_line, vars_bcr)
  call count_vars(vars_bcr, nvars_bcr) ! number of fiels to read
  ! allocate tmp tables
  allocate(time_bcr  (1:nlines))
  allocate(field_bcr (1:nlines, nvars_bcr-1))
  ! read file content
  do i = 1, nlines
    read(201, *) time_bcr(i), field_bcr(i,:)
  end do
  close(unit = 201)
  ! check that given data at cdl is consistent with bc type
  call check_bc_inputs(vars_bcr, bctype_r)
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
! *****************************************************************************
!                                 OUTPUTS
! *****************************************************************************
! *****************************************************************************



! *****************************************************************************
subroutine output_trigger()
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes solution into data files
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: ndisplay, dt, time, dt_const, &
                       write_listing, write_results, nrecord, &
                       nrecord_probes, write_probes
implicit none
real(kind=dp) :: dt_output, rec_time
! -----------------------------------------------------------------------------
! write condition for display
dt_output = ndisplay*dt_const
rec_time = ceiling(time/dt_output)*dt_output
if ((rec_time.ge.time) .and. &
    (rec_time.lt.(time+dt))) then
  write_listing = .True.
else
  write_listing = .False.
endif
! write condition for results
dt_output = nrecord*dt_const
rec_time = ceiling(time/dt_output)*dt_output
if ((rec_time.gt.time) .and. &
    (rec_time.le.(time+dt))) then
  write_results = .True.
else
  write_results = .False.
endif
! write condition for probes and mass balance
if (nrecord_probes == 0) then
  write_probes = .True.
else
  dt_output = nrecord_probes*dt_const
  rec_time = ceiling(time/dt_output)*dt_output
  if ((rec_time.gt.time) .and. &
      (rec_time.le.(time+dt))) then
    write_probes = .True.
  else
    write_probes = .False.
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************




! *****************************************************************************
subroutine display
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes solution into data files
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: niter, dt, time, &
                       write_listing
implicit none
integer       :: jj, hh, mm, ss
! -----------------------------------------------------------------------------
! listing output
if(write_listing.eqv..True.)then
    call convert_time_to_jjhhmmss(time, jj, hh, mm, ss)
    if (jj.lt.1) then
      write(*, '(a,e11.5,a,i9,a,i2.2,":",i2.2,":",i2.2)') &
      'dt (s) = ', dt, &
      ' / niter = ', niter, &
      ' / time (h:m:s) = ', hh, mm, ss
    else
      write(*, '(a,e11.5,a,i9,a,i2.2,":",i2.2,":",i2.2)') &
      'dt (s) = ', dt, &
      ' / niter = ', niter, &
      ' / time (j:h:m) = ', jj, hh, mm
    endif
endif
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine convert_time_to_jjhhmmss(time_in_seconds, jj, hh, mm, ss)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Converts time in seconds to days, hours, minutes, and seconds
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
real(kind=dp), intent(in)  :: time_in_seconds
integer,       intent(out) :: jj,hh, mm, ss
! -----------------------------------------------------------------------------
jj = int(time_in_seconds/86400d0) ! 1 day = 86400 seconds
hh = int(mod(time_in_seconds, 86400d0)/3600d0)
mm = int(mod(time_in_seconds, 3600d0)/60d0)
ss = int(mod(time_in_seconds, 60d0))
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine write_result
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes solutions into data files
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: when, time, dt, dt_const, time_max, tstart, irec, &
                       nx, x, Lx, res_file, vars_outputs, nvars_outputs, &
                       probes, probes_x, pbr_file, vars_probes, nvars_probes, &
                       snaps, snaps_i, snaps_t, dformat, write_results, &
                       probe_buf_n, probe_buf_size, probe_buf_t, probe_buf_f, &
                       nrecord_probes, write_probes, &
                       mass_balance, mb_file, mb_buf_n, mb_buf_size, &
                       mb_buf_t, mb_buf_fh, mb_buf_fq, mb_buf_fqsed, &
                       mb_fh_l, mb_fq_l, mb_fqsed_l, &
                       mb_fh_r, mb_fq_r, mb_fqsed_r, mb_option, exner
implicit none
real(kind=dp)      :: fi (1:nx,1:nvars_outputs)
real(kind=dp)      :: fp (1:nx,1:nvars_probes)
real(kind=dp)      :: fp_interp (1:nvars_probes)
real(kind=dp)      :: alpha
integer            :: i, j, k, l
character(len=4)   :: kchar
character(len=255) :: probe_file
logical            :: found_value
real(kind=dp), pointer :: tmp_buf_t(:)
real(kind=dp), pointer :: tmp_buf_f(:,:,:)
real(kind=dp), pointer :: tmp_mb_t    (:)
real(kind=dp), pointer :: tmp_mb_fh   (:,:)
real(kind=dp), pointer :: tmp_mb_fq   (:,:)
real(kind=dp), pointer :: tmp_mb_fqsed(:,:)
real(kind=dp)          :: mb_row(6)
! -----------------------------------------------------------------------------
! WRITE MAIN RESULT FILE
! fi -> contains the selected variables to write
! -----------------------------------------------------------------------------
! load fields to write
if((when == 0).or.(write_results.eqv..True.).or.(when == 2))then
  do l=1,nvars_outputs
    call prepare_output_var(vars_outputs(l), fi(:,l))
  enddo
endif
! -----------------------------------------------------------------------------
if(when == 0) then
  irec = 0
  call write_sol(x(1:nx), fi, trim(res_file) // "ini.dat")
endif
! -----------------------------------------------------------------------------
if(write_results.eqv..True.)then
  irec = irec + 1
  call write_sol(x(1:nx), fi, &
    & trim(res_file) // &
    & trim(integer_to_string(irec)) // ".dat")
endif
! -----------------------------------------------------------------------------
if(when == 2) then
  irec = irec + 1
  call write_sol(x(1:nx), fi, &
    & trim(res_file) // "fin.dat")
endif
! -----------------------------------------------------------------------------
! WRITE SNAPSHOTS RESULT FILE
! fi -> contains the selected variables to write
! -----------------------------------------------------------------------------
if (snaps_i.le.snaps) then
  if ( (snaps_t(snaps_i).gt.(time-dt)) .and. &
  &     (snaps_t(snaps_i).le.(time)) ) then
    ! get variables
    do l=1,nvars_outputs
      call prepare_output_var(vars_outputs(l), fi(:,l))
    enddo
    ! write
    call write_sol(x(1:nx), fi, &
      & trim(res_file) // "_SNAP" // &
      & trim(integer_to_string(snaps_i)) // ".dat")
    snaps_i = snaps_i + 1
  endif
endif
! -----------------------------------------------------------------------------
! PROBE BUFFER: collect probe data in memory, write to files only at the end
! fp -> contains the probe selected variables to write
! -----------------------------------------------------------------------------
if (probes.gt.0) then
  ! ---------------------------------------------------------------------------
  ! Initialise buffer at the very first call (when == 0)
  ! ---------------------------------------------------------------------------
  if (when == 0) then
    probe_buf_n    = 0
    if (nrecord_probes == 0) then
      probe_buf_size = min(ceiling(time_max / dt_const) + 1, 86400)
    else
      probe_buf_size = min(ceiling(time_max / (dt_const*nrecord_probes)) + 1, 86400)
    endif
    allocate(probe_buf_t(probe_buf_size))
    allocate(probe_buf_f(probe_buf_size, probes, nvars_probes))
  endif
  ! ---------------------------------------------------------------------------
  ! Collect probe values into the buffer
  ! ---------------------------------------------------------------------------
  if (write_probes) then
  ! get field variables on the grid
  do l=1,nvars_probes
    call prepare_output_var(vars_probes(l), fp(:,l))
  enddo
  ! grow buffer if full
  if (probe_buf_n >= probe_buf_size) then
    probe_buf_size = probe_buf_size * 2
    allocate(tmp_buf_t(probe_buf_size))
    allocate(tmp_buf_f(probe_buf_size, probes, nvars_probes))
    tmp_buf_t(1:probe_buf_n)       = probe_buf_t(1:probe_buf_n)
    tmp_buf_f(1:probe_buf_n,:,:)   = probe_buf_f(1:probe_buf_n,:,:)
    deallocate(probe_buf_t)
    deallocate(probe_buf_f)
    probe_buf_t => tmp_buf_t
    probe_buf_f => tmp_buf_f
  endif
  ! append current time step
  probe_buf_n = probe_buf_n + 1
  probe_buf_t(probe_buf_n) = time
  ! interpolate and store each probe
  do k=1,probes
    found_value = .False.
    do i=0,nx
      if ((x(i) <= probes_x(k)).and.(probes_x(k) < x(i+1))) then
        found_value = .True.
        alpha = (x(i+1)-probes_x(k))/(x(i+1)-x(i))
        do j=1,nvars_probes
          fp_interp(j) = alpha*fp(i,j) + (1.d0-alpha)*fp(min(i+1,nx),j)
        enddo
      endif
    enddo
    if (found_value.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR WRITING PROBE FILE, CHECK PROBE COORDINATES'
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    end if
    do j=1,nvars_probes
      probe_buf_f(probe_buf_n, k, j) = fp_interp(j)
    enddo
  enddo
  endif ! write_probes
  ! ---------------------------------------------------------------------------
  ! Flush buffer to files at the end of the computation (when == 2)
  ! ---------------------------------------------------------------------------
  if (when == 2) then
    do k=1,probes
      write(kchar,'(I4.4)') k
      probe_file = trim(pbr_file)//'_'//trim(kchar)//".dat"
      open(unit = 34, file=probe_file)
      write(34, *) '# Probe file'
      write(34, fmt='(a4, a1)', advance='no') 'Time', ','
      do j=1,nvars_probes
        if (j<nvars_probes) then
          write(34, fmt='(a, a1)', advance='no') &
                               trim(adjustl(vars_probes(j))), ','
        else
          write(34, fmt='(a)') trim(adjustl(vars_probes(j)))
        endif
      enddo
      do i=1,probe_buf_n
        write(34, fmt='('//dformat//', A)', advance='no') probe_buf_t(i), ','
        do j=1,nvars_probes
          if (j<nvars_probes) then
            write(34, fmt='('//dformat//', A)', advance='no') &
                                             probe_buf_f(i,k,j), ','
          else
            write(34, fmt='('//dformat//')') probe_buf_f(i,k,j)
          endif
        enddo
      enddo
      close(34)
    enddo
    deallocate(probe_buf_t)
    deallocate(probe_buf_f)
  endif
endif
! -----------------------------------------------------------------------------
! MASS BALANCE BUFFER: collect boundary fluxes in memory, write at the end
! Columns : Time | FH_L | FQ_L | FQSED_L | FH_R | FQ_R | FQSED_R
! FH    : h numerical flux at boundary  [m^2/s]
! FQ    : q numerical flux at boundary  [m^3/s^2]
! FQSED : qsed numerical flux           [m^2/s]  (0 if no sediment)
! _L    : left boundary,  _R : right boundary
! -----------------------------------------------------------------------------
if (mass_balance) then
  ! ---------------------------------------------------------------------------
  ! Initialise buffer at the very first call (when == 0)
  ! ---------------------------------------------------------------------------
  if (when == 0) then
    mb_buf_n    = 0
    if (nrecord_probes == 0) then
      mb_buf_size = min(ceiling(time_max / dt_const) + 1, 86400)
    else
      mb_buf_size = min(ceiling(time_max / (dt_const*nrecord_probes)) + 1, 86400)
    endif
    allocate(mb_buf_t    (mb_buf_size))
    allocate(mb_buf_fh   (mb_buf_size, 2))
    allocate(mb_buf_fq   (mb_buf_size, 2))
    allocate(mb_buf_fqsed(mb_buf_size, 2))
    mb_buf_fh   (:,:) = 0d0
    mb_buf_fq   (:,:) = 0d0
    mb_buf_fqsed(:,:) = 0d0
  endif
  ! ---------------------------------------------------------------------------
  ! Grow buffer if full (double strategy)
  ! ---------------------------------------------------------------------------
  if (write_probes) then
  if (mb_buf_n >= mb_buf_size) then
    mb_buf_size = mb_buf_size * 2
    allocate(tmp_mb_t    (mb_buf_size))
    allocate(tmp_mb_fh   (mb_buf_size, 2))
    allocate(tmp_mb_fq   (mb_buf_size, 2))
    allocate(tmp_mb_fqsed(mb_buf_size, 2))
    tmp_mb_t    (1:mb_buf_n)   = mb_buf_t    (1:mb_buf_n)
    tmp_mb_fh   (1:mb_buf_n,:) = mb_buf_fh   (1:mb_buf_n,:)
    tmp_mb_fq   (1:mb_buf_n,:) = mb_buf_fq   (1:mb_buf_n,:)
    tmp_mb_fqsed(1:mb_buf_n,:) = mb_buf_fqsed(1:mb_buf_n,:)
    deallocate(mb_buf_t)
    deallocate(mb_buf_fh)
    deallocate(mb_buf_fq)
    deallocate(mb_buf_fqsed)
    mb_buf_t     => tmp_mb_t
    mb_buf_fh    => tmp_mb_fh
    mb_buf_fq    => tmp_mb_fq
    mb_buf_fqsed => tmp_mb_fqsed
  endif
  ! ---------------------------------------------------------------------------
  ! Append current time step
  ! ---------------------------------------------------------------------------
  mb_buf_n = mb_buf_n + 1
  ! time before update (fluxes are computed before update)
  mb_buf_t    (mb_buf_n)   = time - dt
  mb_buf_fh   (mb_buf_n,1) = mb_fh_l
  mb_buf_fh   (mb_buf_n,2) = mb_fh_r
  mb_buf_fq   (mb_buf_n,1) = mb_fq_l
  mb_buf_fq   (mb_buf_n,2) = mb_fq_r
  mb_buf_fqsed(mb_buf_n,1) = mb_fqsed_l
  mb_buf_fqsed(mb_buf_n,2) = mb_fqsed_r
  endif ! write_probes
  ! ---------------------------------------------------------------------------
  ! Flush buffer to file at end of computation (when == 2)
  ! ---------------------------------------------------------------------------
  if (when == 2) then
    open(unit = 35, file=trim(mb_file)//".dat")
    if (mb_option == 0) then
      ! -----------------------------------------------------------------------
      ! option 0: per-unit-width fluxes [m^2/s]
      write(35, fmt='(a)') '# Mass balance file (option 0)'
      write(35, fmt='(a)') '# FH  : h numerical flux [m^2/s]'
      write(35, fmt='(a)') '# FQ  : q numerical flux [m^3/s^2]'
      if (exner) write(35, fmt='(a)') '# FQS : qsed numerical flux [m^2/s]'
      write(35, fmt='(a)') '# _L = left boundary, _R = right boundary'
      if (exner) then
        write(35, fmt='(a)') 'Time,FH_L,FQ_L,FQS_L,FH_R,FQ_R,FQS_R'
      else
        write(35, fmt='(a)') 'Time,FH_L,FQ_L,FH_R,FQ_R'
      endif
      do i=1,mb_buf_n
        write(35, fmt='('//dformat//', A)', advance='no') mb_buf_t(i), ','
        if (exner) then
          mb_row(1) = mb_buf_fh   (i,1)
          mb_row(2) = mb_buf_fq   (i,1)
          mb_row(3) = mb_buf_fqsed(i,1)
          mb_row(4) = mb_buf_fh   (i,2)
          mb_row(5) = mb_buf_fq   (i,2)
          mb_row(6) = mb_buf_fqsed(i,2)
          do j=1,6
            if (j<6) then
              write(35, fmt='('//dformat//', A)', advance='no') mb_row(j), ','
            else
              write(35, fmt='('//dformat//')') mb_row(j)
            endif
          enddo
        else
          mb_row(1) = mb_buf_fh(i,1)
          mb_row(2) = mb_buf_fq(i,1)
          mb_row(3) = mb_buf_fh(i,2)
          mb_row(4) = mb_buf_fq(i,2)
          do j=1,4
            if (j<4) then
              write(35, fmt='('//dformat//', A)', advance='no') mb_row(j), ','
            else
              write(35, fmt='('//dformat//')') mb_row(j)
            endif
          enddo
        endif
      enddo
    else
      ! -----------------------------------------------------------------------
      ! option 1: volumic fluxes (x channel width L) [m^3/s]
      write(35, fmt='(a)') '# Mass balance file (option 1)'
      write(35, fmt='(a)') '# QL  : volumic water flux = FH*L [m^3/s]'
      if (exner) write(35, fmt='(a)') '# QSL : volumic solid flux = FQS*L [m^3/s]'
      write(35, fmt='(a)') '# _L = left boundary (L=Lx(1)), _R = right boundary (L=Lx(nx))'
      if (exner) then
        write(35, fmt='(a)') 'Time,QL_L,QSL_L,QL_R,QSL_R'
      else
        write(35, fmt='(a)') 'Time,QL_L,QL_R'
      endif
      do i=1,mb_buf_n
        write(35, fmt='('//dformat//', A)', advance='no') mb_buf_t(i), ','
        if (exner) then
          mb_row(1) = mb_buf_fh   (i,1) * Lx(1)
          mb_row(2) = mb_buf_fqsed(i,1) * Lx(1)
          mb_row(3) = mb_buf_fh   (i,2) * Lx(nx)
          mb_row(4) = mb_buf_fqsed(i,2) * Lx(nx)
          do j=1,4
            if (j<4) then
              write(35, fmt='('//dformat//', A)', advance='no') mb_row(j), ','
            else
              write(35, fmt='('//dformat//')') mb_row(j)
            endif
          enddo
        else
          mb_row(1) = mb_buf_fh(i,1) * Lx(1)
          mb_row(2) = mb_buf_fh(i,2) * Lx(nx)
          do j=1,2
            if (j<2) then
              write(35, fmt='('//dformat//', A)', advance='no') mb_row(j), ','
            else
              write(35, fmt='('//dformat//')') mb_row(j)
            endif
          enddo
        endif
      enddo
    endif
    close(35)
    deallocate(mb_buf_t)
    deallocate(mb_buf_fh)
    deallocate(mb_buf_fq)
    deallocate(mb_buf_fqsed)
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
!******************************************************************************



!******************************************************************************
subroutine prepare_output_var(var_in, f)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! do
! -> check var name
! -> prepare the corresponding variable (compute if not primitive)
! -> load var in temporary array f for printout
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ZB : bottom elevation
!  ZS : water surface elevation
!  H  : water depth
!  U  : water depth-averaged velocity
!  Q  : water depth-averaged discharge
!  QL : water volumic discharge, QL = Q*L
!  L  : rectangular channel width
!  CI : finite volume cell size [m]
!  RH : hydraulic radius
!  FR : Froude number
!  FC : friction coefficient
!  B  : sediment depth (if exner)
!  ZF : hard bottom elevation (if exner)
!  QS : solid discharge (if exner)
!  QSL: solid volumic discharge 
!  TB : bottom shear stress
!  TS : dimensionless bottom shear stress tau* (if exner)
!  A  : wetted area, A=H*L
!  V  : water volume in each cell, V=H*L*ci=A*ci
!  AS : sediment area, AS=(zb-zf)*L
!  VS : sediment volume in each cell, VS=(zb-zf)*L*ci=AS*ci
!  Ti : tracer number i (1<=i<=99)
!  CE : suspension equilibirum concentration (if suspension)
!  WS : settling velocity (if suspension)
!  RO : Rouse number (if suspension)  
!  A1 : auxiliary array 1 (for debug)
!  A2 : auxiliary array 2 (for debug)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: nx, zb, h, u, q, zf, fric, qsed, aux1, aux2, g, eps, &
                       Rh, taub, taus, tau_adim, Lx, exner, ntrac, t, ci, &
                       suspension, Ceq_susp, ws_susp, d50, rho, bed_fric_coef, &
                       trac_names
use module_friction, only: tau, friction_velocity
use module_suspension, only: settling_velocity
implicit none
character(len=3), intent(in)    :: var_in
real(kind=dp),    intent(inout) :: f (1:nx)
real(kind=dp)                   :: ustar
integer                         :: i, k
character(len=3)                :: var
character(len=4)                :: char_trac
logical                         :: main_var_flag
logical                         :: tracer_flag
! -----------------------------------------------------------------------------
var = trim(adjustr(var_in)) ! remove trailing blank space
main_var_flag = .FALSE.
tracer_flag = .TRUE.
! -----------------------------------------------------------------------------
if (var.ne.'') then
  ! READ MAIN VARIABLES
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (trim(var).eq.' ZB') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = zb(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' ZS') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do i=1,nx
      f(i) = zb(i) + h(i)
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  H') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = h(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  U') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = u(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  Q') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = q(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  L') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = Lx(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' CI') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = ci(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' RH') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = Rh(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' FR') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do i=1,nx
      if (h(i).ge.eps) then
        f(i) = u(i)/sqrt(g*h(i))
      else
        f(i) = 0d0
      endif
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' FC') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = fric(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' ZF') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print ZF '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    f = zf(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  B') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print B '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    ! not a primitive variable, compute res
    do i=1,nx
      f(i) = max(0.d0, zb(i)- zf(i))
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' QS') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print QS '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    f = qsed(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'QSL') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print QSL '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    f = qsed(1:nx)*Lx(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' TB') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! not a primitive variable, compute res
    do i=1,nx
      taub(i) = tau(h(i), u(i), Lx(i), 0d0, fric(i)*bed_fric_coef)
    enddo
    f = taub(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' TS') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print TS '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    ! not a primitive variable, compute res
    do i=1,nx
      taub(i) = tau(h(i), u(i), Lx(i), 0d0, fric(i)*bed_fric_coef)
      taus(i) = taub(i)*tau_adim
    enddo
    f = taus(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  A') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! not a primitive variable
    f = h(1:nx)*Lx(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.'  V') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! not a primitive variable
    f = h(1:nx)*Lx(1:nx)*ci(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' QL') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! not a primitive variable
    f = q(1:nx)*Lx(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' AS') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print AS '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    ! not a primitive variable, compute res
    do i=1,nx
      f(i) = max(0.d0, zb(i)- zf(i))*Lx(i)
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' VS') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (exner.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print VS '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    ! not a primitive variable, compute res
    do i=1,nx
      f(i) = max(0.d0, zb(i)- zf(i))*Lx(i)*ci(i)
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(var).eq.' CE') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (suspension.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print CE '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    f = Ceq_susp(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(var).eq.' WS') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (suspension.eqv..False.) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' ERROR : Unable to print WS '
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      stop
    endif
    do i=1,nx
      f(i) = ws_susp
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif (trim(var).eq.' RO') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do i=1,nx
      if (ws_susp.eq.0.d0) then
        ws_susp = settling_velocity(d50)
      endif
      ustar = friction_velocity(h(i), u(i), Lx(i), fric(i)*bed_fric_coef)
      f(i) = ws_susp/(0.4d0*ustar)
    enddo
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' A1') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = aux1(1:nx)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (trim(var).eq.' A2') then
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f = aux2(1:nx)
  else
    main_var_flag = .TRUE.
  endif
  ! ~~~~~~~~~~~~
  ! READ TRACERS
  ! ~~~~~~~~~~~~
  if (ntrac.gt.0) then
    do k=1,ntrac
      char_trac = trac_names(k)
      if ((trim(var).eq." T"//char_trac(4:4)).or. &
      &    (trim(var).eq."T"//char_trac(3:4))) then
        f = t(1:nx, k)
        tracer_flag = .FALSE.
      endif
    enddo
  endif 
  ! ~~~~~
  ! ERROR
  ! ~~~~~
  if ((main_var_flag.eqv..TRUE.).and.(tracer_flag.eqv..TRUE.)) then
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(*,*) ' ERROR : Unkown variable name, see dico VARIABLES '
    write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    stop
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
!******************************************************************************



!******************************************************************************
subroutine write_sol(f0, fi, fichier)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes two vectors in a specified file. File will be owerwritten
! if existing
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: nx, time, vars_outputs, nvars_outputs, dformat 
implicit none
real(kind=dp),   intent(in) :: f0 (1:nx)
real(kind=dp),   intent(in) :: fi (1:nx,1:nvars_outputs)
character(len=*),intent(in) :: fichier
integer                     :: i, j
! -----------------------------------------------------------------------------
open(unit = 24, file = fichier)
! write header with field names
write(24, fmt='('//dformat//', A1, I0)') time, ',', nx
write(24, fmt='(a1, a1)', advance='no') 'X', ','
do j=1,nvars_outputs
  if (j<nvars_outputs) then
    write(24, fmt='(a, a1)' , advance='no') &
                         trim(adjustl(vars_outputs(j))), ','
  else
    write(24, fmt='(a)') trim(adjustl(vars_outputs(j)))
  endif
enddo
! write fields 
do i=1,nx
  write(24, fmt='('//dformat//', A)', advance='no') f0(i), ','
  do j=1,nvars_outputs
    if (j<nvars_outputs) then
      write(24, fmt='('//dformat//', A)', advance='no') fi(i,j), ','
    else
      write(24, fmt='('//dformat//')') fi(i,j)
    endif
  enddo
end do
close(24)
! -----------------------------------------------------------------------------
end subroutine
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                                     UTILS
!******************************************************************************
!******************************************************************************



!******************************************************************************
function integer_to_string(n0) result(chaine)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Convert an integer lower than 9999 into a string
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
!Input :  + no : integer, must be < 10000
!Outup :  + string
integer, intent(in) :: n0
character(len=4)    :: chaine
integer             :: n
! -----------------------------------------------------------------------------
n = n0
chaine(4:4) = achar(mod(n,10)+48)
n = n/10
chaine(3:3) = achar(mod(n,10)+48)
n = n/10
chaine(2:2) = achar(mod(n,10)+48)
n = n/10
chaine(1:1) = achar(mod(n,10)+48)
! -----------------------------------------------------------------------------
end function integer_to_string
!******************************************************************************



!******************************************************************************
subroutine tab55(line)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
character(*) :: line
read(65,'(A)') line
backspace(65)
read(65,'(t55)',advance='no')
end subroutine
!******************************************************************************



!******************************************************************************
subroutine split_string(str, split_str)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Split a header string into tokens.
! If commas are present, splits on commas (preserving empty fields, e.g.
! "ZB,H,U,,,,,," from the VARIABLES keyword).
! Otherwise splits on whitespace (multiple spaces treated as one separator).
! Both input styles are therefore accepted transparently.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
character(len=255), intent(in)              :: str
character(len=3), dimension(:), intent(out) :: split_str
integer                                     :: i, n, slen, start
logical                                     :: use_comma, in_token
! -----------------------------------------------------------------------------
slen      = len_trim(str)
n         = 1
start     = 1
split_str = ''
use_comma = (index(str(1:slen), ',') > 0)
if (use_comma) then
  ! -------------------------------------------------------------------------
  ! Comma-separated: preserves empty fields (e.g. "ZB,H,U,,,,,,")
  ! -------------------------------------------------------------------------
  do i = 1, slen
    if (str(i:i) == ',' .or. i == slen) then
      if (i == slen .and. str(i:i) /= ',') then
        split_str(n) = adjustl(str(start:i))      ! last token
      elseif (i > start) then
        split_str(n) = adjustl(str(start:i-1))    ! normal token
      else
        split_str(n) = ''                         ! empty field
      endif
      n     = n + 1
      start = i + 1
    endif
  enddo
else
  ! -------------------------------------------------------------------------
  ! Whitespace-separated: multiple spaces/tabs treated as one separator
  ! -------------------------------------------------------------------------
  in_token = .False.
  do i = 1, slen
    if (str(i:i) /= ' ' .and. str(i:i) /= char(9)) then
      if (.not. in_token) then
        start    = i
        in_token = .True.
      endif
    else
      if (in_token) then
        split_str(n) = adjustl(str(start:i-1))
        n        = n + 1
        in_token = .False.
      endif
    endif
  enddo
  if (in_token) then
    split_str(n) = adjustl(str(start:slen))       ! last token
  endif
endif
! -----------------------------------------------------------------------------
end subroutine
!******************************************************************************



!******************************************************************************
subroutine count_vars(vars, nvars)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Count the number of char of len(3) different than blank ''
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
character(len=3), dimension(:), intent(in)  :: vars
integer,                        intent(out) :: nvars
integer                                     :: i
character(len=3)                            :: var
! -----------------------------------------------------------------------------
nvars = 0
do i=1,size(vars)
  var = vars(i)
  if (len_trim(var).ne.0) then
    nvars = nvars + 1
  endif
enddo
! -----------------------------------------------------------------------------
end subroutine
!******************************************************************************



end module
!******************************************************************************