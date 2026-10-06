! *****************************************************************************
! 
! -----------------------------------------------------------------------------
!  Project name: Arrakis
!  Description: 1D Finite Volume methods for the Saint-Venant-Exner equations
!  Copyright (C) 2023 Fabien Souille
! -----------------------------------------------------------------------------
!   
!  This program is free: you can redistribute it and/or modify it
!  under the terms of the GNU General Public License published by the 
!  Free Software Foundation, either version 3 of the license, 
!  or (at your option) any later version.
!   
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   
!  See the GNU General Public License for details.
!   
!  You should have received a copy of the GNU General Public License
!  with this program. If not, see <https://www.gnu.org/licenses/>.
! 
! *****************************************************************************
program main
use module_data
use module_alloc,   only: allocation, deallocation
use module_io,      only: read_inputs, display, write_result, output_trigger
use module_checks,  only: check_inputs
use module_init,    only: init_geometry, init_variables
use module_forward, only: forward
! -----------------------------------------------------------------------------
implicit none
! -----------------------------------------------------------------------------

write(*,*) ''
write(*,*) '*****************************************************************'
write(*,*) '*                          Arrakis                              *'
write(*,*) '*****************************************************************'
write(*,*) ''

! -----------------------------------------------------------------------------
! Initialisation
! -----------------------------------------------------------------------------
time = tstart ; when = 0
call read_inputs
call init_geometry
call allocation
call init_variables
call check_inputs
call write_result
! -----------------------------------------------------------------------------

write(*,*) ''
write(*,*) '*****************************************************************'
write(*,*) '*                     starting iterations ...                   *'
write(*,*) '*****************************************************************'
write(*,*) ''

! -----------------------------------------------------------------------------
! Temporal iterations
! -----------------------------------------------------------------------------
niter = 0 ; when = 1
do while (niter<nitermax .and. time<time_max)
  niter = niter+1
  call output_trigger
  call forward
  time = time + dt
  call display
  call write_result
enddo
! -----------------------------------------------------------------------------

write(*,*) ''
write(*,*) '*****************************************************************'
write(*,*) '*                        Finalisation                           *'
write(*,*) '*****************************************************************'
write(*,*) ''

! -----------------------------------------------------------------------------
! Finalisation 
! -----------------------------------------------------------------------------
when = 2
print*, 'niter_fin = ', niter
call write_result
call deallocation
! -----------------------------------------------------------------------------

write(*,*) ''
write(*,*) '*****************************************************************'
write(*,*) '*                            End                                *'
write(*,*) '*****************************************************************'
write(*,*) ''

end program
! *****************************************************************************
