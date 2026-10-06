! *****************************************************************************
module module_linalg
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Allocation of arrays
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: dp
implicit none
! -----------------------------------------------------------------------------



contains



! *****************************************************************************
subroutine invert_2x2_matrix(A, AINV, ok_flag)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute inverse of A and store it in AINV      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps_num
implicit none
! -----------------------------------------------------------------------------
real(kind=dp), dimension(2,2), intent(in)  :: A
real(kind=dp), dimension(2,2), intent(out) :: AINV
logical, intent(out)                       :: ok_flag
real(kind=dp)                              :: det
! -----------------------------------------------------------------------------
det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
if (abs(det) .le. eps_num) then
  AINV = 0.0D0
  ok_flag = .false.
  return
end if
AINV(1,1) = A(2,2) / det
AINV(1,2) =-A(1,2) / det
AINV(2,1) =-A(2,1) / det
AINV(2,2) = A(1,1) / det
ok_flag = .true.
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine invert_3x3_matrix(A, AINV, ok_flag)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute inverse of A and store it in AINV      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps_num
implicit none
! -----------------------------------------------------------------------------
real(kind=dp), dimension(3,3), intent(in)  :: A
real(kind=dp), dimension(3,3), intent(out) :: AINV
logical, intent(out)                       :: ok_flag
real(kind=dp)                              :: det
real(kind=dp), dimension(3,3)              :: cofactor
! -----------------------------------------------------------------------------
det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3)  &
    + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
if (abs(det) .le. eps_num) then
  AINV = 0.0D0
  ok_flag = .false.
  return
end if
cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
AINV = transpose(cofactor) / det
ok_flag = .true.
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine invert_matrix(n, A, AINV, ok_flag)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute inverse of A with Lapack LU decomposition
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use module_data, only: eps_num
implicit none
! -----------------------------------------------------------------------------
integer, intent(in)                        :: n
real(kind=dp), dimension(n,n), intent(in)  :: A
real(kind=dp), dimension(n,n), intent(out) :: AINV
logical, intent(out)                       :: ok_flag
integer                                    :: info
integer, dimension(n)                      :: ipiv
real(kind=dp), dimension(n)                :: work
! -----------------------------------------------------------------------------
AINV = A
call DGETRF(n, n, AINV, n, ipiv, info)
if (info .ne. 0) then
  AINV = 0.0D0
  ok_flag = .false.
  return
end if
call DGETRI(n, AINV, n, ipiv, work, N, info)
if (info .ne. 0) then
  AINV = 0.0D0
  ok_flag = .false.
  return
end if
ok_flag = .true.
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



! *****************************************************************************
subroutine matrix_product(n, A, B, C)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Compute product of nxn matrices: C = A x B 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
! -----------------------------------------------------------------------------
integer, intent(in)                        :: n
real(kind=dp), dimension(n,n), intent(in)  :: A, B
real(kind=dp), dimension(n,n), intent(out) :: C
integer                                    :: i,j,k
! -----------------------------------------------------------------------------
do i=1,n
   do j=1,n
      C(i,j) = 0d0
      do k=1,n
        C(i,j) = C(i,j) + A(i,k)*B(k,j)
      enddo
   enddo
enddo
! -----------------------------------------------------------------------------
end subroutine
! *****************************************************************************



end module
! *****************************************************************************