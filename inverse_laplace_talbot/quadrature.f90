!> @file quadrature.f90
!> @author Mainak
!> @brief contains all quadrature routines
module quadratures
contains
!! ---------   Variable declaration block ------------
subroutine simpson_onethird(integrand, spacing, ubound, lbound, status, value)
implicit none
double precision, dimension(:),intent(in):: integrand
double precision, intent(in)             :: spacing
double precision, intent(in)             :: lbound
double precision, intent(in)             :: ubound
double precision, intent(out)            :: value
logical, intent(out)                     :: status
!! ---------------------------------------------------
integer:: Narray
integer:: kidx
!-----------------------------------------------------
!
!-----------------------------------------------------
!! ---------------- Initialization block
   Narray = size(integrand)
   value = 0.d0
   status = .false.
! ---------- Formula to implement --------------------
! -------   \int_{a}^{b} f(x)dx = \frac{h}{3}
!!            [f(x_0) + 4\sum_{k=1}^{N/2} f(x_{2k-1})
!!            + 2\sum_{k=1}^{N/2-1} f(x_{2k}) + f(x_N)]
!!           x_0 = a; x_N = b
!! ---------------------------------------------------
!    if (abs(Narray*spacing - (ubound-lbound)) >= tolerance)then
  !    call assert()
!    end if

do kidx = 2, Narray-1
  if ((kidx - 2*(kidx/2)) == 0)then
! kidx even
    value = value + 4.d0*integrand(kidx)
  else if ((kidx - 2*(kidx/2)) == 1)then
! kidx odd
    value = value + 2.d0*integrand(kidx)
  else
! something wrong
   stop "Problem with simpson counter"
  end if
end do
   value = value + (integrand(1)+integrand(Narray))

 value = value*(spacing/ 3.d0)

end subroutine simpson_onethird
!!!!
subroutine test_simpson()
implicit none

end subroutine test_simpson
!!
end module quadratures
